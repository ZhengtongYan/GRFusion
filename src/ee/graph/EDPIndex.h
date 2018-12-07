#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include "GraphView.h"
#include <unordered_map>

using namespace std;

namespace voltdb {

class PartitionVertex;
class EDPIndex;

class PartitionEdge {
	friend class Partition;
public:
	PartitionEdge(Edge* edge, PartitionVertex* fromVertex, PartitionVertex* toVertex, int weightField, int labelField);
	~PartitionEdge(void);

	Edge* getEdge() { return edge; }

	int getID() { return eID; }
	PartitionVertex* getFrom() { return fromV; }
	PartitionVertex* getTo() { return toV; }

	double getWeight() { return weight; }
	int getLabel() { return label; }

private:
	Edge* edge;

	int eID;
	PartitionVertex* fromV;
	PartitionVertex* toV;

	double weight;
	int label;
};

class Partition;

class PartitionVertex {
	friend class Partition;
public:
	PartitionVertex(Vertex * vertex);
	~PartitionVertex(void);

	Vertex* getVertex() { return vertex; }

	int getID() { return vID; }
	bool isBridge() { return bridge; }
	bool isBridgeBackwards() { return bridgeBackwards; }
	int getComponentID() { return componentId; }
	Partition* getPartition() { return partition; }

	map<int, PartitionEdge*>* getInEdges() { return inEdges; }
	map<int, PartitionEdge*>* getOutEdges() { return outEdges; }

	map<int, Partition*>* getOtherHostsOut() { return otherHosts_out; }
	map<int, Partition*>* getOtherHostsIn() { return otherHosts_in; }

	void setPartition(Partition * p) { partition = p; }
	void setBridge(bool b) { bridge = b; }
	void setBridgeBackwards(bool b) { bridgeBackwards = b; }
	void setComponentId(int cId) { componentId = cId; }
	//Adds an edge to this vertex
	void addEdge(PartitionEdge* edge, bool outEdge);



private:
	Vertex* vertex;

	int vID;
	bool bridge;
	bool bridgeBackwards;
	int componentId;

	Partition* partition;

	map<int, PartitionEdge*> inEdges;
	map<int, PartitionEdge*> outEdges;

	map<int, Partition*> otherHosts_out;
	map<int, Partition*> otherHosts_in;
};

class MonochromePath{

public:
	MonochromePath(int src, int dest, int cost);
protected:
	int src_id;
	int dest_id;
	int cost;

	vector<PartitionEdge*> pathEdges;
};

class PathCacheEntry {
	friend class Partition;
public:
	PathCacheEntry(int label, int src, int dest, int cost, vector<PartitionEdge*>* pathEdges);
	~PathCacheEntry();

	int getLabel() { return label; }
	int getSource() { return src; }
	int getDestination() {return dest; }
	int getCost() { return cost; }

	vector<PartitionEdge*>* getPath() { return &pathEdges; }

	//Returns the size of the cache entry in bytes
	uint64_t getEntrySize();

	int getLRU() { return lru_access; }
	void setLRU(int lru) { lru_access = lru; }

private:
	int label;
	int src;
	int dest;
	int cost;

	vector<PartitionEdge*> pathEdges;

	int lru_access;
	int timestamp;
};

class PathCacheComparator {
	friend class PathCacheEntry;
public:
	PathCacheComparator();
	bool operator() (PathCacheEntry& lhs, PathCacheEntry&rhs) const
	  {
	    return lhs.lru_access > rhs.lru_access;
	  }
};

typedef std::priority_queue<PathCacheEntry*,std::vector<PathCacheEntry*>,PathCacheComparator> pathcachepq_type;

//Data structure to track last-updated-timestamps of each SCC in a partition
class PartitionConnectedComponents{
public:
	PartitionConnectedComponents(int initialNumOfComponents, Partition* partition);
	void addConnectedComponent(int componentId);
	void removeConnectedComponent(int componentId);
	void advanceComponentTimestamp(int componentId);
	//Merge component2 into component1
	void mergeConnectedComponents(int componentId1, int componentId2);

	int getLastTimestamp() { return lastTimestamp; }
	int advanceTimestamp() { return ++lastTimestamp; }

	bool inSameComponent(PartitionVertex* v1, PartitionVertex* v2) { return v1->getComponentID() != -1 && (v1->getComponentID() == v2->getComponentID()); }

	int getConnectedComponentCount() { return numOfConnectedComponents; }

private:
	Partition* partition;
	map<int, int> componentLastUpdateTimestamp;

	int lastTimestamp = -1;
	int numOfConnectedComponents = 0;
};

class DistFromSource{
	friend class DistFromSourceComparator;
	friend class Partition;
public:
	DistFromSource(int vId, double dist);
private:
	int vertexId;
	double distance;
};

class DistFromSourceComparator {
	friend class PathCacheEntry;
public:
	DistFromSourceComparator();
	bool operator() (DistFromSource& lhs, DistFromSource&rhs) const
	  {
	    return lhs.distance > rhs.distance;
	  }
};

typedef std::priority_queue<DistFromSource*,std::vector<DistFromSource*>,DistFromSourceComparator> distfromsourcepq_type;

class Partition{
friend class PartitionConnectedComponents;
public:
	Partition(EDPIndex* index, int label);
	~Partition();

	int getLabel() { return label; }

	map<int, PartitionVertex>* getVertexes() { return p_vertexes; }
	map<int, PartitionEdge>* getEdges() { return p_edges; }

	PartitionVertex* getVertex(int vId);
	PartitionEdge* getEdge(int eId);

	PartitionVertex* addVertex(int id);
	PartitionEdge* addEdge(int id, int from, int to, int weightField, int labelField);

	PartitionConnectedComponents* getConnectedComponents() { return &connectedComponents; }
	float getEdgeWeight(int from, int to);

private:
	EDPIndex* index;
	int label;
	//Map of VertexID to PartitionVertex object in this Partition
	map<int, PartitionVertex* > p_vertexes;
	map<int, PartitionEdge* > p_edges;

	//Index of PartitionVertex by ConnectedComponentId
	map<int, vector<PartitionVertex*>> cc_vertexes;

	//Strongly Connected Components calculation
	PartitionConnectedComponents connectedComponents;
	int numOfComponents;
	void buildCC();
	void CC_dfs(PartitionVertex* startVertex, vector<bool>& isVisited, map<int, int>* gToLId, int& curCompId);
	float dijkstra(PartitionVertex* startVertex, PartitionVertex* endVertex);

};

class EDPIndex{
	friend class GraphView;
	friend class Partition;

public:
	typedef pair<double, pair<int, int> > PQEntryWithLength;
	typedef pair<int, int > PQEntry;

	EDPIndex();
	~EDPIndex();
	static EDPIndex buildIndex(GraphView* g, int numLabels, int labelColumn, int costColumn, uint32_t cacheSize);

	PartitionVertex* addVertex(Partition * p, Vertex * v);
	PartitionEdge* addEdge(Edge * e);

	Partition* getPartition(int label);

	PartitionVertex* findVertex(int v, vector<int>* labels);

	//Queries
	void SP_TopK(int src, int dest, int k, vector<int>* labels);


protected:

private:
	GraphView* graph;
	vector<Partition*> partitions;

	//Direct references from vertex ID to vertex, outside of partitions
	map<int, PartitionVertex*> i_vertexes;

	//Column index of label for edges
	int labelColumn;
	//Column index of cost for edges
	int costColumn;

	//Cache fields
	uint64_t maxCacheSize;
	uint64_t cacheSize;
	//Cache. Index of vector is partition label, second key is source node, third key is destination node
	uint64_t lru_clock;
	vector<map<int, map<int, PathCacheEntry>>> cache;
	pathcachepq_type lru;

	//Adds a cache entry, replacing using LRU if necessary
	void insertPath(int label, int src, int dest, double cost, vector<PartitionEdge*>* pathEdges);
	PathCacheEntry* getPath(int label, int src, int dest);
};

}

#endif

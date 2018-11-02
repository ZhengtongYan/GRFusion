#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include "GraphView.h"
#include <unordered_map>

using namespace std;

namespace voltdb {

class PartitionVertex;
class EDPIndex;

class PartitionEdge {
public:
	PartitionEdge(Edge* edge, PartitionVertex* fromVertex, PartitionVertex* toVertex, int weightField, int labelField, bool directed);
	~PartitionEdge(void);

	Edge* getEdge();

	int getID();
	PartitionVertex* getFrom();
	PartitionVertex* getTo();

	double getWeight();
	int getLabel();

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
public:
	PartitionVertex(Vertex * vertex);
	~PartitionVertex(void);

	Vertex* getVertex();

	int getID();
	bool isBridge();
	int getComponentID();
	Partition* getPartition();

	map<int, PartitionEdge*>* getInEdges();
	map<int, PartitionEdge*>* getOutEdges();

	vector<Partition*>* getOtherHosts();

	void setPartition(Partition * p) { partition = p; }
	//Adds an edge to this vertex and updates OtherHostsList if applicable. Partition must be set first
	void addEdge(PartitionEdge* edge, EDPIndex* edp);

private:
	Vertex* vertex;

	int vID;
	bool bridge;
	int componentID;

	Partition* partition;

	map<int, PartitionEdge*> inEdges;
	map<int, PartitionEdge*> outEdges;

	map<int, Partition*> otherHosts;
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
public:
	PathCacheEntry(int src, int dest, int cost, vector<PartitionEdge*>* pathEdges);
	~PathCacheEntry();

	int getSource();
	int getDestination();
	int getCost();

	vector<PartitionEdge*> getPath();

	//Returns the size of the cache entry in bytes
	int getEntrySize();

private:
	int src;
	int dest;
	int cost;

	vector<PartitionEdge*> pathEdges;
};

class Partition{

public:
	Partition(int label);
	~Partition();

	int getLabel();

	map<int, PartitionVertex>* getVertexes();
	map<int, PartitionEdge>* getEdges();

	PartitionVertex* addVertex(int id);
	PartitionEdge* addEdge(int id, int from, int to, float weight, int label, bool addVertexIfNotFound);
private:
	int label;
	//Map of VertexID to PartitionVertex object in this Partition
	map<int, PartitionVertex > p_vertexes;
	map<int, PartitionEdge > p_edges;
};

class EDPIndex{
	friend class GraphView;
	friend class Partition;

public:
	EDPIndex();
	static EDPIndex buildIndex(GraphView* g, int numLabels, int labelColumn, uint32_t cacheSize);

	void addVertex(Vertex * v);
	void addEdge(Edge * e);

	Partition* getPartition(int label);

protected:

private:
	GraphView* graph;
	Partition* partitions;

	//Direct references from vertex ID to vertex, outside of partitions
	map<int, PartitionVertex*> i_vertexes;

	uint64_t cacheSize;

};

}

#endif

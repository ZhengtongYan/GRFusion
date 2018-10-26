#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include "GraphView.h"
#include <unordered_map>

using namespace std;

namespace voltdb {

class Partition;

class PartitionEdge {
public:
	PartitionEdge(void);
	~PartitionEdge(void);

	Edge* getEdge();

	int getID();
	int getFromID();
	int getToID();

	float getWeight();
	int getLabel();
	bool isDirected();

private:
	Edge* edge;

	int eID;
	int fromID;
	int toID;

	float weight;
	int label;
	bool directed;
};

class PartitionVertex {
public:
	PartitionVertex(void);
	~PartitionVertex(void);

	Vertex* getVertex();

	int getID();
	bool isBridge();
	int getComponentID();

	vector<PartitionEdge*>* getInEdges();
	vector<PartitionEdge*>* getOutEdges();

	vector<Partition*>* getOtherHosts();

	//Adds an edge to this vertex and updates OtherHostsList if applicable
	void addEdge(PartitionEdge* edge);

private:
	Vertex* vertex;

	int vID;
	bool bridge;
	int componentID;

	vector<PartitionEdge*> inEdges;
	vector<PartitionEdge*> outEdges;

	vector<Partition*> otherHosts;
};

class MonochromePath{

public:
	MonochromePath(int src, int dest, int weight);
protected:
	int src_id;
	int dest_id;
	int weight;

	vector<PartitionEdge*> pathEdges;
};

class PathCacheEntry {
public:
	PathCacheEntry(int src, int dest, int cost, vector<PartitionEdge*>* pathEdges);
	~PathCacheEntry();

	int getSource();
	int getDestination();
	int getCost();

	vector<PartitionEdge*>* getPath();

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

	vector<map<int, PartitionVertex>>* getPartitions();
	map<int, PartitionEdge>* getEdges();

	PartitionVertex* addVertex(int id);
	PartitionEdge* addEdge(int id, int from, int to, float weight, int label, bool addVertexIfNotFound);
private:
	int label;
	//Map of VertexID to PartitionVertex object in this Partition
	map<int, PartitionVertex > p_vertexes;
	map<int, PartitionEdge > m_edges;
};

class EDPIndex{
	friend class GraphView;
	friend class Partition;

public:

protected:
	GraphView graph;
private:

};

}

#endif

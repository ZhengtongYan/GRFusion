#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include "GraphView.h"

using namespace std;

namespace voltdb {

class MonochromePath{
	friend class Vertex;

public:
	MonochromePath(int src, int dest, int weight)
protected:
	int src_id;
	int dest_id;
	int weight;

};

class Partition{
	friend class Vertex;
	friend class Edge;

public:

protected:
	int label;
	std::map<int, Vertex* > m_vertexes;
	std::map<int, Edge* > m_edges;
private:
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

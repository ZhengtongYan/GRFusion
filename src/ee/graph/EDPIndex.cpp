#include "EDPIndex.h"
#include "Edge.h"
#include "common/ValuePeeker.hpp"


using namespace std;

namespace voltdb
{

	//PartitionEdge member definitions

	PartitionEdge::PartitionEdge(Edge* edge, PartitionVertex* fromVertex, PartitionVertex* toVertex, int weightField, int labelField, bool directed)
	{
		this->edge = edge;

		eID = edge->getId();
		fromV = fromVertex;
		toV = toVertex;

		TableTuple edgeTuple = edge->getGraphView()->getEdgeTuple(edge->getTupleData());
		weight = ValuePeeker::peekDouble(edgeTuple.getNValue(weightField));
		label = ValuePeeker::peekInteger(edgeTuple.getNValue(labelField));
	}

	Edge * PartitionEdge::getEdge() { return edge; }

	int PartitionEdge::getID() { return eID; }
	PartitionVertex* PartitionEdge::getFrom() { return fromV; }
	PartitionVertex* PartitionEdge::getTo() { return toV; }
	int PartitionEdge::getLabel() { return label; }
	double PartitionEdge::getWeight() { return weight; }

	//PartitionVertex member definitions

	PartitionVertex::PartitionVertex(Vertex* vertex)
	{
		this->vertex = vertex;
		vID = vertex->getId();

		componentID = -1;
		bridge = false;
		partition = NULL;
	}

	Vertex * PartitionVertex::getVertex() { return vertex; }

	int PartitionVertex::getID() { return vID; }
	bool PartitionVertex::isBridge() { return bridge; }
	int PartitionVertex::getComponentID() { return componentID; }
	Partition* PartitionVertex::getPartition() { return partition; }

	map<int,PartitionEdge*>* PartitionVertex::getInEdges() { return inEdges; }
	map<int,PartitionEdge*>* PartitionVertex::getOutEdges() { return outEdges; }

	vector<Partition*>* PartitionVertex::getOtherHosts() { return otherHosts; }

	void PartitionVertex::addEdge(PartitionEdge * edge, EDPIndex* edp)
	{
		assert(partition != NULL);
		assert(edge->getFrom() == this);

		outEdges.insert(std::pair<int,PartitionEdge*>(edge->getID(),edge));

		//Insert into in edges in To if not already exists
		edge->getTo()->getInEdges()->insert(std::pair<int,PartitionEdge*>(edge->getID(),edge));

		//Check if the addition of this edge has added to this vertex's OtherHosts list
		otherHosts.insert(std::pair<int,Partition*>(edge->getLabel(),edp->getPartition(edge->getLabel())));
		bridge = otherHosts.empty();
	}

	//Partition member definitions

	Partition::Partition(int label)
	{
		this->label = label;
	}

	Partition::~Partition()
	{
		delete p_vertexes;
		delete p_edges;
	}

	int Partition::getLabel() { return label; }
	map<int, PartitionVertex>* Partition::getVertexes() { return p_vertexes; }
	map<int, PartitionEdge>* Partition::getEdges() { return p_edges; }

	//EDPIndex member definitions

	EDPIndex::EDPIndex()
	{
		graph = NULL;
		partitions = NULL;
		cacheSize = 0;
	}

	static EDPIndex EDPIndex::buildIndex(GraphView* g, int numLabels, int labelColumn, uint32_t cacheSize)
	{
		this->cacheSize = cacheSize;
		graph = g;
		partitions = new Partition[numLabels];


	}

}

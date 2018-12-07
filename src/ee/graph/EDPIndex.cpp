#include "EDPIndex.h"
#include "Edge.h"
#include "Vertex.h"
#include "common/ValuePeeker.hpp"
#include <algorithm>


using namespace std;

namespace voltdb
{

	//PartitionEdge member definitions

	PartitionEdge::PartitionEdge(Edge* edge, PartitionVertex* fromVertex, PartitionVertex* toVertex, int weightField, int labelField)
	{
		this->edge = edge;

		eID = edge->getId();
		fromV = fromVertex;
		toV = toVertex;

		TableTuple edgeTuple = edge->getGraphView()->getEdgeTuple(edge->getTupleData());
		weight = ValuePeeker::peekDouble(edgeTuple.getNValue(weightField));
		label = ValuePeeker::peekInteger(edgeTuple.getNValue(labelField));
	}

	//PartitionVertex member definitions

	PartitionVertex::PartitionVertex(Vertex* vertex)
	{
		this->vertex = vertex;
		vID = vertex->getId();

		componentId = -1;
		bridge = false;
		bridgeBackwards = false;
		partition = NULL;
	}

	void PartitionVertex::addEdge(PartitionEdge * edge, bool outEdge)
	{
		assert(partition != NULL);
		assert((outEdge && edge->getFrom() == this) || (!outEdge && edge->getTo() == this));

		if(outEdge)
			outEdges.insert(std::pair<int,PartitionEdge*>(edge->getID(),edge));
		else
			inEdges.insert(std::pair<int,PartitionEdge*>(edge->getID(),edge));
	}

	//Partition member definitions

	Partition::Partition(EDPIndex* index, int label)
	{
		this->index = index;
		this->label = label;
		this->numOfComponents = 0;
	}

	Partition::~Partition()
	{
		delete p_vertexes;
		delete p_edges;
	}


	PartitionVertex* Partition::addVertex(int vId)
	{
		PartitionVertex* pv = new PartitionVertex(index->graph->getVertex(vId));
		p_vertexes[vId] = pv;
		return pv;
	}

	PartitionVertex* Partition::getVertex(int vId)
	{
		PartitionVertex *pv;
		try
		{
			pv = &(getVertexes()->at(vId));
		} catch(const std::out_of_range& oor){}
		return pv;
	}

	PartitionEdge* Partition::addEdge(int id, int from, int to, int weightField, int labelField)
	{
		PartitionEdge* pe = new PartitionEdge(index->graph->getEdge(id), index->graph->getVertex(from), index->graph->getVertex(to),
				weightField, labelField);

		p_edges[id] = pe;
		return pe;
	}

	PartitionEdge* Partition::getEdge(int eId)
	{
		PartitionEdge *pe;
		try
		{
			pe = &(getEdges()->at(eId));
		} catch(const std::out_of_range& oor){}
		return pe;
	}

	float Partition::getEdgeWeight(int from, int to)
	{
		//Probe cache for
		PathCacheEntry *pce = NULL;
		try
		{
			pce = &(index->cache.at(label).at(from).at(to));
		} catch(const std::out_of_range& oor){}
		if(pce == NULL || (pce != NULL && pce->timestamp < connectedComponents.getLastTimestamp()))
		{

		}
	}

	void Partition::buildCC()
	{
		int numVertexes = p_vertexes.size();
		vector<bool>isVisited = vector<bool>(numVertexes, false);

		map<int, int> globalToLocalId;
		int nextLocalId = 0;

		numOfComponents = 0;
		int currentComponentId = 0;

		for(map<int, PartitionVertex>::iterator it = p_vertexes.begin(); it != p_vertexes.end(); it++)
		{
			it->second.setComponentId(-1);
			globalToLocalId.insert(it->first, nextLocalId++);
		}
		for(map<int, PartitionVertex>::iterator it = p_vertexes.begin(); it != p_vertexes.end(); it++)
		{
			CC_dfs(&(it->second), isVisited, &globalToLocalId, currentComponentId);
		}
		connectedComponents = PartitionConnectedComponents(numOfComponents, this);
	}

	void Partition::CC_dfs(PartitionVertex* startVertex, vector<bool>& isVisited, map<int, int>* gToLId, int& curCompId)
	{
		vector<PartitionVertex*> stack;
		stack.push_back(startVertex);
		PartitionVertex *currentVertex = NULL, *reachableVertex = NULL;

		while(!stack.empty())
		{
			currentVertex = stack.pop_back();
			currentVertex->componentId = curCompId;
			isVisited[gToLId[currentVertex->getID()]] = true;
			//Iterate over out edges
			for(map<int,PartitionEdge*>::iterator out = currentVertex->outEdges.begin(); out != currentVertex->outEdges.end(); out++)
			{
				PartitionVertex* rVertex = out->second->toV;
				if(!isVisited[gToLId[rVertex->getID()]])
				{
					stack.push_back(rVertex);
				}
			}
			//Iterate over in edges
			for(map<int,PartitionEdge*>::iterator in = currentVertex->inEdges.begin(); in != currentVertex->outEdges.end(); in++)
			{
				PartitionVertex* rVertex = in->second->fromV;
				if(!isVisited[gToLId[rVertex->getID()]])
				{
					stack.push_back(rVertex);
				}
			}
		}
		numOfComponents = ++curCompId;
	}

	double Partition::dijkstra(PartitionVertex* startVertex, PartitionVertex* endVertex)
	{
		float distance = -1;
		DistFromSource* dist = new DistFromSource(startVertex->vID, 0), *toDist;
		map<int, DistFromSource*> distMap = new map<int, DistFromSource>;
		distMap[startVertex->vID] = *dist;
		distfromsourcepq_type pq;

		PartitionVertex* u, v;

		pq.push(dist);
		while(!pq.empty())
		{
			dist = pq.pop();
			u = getVertex(dist->vertexId);
			//Explore only direct monoedges
			for(auto e = u->outEdges.begin(); e != u->outEdges.end(); e++)
			{
				v = e->second->toV;
				toDist = NULL;
				try
				{
					toDist = distMap.at(v.vID);
				} catch(const std::out_of_range& oor){}

				if(toDist == NULL)
				{
					toDist = new DistFromSource(v.vID, dist->distance + e->second->weight);
					pq.push(toDist);
					distMap[v.vID] = toDist;
				}
				else if(toDist->distance > dist->distance + e->second->weight)
				{
					toDist->distance = dist->distance + e->second->weight;
				}
			}
		}

		delete distMap;

	}

	DistFromSource::DistFromSource(int vId, double dist)
	{
		vertexId = vId;
		distance = dist;
	}

	//PathCacheEntry member definitions

	PathCacheEntry::PathCacheEntry(int label, int src, int dest, int cost, vector<PartitionEdge*>* pathEdges)
	{
		this->label = label;
		this->src = src;
		this->dest = dest;
		this->cost = cost;
		this->pathEdges = *(new vector<PartitionEdge*>);
		lru_access = 0;
		this->timestamp = INT_MAX;
		for(vector<PartitionEdge*>::iterator peI=pathEdges->begin(); peI!=pathEdges->end(); peI++)
		{
			this->pathEdges.push_back(*peI);
		}
	}

	uint64_t PathCacheEntry::getEntrySize()
	{
		return (sizeof(int) * 5) + sizeof(vector<PartitionEdge*>*)+ sizeof(PartitionEdge*) * pathEdges.size();
	}

	//PartitionConnectedComponents member definitions

	PartitionConnectedComponents::PartitionConnectedComponents(int initialNumOfComponents, Partition* partition)
	{
		numOfConnectedComponents = initialNumOfComponents;
		for(int i = 0; i < numOfConnectedComponents; i++)
		{
			componentLastUpdateTimestamp.insert(i,0);
		}
		lastTimestamp = 0;
		this->partition = partition;
	}

	void PartitionConnectedComponents::addConnectedComponent(int componentId)
	{
		componentLastUpdateTimestamp.insert(componentId, getLastTimestamp()+1);
		advanceTimestamp();
	}

	void PartitionConnectedComponents::removeConnectedComponent(int componentId)
	{
		map<int,int>::iterator it = componentLastUpdateTimestamp.find(componentId);
		if(it != componentLastUpdateTimestamp.end())
			componentLastUpdateTimestamp.erase(it);
	}

	void PartitionConnectedComponents::advanceComponentTimestamp(int componentId)
	{
		componentLastUpdateTimestamp.at(componentId) = advanceTimestamp();
	}

	void PartitionConnectedComponents::mergeConnectedComponents(int componentId1, int componentId2)
	{
		map<int,int>::iterator it1, it2;
		it1 = componentLastUpdateTimestamp.find(componentId1);
		it2 = componentLastUpdateTimestamp.find(componentId2);
		assert(it1 != componentLastUpdateTimestamp.end() && it2 != componentLastUpdateTimestamp.end());

		componentLastUpdateTimestamp.erase(it2);
		it1->second = advanceTimestamp();

		vector<PartitionVertex*>* cc2_vertexes = partition->cc_vertexes[componentId2];

		//Move vertexes from CC2 to CC1
		for(vector<PartitionVertex*>::iterator it = cc2_vertexes->begin(); it != partition->cc_vertexes.end(); it++)
		{
			(*it)->setComponentId(componentId1);
		}

		partition->cc_vertexes.erase(componentId2);
	}

	//EDPIndex member definitions

	EDPIndex::EDPIndex()
	{
		graph = NULL;
		partitions = NULL;
		maxCacheSize = 0;
		cacheSize = 0;
		cache = NULL;
		labelColumn = -1;
		costColumn = -1;
		lru_clock = 0;
	}

	EDPIndex::~EDPIndex()
	{
		delete partitions;
		delete cache;

	}

	PartitionEdge* EDPIndex::addEdge(Edge* e)
	{
		int edgeLabel = ValuePeeker::peekInteger(graph->getEdgeTuple(e->getId())->getNValue(labelColumn));
		int edgeCost = ValuePeeker::peekDouble(graph->getEdgeTuple(e->getId())->getNValue(costColumn));
		Partition * p = &(partitions[edgeLabel]);
		Vertex *vFrom = e->getStartVertex(), *vTo = e->getEndVertex();
		PartitionEdge* pEdge = p->addEdge(e->getId(), e->getStartVertexId(), e->getEndVertexId(), edgeCost, edgeLabel);
		//Check if vertices already exist in partition
		PartitionVertex *pvFrom = p->getVertex(vFrom->getId()), *pvTo = p->getVertex(vTo->getId());

		//If vertices don't exist in partition yet, add them
		if(pvFrom == NULL)
			pvFrom = addVertex(p, vFrom);
		if(pvTo == NULL)
			pvTo = addVertex(p, vTo);

		//Add edge to vertices
		pvFrom->addEdge(pEdge, true);
		pvTo->addEdge(pEdge, false);

		return pEdge;
	}

	PartitionVertex* EDPIndex::addVertex(Partition* p, Vertex* v)
	{
		PartitionVertex *pvFrom = p->addVertex(v->getId());
		//Build OtherHosts list
		map<int, Partition*>* otherHostsOut = pvFrom->getOtherHostsOut();
		map<int, Partition*>* otherHostsIn = pvFrom->getOtherHostsIn();
		//Build forwards OtherHosts list with out edges
		int fanOut = v->fanOut();
		Edge *outEdge;
		for(int i = 0; i < fanOut; i++)
		{
			int outEdgeId = v->getOutEdgeId(i);
			int outEdgeLabel = ValuePeeker::peekInteger(graph->getEdgeTuple(outEdgeId)->getNValue(labelColumn));
			otherHostsOut->insert(std::pair<int, Partition*>(outEdgeLabel, &(partitions[outEdgeLabel])));
		}
		//Build backwards OtherHosts list with in edges
		int fanIn = v->fanIn();
		Edge *inEdge;
		for(int i = 0; i < fanOut; i++)
		{
			int inEdgeId = v->getInEdgeId(i);
			int inEdgeLabel = ValuePeeker::peekInteger(graph->getEdgeTuple(inEdgeId)->getNValue(labelColumn));
			otherHostsIn->insert(std::pair<int, Partition*>(inEdgeLabel, &(partitions[inEdgeLabel])));
		}
		//Set bridge flags for forwards/backwards
		pvFrom->setBridge(!otherHostsOut->empty());
		pvFrom->setBridgeBackwards(!otherHostsIn->empty());

		return pvFrom;
	}

	static EDPIndex EDPIndex::buildIndex(GraphView* g, int numLabels, int labelColumn, int costColumn, uint32_t cacheSize)
	{
		EDPIndex* edp = new EDPIndex();
		edp->maxCacheSize = cacheSize;
		edp->labelColumn = labelColumn;
		edp->costColumn = costColumn;
		edp->graph = g;
		edp->partitions = new vector<Partition>(numLabels, NULL);
		for(int i = 0; i < numLabels; i++)
		{
			edp->partitions[i] = new Partition(edp, i);
		}
		edp->cache = new vector<map<int, map<int, PathCacheEntry>>>(numLabels, new map<int, map<int, PathCacheEntry>>);
		std::map<int, Edge*>* edges = &(g->m_edges);
		std::map<int, Vertex*>* vertices = &(g->m_vertexes);
		Partition p;

		bool isBridge;

		for(std::map<int, Edge*>::iterator eIt=edges->begin(); eIt!=edges->end(); eIt++)
		{
			edp->addEdge(eIt->second);
		}

		return edp;
	}

	void EDPIndex::insertPath(int label, int src, int dest, double cost, vector<PartitionEdge*>* pathEdges)
	{
		PathCacheEntry pce = new PathCacheEntry(label, src, dest, cost, pathEdges);
		uint64_t pce_size = pce.getEntrySize();
		if(maxCacheSize - cacheSize >= pce_size)
		{
			//We have enough space in cache
			cacheSize+=pce_size;
			cache[label][src][dest] = pce;
			lru.push(&pce);
		}
		else
		{
			PathCacheEntry* lru_pce;
			uint64_t lru_pce_size;
			//Remove cache entries by LRU while there's not enough space in the cache
			while(maxCacheSize - cacheSize < pce_size)
			{
				lru_pce = lru.pop();
				int llabel = lru_pce->getLabel(), lsrc = lru_pce->getSource(), ldest = lru_pce->getDestination();
				lru_pce_size = lru_pce->getEntrySize();
				cacheSize -= lru_pce_size;
				delete cache[llabel][lsrc][ldest];
				cache[llabel][lsrc].erase(ldest);
				if(cache[llabel][lsrc].empty())
					cache[llabel].erase(lsrc);
			}
		}

		cache[label][src][dest].setLRU(lru_clock++);
	}

	PathCacheEntry* EDPIndex::getPath(int label, int src, int dest)
	{
		PathCacheEntry* pce = NULL;
		try
		{
			pce = cache.at(label).at(src).at(dest);
			pce->setLRU(lru_clock++);
		}catch(const std::out_of_range& oor){}
		return pce;
	}

	PartitionVertex* EDPIndex::findVertex(int v, vector<int>* labels)
	{
		PartitionVertex* pv = NULL;
		for(vector<int>::iterator lI=labels->begin(); lI!=labels->end();lI++)
		{
			pv = partitions[(*lI)].getVertex(v);
			if(pv != NULL)
			{
				break;
			}
		}
		return pv;
	}

	//EDPIndex Query Functions

	void EDPIndex::SP_TopK(int src, int dest, int k, vector<int>* labels)
	{
		double minCost = DBL_MAX;
		int foundPathsSoFar = 0;
		//PQEntryWithLength.first is the cost, PQEntryWithLength.second.first is the vertexId, PQEntryWithLength.second.second is the path length
		priority_queue<PQEntryWithLength, vector<PQEntryWithLength>, std::greater<PQEntryWithLength>> pq;
		Partition* currentPartition = NULL;
		PartitionVertex *pv_src = findVertex(src, labels, *currentPartition), *pv_dest;
		if(pv_src != NULL && pv_dest != NULL)
			pq.push(make_pair(0, make_pair(src, 0))); //zero cost to reach src from itself
		PartitionVertex* pv = NULL;
		PartitionEdge* pe = NULL;
		int currVId, fanOut = -1, candVertexId = -1;

		//upper-bound to avoid loops (considering an average fan-out of 10
		int maxPQOperations = graph->numOfVertexes() * 10;
		int iterationNum = 0;
		TableTuple* edgeTuple;
		int length;

		while(!pq.empty())
		{
			//select next vertex to explore
			currVId = ((pair<int, int >)( pq.top().second)).first;
			length = ((pair<int, int >)( pq.top().second)).second;
			minCost = pq.top().first;
			if(currVId == dest)
			{
				foundPathsSoFar++;

				//add a tuple here
				TableTuple temp_tuple = graph->m_pathTable->tempTuple();
				//start vertex, end vertex, length, cost, path
				temp_tuple.setNValue(0, ValueFactory::getIntegerValue(src));
				temp_tuple.setNValue(1, ValueFactory::getIntegerValue(dest));
				temp_tuple.setNValue(2, ValueFactory::getIntegerValue(length));
				temp_tuple.setNValue(3, ValueFactory::getDoubleValue(minCost));
				//temp_tuple.setNValue(4, ValueFactory::getStringValue("Test", NULL) );
				if(graph->m_pathTable->activeTupleCount() <= 1000)
				graph->m_pathTable->insertTempTuple(temp_tuple);
			}
			iterationNum++;

			if(foundPathsSoFar == k || iterationNum == maxPQOperations)
			{
				break;
			}
			pq.pop();
			pv_dest = currentPartition->getVertex(dest);
			//explore the outgoing vertexes
			//First, check if the destination is in the same connected component. If so, just search for the shortest path
			pv = findVertex(currVId, labels);
			if(pv != NULL && currentPartition->)
			fanOut = pv->getVertex()->fanOut();
			double edgeCost = 1;
			for(int i = 0; i < fanOut; i++)
			{
				pe = pv->getOutEdges()->at(i);
				if(!std::find(labels->begin(), labels->end(), pe->getLabel())) continue;
				candVertexId = pe->getTo()->getID();

				if (spColumnIndexInEdgesTable >= 0)
				{
					edgeTuple = this->getEdgeTuple(e->getTupleData());
					edgeCost = ValuePeeker::peekDouble(edgeTuple->getNValue(spColumnIndexInEdgesTable));
				}


				//these lines are commented to allow top k search
				//if ( (costMap.find(candVertexId) == costMap.end()) ||
				//	 (costMap[candVertexId] > minCost + 1) )
				{
					//costMap[candVertexId] = minCost + 1;
					pq.push(make_pair(minCost + edgeCost, make_pair(candVertexId, length+1)));
				}

			}
		}

		std::stringstream paramsToPrint;
		paramsToPrint << "TopK SP: from = " << src << ", to = " << dest
				<< ", k = " << k
				<< ", foundPaths = " << foundPathsSoFar
				<< ", numOfRowsAdded = " << m_pathTable->activeTupleCount();
		LogManager::GLog("GraphView", "TopK_SP", 334, paramsToPrint.str());
	}

}

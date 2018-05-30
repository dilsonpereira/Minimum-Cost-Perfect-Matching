#include "Graph.h"

Graph::Graph(int n, list< pair<int, int> > edges):
	n(n),
	m(edges.size()),
	adjMat(n, vector<bool>(n, false)),
	adjList(n),
	edges(),
	edgeIndex(n, vector<int>(n, -1))
{
	int index = 0;
	for(list< pair<int, int> >::iterator it = edges.begin(); it != edges.end(); it++)
	{
		int u = (*it).first;
		int v = (*it).second;

		AddEdge(u, v);
	}
}

pair<int, int> Graph::GetEdge(int e)
{
	if(e > edges.size())
		throw "Error: edge does not exist";

	return edges[e];
}

int Graph::GetEdgeIndex(int u, int v)
{
	if( u > n ||
		v > n )
		throw "Error: vertex does not exist";

	if(edgeIndex[u][v] == -1)
		throw "Error: edge does not exist";

	return edgeIndex[u][v];
}

void Graph::AddVertex()
{
	for(int i = 0; i < n; i++)
	{
		adjMat[i].push_back(false);
		edgeIndex[i].push_back(-1);
	}
	n++;
	adjMat.push_back( vector<bool>(n, false) );
	edgeIndex.push_back( vector<int>(n, -1) );
	adjList.push_back( list<int>() );
}

void Graph::AddEdge(int u, int v)
{
	if( u > n ||
		v > n )
		throw "Error: vertex does not exist";

	if(adjMat[u][v]) return;

	adjMat[u][v] = adjMat[v][u] = true;
	adjList[u].push_back(v);
	adjList[v].push_back(u);

	edgeIndex[u][v] = edgeIndex[v][u] = m++;
}

const list<int> & Graph::GetAdjList(int v)
{
	if(v > n)
		throw "Error: vertex does not exist";

	return adjList[v];
}

bool Graph::Adj(int u, int v)
{
	if( u > n ||
		v > n )
		throw "Error: vertex does not exist";

	return adjMat[u][v];
}


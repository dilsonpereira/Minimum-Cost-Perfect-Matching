#pragma once

#include <list>
#include <vector>


namespace mincostmatching
{

class Graph
{
public:
	//n is the number of vertices
	//edges is a list of pairs representing the edges (default = empty list)
	Graph(int n, const std::list< std::pair<int, int> > & edges = std::list< std::pair<int, int> >());

	//Default constructor creates an empty graph
	Graph(): n(0), m(0) {};

	//Returns the number of vertices
	int GetNumVertices() const { return n; };
	//Returns the number of edges
	int GetNumEdges() const { return m; };

	//Given the edge's index, returns its endpoints as a pair
	std::pair<int, int> GetEdge(int e) const;
	//Given the endpoints, returns the index
	int GetEdgeIndex(int u, int v) const;

	//Adds a new vertex to the graph
	void AddVertex();
	//Adds a new edge to the graph
	void AddEdge(int u, int v);

	//Returns the adjacency list of a vertex
	const std::list<int> & AdjList(int v) const;

	//Returns the graph's adjacency matrix
	const std::vector< std::vector<bool> > & AdjMat() const;
private:
	//Number of vertices
	int n;
	//Number of edges
	int m;

	//Adjacency matrix
	std::vector< std::vector<bool> > adjMat;

	//Adjacency lists
	std::vector< std::list<int> > adjList;

	//Array of edges
	std::vector< std::pair<int, int> > edges;

	//Indices of the edges
	std::vector< std::vector<int> > edgeIndex;
};

}

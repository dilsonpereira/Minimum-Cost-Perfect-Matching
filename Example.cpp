#include "Matching.h"
#include "Graph.h"
#include <cstdlib>
#include <iostream>
using namespace std;

Graph CreateRandomGraph()
{
	//Please see Graph.h for a description of the interface
	
	int n = 100;

	Graph G(n);

	for(int i = 0; i < n; i++)
		for(int j = i+1; j < n; j++)
			if(rand()%10 == 0)
				G.AddEdge(i, j);
	
	return G;
}

void MinimumCostPerfectMatchingExample(Graph & G)
{
	//to get the edge given the index:
	//pair<int, int> = G.GetEdge(i);
	//int u = pair.first, v = pair.second
	
	//to get the index given the edge:
	//int i = G.GetEdgeIndex(u, v);

	vector<double> cost(G.GetNumEdges());
	//cost[i] is the cost of the edge with index i
	for(int i = 0; i < G.GetNumEdges(); i++)
		cost[i] = rand()%1000;

	//Create a Matching instance passing the graph
	Matching M(G);

	//Pass the costs to solve the problem
	pair< list<int>, double > p;
	p = M.SolveMinimumCostPerfectMatching(cost);

	list<int> matching = p.first;
	double obj = p.second;

	cout << "Optimal matching cost: " << obj << endl;
	/*
	cout << "Edges in the matching:" << endl;
	for(list<int>::iterator it = matching.begin(); it != matching.end(); it++)
	{
		pair<int, int> e = G.GetEdge( *it );

		cout << e.first << " " << e.second << endl;
	}
	*/
}

void MaximumMatchingExample(Graph & G)
{
	Matching M(G);

	list<int> matching;
	matching = M.SolveMaximumMatching();

	cout << "Number of edges in the maximum matching: " << matching.size() << endl;
	/*
	cout << "Edges in the matching:" << endl;
	for(list<int>::iterator it = matching.begin(); it != matching.end(); it++)
	{
		pair<int, int> e = G.GetEdge( *it );

		cout << e.first << " " << e.second << endl;
	}
	*/
}

int main(int argc, char* argv[])
{
	//Class Graph will throw const char * exceptions
	try
	{
		int x;
		cout << "Type random seed: ";
		cin >> x;
		srand( x );

		Graph G = CreateRandomGraph();

		MinimumCostPerfectMatchingExample(G);
		//MaximumMatchingExample(G);
	}
	catch(const char * error)
	{
		cout << error << endl;
	}

	return 0;
}




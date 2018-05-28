#pragma once

#include "BinaryHeap.h"
#include <cstring>
#include <list>
#include <vector>
#include <iostream>
using namespace std;

#define EVEN 2
#define ODD 1
#define UNLABELED 0

class Matching
{
public:
	//Parametric constructor receives the number of vertices
	Matching(int n);
	~Matching();

	//Adds a new edge to the graph, without a cost
	//Use this with SolvePerfectMatching()
	void AddEdge(int u, int v);
	//Adds a new edge with a cost
	void AddEdge(int u, int v, double c);
	//Changes the cost
	void SetCost(int u, int v, double c);

	void SolveMinimumCostPerfectMatching();
	void SolvePerfectMatching();

	//Returns 1 if {u,v} is in the matching, 0 otherwise
	int IsInMatching(int u, int v);

	double getObj();

private:
	void Grow();
	void Expand(int u, bool expandBlocked);
	void Augment(int u, int v);
	void Reset();
	int Blossom(int u, int v);
	void UpdateDualCosts();
	void Clear();
	void DestroyBlossom(int t);
	void Heuristic();
	void PositiveCosts();
	void RetrieveMatching();

	int GetFreeBlossomIndex();
	int AddFreeBlossomIndex(int i);
	void ClearBlossomIndices();

	bool IsEdgeBlocked(int u, int v);
	bool IsEdgeBlocked(int e);
	bool IsAdjacent(int u, int v);

	list<int> free;//List of free blossom indices

	vector<int> outer;//outer[v] gives the index of the blossom that contains v but is not contained in any other blossom (default is outer[v] = v)
	vector< list<int> > deep;//deep[v] is a list of all the original vertices contained inside v
	int **shallow;//shallow[v] is a list of the vertices immediately contained inside v, the list has the exact order of the odd circuit of the blossom
	int *sizeShallow;
	vector<int> tip;//tip of the blossom 	
	vector<bool> active;

	vector<int> type;//Even, odd, neither (2, 1, 0)
	vector<int> forest;//forest[v] gives the father of v in the alternating forest
	vector<int> root;//root[v] gives the root of the alternating forest 

	vector<bool> blocked;//A blossom can be blocked, this means that it behaves as if it were an original vertex and cannot be expanded
	vector<double> dual;//dual multipliers associated to the blossoms, if dual[v] > 0, the blossom is blocked and full
	double *slack;//slack associated to each edge, if slack[e] > 0, the edge cannot be used
	vector<int> mate;//mate[v] gives the mate of v

	int *matching;
	double obj;

	int m, n;
	int sizeE;
	int *E1, *E2;//Lists of edges (endpoints)
	double *cost;
	double minEdge;

	int perfect;

	int **AdjMat;
	int **I;
	double **C;
	int **AdjList;
	int *sizeAdjList;

	list<int> BFSList;
	vector<int> visited;
};


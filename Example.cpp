#include "Matching.h"
#include <cstdlib>
#include <stdio.h>
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	//Number of vertices
	//IT MUST BE EVEN!!!
	int n = 6;

	//Create a Matching instance passing the number of vertices
	Matching *M = new Matching(n);

	//Add the edges
	//M->AddEdge(i, j, cost);
	M->AddEdge(0, 1, 4.0);
	M->AddEdge(0, 3, 3.0);
	M->AddEdge(1, 3, 9.0);
	M->AddEdge(1, 4, 12.0);
	M->AddEdge(2, 3, 2.0);
	M->AddEdge(2, 4, 0.0);
	M->AddEdge(2, 5, 8.0);
	M->AddEdge(4, 5, 1.0);

	//The cost can be changed with
	//M->setCost(v1, v2, cost);
	
	//To solve a unweighted perfect matching problem
	//M->SolvePerfectMatching();
	//To solve a weighted perfect matching problem
	M->SolveMinimumCostPerfectMatching();

	//To check whether an edge is in the matching
	//M->IsInMatching(i, j);
	if(M->IsInMatching(0, 1)) printf("{%d, %d} ", 0,1);
	if(M->IsInMatching(0, 3)) printf("{%d, %d} ", 0,3);
	if(M->IsInMatching(1, 3)) printf("{%d, %d} ", 1,3);
	if(M->IsInMatching(1, 4)) printf("{%d, %d} ", 1,4);
	if(M->IsInMatching(2, 3)) printf("{%d, %d} ", 2,3);
	if(M->IsInMatching(2, 4)) printf("{%d, %d} ", 2,4);
	if(M->IsInMatching(2, 5)) printf("{%d, %d} ", 2,5);
	if(M->IsInMatching(4, 5)) printf("{%d, %d} ", 4,5);
	printf("\n");

	//To get optimal cost
	//M->getObj()
	printf("%.2f\n", M->getObj());

	delete M;

	return 0;
}

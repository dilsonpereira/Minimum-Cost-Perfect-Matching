#pragma once

#include "Globals.h"

class BinaryHeap
{
public:
	//The constructor receives the maximum size of the structure
	//It will not resize itself!!!
	BinaryHeap(int maxSize);
	~BinaryHeap();

	//Inserts key k and satellite information s in the heap
	void Insert(double k, int s);

	//Deletes the element with minimum key and returns its satellite information
	int DeleteMin();

	//Resets the structure
	void Clear();

private:
	double * key; 
	int * satellite; 

	int size;
};




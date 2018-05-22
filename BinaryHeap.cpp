#include "BinaryHeap.h"

BinaryHeap::BinaryHeap(int maxSize)
{
	key = new double[maxSize+1];
	satellite = new int[maxSize+1];
	size = 0;
}

BinaryHeap::~BinaryHeap()
{
	delete [] key;
	delete [] satellite;
}

void BinaryHeap::Clear()
{
	size = 0;
}

void BinaryHeap::Insert(double k, int s)
{
	key[s] = k;
	int i;
	for(i = ++size; i/2 > 0 && GREATER(key[i/2], k); i /= 2)
	{
		satellite[i] = satellite[i/2];
		key[i] = key[i/2];
	}
	satellite[i] = s;
	key[i] = k;
}

int BinaryHeap::DeleteMin()
{
	if(size == 0) return -1;

	int min = satellite[1];
	int slast = satellite[size];
	double klast = key[size--];

	int child;
	int i;
	for(i = 1, child = 2; child  <= size; i = child, child *= 2)
	{
		if(child < size && GREATER(key[child], key[child+1]))
			child++;

		if(GREATER(klast, key[child]))
		{
			satellite[i] = satellite[child];
			key[i] = key[child];
		}
		else
			break;
	}
	satellite[i] = slast;
	key[i] = klast;

	return min;
}





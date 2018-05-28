#include "Matching.h"

Matching::Matching(int n):
	n(n),
	outer(2*n),
	tip(2*n),
	active(2*n),
	deep(2*n),
	shallow(2*n),
	blocked(2*n),
	type(2*n),
	mate(2*n),
	dual(2*n),
	forest(2*n),
	root(2*n),
	visited(2*n)
{
	AdjMat = new int*[n];
	AdjList = new int*[n];
	sizeAdjList = new int[n];
	I = new int*[n];
	C = new double*[n];

	for(int i = 0; i < 2*n; i++)
	{
		if(i < n)
		{
			sizeAdjList[i] = 0;
			AdjMat[i] = new int[n];
			AdjList[i] = new int[n];
			I[i] = new int[n];
			C[i] = new double[n];
			for(int j = 0; j < n; j++)
			{
				I[i][j] = -1;
				C[i][j] = 0;
				AdjMat[i][j] = 0;
			}
		}
	}

	E1 = new int[(n*(n-1))/2];
	E2 = new int[(n*(n-1))/2];
	matching = new int[(n*(n-1))/2];
	slack = new double[(n*(n-1))/2];
	cost = new double[(n*(n-1))/2];
	
	for(int i = 0; i < (n*(n-1))/2; i++)
	{
		cost[i] = 0;
	}
	m = 0;	
}

Matching::~Matching()
{
	for(int i = 0; i < 2*n; i++)
	{
		if(i < n)
		{
			delete [] AdjList[i];
			delete [] AdjMat[i];
			delete [] C[i];		
			delete [] I[i];
		}
	}
	delete [] AdjMat;
	delete [] AdjList;
	delete [] sizeAdjList;
	delete [] C;
	delete [] I;

	delete [] E1;
	delete [] E2;
	delete [] slack;
	delete [] cost;
	delete [] matching;
}

void Matching::AddEdge(int u, int v)
{
	E1[m] = u;
	E2[m] = v;	
	AdjMat[u][v] = AdjMat[v][u] = 1;
	I[u][v] = m; 
	I[v][u] = m++;
	AdjList[u][sizeAdjList[u]++]=v;
	AdjList[v][sizeAdjList[v]++]=u; 
}
void Matching::AddEdge(int u, int v, double c)
{	
	cost[m] = c;
	E1[m] = u;
	E2[m] = v;
	AdjMat[u][v] = AdjMat[v][u] = 1;
	C[u][v] = C[v][u] = c;
	I[u][v] = m; 
	I[v][u] = m++; 
	AdjList[u][sizeAdjList[u]++]=v;
	AdjList[v][sizeAdjList[v]++]=u;
}

int Matching::IsInMatching(int u, int v)
{
	return matching[I[u][v]];
}

//Grows an alternating forest
void Matching::Grow()
{
	Reset();

	while(!BFSList.empty())
	{
		int w = outer[BFSList.front()];
		BFSList.pop_front();

		for(list<int>::iterator it = deep[w].begin(); it != deep[w].end(); it++)
		{
			int u = *it;

			int cont = false;
			for(int i = 0; i < sizeAdjList[u]; i++)
			{
				int v = AdjList[u][i];

				if(IsEdgeBlocked(u, v)) continue;

				//u is even and v is odd
				if(type[outer[v]] == ODD) continue;	

				//if v is unlabeled
				if(type[outer[v]] != EVEN)
				{
					//We grow the alternating forest
					int vm = mate[outer[v]];

					forest[outer[v]] = u;
					type[outer[v]] = ODD;
					root[outer[v]] = root[outer[u]];
					forest[outer[vm]] = v;
					type[outer[vm]] = EVEN;
					root[outer[vm]] = root[outer[u]];

					if(!visited[outer[vm]])
					{
						BFSList.push_back(vm);
						visited[outer[vm]] = true;
					}
				}
				//If v is even and u and v are on different trees
				//we found an augmenting path
				else if(root[outer[v]] != root[outer[u]])
				{
					Augment(u,v);
					Reset();

					cont = true;
					break;
				}
				//If u and v are even and on the same tree
				//we found a blossom
				else if(outer[u] != outer[v])
				{
					int b = Blossom(u,v);

					BFSList.push_front(b);
					visited[b] = true;

					cont = true;
					break;
				} 
			}
			if(cont) break;
		}
	}
	
	//Check if the matching is perfect
	perfect = true;
	for(int i = 0; i < n; i++)
		if(mate[outer[i]] == -1)
			perfect = false;
}

bool Matching::IsAdjacent(int u, int v)
{
	return (AdjMat[u][v] and not IsEdgeBlocked(u, v));
}

bool Matching::IsEdgeBlocked(int u, int v)
{
	return GREATER(slack[ I[u][v] ], 0);
}

bool Matching::IsEdgeBlocked(int e)
{
	return GREATER(slack[e], 0);
}

void Matching::Heuristic()
{
	vector<int> degree(n, 0);
	BinaryHeap B(n);

	for(int i = 0; i < m; i++)
	{
		if(IsEdgeBlocked(i)) continue;

		int u = E1[i];
		int v = E2[i];

		degree[u]++;
		degree[v]++;
	}

	for(int i = 0; i < n; i++)
		B.Insert(degree[i], i);

	int u =	B.DeleteMin();
	while(u != -1)
	{
		if(mate[outer[u]] == -1)
		{
			int min = -1;
			for(int i = 0; i < sizeAdjList[u]; i++)
			{
				int v = AdjList[u][i];
				if(IsEdgeBlocked(u, v) ||
					(outer[u] == outer[v]) ||
					(mate[outer[v]] != -1) )
					continue;

				if(min == -1 || degree[v] < degree[min])
					min = v;	
			}
			if(min != -1)
			{
				mate[outer[u]] = min;
				mate[outer[min]] = u;
			}
		}

		u = B.DeleteMin();
	}
}

//Destroys a blossom recursively
void Matching::DestroyBlossom(int t)
{
	if((t < n) ||
		(blocked[t] && GREATER(dual[t], 0))) return;

	for(list<int>::iterator it = shallow[t].begin(); it != shallow[t].end(); it++)
	{
		int s = *it;
		outer[s] = s;
		for(list<int>::iterator jt = deep[s].begin(); jt != deep[s].end(); jt++)
			outer[*jt] = s;	

		DestroyBlossom(s);
	}

	active[t] = false;
	blocked[t] = false;
	AddFreeBlossomIndex(t);
	mate[t] = -1;
}

void Matching::Expand(int u, bool expandBlocked = false)
{
	int v = mate[u];

	int index = m;
	int p, q;
	//Find the regular edge {p,q} of minimum index I[p,q] connecting u and its mate
	for(list<int>::iterator it = deep[u].begin(); it != deep[u].end(); it++)
	{	
		int di = *it;
		for(list<int>::iterator jt = deep[v].begin(); jt != deep[v].end(); jt++)
		{
			int dj = *jt;
			if(IsAdjacent(di, dj) && I[di][dj] < index)
			{
				index = I[di][dj];
				p = di;
				q = dj;
			}
		}
	}
	
	mate[u] = q;
	//If u is a regular vertex, we are done
	if(u < n || (blocked[u] and not expandBlocked)) return;

	bool found = false;
	//Find the position t of the new tip of the blossom
	for(list<int>::iterator it = shallow[u].begin(); it != shallow[u].end() and not found; )
	{
		int si = *it;
		for(list<int>::iterator jt = deep[si].begin(); jt != deep[si].end() and not found; jt++)
		{
			if(*jt == p )
				found = true;
		}
		it++;
		if(not found)
		{
			shallow[u].push_back(si);
			shallow[u].pop_front();
		}
	}

	list<int>::iterator it = shallow[u].begin();
	//Adjust the mate of the tip
	mate[*it] = mate[u];
	it++;
	//
	//Now we go through the odd circuit adjusting the new mates
	while(it != shallow[u].end())
	{
		list<int>::iterator itnext = it;
		itnext++;
		mate[*it] = *itnext;
		mate[*itnext] = *it;
		itnext++;
		it = itnext;
	}

	//We update the sets blossom, shallow, and outer since this blossom is being deactivated
	for(list<int>::iterator it = shallow[u].begin(); it != shallow[u].end(); it++)
	{
		int s = *it;
		outer[s] = s;
		for(list<int>::iterator jt = deep[s].begin(); jt != deep[s].end(); jt++)
			outer[*jt] = s;	

		Expand(s, expandBlocked);
	}
	active[u] = false;
	AddFreeBlossomIndex(u);
}

//Augment the path root[u], ..., u, v, ..., root[v]
void Matching::Augment(int u, int v)
{
	//We go from u and v to its respective roots, alternating the matching
	int p = outer[u];
	int q = outer[v];
	int fp = forest[p];
	mate[p] = q;
	mate[q] = p;
	Expand(p);
	Expand(q);
	while(fp != -1)
	{
		q = outer[forest[p]];
		p = outer[forest[q]];
		fp = forest[p];

		mate[p] = q;
		mate[q] = p;
		Expand(p);
		Expand(q);
	}

	p = outer[v];
	fp = forest[p];
	while(fp != -1)
	{
		q = outer[forest[p]];
		p = outer[forest[q]];
		fp = forest[p];

		mate[p] = q;
		mate[q] = p;
		Expand(p);
		Expand(q);
	}
}

void Matching::Reset()
{
	for(int i = 0; i < 2*n; i++)
	{
		forest[i] = -1;
		root[i] = i;

		if(i >= n && active[i] && outer[i] == i)
			DestroyBlossom(i);
	}

	visited.assign(2*n, 0);
	BFSList.clear();
	for(int i = 0; i < n; i++)
	{
		if(mate[outer[i]] == -1)
		{
			type[outer[i]] = 2;
			if(!visited[outer[i]])
				BFSList.push_back(i);
			visited[outer[i]] = true;
		}
		else type[outer[i]] = 0;
	}
}

int Matching::GetFreeBlossomIndex()
{
	int i = free.back();
	free.pop_back();
	return i;
}

int Matching::AddFreeBlossomIndex(int i)
{
	free.push_back(i);
}

void Matching::ClearBlossomIndices()
{
	free.clear();
	for(int i = n; i < 2*n; i++)
		AddFreeBlossomIndex(i);
}

//Contracts the blossom w, ..., u, v, ..., w, where w is the first vertex that appears in the paths from u and v to their respective roots
int Matching::Blossom(int u, int v)
{
	int t = GetFreeBlossomIndex();

	vector<bool> isInPath(2*n, false);

	//Find the tip of the blossom
	int u_ = u; 
	while(u_ != -1)
	{
		isInPath[outer[u_]] = true;

		u_ = forest[outer[u_]];
	}

	int v_ = outer[v];
	while(not isInPath[v_])
		v_ = outer[forest[v_]];
	tip[t] = v_;

	//Find the odd circuit, update shallow, outer, blossom and deep
	//First we construct the set shallow (the odd circuit)
	list<int> circuit;
	u_ = outer[u];
	circuit.push_front(u_);
	while(u_ != tip[t])
	{
		u_ = outer[forest[u_]];
		circuit.push_front(u_);
	}

	shallow[t].clear();
	deep[t].clear();
	for(list<int>::iterator it = circuit.begin(); it != circuit.end(); it++)
	{
		shallow[t].push_back(*it);
	}

	v_ = outer[v];
	while(v_ != tip[t])
	{
		shallow[t].push_back(v_);
		v_ = outer[forest[v_]];
	}

	//Now we construct deep and update outer
	for(list<int>::iterator it = shallow[t].begin(); it != shallow[t].end(); it++)
	{
		u_ = *it;
		outer[u_] = t;
		for(list<int>::iterator jt = deep[u_].begin(); jt != deep[u_].end(); jt++)
		{
			deep[t].push_back(*jt);
			outer[*jt] = t;
		}
	}

	forest[t] = forest[tip[t]];
	type[t] = EVEN;
	root[t] = root[tip[t]];
	active[t] = true;
	outer[t] = t;
	mate[t] = mate[tip[t]];

	return t;
}

void Matching::UpdateDualCosts()
{
	double e1, e2, e3;
	int inite1 = false, inite2 = false, inite3 = false;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i], v = E2[i];
		if( (type[outer[u]] == EVEN && type[outer[v]] == UNLABELED) || (type[outer[v]] == EVEN && type[outer[u]] == UNLABELED) )
		{
			if(!inite1 || GREATER(e1, slack[i]))
			{
				e1 = slack[i];
				inite1 = true;
			}
		}
		else if( (outer[u] != outer[v]) && type[outer[u]] == EVEN && type[outer[v]] == EVEN )
		{
			if(!inite2 || GREATER(e2, slack[i]))
			{
				e2 = slack[i];
				inite2 = true;
			}
		}
	}
	for(int i = n; i < 2*n; i++)
	{
		if(active[i] && i == outer[i] && type[outer[i]] == ODD && (!inite3 || GREATER(e3, dual[i])))
		{
			e3 = dual[i]; 
			inite3 = true;
		}	
	}
	double e;
	if(inite1) e = e1;
	else if(inite2) e = e2;
	else if(inite3) e = e3;

	if(GREATER(e, e2/2.0) && inite2)
		e = e2/2.0;
	if(GREATER(e, e3) && inite3)
		e = e3;
	 
	for(int i = 0; i < 2*n; i++)
	{
		if(i != outer[i]) continue;

		if(active[i] && type[outer[i]] == EVEN)	
		{
			dual[i] += e; 
		}
		else if(active[i] && type[outer[i]] == ODD)
		{
			dual[i] -= e; 
		}
	}

	sizeE = 0;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i], v = E2[i];			

		if(outer[u] != outer[v])
		{	
			if(type[outer[u]] == EVEN && type[outer[v]] == EVEN)
				slack[i] -= 2.0*e;
			else if(type[outer[u]] == ODD && type[outer[v]] == ODD)
				slack[i] += 2.0*e;
			else if( (type[outer[v]] == UNLABELED && type[outer[u]] == EVEN) || (type[outer[u]] == UNLABELED && type[outer[v]] == EVEN) )
				slack[i] -= e;
			else if( (type[outer[v]] == UNLABELED && type[outer[u]] == ODD) || (type[outer[u]] == UNLABELED && type[outer[v]] == ODD) )
				slack[i] += e;
		}
		
		//if(!GREATER(slack[i], 0))
		//{ 
		//	E[sizeE++] = i;
		//	AdjList[u][sizeAdjList[u]++] = v;
		//	AdjList[v][sizeAdjList[v]++] = u;
		////}
	}
	for(int i = n; i < 2*n; i++)
	{
		if(GREATER(dual[i], 0))
		{
			//if(blocked[i] == false)
			//{
				//the blossom is becoming blocked
			//	mate[i] = mate[tip[i]];	
			//}

			blocked[i] = true;
		}
		else if(active[i] && blocked[i])
		{
			//The blossom is becoming unblocked
			if(mate[i] == -1)
			{
				DestroyBlossom(i);
			}
			else
			{
				blocked[i] = false;
				//int q = mate[i];
				//int p = mate[outer[q]];
				//Expand2(i, p, q);
				Expand(i);
			}
		}
	}	
}

void Matching::SolveMinimumCostPerfectMatching()
{
	SolvePerfectMatching();
	if(!perfect)
		throw("Error: The graph does not have a perfect matching");

	Clear();
	PositiveCosts();

	//Initialize slacks (reduced costs for the edges)
	for(int i = 0; i < m; i++)
		slack[i] = cost[i];

	//Run an heuristic maximum matching algorithm
	Heuristic();
	//Grow a hungarian forest
	Grow();

	//If the matching on the compressed graph is perfect, we are done
	while(not perfect)
	{
		UpdateDualCosts();
		//Set up the algorithm for a new grow step
		Reset();
		Heuristic();
		Grow();
	}

	RetrieveMatching();
}

void Matching::PositiveCosts()
{
	minEdge = 0;
	for(int i = 0; i < m ;i++)
		if(GREATER(minEdge - cost[i], 0)) 
			minEdge = cost[i];

	for(int i = 0; i < m; i++)
		cost[i] -= minEdge;
}

void Matching::SolvePerfectMatching()
{
	Clear();
	Grow();
	RetrieveMatching();
}

//Sets up the algorithm for a new run
void Matching::Clear()
{
	ClearBlossomIndices();

	for(int i = 0; i < 2*n; i++)
	{
		outer[i] = i;
		deep[i].clear();
		if(i<n)
			deep[i].push_back(i);
		shallow[i].clear();
		if(i < n) active[i] = true;
		else active[i] = false;
	
		type[i] = 0;
		forest[i] = -1;
		root[i] = i;

		blocked[i] = false;
		dual[i] = 0;
		mate[i] = -1;
		tip[i] = i;
	}
	for(int i = 0; i < (n*(n-1))/2; i++)
	{
		slack[i] = 0;
		matching[i] = 0;
	}
	
	sizeE = 0;
}

void Matching::SetCost(int u, int v, double c)
{
	C[u][v] = C[v][u] = c;
	cost[I[u][v]] = c;			
}

double Matching::getObj()
{
	return obj;
}

void Matching::RetrieveMatching()
{
	for(int i = 0; i < 2*n; i++)
		if(active[i] && outer[i] == i)
			Expand(i, true);

	obj = 0;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i];
		int v = E2[i];

		cost[i] += minEdge;

		if(mate[u] == v)
		{
			matching[i] = 1;
			obj += cost[i];
		}
		else
		{
			matching[i] = 0;
		}
	}

	double dualObj = 0;
	for(int i = 0; i < 2*n; i++)
	{
		if(i < n) dualObj += dual[i];
		else if(blocked[i]) dualObj += dual[i];	
	}
}

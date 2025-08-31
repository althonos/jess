// ==================================================================
// Scanner.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type Scanner (the main Jess query object).
// ==================================================================

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "Scanner.h"
#include "KdTree.h"
#include "Region.h"
#include "Annulus.h"
#include "Join.h"
#include "CandidateSet.h"

// ==================================================================
// type ScannerData
//===================================================================

struct _ScannerData
{
	CandidateSet** candidates;
	int* index;
	Atom** atom;
	Join** regions;
	bool *active;
	KdTree **trees;
	KdTreeQuery **queries;
	double* weights;
	int* order;
	size_t count;
	size_t capacity;
};

// ==================================================================
// type Scanner
//===================================================================
// template				The template object
// set[k]				Set of candidates for atom k
// tree[k]				Tree of candidate positions for atom k
// query[k]				Current query state for tree k
// index[k]				Index of result[k] in set[k].
// result[k]			kth atom of current result set
// region[i]			Temporary region pointer
// count				= template->count(template)
// threshold			The global distance cutoff  
// max_total_threshold		The maximum value the distance cutoff can take
// 				after adding the global and single-residue
// 				distance cutoff
//===================================================================

struct _Scanner
{
	Template *template;
	CandidateSet **set;
	KdTree **tree;
	KdTreeQuery **query;
	bool *active;
	int *index;
	Atom **atom;
	Join **regions;
	double *weights;
	int* order;
	int count;
	double threshold;
	double max_total_threshold;
};

// ==================================================================
// Methods of type Scanner
// ==================================================================

#ifdef HAVE_THREADLOCALSTORAGE
static __thread CandidateSet **candidates;
static int _CandidateSet_size_compare(const void* a, const void* b)
{
	int x = *((int*)a);
	int y = *((int*)b);
	int diff = candidates[x]->count - candidates[y]->count;
	return (diff == 0) ? x-y : diff;
}
#endif

#ifdef HAVE_GNU_QSORT_R
static int _CandidateSet_size_compare_r(const void* a, const void* b, void* data)
#else
static int _CandidateSet_size_compare_r(void* data, const void* a, const void* b)
#endif
{
	CandidateSet** candidates = (CandidateSet**) data;
	int x = *((int*)a);
	int y = *((int*)b);
	int diff = candidates[x]->count - candidates[y]->count;
	return (diff == 0) ? x-y : diff;
}

Scanner *Scanner_create(Molecule *M, Template *T, ScannerData* D, double r, double s, bool reorder)
{
	return Scanner_reuse(NULL,M,T,D,r,s,reorder);
}

Scanner *Scanner_reuse(Scanner *S, Molecule *M, Template *T, ScannerData* D, double r, double s, bool reorder)
{
	int k,n=T->count(T);
	int m;

	if(!D) return NULL;
	if(!S)
	{
		S=(Scanner*) malloc(sizeof(Scanner));
		if(!S) return NULL;
	}

	if(ScannerData_resize(D,n)!=0) 
		return NULL;

	S->set=D->candidates;
	S->tree=D->trees;
	S->query=D->queries;
	S->index=D->index;
	S->atom=D->atom;
	S->regions=D->regions;
	S->active=D->active;
	S->weights=D->weights;
	S->order=D->order;

	S->template=T;
	S->threshold=r;
	S->max_total_threshold=s;
	S->count=n;

	for(k=0; k<n; k++)
	{
		S->atom[k]=NULL;
		S->order[k]=k;
		S->index[k]=-1;
		S->active[k]=false;
		S->weights[k]=S->template->distWeight(S->template, k);

		T->candidates(T,M,k, &S->set[k]);
		if(S->set[k]->count==0)
		{
			Scanner_free(S);
			return NULL;
		}

		S->tree[k]=KdTree_reuse(S->tree[k],S->set[k]->coord,S->set[k]->count,3);
		if(!S->tree[k])
		{
			Scanner_free(S);
			return NULL;
		}

	}

#if defined(HAVE_GNU_QSORT_R)
	if(reorder) qsort_r(S->order,n,sizeof(int),_CandidateSet_size_compare_r,S->set);
#elif defined(HAVE_APPLE_QSORT_R)
	if(reorder) qsort_r(S->order,n,sizeof(int),S->set,_CandidateSet_size_compare_r);
#elif defined(HAVE_WIN32_QSORT_S)
	if(reorder) qsort_s(S->order,n,sizeof(int),_CandidateSet_size_compare_r,S->set);
#elif defined(HAVE_THREADLOCALSTORAGE)
	if(reorder)
	{
		candidates = S->set;
		qsort(S->order,n,sizeof(int),_CandidateSet_size_compare);
	}
#else
#warning "No re-entrant `qsort` implementation or thread-local storage, no iteration reordering will be done."
#endif	

	if(S->count>0 && S->set[S->order[0]]->count>0)
	{
		S->index[S->order[0]]=0;
		S->atom[S->order[0]]=S->set[S->order[0]]->atom[0];
	}

	return S;
}

void Scanner_free(Scanner *S)
{
	if(S)
	{
		free(S);
	}
}

Atom **Scanner_next(Scanner *S, int ignore_chain)
{
	int j,k;
	double min,max;
	double dynamic_threshold = S->threshold;

	k=S->count-1;

	// Attempt to find the next query result.

	while(k>=0)
	{
		// If k==S->count, we have a hit!

		if(k==S->count) break;

		// If k==0 we must find the next left-most
		// atom of the query...

		if(k==0)
		{
			S->index[S->order[0]]++;

			if(S->index[S->order[0]]>=S->set[S->order[0]]->count)
			{
				// End of query...
				k=-1;
			}
			else
			{
				S->atom[S->order[0]]=S->set[S->order[0]]->atom[S->index[S->order[0]]];
				k++;
			}

			continue;
		}

		// So k>0. If there is an active query for this
		// set then query it now...

		if(S->active[S->order[k]])
		{
			S->index[S->order[k]]=KdTreeQuery_next(S->query[S->order[k]]);
			if(S->index[S->order[k]]<0)
			{
				// The query ended. So we need to destroy
				// this query, then drop down a level...

				S->active[S->order[k]]=false;
				S->atom[S->order[k]]=NULL;
				k--;
			}
			else
			{
				// The query was successful(?) Remember the
				// atom, check n-ary constraints and continue
				// up...

				S->atom[S->order[k]]=S->set[S->order[k]]->atom[S->index[S->order[k]]];
				if(S->template->check(S->template,S->atom,S->order,k+1,ignore_chain))
				{
					k++;
				}
			}

			continue;
		}

		// There is no active query for set k. If there is
		// no active query result for set k-1 then we need
		// to drop down again...

		if(S->index[S->order[k-1]]<0)
		{
			k--;
			continue;
		}

		// So, there is an active query result at k-1 and
		// no active query at k; create a new query at
		// index k and try again (with the same k)

		for(j=0; j<k; j++)
		{
			S->template->range(S->template,S->order[j],S->order[k],&min,&max);

			dynamic_threshold = S->threshold + S->weights[S->order[j]] + S->weights[S->order[k]];
			// Limit threshold to a hard cutoff so execution does not suffer
			if(dynamic_threshold > S->max_total_threshold){
				dynamic_threshold = S->max_total_threshold;
			}
			min -= dynamic_threshold;
			max += dynamic_threshold;
			if(min<0.5) min=0.5;

			assert(S->atom[S->order[j]]);
			S->regions[k]->R[j]=Annulus_reuse(S->regions[k]->R[j],S->atom[S->order[j]]->x,min,max,3);
		}

		S->active[S->order[k]]=true;
		S->query[S->order[k]]=KdTreeQuery_reuse(S->query[S->order[k]],S->tree[S->order[k]],S->regions[k]);
		if(!S->query[S->order[k]]) return NULL;
	}

	// If k<0 there is no more!

	if(k<0) return NULL;

	// Otherwise, the atoms are listed in S->atom.

	return S->atom;
}

// ==================================================================
// Methods of type ScannerData
// ==================================================================

ScannerData *ScannerData_create()
{
	ScannerData* D = malloc(sizeof(ScannerData));
	if(!D) return NULL;
	D->count=0;
	D->capacity=0;
	D->candidates=NULL;
	D->index=NULL;
	D->atom=NULL;
	D->regions=NULL;
	D->active=NULL;
	D->queries=NULL;
	D->trees=NULL;
	D->weights=NULL;
	D->order=NULL;
	return D;
}

void ScannerData_free(ScannerData* D)
{
	if(D)
	{
		if(D->candidates) 
		{
			for(int i=0;i<D->capacity;i++) CandidateSet_free(D->candidates[i]);
			free(D->candidates);
		}
		if(D->regions)
		{
			for(int i=0;i<D->capacity;i++) Join_free(D->regions[i]);
			free(D->regions);
		}
		if(D->queries)
		{
			for(int i=0;i<D->capacity;i++) KdTreeQuery_free(D->queries[i]);
			free(D->queries);
		}
		if(D->trees)
		{
			for(int i=0;i<D->capacity;i++) KdTree_free(D->trees[i]);
			free(D->trees);
		}
		if(D->index) free(D->index);
		if(D->atom) free(D->atom);
		if(D->active) free(D->active);
		if(D->weights) free(D->weights);
		if(D->order) free(D->order);
		free(D);
	}
}

int ScannerData_resize(ScannerData* D, int n)
{
	if(D->capacity<n)
	{
		D->candidates=(CandidateSet**)realloc(D->candidates,n*sizeof(CandidateSet*));
		if(!D->candidates) return -1;
		for(int i=D->capacity;i<n;i++)
		{
			D->candidates[i]=CandidateSet_create();
			if(!D->candidates[i]) return -1;
		}

		D->regions=(Join**)realloc(D->regions,n*sizeof(Join*));
		if(!D->regions) return -1;
		for(int i=D->capacity;i<n;i++)
		{
			D->regions[i]=Join_allocate(i,innerJoin);
			if(!D->regions[i]) return -1;
		}

		D->trees=(KdTree**)realloc(D->trees,n*sizeof(KdTree*));
		if(!D->trees) return -1;
		for(int i=D->capacity;i<n;i++) D->trees[i]=NULL;

		D->queries=(KdTreeQuery**)realloc(D->queries,n*sizeof(KdTreeQuery*));
		if(!D->queries) return -1;
		for(int i=D->capacity;i<n;i++) D->queries[i]=NULL;

		D->index=(int*)realloc(D->index,n*sizeof(int));
		if(!D->index) return -1;

		D->atom=(Atom**)realloc(D->atom,n*sizeof(Atom*));
		if(!D->atom) return -1;

		D->active=(bool*)realloc(D->active,n*sizeof(bool));
		if(!D->active) return -1;

		D->weights=(double*)realloc(D->weights,n*sizeof(double));
		if(!D->weights) return -1;

		D->order=(int*)realloc(D->order,n*sizeof(int));
		if(!D->order) return -1;

		D->capacity=n;
	}

	D->count=n;
	return 0;
}

// ==================================================================

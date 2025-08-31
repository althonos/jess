// ==================================================================
// KdTree.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type KdTree, related types and methods.
// ==================================================================

#include "KdTree.h"
#include "qselect.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#define NO_NODE	SIZE_MAX
#define MAX_DIM 3

// ==================================================================
// Forward declarations of local types
// ==================================================================
// KdTreeNode			One node of a KdTree
// ==================================================================

typedef size_t index_t; 
typedef struct _KdTreeNode KdTreeNode;

// ==================================================================
// Local type KdTreeNode
// ==================================================================
// type					One of 0,1,2 for a node or -1 for a leaf
// index				Index of splitting coordinate or point
// left,right			Branches of tree (unless this is a leaf)
// min,max				The box region encompassed by this node
// depth				The depth of the tree below this node
// ==================================================================

struct _KdTreeNode
{
	index_t left;
	index_t right;
	double min[MAX_DIM];
	double max[MAX_DIM];
	int type;
	int index;
	int depth;
};

// ==================================================================
// Declaration of methods of local type KdTreeNode
// ==================================================================
// create(...)			Creates a new node and its descendants
// free(N)				Frees a node and all its descendants
// ==================================================================

static index_t KdTreeNode_create(KdTree*,int*,int,int,double**,int);
// static void KdTreeNode_free(KdTreeNode*);

// ==================================================================
// type KdTree
// ==================================================================
// root					The root of the tree
// dim					Dimension of the points
// node					The nodes stored in an array
// capacity				The capacity of the node array
// count				The number of element in the node array
// dim					The dimension of the tree (at most MAX_DIM)
// ==================================================================

struct _KdTree
{
	KdTreeNode *nodes;
	size_t count;
	size_t capacity;
	index_t root;
	int dim;
};

// ==================================================================
// type KdTreeQuery
// ==================================================================
// tree					The tree which the query relates to
// region				The region being queried (see Region.h)
// count				Number of nodes on the stack
// stack				The stack of nodes in the query
// ==================================================================

struct _KdTreeQuery
{
	KdTree *tree;
	Join *region;
	Box box;
	int count;
	int maxdepth;
	index_t stack[0];
};


// ==================================================================
// Declaration of private methods/members of type KdTree
// ==================================================================
// compare(pa,pb)		Used during node creation (qsort)
// data,index			Static globals (see KdTreeNode_create)
// ==================================================================

static int KdTree_compare(const void*, const void*);
static int KdTree_compare_r(const void*, const void*, void*);

#ifdef HAVE_THREADLOCALSTORAGE
static __thread double **KdTree_data;
static __thread int KdTree_index;
#else
static double **KdTree_data;
static int KdTree_index;
#endif

struct _KdTreeCompareData
{
      double** data;
      int index;
};

// ==================================================================
// Public methods of type KdTree
// ==================================================================

KdTree *KdTree_create(double **u, int n, int d)
{
	KdTree *K;
	int i,j;
	int *tmp;

	if(n<1 || d<1 || !u) return NULL;
	if(d>MAX_DIM) return NULL;

	// 1. Create memory for the object.

	K = (KdTree*)malloc(sizeof(KdTree));
	if(!K) return NULL;

	K->root=NO_NODE;
	K->nodes=NULL;
	K->capacity=0;
	K->count=0;

	return KdTree_reuse(K,u,n,d);
}

KdTree *KdTree_reuse(KdTree *K, double **u, int n, int d)
{
	int i,j;
	int *tmp;

	if(!K) return KdTree_create(u,n,d);

	if((n<1 || d<1 || !u) || (d>MAX_DIM)) {
		KdTree_free(K);
		return NULL;
	}

	K->dim=d;
	K->count=0;
	K->root=NO_NODE;

	// 3a. Create a temporary array to hold indices

	tmp = (int*)malloc(n*sizeof(int));
	for(i=0; i<n; i++) tmp[i]=i;

	// 3b. Create the tree recursively. This takes time
	// of order at most n.log(n)^2, assuming that qsort
	// always manages n.log(n) and d is constant.

	K->root = KdTreeNode_create(K,tmp,n,0,u,d);
	free(tmp);

	if(K->root == NO_NODE) 
	{
		KdTree_free(K);
		return NULL;
	}

	// 4. Return the result!

	return K;
}

void KdTree_free(KdTree *K)
{
	int i;

	if(K)
	{
		if(K->nodes) free(K->nodes);
		free(K);
	}
}

KdTreeQuery *KdTree_query(KdTree *K, Join *J)
{
	return KdTreeQuery_reuse(NULL,K,J);
}

// ==================================================================
// Methods of type KdTreeQuery
// ==================================================================

KdTreeQuery *KdTreeQuery_reuse(KdTreeQuery *Q, KdTree *K, Join *J)
{
	int rq;

	if(!(Q) || (Q->maxdepth < K->nodes[K->root].depth))
	{
		rq = sizeof(KdTreeQuery)+K->nodes[K->root].depth*sizeof(index_t);
		Q=(KdTreeQuery*)realloc(Q,rq);
		if(!Q) return NULL;
		Q->maxdepth=K->nodes[K->root].depth;
	}

	Q->tree=K;
	Q->region=J;
	Q->count=1;
	Q->stack[0]=K->root;

	// NB: compute a bounding box around the `Join` so we can use
	//	   the bounding box to compute intersections in `KdTreeQuery_next`
	//	   instead of computing the individual `Annulus` intersections.

	Q->box.dim = K->dim;
	Join_computeBox(Q->region,&Q->box);

	return Q;
}

int KdTreeQuery_next(KdTreeQuery *Q)
{
	KdTreeNode *N;
	Join *J = Q->region;
	Box *B = &Q->box;
	index_t *stack=&(Q->stack[0]);
	int dim = Q->tree->dim;
	int *count = &(Q->count);

	// Until the stack is empty (or we return inside
	// the while loop...

	while(*count>0)
	{
		// Pull the top node off the stack.

		N = &Q->tree->nodes[stack[--(*count)]];

		// If the node is a leaf we simply test it and
		// remove it from the stack. If the point is in
		// the query region, return it; otherwise continue
		// with the rest of the stack.

		if(N->type<0)
		{
			if(_Join_po(J,N->min,dim))
			{
				return N->index;
			}
			else
			{
				continue;
			}
		}

		// So the node is internal (ie not a leaf).
		// If the query region does not intersect
		// the node's region then we can remove it
		// and continue with the rest of the stack.

		if(J->type==innerJoin)
		{
			if(!_Box_ro(B,N->min,N->max,dim)) continue;
		}
		else
		{
			if(!_Join_ro(J,N->min,N->max,dim)) continue;
		}

		// The query region *does* intersect the node's
		// region. So now we must place the child nodes
		// onto the stack.

		stack[(*count)++]=N->left;
		stack[(*count)++]=N->right;
	}

	return -1;
}

void KdTreeQuery_free(KdTreeQuery *Q)
{
	if(Q)
	{
		free(Q);
	}
}

// ==================================================================
// Private methods of type KdTree
// ==================================================================

static int KdTree_compare(const void *pa, const void *pb)
{
	const int a = *((const int*)pa);
	const int b = *((const int*)pb);
	double c = KdTree_data[a][KdTree_index];
	double d = KdTree_data[b][KdTree_index];

	return c<d ? -1 : c>d ? 1 : 0;
}

static int KdTree_compare_r(const void* pa, const void* pb, void* data)
{
	struct _KdTreeCompareData _data  = *((struct _KdTreeCompareData*) data);
	const int a = *((const int*)pa);
	const int b = *((const int*)pb);
	double c = _data.data[a][_data.index];
	double d = _data.data[b][_data.index];

	return c<d ? -1 : c>d ? 1 : 0;
}

// ==================================================================
// Methods of local type KdTreeNode
// ==================================================================

#ifndef Jess_min
#define Jess_min(x,y) (x<y ? x:y)
#endif

#ifndef Jess_max
#define Jess_max(x,y) (x>y ? x:y)
#endif

static index_t KdTreeNode_create(KdTree *K, int *idx, int n, int type,double **u,int dim)
{
	KdTreeNode *N;
	int split;
	int i,rq;
	index_t k, left, right;

	// 1. The really easy case. If n is 0 do nothing!

	if(n<=0) return NO_NODE;

	// 1.5. We'll need to create a node in all other cases.

	if(K->count>=K->capacity) {
		K->capacity = K->capacity + (K->capacity >> 3) + 6;
		K->nodes = realloc(K->nodes, K->capacity*sizeof(KdTreeNode));
		if(!K->nodes) return NO_NODE;
	}

	k=K->count;
	N=&K->nodes[k];
	K->count++;
	N->left = NO_NODE;
	N->right = NO_NODE;

	// 2. The easy case. If n is 1, create a leaf.

	if(n==1)
	{
		N->type=-1;
		N->index=idx[0];
		N->depth=1;
		N->left = NO_NODE;
		N->right = NO_NODE;
		memcpy(N->min,u[idx[0]],sizeof(double)*dim);
		memcpy(N->max,u[idx[0]],sizeof(double)*dim);
		return k;
	}
 
	// 2.5. Now we need to order the indices by coordinate
	// numbered type.

	struct _KdTreeCompareData _data = { u, type };
	split = qselect_r(idx, n, sizeof(int), n/2, KdTree_compare_r, &_data);

	// 3. The recursive case. Find [n/2] and split the array into
	// two pieces. Create a node whose splitting value is the median.
	// But make sure that if there are several entries with the same
	// coordinate that we take the right-most.

	N->type=type;
	N->index=idx[split-1];

	// Now create the left and right branches of the node.

	type = (type+1)%dim;
	left = KdTreeNode_create(K, idx,split,type,u,dim);
	right = KdTreeNode_create(K, &idx[split],n-split,type,u,dim);
	if((left==NO_NODE) || (right==NO_NODE)) return NO_NODE;

	// DANGER: we need to update the pointer `N` because the memory for
	//		   storing nodes may have been reallocated in the recursive
	// 		   `KdTreeNode_create` calls. Then it's safe again to access
	//		   and update it.

	N=&K->nodes[k];
	N->left=left;
	N->right=right;

	// Compute max,min and depth...
	N->depth=Jess_max(K->nodes[left].depth,K->nodes[right].depth)+1;

	for(i=0; i<dim; i++)
	{
		N->min[i]=Jess_min(K->nodes[left].min[i],K->nodes[right].min[i]);
		N->max[i]=Jess_max(K->nodes[left].max[i],K->nodes[right].max[i]);
	}

	// We're done...
	return k;
}

#undef Jess_min
#undef Jess_max

// ==================================================================


// ==================================================================
// KdTree.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type KdTree, related types and methods.
// ==================================================================

#include "KdTree.h"
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
// Local functions
// ==================================================================

#define min(x,y) (x<y ? x:y)
#define max(x,y) (x>y ? x:y)

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
	memset(K->nodes,0,sizeof(KdTreeNode)*K->capacity);

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
	KdTreeQuery *Q;
	int rq;

	rq = sizeof(KdTreeQuery)+K->nodes[K->root].depth*sizeof(index_t);

	Q = (KdTreeQuery*)malloc(rq);
	Q->tree=K;
	Q->region=J;
	Q->count=1;
	Q->maxdepth=(K->root==NO_NODE) ? 0 : K->nodes[K->root].depth;
	Q->stack[0]=K->root;

	return Q;
}

// ==================================================================
// Methods of type KdTreeQuery
// ==================================================================

KdTreeQuery *KdTreeQuery_reuse(KdTreeQuery *Q, KdTree *K, Join *J)
{
	int rq;

	if((!Q)) return KdTree_query(K,J);

	if((Q->maxdepth < K->nodes[K->root].depth))
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

	return Q;
}

int KdTreeQuery_next(KdTreeQuery *Q)
{
	KdTreeNode *N;
	Join *J = Q->region;
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

		if(!_Join_ro(J,N->min,N->max,dim))
		{
			continue;
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

typedef int(*compare_t)(const void*, const void*, void*);


static inline void
_qselect_memswap(void *restrict p1, void *restrict p2, size_t n)
{
  while (n > 0)
    {
      unsigned char t = ((unsigned char *)p1)[--n];
      ((unsigned char *)p1)[n] = ((unsigned char *)p2)[n];
      ((unsigned char *)p2)[n] = t;
    }
}

static size_t _qselect_partition(void* base, size_t size, size_t left, size_t right, size_t pivot_index, compare_t compare, void* arg)
{
	size_t store_index = left;
	_qselect_memswap(base+right*size, base+pivot_index*size, size);
	for(size_t i=left; i<right; i++)
	{
		if (compare(base+i*size, base+right*size, arg) <= 0)
		{
			_qselect_memswap(base+i*size, base+store_index*size, size);
			store_index += 1; 
		}
	}
	_qselect_memswap(base+store_index*size, base+right*size, size);
	return store_index;
}

static size_t _qselect_median3(void* base, size_t size, size_t left, size_t right, compare_t compare, void* arg)
{
	size_t mid = (left + right + 1) / 2;
	if ((compare(base+left*size, base+mid*size, arg) > 0) != (compare(base+left*size, base+right*size, arg) > 0))
		return left;
	else if ((compare(base+mid*size, base+left*size, arg) < 0) != (compare(base+mid*size, base+right*size, arg) < 0))
		return mid;
	else
		return right;
}

static size_t qselect_r(void* base, size_t n, size_t size, size_t k, compare_t compare, void* arg) 
{
	size_t pivot_index;
	size_t left = 0;
	size_t right = n-1;

	while(1) {
		if (left == right)
			return left;

		pivot_index = _qselect_median3(base, size, left, right, compare, arg); 
		pivot_index = _qselect_partition(base, size, left, right, pivot_index, compare, arg);

		if(pivot_index == k) 
			return pivot_index;
		else if (k < pivot_index)
			right = pivot_index - 1;
		else
			left = pivot_index + 1;
	}
}


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
		K->capacity = (K->capacity == 0) ? 32 : K->capacity*2;
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
	N->depth = max(K->nodes[left].depth,K->nodes[right].depth)+1;

	for(i=0; i<dim; i++)
	{
		N->min[i]=min(K->nodes[left].min[i],K->nodes[right].min[i]);
		N->max[i]=max(K->nodes[left].max[i],K->nodes[right].max[i]);
	}

	// We're done...
	return k;
}

#undef min
#undef max

// ==================================================================


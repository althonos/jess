// ==================================================================
// Join.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of region type Join. This joins two regions as a union
// or intersection. Caveat: the intersection oracle is sub-optimal
// but I don't know how to fix it!
// ==================================================================

#ifndef JOIN_H
#define JOIN_H

#include "Annulus.h"
#include "Region.h"

// ==================================================================
// Join types
// ==================================================================
// innerJoin				An intersection of two regions
// outerJoin				The union of two regions
// ==================================================================

typedef enum {innerJoin,outerJoin} JoinType;

// ==================================================================
// type Join
// ==================================================================
// count				The number of regions in the join
// R[k]					The kth region in the join
// ==================================================================

typedef struct _Join Join;

struct _Join
{
	int count;
	Region *R[0];
};

// ==================================================================
// Construction methods for region type Join
// ==================================================================
// create(R,n,type)			Create inner/outer join on R[0...n-1]
// free(J)					Frees join AND nested regions (J->free)
// ==================================================================

extern Region *Join_create(Region**,int,JoinType);
extern void Join_free(Region*);

// ==================================================================
// Oracles
// ==================================================================

int Join_oro(Region *R,double *min,double *max,int dim);
int Join_iro(Region *R,double *min,double *max,int dim);
int Join_opo(Region *R,double *x,int dim);
int Join_ipo(Region *R,double *x,int dim);

// ==================================================================
// The oracles
// ==================================================================

static inline int _Join_oro(Join *J,double *min,double *max,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) &(J->R[k])[1];
		if(_Annulus_ro(A,min,max,dim))
		{
			return 1;
		}
	}

	return 0;
}

static inline int _Join_iro(Join *J,double *min,double *max,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) &(J->R[k])[1];
		if(!(_Annulus_ro(A,min,max,dim)))
		{
			return 0;
		}
	}

	return 1;
}

static inline int _Join_opo(Join *J,double *x,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) &(J->R[k])[1];
		if(_Annulus_po(A,x,dim))
		{
			return 1;
		}
	}

	return 0;
}

static inline int _Join_ipo(Join *J,double *x,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) &(J->R[k])[1];
		if(!(_Annulus_po(A,x,dim)))
		{
			return 0;
		}
	}

	return 1;
}

// ==================================================================

#endif


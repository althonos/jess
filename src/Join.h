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
#include "Box.h"
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
// type					The type of the join (inner or outer)
// R[k]					The kth region in the join
// Box					A bounding box around the join region.
// ==================================================================

typedef struct _Join Join;

struct _Join
{
	int count;
    JoinType type;
	Annulus *R[0];
};

// ==================================================================
// Construction methods for region type Join
// ==================================================================
// create(R,n,type)			Create inner/outer join on R[0...n-1]
// free(J)					Frees join AND nested regions (J->free)
// ==================================================================

extern Join *Join_allocate(int,JoinType);
extern Join *Join_create(const Annulus**,int,JoinType);
extern void Join_computeBox(const Join*,Box*);
extern void Join_free(Join*);

// ==================================================================
// Oracles
// ==================================================================

// int Join_ro(Region *R,double *min,double *max,int dim);
// int Join_po(Region *R,double *min,double *max,int dim);

// int Join_oro(Region *R,double *min,double *max,int dim);
// int Join_iro(Region *R,double *min,double *max,int dim);
// int Join_opo(Region *R,double *x,int dim);
// int Join_ipo(Region *R,double *x,int dim);

// ==================================================================
// The oracles
// ==================================================================

static inline int _Join_oro(const Join *J,const double* restrict min, const double* restrict max,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) J->R[k];
		if(_Annulus_ro(A,min,max,dim))
		{
			return 1;
		}
	}

	return 0;
}

static inline int _Join_iro(const Join *J,const double* restrict min,const double* restrict max,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) J->R[k];
		if(!(_Annulus_ro(A,min,max,dim)))
		{
			return 0;
		}
	}

	return 1;
}

static inline int _Join_opo(const Join *J,const double *x,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) J->R[k];
		if(_Annulus_po(A,x,dim))
		{
			return 1;
		}
	}

	return 0;
}

static inline int _Join_ipo(const Join *J,const double *x,int dim)
{
	Annulus* A;
	int k;

	for(k=0; k<J->count; k++)
	{
		A=(Annulus*) J->R[k];
		if(!(_Annulus_po(A,x,dim)))
		{
			return 0;
		}
	}

	return 1;
}

static inline int _Join_ro(const Join *J,const double* restrict min, const double* restrict max, int dim)
{
    return (J->type == innerJoin) ? _Join_iro(J,min,max,dim) : _Join_oro(J,min,max,dim);
}

static inline int _Join_po(const Join *J,const double *x, int dim)
{
    return (J->type == innerJoin) ? _Join_ipo(J,x,dim) : _Join_opo(J,x,dim);
}

// ==================================================================

#endif


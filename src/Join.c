// ==================================================================
// Join.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Join region code.
// ==================================================================

#include "Join.h"
#include "Box.h"
#include "Annulus.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ==================================================================
// Oracles for type Join
// ==================================================================

int Join_oro(Region *R,double *min,double *max,int dim) 
{
	Join *J=(Join*)&R[1];
	return _Join_oro(J,min,max,dim);
}

int Join_iro(Region *R,double *min,double *max,int dim) 
{
	Join *J=(Join*)&R[1];
	return _Join_iro(J,min,max,dim);
}

int Join_opo(Region *R,double *x,int dim)
{
	Join *J=(Join*)&R[1];
	return _Join_opo(J,x,dim);
}

int Join_ipo(Region *R,double *x,int dim)
{
	Join *J=(Join*)&R[1];
	return _Join_ipo(J,x,dim);
}

// ==================================================================
// Construction and destruction
// ==================================================================

Join *Join_allocate(int count,JoinType type)
{
	Join *J;
	int rq;

	rq = sizeof(Join)+sizeof(Annulus*)*count;
	J = (Join*)malloc(rq);
	if(!J) return NULL;

	J->count=count;
	J->type=type;
	for(int i=0;i<count;i++) J->R[i]=NULL;

	return J;
}

Join *Join_create(const Annulus **S,int count,JoinType type)
{
	Join *J = Join_allocate(count,type);
	if(!J) return NULL;
	memcpy(J->R,S,sizeof(Annulus*)*count);

	return J;
}

#ifndef Jess_min
#define Jess_min(x,y) (x<y ? x:y)
#endif

#ifndef Jess_max
#define Jess_max(x,y) (x>y ? x:y)
#endif

void Join_computeBox(const Join* J, Box* B)
{
	int i;
	int r;
	int dim;

	assert(J);
	assert(B);
	
	if(J->count == 0)
	{
		for(i=0;i<MAX_BOX_DIM;i++)
		{
			B->min[i] = 0.0;
			B->max[i] = 0.0;
		}
	}
	else
	{
		dim = J->R[0]->dim;
		for(i=0;i<dim;i++)
		{
			B->min[i] = -INFINITY;
			B->max[i] =  INFINITY;
			for(r=0;r<J->count;r++)
			{
				B->min[i] = Jess_max(B->min[i], J->R[r]->minBox[i]);
				B->max[i] = Jess_min(B->max[i], J->R[r]->maxBox[i]);
			}
		}
	}
}

#undef Jess_min
#undef Jess_max

void Join_free(Join *J)
{
	int k;

	if(J)
	{
		for(k=0; k<J->count; k++)
		{
			if(J->R[k]) Annulus_free(J->R[k]);
		}
		free(J);
	}
}

// ==================================================================



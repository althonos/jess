// ==================================================================
// Join.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Join region code.
// ==================================================================

#include "Join.h"
#include "Annulus.h"
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

Join *Join_create(Annulus **S,int count,JoinType type)
{
	Join *J;
	int rq;

	rq = sizeof(Join)+sizeof(Annulus*)*count;
	J = (Join*)calloc(1,rq);
	if(!J) return NULL;

	J->count=count;
	J->type=type;
	memcpy(J->R,S,sizeof(Annulus*)*count);

	return J;
}

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



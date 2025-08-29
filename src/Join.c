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

Region *Join_create(Region **S,int count,JoinType type)
{
	Region *R;
	Join *J;
	int rq;

	rq = sizeof(Region)+sizeof(Join)+sizeof(Region*)*count;
	R = (Region*)calloc(1,rq);
	J = (Join*)&R[1];
	memcpy(J->R,S,sizeof(Region*)*count);

	R->free=Join_free;
	J->count=count;

	if(type==innerJoin)
	{
		R->intersectionQ = Join_iro;
		R->inclusionQ = Join_ipo;
	}
	else // if type==outerJoin
	{
		R->intersectionQ = Join_oro;
		R->inclusionQ = Join_opo;
	}

	return R;
}

void Join_free(Region *R)
{
	Join *J;
	int k;

	if(R)
	{
		J = (Join*)&R[1];

		for(k=0; k<J->count; k++)
		{
			if(J->R[k]) J->R[k]->free(J->R[k]);
		}

		free(R);
	}
}

// ==================================================================



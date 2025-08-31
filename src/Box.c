// ==================================================================
// Box.c
// Copyright (c) Martin Larralde, 2025
// ==================================================================
// Implementation of the Box region and its oracles
// ==================================================================

#include "Box.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ==================================================================
// Oracles for type Annulus
// ==================================================================

int Box_po(Region *vA, double *x, int d)
{
	Box *B=(Box*)&vA[1];
	return _Box_po(B, x, d);
}

int Box_ro(Region *vA, double *minBox, double *maxBox, int d)
{
	Box *B=(Box*)&vA[1];
	return _Box_ro(B, minBox, maxBox, d);
}

// ==================================================================
// Methods of for regions of type Annulus
// ==================================================================

Region *Box_create(double* a, double* b, int d)
{
	Box *B;
	Region *R;
	int rq;
	double tmp;

    if(d>MAX_BOX_DIM)
        return NULL;

	rq = sizeof(Region)+sizeof(Box);
	R = (Region*)calloc(1,rq);
	B = (Box*)&R[1];
	R->intersectionQ=Box_ro;
	R->inclusionQ=Box_po;
	R->free=Box_free;

	memcpy(B->min,a,sizeof(double)*d);
	memcpy(B->max,b,sizeof(double)*d);
	B->dim=d;

	return R;
}

void Box_free(Region *R)
{
	if(R) free(R);
}

// ==================================================================


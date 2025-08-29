// ==================================================================
// Annulus.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Annulus region and its oracles
// ==================================================================

#include "Annulus.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ==================================================================
// Oracles for type Annulus
// ==================================================================

int Annulus_po(Region *vA, double *x, int d)
{
	Annulus *A=(Annulus*)&vA[1];
	return _Annulus_po(A, x, d);
}

int Annulus_ro(Region *vA, double *minBox, double *maxBox, int d)
{
	Annulus *A=(Annulus*)&vA[1];
	return _Annulus_ro(A, minBox, maxBox, d);
}

// ==================================================================
// Methods of for regions of type Annulus
// ==================================================================

Region *Annulus_create(double *u, double a, double b, int d)
{
	Annulus *A;
	Region *R;
	int rq;
	double tmp;

	rq = sizeof(Region)+sizeof(Annulus)+d*sizeof(double);
	R = (Region*)calloc(1,rq);
	A = (Annulus*)&R[1];
	R->intersectionQ=Annulus_ro;
	R->inclusionQ=Annulus_po;
	R->free=Annulus_free;

	if(b<a)
	{
		tmp=b;
		b=a;
		a=tmp;
	}

	if(a<0.0) a=0.0;
	if(b<0.0) b=0.0;

	A->centre=(double*)&A[1];
	memcpy(A->centre,u,sizeof(double)*d);
	A->min=a*a;
	A->max=b*b;
	A->dim=d;

	return R;
}

void Annulus_free(Region *R)
{
	if(R) free(R);
}

// ==================================================================


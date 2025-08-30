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

Annulus *Annulus_create(double *u, double a, double b, int d)
{
	Annulus *A;
	Region *R;
	int rq;
	double tmp;

	if(d!=ANNULUS_DIM) return NULL;
	A = (Annulus*)malloc(sizeof(Annulus));

	if(b<a)
	{
		tmp=b;
		b=a;
		a=tmp;
	}

	if(a<0.0) a=0.0;
	if(b<0.0) b=0.0;

	memcpy(A->centre,u,sizeof(double)*ANNULUS_DIM);
	A->min=a*a;
	A->max=b*b;
	A->dim=d;

	return A;
}

Annulus *Annulus_reuse(Annulus *A, double *u, double a, double b, int d)
{
	double tmp;

	if(!A) return Annulus_create(u,a,b,d);
	if(d!=ANNULUS_DIM) return NULL;

	if(b<a)
	{
		tmp=b;
		b=a;
		a=tmp;
	}

	if(a<0.0) a=0.0;
	if(b<0.0) b=0.0;

	memcpy(A->centre,u,sizeof(double)*ANNULUS_DIM);
	A->min=a*a;
	A->max=b*b;
	A->dim=d;

	return A;
}

void Annulus_free(Annulus *A)
{
	if(A) free(A);
}

// ==================================================================


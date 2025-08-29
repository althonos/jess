// ==================================================================
// Annulus.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of Annulus construction and destruction.
// ==================================================================

#ifndef ANNULUS_H
#define ANNULUS_H

#include "Region.h"
#include <math.h>

#define ANNULUS_DIM 3

// ==================================================================
// type Annulus
// ==================================================================
// centre				The center of the annulus
// min,max				Limits of the radii
// dim					The dimension of the space
// a,b					Inner and outer radii
// ==================================================================

struct _Annulus
{
	double centre[ANNULUS_DIM];
	double min;
	double max;
	int dim;
};

typedef struct _Annulus Annulus;

// ==================================================================
// Methods for Annulus manipulation
// ==================================================================
// create(u,a,b,d)		Make region {x in R^d : a <= |x-u| <= b }.
// free(A)				Free the region given (or use R->free)
// ==================================================================

extern Annulus *Annulus_create(double*,double,double,int);
extern void Annulus_free(Annulus*);

// ==================================================================
// Local "functions"
// ==================================================================

#define min(x,y) (x<y ? x:y)
#define max(x,y) (x>y ? x:y)

// ==================================================================
// Oracles
// ==================================================================

int Annulus_po(Region *vA, double *x, int d);
int Annulus_ro(Region *vA, double *minBox, double *maxBox, int d);

static inline int _Annulus_po(Annulus *A, double *x, int d)
{
	double tmp,sum;
	int i;

	// Does x lie within annulus A?

	if(A->dim!=d) return 0;

	for(sum=0.0,i=0; i<ANNULUS_DIM; i++)
	{
		tmp = A->centre[i]-x[i];
		sum += tmp*tmp;
	}

	return sum<A->min || sum>A->max ? 0:1;
}

static inline int _Annulus_ro(Annulus *A, double *minBox, double *maxBox, int d)
{
	double minSum;
	double maxSum;
 	double t1,t2;
	double t3,t4;
	int i;

	if(d!=A->dim) return 0;

	// Does the box region [minBox,maxBox] intersect the annulus A?

	minSum=0.0;
	maxSum=0.0;
	for(i=0; i<ANNULUS_DIM; i++)
	{
		t1 = A->centre[i]-minBox[i];
		t2 = A->centre[i]-maxBox[i];
		t1 *= t1;
		t2 *= t2;

		if(minBox[i]>A->centre[i] || maxBox[i]<A->centre[i])
		{
			minSum += min(t1,t2);
		}

		maxSum += max(t1,t2);
	}

	return minSum>A->max || maxSum<A->min ? 0:1;
}

// ==================================================================

#undef min
#undef max
#endif


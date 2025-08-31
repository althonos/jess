// ==================================================================
// Box.h
// Copyright (c) Martin Larralde, 2025
// ==================================================================
// Declaration of Box construction and destruction.
// ==================================================================

#ifndef BOX_H
#define BOX_H

#include "Region.h"
#include <math.h>

#define MAX_BOX_DIM 3

// ==================================================================
// type Annulus
// ==================================================================
// centre				The center of the annulus
// min,max				Limits of the radii
// dim					The dimension of the space
// ==================================================================

struct _Box
{
	double min[MAX_BOX_DIM];
	double max[MAX_BOX_DIM];
	int dim;
};

typedef struct _Box Box;

// ==================================================================
// Methods for Annulus manipulation
// ==================================================================
// create(a,b,d)		Make region {x in R^d : a_i <= x_i <= b_i }.
// free(A)				Free the region given (or use R->free)
// ==================================================================

extern Region *Box_create(double*,double*,int);
extern void Box_free(Region*);

// ==================================================================
// Local "functions"
// ==================================================================

#define min(x,y) (x<y ? x:y)
#define max(x,y) (x>y ? x:y)

// ==================================================================
// Oracles
// ==================================================================

int Box_po(Region *vA, double *x, int d);
int Box_ro(Region *vA, double *minBox, double *maxBox, int d);

static inline int _Box_po(Box *B, double *x, int d)
{
	int i;

	// Does x lie within box B?

	if(B->dim!=d) return 0;

    for(i=0;i<d;i++)
        if(!( (x[i] >= B->min[i]) && (x[i] <= B->max[i]) ))
            return 0;

    return 1;
}

static inline int _Box_ro(Box *B, double *minBox, double *maxBox, int d)
{
	int i;

	if(d!=B->dim) return 0;

	// Does the box region [minBox,maxBox] intersect the box B?

	for(i=0; i<d; i++)
        if(! ((minBox[i]<=B->max[d]) && (B->min[d]<=maxBox[d])) )
            return 1;

	return 1;
}

// ==================================================================

#undef min
#undef max
#endif


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

int Box_po(Region *vA, double *x)
{
	Box *B=(Box*)&vA[1];
	return _Box_po(B, x);
}

int Box_ro(Region *vA, double *minBox, double *maxBox)
{
	Box *B=(Box*)&vA[1];
	return _Box_ro(B, minBox, maxBox);
}

// ==================================================================
// Methods of for regions of type Annulus
// ==================================================================

Box *Box_create(const double* restrict a, const double* restrict b)
{
	Box *B;
	double tmp;

	B = (Box*)malloc(sizeof(Box));
	memcpy(B->min,a,sizeof(double)*BOX_DIM);
	memcpy(B->max,b,sizeof(double)*BOX_DIM);

	return B;
}

void Box_free(Box *B)
{
	if(B) free(B);
}

// ==================================================================


// ==================================================================
// TessTemplate.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of TessTemplate creation.
// ==================================================================

#ifndef TESSTEMPLATE_H
#define TESSTEMPLATE_H

#include "Template.h"
#include "TessAtom.h"
#include <stdio.h>

// ==================================================================
// Forward declarations
// ==================================================================
// TessTemplate			The implementation of type TessTemplate
// ==================================================================

typedef struct _TessTemplate TessTemplate;

// ==================================================================
// Implementation of type TessTemplate
// ==================================================================
// count				Number of atoms in the template
// atom[k]				Ptr to TessAtom k
// distance[i][j]		Distance between atoms i and j
// symbol				The name of the template
// dim					The dimension of the template (# residues)
// ==================================================================

struct _TessTemplate
{
	int count;
	TessAtom **atom;
	double **distance;
	char *symbol;
	int dim;
};

// ==================================================================
// Oracles of type TessTemplate
// ==================================================================

extern int TessTemplate_count(const Template *T);
extern int TessTemplate_match(const Template *T,int k,const Atom *A);
extern int TessTemplate_range(const Template *T,int i,int j,double *a,double *b);
extern const double *TessTemplate_position(const Template *T, int k);
extern double TessTemplate_distWeight(const Template *T, int k);
extern int TessTemplate_check(const Template *T, Atom **A, int k, int ignore_chain);
extern const char *TessTemplate_name(const Template *T);
extern double TessTemplate_logE(const Template *T,double rmsd, int n);

// ==================================================================
// Private methods of type TessTemplate
// ==================================================================

extern void TessTemplate_free(Template *T);

// ==================================================================
// Creation of a TessTemplate object
// ==================================================================
// create(file,s)			Parse file & create template
// ==================================================================

extern Template *TessTemplate_create(FILE*,const char*);
extern Template *TessTemplate_copy(const Template*);

// ==================================================================

#endif


// ==================================================================
// Template.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of abstract interface Template and related types.
// ==================================================================

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "Atom.h"
#include "CandidateSet.h"

// ==================================================================
// Forward declarations
// ==================================================================
// Template				The template interface
// ==================================================================

typedef struct _Template Template;

// ==================================================================
// type Template
// ==================================================================
// free(T)				Free memory associated with template T
// count(T)				Return number of atoms in template
// match(T,k,A)			True if A matches atom k of T
// range(T,i,j,a,b)		[*a,*b] <- range of |atom i - atom j|
// check(T,A,k,ignore_chain)	Check n-ary rules on atom k-1 and 0,...,k-2
// candidates(T,M,k)	Generate a set of candidate atoms of M matching atom k of T
// position(T,i)		Position of atom i (example position)
// name(T)			Returns the symbolic name for the template
// logE(T,x,n)			Provide an estimate of logE for a hit
// ==================================================================

struct _Template
{
	void (*free)(Template*);
	int (*count)(const Template*);
	int (*match)(const Template*,int,const Atom*);
	int (*range)(const Template*,int,int,double*,double*);
	int (*check)(const Template*,Atom**,int*,int,int);
	void (*candidates)(const Template*, const Molecule*, int, CandidateSet**);
	const double *(*position)(const Template*,int);
	const char *(*name)(const Template*);
	double (*logE)(const Template*,double,int);
	double (*distWeight)(const Template*,int);
	Template* (*copy)(const Template*);
};

// ==================================================================

#endif



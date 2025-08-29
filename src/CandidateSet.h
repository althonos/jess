// ==================================================================
// CandidateSet.h
// Copyright (c) Jonathan Barker, 2002
// Copyright (c) Martin Larralde, 2025
// ==================================================================
// Declaration of type CandidateSet.
// ==================================================================

#ifndef CANDIDATESET_H
#define CANDIDATESET_H

#include "Atom.h"
#include "Molecule.h"

// ==================================================================
// Type CandidateSet
// ==================================================================
// count				Number of atoms in the set
// atom[k]				Points to ATOM record for kth candidate
// coord[k]				Points to coordinates for kth candidate
// ==================================================================

struct _CandidateSet
{
	size_t count;
	size_t capacity;
	Atom **atom;
	double **coord;
};

typedef struct _CandidateSet CandidateSet;

// ==================================================================
// Declaration of methods of local type CandidateSet
// ==================================================================
// create(M,T,k)		Create from molecule M, atom k of T
// free(S)				Free candidate set
// ==================================================================

CandidateSet *CandidateSet_create();
CandidateSet *CandidateSet_reuse(CandidateSet *S, const Molecule *M);
void CandidateSet_addAtom(CandidateSet*, Atom*);
void CandidateSet_recordCoordinates(CandidateSet*);
void CandidateSet_free(CandidateSet*);

// ==================================================================
// Type CandidateSetArray
// ==================================================================
// count				Number of atoms in the set
// atom[k]				Points to ATOM record for kth candidate
// coord[k]				Points to coordinates for kth candidate
// ==================================================================

struct _CandidateSetArray
{
	size_t count;
	size_t capacity;
	CandidateSet **items;
};

typedef struct _CandidateSetArray CandidateSetArray;

// ==================================================================
// Declaration of methods of local type CandidateSetArray
// ==================================================================
// create()				Create new candidate set array
// free(C)				Free candidate set array
// get(C,M,k)		    Get a new candidate set for M at index k
// ==================================================================

CandidateSetArray *CandidateSetArray_create();
void CandidateSetArray_free(CandidateSetArray*);
CandidateSet* CandidateSetArray_get(CandidateSetArray *C, int k);

#endif
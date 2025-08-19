// ==================================================================
// CandidateSet.h
// Copyright (c) Jonathan Barker, 2002
// Copyright (c) Martin Larralde, 2025
// ==================================================================
// Implementation of type CandidateSet.
// ==================================================================

#include <stdlib.h>

#include "CandidateSet.h"
#include "Molecule.h"
#include "Template.h"

// ==================================================================
// Methods of local type CandidateSet
// ==================================================================

CandidateSet *CandidateSet_create(const Molecule *M)
{
	CandidateSet *S;
	Atom *A;
	int n = Molecule_count(M);
	int m;

	S = (CandidateSet*)calloc(1,sizeof(CandidateSet));
    if(!S) return NULL;

	S->atom=(Atom**)calloc(n,sizeof(Atom*));
	S->coord = NULL;
	S->count = 0;
    
    return S;
}

void CandidateSet_free(CandidateSet *S)
{
	if(S)
	{
		if(S->atom) free(S->atom);
		if(S->coord) free(S->coord);
		free(S);
	}
}

// ==================================================================

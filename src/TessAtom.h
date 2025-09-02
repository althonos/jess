// ==================================================================
// TessAtom.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type TessAtom and its methods
// ==================================================================

#ifndef TESSATOM_H
#define TESSATOM_H

#include "Atom.h"
#include <stdlib.h>

// ==================================================================
// Forward declarations
// ==================================================================
// TessAtom				A template for a single atom
// ==================================================================

typedef struct _TessAtom TessAtom;

// ==================================================================
// Type TessAtom
// ==================================================================
// code					The code describing the match mode
// resSeq				The residue sequence number
// nameCount			Number of atom names present
// resNameCount			Number of residue names present
// chainID				The chain field
// name[k]				kth atom name alternate
// resName[k]			kth residue name alternate
// pos[k]				kth coordinate of position of atom
// distWeight[k]			kth atom distance threshold modifier (weight)
// ==================================================================

struct _TessAtom
{
	int code;
	int resSeq;
	int nameCount;
	int resNameCount;
	char chainID1; //Riziotis edit
	char chainID2;
	//char chainID; 
	char **name;
	char **resName;
	double pos[3];
	double distWeight;
};

// ==================================================================
// Methods of type TessAtom
// ==================================================================
// create(s)			Create from TESS template record
// free(J)				Free memory associated with J
// position(J)			Return coordinates of J
// match(J,A)			True if A matches J
// resSeq(A)			Return resSeq field of A
// chainID(A)			Return the chain ID of A
// ==================================================================

extern TessAtom *TessAtom_create(const char*);
extern TessAtom* TessAtom_copy(const TessAtom*);
extern int TessAtom_match(const TessAtom*,const Atom*);
extern double TessAtom_distance(const TessAtom*, const TessAtom*);

static inline void TessAtom_free(TessAtom *A)
{
	if(A) free(A);
}

static inline const double *TessAtom_position(const TessAtom *A)
{
	return A->pos;
}

#ifndef TessAtom_toupper
#define TessAtom_toupper(c) ((( c >= 'a' ) && (c <= 'z')) ? c & (~0x20) : c);
#endif

static int TessAtom_compareName(const char* restrict a, const char* restrict b)
{
	for(int i=0; i<4;i++)
	{
		char ca = TessAtom_toupper(a[i]);
		char cb = TessAtom_toupper(b[i]);
		int cmp = ca - cb;
		if(cmp != 0) return cmp;
	}
	return 0;
}


static inline int TessAtom_resSeq(const TessAtom *A)
{
	return A->resSeq;
}

//Riziotis edit

static inline char TessAtom_chainID1(const TessAtom *A)
{
	return A->chainID1;
}

static inline char TessAtom_chainID2(const TessAtom *A)
{
	return A->chainID2;
}

static inline double TessAtom_distWeight(const TessAtom *A)
{
	return A->distWeight;
}

//extern char TessAtom_chainID(const TessAtom*);

static inline const char* TessAtom_resName(const TessAtom* A, int k)
{
	return A->resName[k];
}

static inline int TessAtom_resNameCount(const TessAtom* A)
{
	return A->resNameCount;
}


static inline int TessAtom_code(const TessAtom *A)
{
	return A->code;
}



// ==================================================================

#endif


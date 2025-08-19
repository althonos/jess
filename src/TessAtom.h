// ==================================================================
// TessAtom.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of type TessAtom and its methods
// ==================================================================

#ifndef TESSATOM_H
#define TESSATOM_H

#include "Atom.h"

// ==================================================================
// Forward declarations
// ==================================================================
// TessAtom				A template for a single atom
// ==================================================================

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

typedef struct _TessAtom TessAtom;

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
extern void TessAtom_free(TessAtom*);
extern const double *TessAtom_position(const TessAtom*);
extern int TessAtom_match(const TessAtom*,const Atom*);
extern int TessAtom_resSeq(const TessAtom*);
//Riziotis edit
extern char TessAtom_chainID1(const TessAtom*);
extern char TessAtom_chainID2(const TessAtom*);
extern double TessAtom_distWeight(const TessAtom*);
//extern char TessAtom_chainID(const TessAtom*);

// ==================================================================

#endif


// ==================================================================
// TessTemplate.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of TessTemplate creation and oracles.
// ==================================================================

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include "TessTemplate.h"
#include "TessAtom.h"
#include "Annulus.h"
#include "Join.h"
#include "CandidateSet.h"

#if HAVE_ALLOCA
#include <alloca.h>
#endif

// ==================================================================
// Forward declaration of local types
// ==================================================================
// Node					A node in a linked list of TessAtoms
// TessTemplate			The implementation of type TessTemplate
// ==================================================================

typedef struct _Node Node;

// ==================================================================
// Local type Node
// ==================================================================

struct _Node
{
	TessAtom *atom;
	Node *succ;
};

// ==================================================================
// Oracles of type TessTemplate
// ==================================================================

int TessTemplate_count(const Template *T)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return J->count;
}

int TessTemplate_match(const Template *T,int k,const Atom *A)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return TessAtom_match(J->atom[k],A);
}

int TessTemplate_range(const Template *T,int i,int j,double *a,double *b)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	*a = *b = J->distance[i][j];
	return 1;
}

const double *TessTemplate_position(const Template *T, int k)
{
	const TessTemplate *J=(const TessTemplate*)&T[1];
	return TessAtom_position(J->atom[k]);
}

double TessTemplate_distWeight(const Template *T, int k)
{
	const TessTemplate *J=(const TessTemplate*)&T[1];
	return TessAtom_distWeight(J->atom[k]);
}

int TessTemplate_check(const Template *T, Atom **A, int *order, int k, IgnoreType ignore_chain)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	int i;
	int a,b,c,d,e,f;

	for(i=0; i<k-1; i++)
	{
		// Compare chain ids

		//Riziotis edit
		c = A[order[i]]->chainID2-A[order[k-1]]->chainID2;
		d = TessAtom_chainID2(J->atom[order[i]])-TessAtom_chainID2(J->atom[order[k-1]]);
		//c = A[i]->chainID-A[k-1]->chainID;
		//d = TessAtom_chainID(J->atom[i])-TessAtom_chainID(J->atom[k-1]);

		if(ignore_chain==ignoreNone)
		{
			if(c==0 && d!=0) return 0;
			if(c!=0 && d==0) return 0;
			if(c!=0)
			{
				continue;
			}
		}

		// Compare residue sequence numbers

		a = A[order[i]]->resSeq-A[order[k-1]]->resSeq;
		b = TessAtom_resSeq(J->atom[order[i]])-TessAtom_resSeq(J->atom[order[k-1]]);

		if(a==0 && b!=0) return 0;
		if(a!=0 && b==0) return 0;

		if(ignore_chain!=ignoreAtoms)
		{
			if(a==0 && c!=0) return 0;
		}
	}

	return 1;
}

const char *TessTemplate_name(const Template *T)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	return J->symbol;
}

double TessTemplate_logE(const Template *T,double rmsd, int n)
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	const double logA = -8.5;
	const double alpha = 2.5;
	const double beta = -0.7;

	// Approximate number of hits per molecule atom is:
	// A * rmsd^alpha * J->dim^beta (empirically). This gives
	// us log(expected number of hits in molecule of size n)

	return logA + alpha*rmsd + beta*(double)J->dim + log((double)n);
}

void TessTemplate_candidates(const Template *T, const Molecule *M, int k, CandidateSet** S) 
{
	const TessTemplate *J = (const TessTemplate*)&T[1];
	Atom *A;
	Atom **it;
	int n = Molecule_count(M);
	int m;
	int i;
	int j;
	int code;

	if(!(*S = CandidateSet_reuse(*S, n)))
		return;

	// The following match codes require a match on residue name, so we
	// can use the residue name index to iterate only on atoms from 
	// (one of) the required residue(s).

	code = TessAtom_code(J->atom[k]);
	if((code >= -1) && (code <= 8))
	{
		// WARNING: A template atom may have several residue names that 
		// are actually equal to each other (because TessAtom_parse does not
		// deduplicate), so to avoid the same atom from being selected more
		// than once, we use an array to remember which of the residue
		// names we have already processed.
#ifdef HAVE_ALLOCA
		char* done = (char*) alloca(M->index->n*sizeof(char));
		memset(done, 0, M->index->n*sizeof(char));
#else
		char* done = (char*) calloc(M->index->n, sizeof(char));
#endif
		for (i=0; i<TessAtom_resNameCount(J->atom[k]); i++) {
			const char* resName = TessAtom_resName(J->atom[k], i);
			int	j = ResIndex_find(M->index, resName);
			if((j == -1) || (done[j])) continue;
			done[j] = 1;
			for (it = ResIndex_values(M->index, j); *it != NULL; it++) {
				A = (*it);
				if(TessAtom_match(J->atom[k],A)) CandidateSet_addAtom(*S, A);
			}
		}
#ifndef HAVE_ALLOCA
		free(done);
#endif
	}
	else
	{
		// For remaining match codes, a match on residue name is not required,
		// so we just fallback to the original implementation.
		for (m=0; m<n; m++) {
			A = (Atom*)Molecule_atom(M,m);
			if(TessAtom_match(J->atom[k],A)) CandidateSet_addAtom(*S, A);
		}
	}

	if ((*S)->count > 0)
		CandidateSet_recordCoordinates(*S);

	// return S;
}

// ==================================================================
// Private methods of type TessTemplate
// ==================================================================

void TessTemplate_free(Template *T)
{
	TessTemplate *J;
	int i;

	if(T)
	{
		J=(TessTemplate*)&T[1];
		if(J->symbol) free(J->symbol);

		for(i=0; i<J->count; i++)
		{
			TessAtom_free(J->atom[i]);
		}

		free(T);
	}
}

// ==================================================================
// Creation of type TessTemplate
// ==================================================================

Template *TessTemplate_create(FILE *file,const char *sym)
{
	Template *T;
	TessTemplate *J;
	TessAtom *A;
	Node *head,*n;
	const double *x;
	const double *y;
	double tmp;
	int rq,i,j,k,count;
	char buf[0x100];

	// head is the head of a list of TessAtoms which
	// we create by parsing file. count is the number
	// of atom templates encountered.

	head=NULL;
	count=0;

	// Loop through the entire stream, line by line.

	while(fgets(buf,0x100,file))
	{
		if(strncmp(buf,"ATOM",4)==0 || strncmp(buf,"HETATM",6)==0)
		{
			// Parse a single TESS/Jess atom record...

			A = TessAtom_create(buf);

			// If some errors occur parsing the template
			// free all memory and return NULL.

			if(!A)
			{
				while(head)
				{
					n=head->succ;
					TessAtom_free(head->atom);
					free(head);
					head=n;
				}

				return NULL;
			}

			// No errors thus far - add the atom template
			// to the list we're creating...

			n = (Node*)malloc(sizeof(Node));
			n->succ=head;
			head=n;
			n->atom=A;
			count++;
		}
	}

	if(count==0) return NULL;

	// Set up the memory for the template. This is a
	// Template, followed by a TessTemplate followed
	// by various arrays needed for the TessTemplate.
	// To allow easy destruction it makes sense to
	// allocate it all in one big chunk, then set
	// pointers later...

	rq = sizeof(Template)+sizeof(TessTemplate);
	rq += count*sizeof(TessAtom*);
	rq += count*count*sizeof(double);
	rq += count*sizeof(double*);

	T = (Template*)calloc(1,rq);
	J = (TessTemplate*)&T[1];
	J->atom=(TessAtom**)&J[1];
	J->distance=(double**)&J->atom[count];

	J->distance[0]=(double*)&J->distance[count];
	for(i=1; i<count; i++)
	{
		J->distance[i]=(double*)&J->distance[i-1][count];
	}

	// Set up the method pointers in the Template
	// structure at the head of the object...

	T->free=TessTemplate_free;
	T->match=TessTemplate_match;
	T->position=TessTemplate_position;
	T->count=TessTemplate_count;
	T->range=TessTemplate_range;
	T->check=TessTemplate_check;
	T->name=TessTemplate_name;
	T->logE=TessTemplate_logE;
	T->distWeight=TessTemplate_distWeight;
	T->candidates=TessTemplate_candidates;

	// Set up the data fields

	J->symbol=strdup(sym);
	J->count=count;

	// Create pointers to the atoms found in the
	// template to allow access by atom index rather
	// than iteration through a list.

	for(i=0; i<count; i++)
	{
		n=head->succ;
		J->atom[count-i-1]=head->atom;
		free(head);
		head=n;
	}

	// Compute the distances in the template.

	for(i=0; i<count; i++)
	{
		J->distance[i][i]=0.0;
		x=TessAtom_position(J->atom[i]);

		for(j=i+1; j<count; j++)
		{
			J->distance[i][j]=0.0;
			y=TessAtom_position(J->atom[j]);

			for(k=0; k<3; k++)
			{
				tmp = x[k]-y[k];
				J->distance[i][j] += tmp*tmp;
			}

			J->distance[i][j]=sqrt(J->distance[i][j]);
			J->distance[j][i]=J->distance[i][j];
		}
	}

	// Compute the dimension of the template...

	for(J->dim=0,i=0; i<count; i++)
	{
		// Have we already seen this residue sequence number?

		for(j=0; j<i; j++)
		{
			if(TessAtom_resSeq(J->atom[i])==TessAtom_resSeq(J->atom[j]))
			{
				j=i+10;
			}
		}

		if(j<=i)
		{
			J->dim++;
		}
	}

	return T;
}

Template* TessTemplate_copy(const Template *T)
{
	int count;
	int rq;
	int i,j;
	Template *T2;
	TessTemplate *J2;
	const TessTemplate *J;

	J = (const TessTemplate*)&T[1];
	count = J->count;

	// Allocate memory for the copy

	rq = sizeof(Template)+sizeof(TessTemplate);
	rq += count*sizeof(TessAtom*);
	rq += count*count*sizeof(double);
	rq += count*sizeof(double*);

	T2 = (Template*)calloc(1,rq);
	J2 = (TessTemplate*)&T2[1];

	// Set up the method pointers

	T2->free=TessTemplate_free;
	T2->match=TessTemplate_match;
	T2->position=TessTemplate_position;
	T2->count=TessTemplate_count;
	T2->range=TessTemplate_range;
	T2->check=TessTemplate_check;
	T2->candidates=TessTemplate_candidates;
	T2->name=TessTemplate_name;
	T2->logE=TessTemplate_logE;
	T2->distWeight=TessTemplate_distWeight;
	T2->copy=TessTemplate_copy;

	// Copy atoms

	J2->atom=(TessAtom**)&J2[1];
	for(i=0; i<count; i++)
	{
		J2->atom[i]=TessAtom_copy(J->atom[i]); 
		if (J2->atom[i] == NULL) {
			free(T2);
			return NULL;
		}
	}

	// Copy distances

	J2->distance=(double**)&J2->atom[count];
	J2->distance[0]=(double*)&J2->distance[count];
	for(i=1; i<count; i++)
	{
		J2->distance[i]=(double*)&J2->distance[i-1][count];
	}

	for (i = 0; i < count; i++)
	{
		for (j = 0; j < count; j++)
		{
			J2->distance[i][j] = J->distance[i][j];
		}
	}

	// Copy data fields

	J2->symbol=(J->symbol==NULL)?NULL:strdup(J->symbol);
	J2->count=J->count;
	J2->dim=J->dim;

	return T2;
}

// ==================================================================

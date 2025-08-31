// ================================================================== Jess.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of types Jess and JessQuery.
// ==================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include "Jess.h"
#include "Molecule.h"
#include "Scanner.h"
#include "TessTemplate.h"
#include "Super.h"


// ==================================================================
// Forward declarations of local types
// ==================================================================
// Node					A linked list of templates
// ==================================================================

typedef struct _Node Node;

// ==================================================================
// type Jess
// ==================================================================
// head					The head of the list of nodes
// ==================================================================

struct _Jess
{
	Node *head;
};

// ==================================================================
// type JessQuery
// ==================================================================
// node					The current node
// scanner				The current scanner
// super				The current superposition
// reverseQ				True if superposition is reversed
// molecule				The molecule being scanned
// atoms				Array of Atoms which are hit
// threshold			The distance threshold
// candidates			CandidateSetArray to recycle between scanners
// reorder				Whether to use atom reordering while scanning
// ==================================================================

struct _JessQuery
{
	Node *node;
	Scanner *scanner;
	bool scan;
	Superposition *super;
	int reverseQ;
	Molecule *molecule;
	Atom **atoms;
	double threshold;
	double max_total_threshold;
	ScannerData* scanner_data;
	bool reorder;
};

// ==================================================================
// type Node
// ==================================================================
// template				The template at this node
// next					The next node in the list
// ==================================================================

struct _Node
{
	Template *template;
	Node *next;
};

// ==================================================================
// Methods of type Jess
// ==================================================================

Jess *Jess_create(void)
{
	return (Jess*)calloc(1,sizeof(Jess));
}

void Jess_free(Jess *J)
{
	Node *n;
	Template *T;

	if(J)
	{
		while(J->head)
		{
			n=J->head->next;
			T=J->head->template;
			if(T) T->free(T);
			free(J->head);
			J->head=n;
		}
		free(J);
	}
}

void Jess_addTemplate(Jess *J, Template *T)
{
	Node *n;

	n=(Node*)malloc(sizeof(Node));
	n->template=T;
	n->next=J->head;
	J->head=n;
}

JessQuery *Jess_query(Jess *J, Molecule *M,double t,double s,bool reorder)
{
	JessQuery *Q;

	Q = (JessQuery*)calloc(1,sizeof(JessQuery));
	if(!Q)
		return NULL;

	Q->node=J->head;
	Q->molecule=M;
	Q->threshold=t;
	Q->max_total_threshold=s;
	Q->scan=false;
	Q->reorder=reorder;
	Q->scanner=NULL;
	Q->scanner_data=ScannerData_create();
	if(!Q->scanner_data)
	{
		JessQuery_free(Q);
		return NULL;
	}

	return Q;
}

// ==================================================================
// Methods of type JessQuery
// ==================================================================

void JessQuery_free(JessQuery *Q)
{
	if(Q)
	{
		Scanner_free(Q->scanner);
		Superposition_free(Q->super);
		ScannerData_free(Q->scanner_data);
		free(Q);
	}
}

Template *JessQuery_template(const JessQuery *Q)
{
	if(!Q->node) return NULL;
	return Q->node->template;
}

const Molecule *JessQuery_molecule(const JessQuery *Q)
{
	return Q->molecule;
}

Atom **JessQuery_atoms(const JessQuery *Q)
{
	if(!Q->atoms) return NULL;
	return Q->atoms;
}

Superposition *JessQuery_superposition(const JessQuery *Q)
{
	int i;
	int count;
	Template *T;
	Atom **A;
	Superposition* super;

	A = Q->atoms;
	T = Q->node->template;
	count = T->count(T);
	super=Superposition_create();

	for(i=0; i<count; i++)
	{
		Superposition_align(super,A[i]->x,T->position(T,i));
	}

	return super;
}

int JessQuery_next(JessQuery *Q, int ignore_chain)
{
	Template *T;
	Atom **A;

	while(Q->node)
	{
		Superposition_free(Q->super);
		Q->super=NULL;

		if(!Q->scan)
		{
			Q->scanner=Scanner_reuse(
				Q->scanner,
				Q->molecule,
				Q->node->template,
				Q->scanner_data,
				Q->threshold,
				Q->max_total_threshold,
				Q->reorder
				);

			if(!Q->scanner)
			{	
				Q->node=Q->node->next;
				continue;
			}
			else
			{
				Q->scan=true;
			}
		}

		if((A=Scanner_next(Q->scanner, ignore_chain)))
		{
			Q->atoms=A;
			T = Q->node->template;

			return 1;
		}

		// Scanner_free(Q->scanner);
		// Q->scanner=NULL;
		Q->scan=false;

		Superposition_free(Q->super);
		Q->super=NULL;

		Q->atoms=NULL;
		Q->node=Q->node->next;
	}

	// All done...

	return 0;
}

// ==================================================================


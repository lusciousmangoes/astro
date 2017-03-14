/*------------------------------------------------------
 * data type definitions and access macros
 * ----------------------------------------------------*/
#ifndef _dnc_h_
#define _dnc_h_

/* set the precision */


/* standard includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>

/* vector defs and math operation macros */
#include "vectors.h"

/* node type constants
   bit 0: cell=0/particle=1
   bit 1: massive=0/massless=1 */
#define PARTICLE 1
#define UPDATE_FORCES 4

/* cell node */
#define CELL 0

/*--------------------------------------------------------
   macros to fetch details for both cells and particles */
#define TYPE(n) ((n)->type)
#define POS(n)  ((n)->x)
#define MASS(n) ((n)->mass)
#define NEXT(n) ((n)->NextNode)

/*------------------------------------------------------
   cell-only macros */
#define QUAD(cn) ((cn)->d.c->Q)
#define D1(cn)   ((cn)->a)
#define D2(cn)   ((cn)->d.c->D2)
#define D3(cn)   ((cn)->d.c->D3)
#define CHILD(cn) ((cn)->d.c->subNode)
#define RMAX(cn) ((cn)->d.c->rmax)
#define THETA(cn) ((cn)->d.c->theta)

/*------------------------------------------------------
 particle-only macros */
#define VEL(pn)  ((pn)->d.p->v)
#define ACC(pn)  ((pn)->a)
#define ID(pn)   ((pn)->d.p->id)

/*-------------------------------------------------------
   general node checking macros */
#define IS_PARTICLE(n) (TYPE(n) & PARTICLE)
#define IS_CELL(n) (!IS_PARTICLE(n))
#define UPDATE(n) (TYPE(n) & UPDATE_FORCES)

#define UPDATE_FORCE(n) TYPE(n) |= UPDATE_FORCES

/* --------------------------------------------------------
   particle specific information
   (you can add lots of stuff here without any
    performance impact since it is never moved
    in memory, though macros are recommended to
    access it -- see above)  */
typedef struct {
	float v[3];		/* particle velocity */
	int id;			/* particle id */
} Particle;

/*--------------------------------------------------------
   cell specific information */
typedef struct {
	float Q[6];
	float D2[6];
	float D3[10];
	float theta;
	float rmax;
	struct _Node *subNode;
} Cell;

/*--------------------------------------------------------
 information common to particles and cells */
typedef struct _Node {
	int type;		/* specifies the properties of the node */
	float mass;		/* mass */
	float x[3];		/* position (or COM) */
	float a[3];		/* acceleration (and D1) */
	struct _Node *NextNode;
	union {
		Particle *p;
		Cell *c;
	} d;
} Node;

/*--------------------------------------------------------
   Box structure */
typedef struct {
	float min[3];
	float max[3];
} Box;

#endif

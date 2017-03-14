#include "dnc.h"
#include "build_tree.h"

void calc_interaction(Node *nodeA, Node *nodeB);
void calc_coefficients(Node *nodeA, Node *nodeB, float *R, const float d);
void calc_coefficients_pp(Node *nodeA, Node *nodeB, const float *R);
void calc_coefficients_tp(Node *nodeA, Node *nodeB, const float *R);
void calc_coefficients_pc(Node *nodeA, Node *nodeB, const float *R, const float *d);
void calc_coefficients_tpc(Node *nodeA, Node *nodeB, const float *R, const float *Dn);
void calc_coefficients_cc(Node *nodeA, Node *nodeB, const float *R, const float *d);
void calc_interaction_threaded(Node *root);
void merge_down(Node *n);
void *process_list(void *ptr);

/* local variables  */
float reps2, reps3, reps5, reps7;
extern int bucketSize;  /* delcared in dnc.c */

/*------------------------------------------------
 * get_bounds: get the boundaries for a set of
 *             N particles
 *-----------------------------------------------*/
Box get_bounds(Node *particles, int N)
{
	int i, j;
	Box b;
	Node *p;
	float maxW;

	p = particles;
	VECTOR_CLEAR(b.min);
	VECTOR_CLEAR(b.max);
	for (i = 0; i < N; i++) {
		for (j = 0; j < 3; j++) {
			if (POS(p)[j] < b.min[j]) {
				b.min[j] = POS(p)[j];
			} else if (POS(p)[j] > b.max[j]) {
				b.max[j] = POS(p)[j];
			}
		}
		p++;
	}
	
	/* make bounds cubic */
	/* force the centre of the root node to be at the origin */
	maxW = 0;
	for (i = 0; i < 3; i++) {
		if (fabs(b.min[i]) > maxW)
			maxW = fabs(b.min[i]);
		if (fabs(b.max[i]) > maxW)
			maxW = fabs(b.max[i]);
	}
	maxW = pow(2, floor(log(maxW) / log(2)) + 1);

	for (i = 0; i < 3; i++) {
		b.min[i] = -maxW;
		b.max[i] = +maxW;
	}

	return b;
}

/*------------------------------
 *  calculate forces
 *------------------------------*/
void calc_forces(Node *particles, int N, float eps)
{
	Node rootNode;
	Cell rootCell;
	Box bounds;

	bounds = get_bounds(particles, N);
	
	rootNode.d.c = &rootCell;
	buildTree(&rootNode, &bounds, particles, N, bucketSize);

	reps2 = 1.0 / (eps * eps);
	reps3 = reps2 / eps;
	reps5 = reps3 * reps2;
	reps7 = reps5 * reps2;

	calc_interaction(&rootNode, &rootNode);

	merge_down(&rootNode);

	freeTree(&rootNode);

	return;
}

/*----------------------------------------------------------------------
 * calculate the accelerations of nodeA and nodeB on each other
 *----------------------------------------------------------------------*/
void calc_interaction(Node *nodeA, Node *nodeB)
{
	int well_separated;
	float delta[3], d2, sep;
	Node *n1, *n2;

	/* if neither node needs to be updated, just return */
	if (!(UPDATE(nodeA) || UPDATE(nodeB)))
		return;

	/* see if nodes are equal */
	if (nodeA == nodeB) {
		/* if this world has subworlds, do their interactions */
		if (IS_CELL(nodeA)) {
			n1 = CHILD(nodeA);
			while (n1 != NULL) {
				n2 = n1;
				while (n2 != NULL) {
					calc_interaction(n1, n2);
					n2 = NEXT(n2);
				}
				n1 = NEXT(n1);
			}
		}
		return;
	}

	/* since worlds are different, check their separation */
	VECTOR_SUB(delta, POS(nodeB), POS(nodeA));
	VECTOR_MAG2(d2, delta);
	sep = (IS_CELL(nodeA) ? RMAX(nodeA) : 0) + (IS_CELL(nodeB) ? RMAX(nodeB) : 0);
	sep *= sep;
	well_separated = (d2 >= sep);

	if (well_separated) {
		/* calc taylor coeffs */
		/* potential at A due to B */
		calc_coefficients(nodeA, nodeB, delta, d2);
	}
	/* otherwise, split larger node (note: if one is a particle, it will
	   be smaller since rmax = 0) */
	else {
		/* store the larger world in n1, smaller (or particle) in n2 */
		if (IS_PARTICLE(nodeA)) {
			n1 = nodeB;
			n2 = nodeA;
		} else if (IS_PARTICLE(nodeB)) {
			n1 = nodeA;
			n2 = nodeB;
		} else if (RMAX(nodeA) > RMAX(nodeB)) {
			n1 = nodeA;
			n2 = nodeB;
		} else {
			n1 = nodeB;
			n2 = nodeA;
		}

		n1 = CHILD(n1);
		while (n1 != NULL) {
			calc_interaction(n1, n2);
			n1 = NEXT(n1);
		}
	}

}

void calc_coefficients(Node *nodeA, Node *nodeB, float *R, const float d2)
{
	float Dn[4];
	float rr, rr2, rr3, rr5, rr7, f2;

	f2 = d2 * reps2;

	if (f2 < 1) {

		/* softening using F2 kernel -- this D's are continuous
		   (but not D3 derivative) */
		Dn[1] = -0.125 * (15. * f2 * f2 - 42. * f2 + 35.) * reps3;
		Dn[2] = 1.5 * (7. - 5. * f2) * reps5;
		Dn[3] = -15. * reps7;

	} else {
		/* unsoftened values */
		RECIP_ROOT(rr, d2);
		rr2 = rr * rr;
		rr3 = rr2 * rr;
		rr5 = rr3 * rr2;
		rr7 = rr5 * rr2;

		Dn[1] = -rr3;
		Dn[2] = 3.0 * rr5;
		Dn[3] = -15.0 * rr7;
	}

	/* cell/cell interaction */
	if (IS_CELL(nodeA) && IS_CELL(nodeB))
		calc_coefficients_cc(nodeA, nodeB, R, Dn);

	/* cell/particle interaction */
	else if (IS_CELL(nodeA)) {
		/* switching order of nodes so switch sign  */
		VECTOR_MUL(R, -1);
		calc_coefficients_pc(nodeB, nodeA, R, Dn);
	}
	
	/* particle/cell interaction */
	else if (IS_CELL(nodeB))
		calc_coefficients_pc(nodeA, nodeB, R, Dn);

	/* particle/particle interaction */
	else {
		VECTOR_MUL(R, Dn[1]);
		calc_coefficients_pp(nodeA, nodeB, R);
	}
	return;
}

/* calculate coefficients for a particle/particle interaction */
void calc_coefficients_pp(Node *nodeA, Node *nodeB, const float *R)
{
	VECTOR_MUL_ADD(D1(nodeA), MASS(nodeB), R);
	VECTOR_MUL_ADD(D1(nodeB), -MASS(nodeA), R);
}

/* calculate coefficients for a particle/cell interaction
   nodeA is the particle, nodeB is the cell */
void calc_coefficients_pc(Node *nodeA, Node *nodeB, const float *R, const float *Dn)
{
	float *Q;
	float D[10], Dtwo[6];
	float D1Q[3];
	float Rxx, Rxy, Rxz, Ryy, Ryz, Rzz;
	float RxxD, RyyD, RzzD;

	Rxx = R[X] * R[X];
	Rxy = R[X] * R[Y];
	Rxz = R[X] * R[Z];
	Ryy = R[Y] * R[Y];
	Ryz = R[Y] * R[Z];
	Rzz = R[Z] * R[Z];

	/* calculate Bij */
	Dtwo[XX] = Dn[1] + Rxx * Dn[2];	/* D2xx */
	Dtwo[XY] = Rxy * Dn[2];	/* D2xy */
	Dtwo[XZ] = Rxz * Dn[2];	/* D2xz */
	Dtwo[YY] = Dn[1] + Ryy * Dn[2];	/* D2yy */
	Dtwo[YZ] = Ryz * Dn[2];	/* D2yz */
	Dtwo[ZZ] = Dn[1] + Rzz * Dn[2];	/* D2zz */
	/* particles do not have a D2 term but cells do */

	RxxD = Rxx * Dn[3] + Dn[2];
	RyyD = Ryy * Dn[3] + Dn[2];
	RzzD = Rzz * Dn[3] + Dn[2];
	
	/* calculate Cijk */
	D[XXX] = R[X] * (2.0 * Dn[2] + RxxD);	/* D3xxx */
	D[XXY] = R[Y] * RxxD;	/* D3xxy */
	D[XXZ] = R[Z] * RxxD;	/* D3xxz */
	D[XYY] = R[X] * RyyD;	/* D3xyy */
	D[XYZ] = R[X] * Ryz * Dn[3];	/* D3xyz */
	D[XZZ] = R[X] * RzzD;	/* D3xzz */
	D[YYY] = R[Y] * (2.0 * Dn[2] + RyyD);	/* D3yyy */
	D[YYZ] = R[Z] * RyyD;	/* D3yyz */
	D[YZZ] = R[Y] * RzzD;	/* D3yzz */
	D[ZZZ] = R[Z] * (2.0 * Dn[2] + RzzD);	/* D3zzz */
	/* particles do not need a D3 term but cells do */

	/* calculate Ai for worldA */
	Q = QUAD(nodeB);
	D1Q[X] = R[X] * Dn[1] + 0.5 * (Q[XX] * D[XXX] + Q[YY] * D[XYY] + Q[ZZ] * D[XZZ])
	         + Q[XY] * D[XXY] + Q[XZ] * D[XXZ] + Q[YZ] * D[XYZ];	
	D1Q[Y] = R[Y] * Dn[1] + 0.5 * (Q[XX] * D[YXX] + Q[YY] * D[YYY] + Q[ZZ] * D[YZZ])
	         + Q[XY] * D[YXY] + Q[XZ] * D[YXZ] + Q[YZ] * D[YYZ];	
	D1Q[Z] = R[Z] * Dn[1] + 0.5 * (Q[XX] * D[ZXX] + Q[YY] * D[ZYY] + Q[ZZ] * D[ZZZ])
	         + Q[XY] * D[ZXY] + Q[XZ] * D[ZXZ] + Q[YZ] * D[ZYZ];

	VECTOR_MUL_ADD(D1(nodeA), MASS(nodeB), D1Q);

	/* calculate Ai for worldB (Q = 0 for a particle) */
	VECTOR_MUL_ADD(D1(nodeB), -Dn[1] * MASS(nodeA), R);
	VECTOR6_MUL_ADD(D2(nodeB), -MASS(nodeA), Dtwo);
	VECTOR10_MUL_ADD(D3(nodeB), -MASS(nodeA), D);
}

/* calculate coefficients for a cell/cell interaction */
void calc_coefficients_cc(Node *nodeA, Node *nodeB, const float *R, const float *Dn)
{
	float D[10], Dtwo[6];
	float D1Q[3];
	float Rxx, Rxy, Rxz, Ryy, Ryz, Rzz;
	float RxxD, RyyD, RzzD;
	float twoD2;

	Rxx = R[X] * R[X];
	Rxy = R[X] * R[Y];
	Rxz = R[X] * R[Z];
	Ryy = R[Y] * R[Y];
	Ryz = R[Y] * R[Z];
	Rzz = R[Z] * R[Z];

	/* calculate Bij */
	Dtwo[XX] = Dn[1] + Rxx * Dn[2];
	Dtwo[XY] = Rxy * Dn[2];
	Dtwo[XZ] = Rxz * Dn[2];
	Dtwo[YY] = Dn[1] + Ryy * Dn[2];
	Dtwo[YZ] = Ryz * Dn[2];
	Dtwo[ZZ] = Dn[1] + Rzz * Dn[2];

	RxxD = Rxx * Dn[3] + Dn[2];
	RyyD = Ryy * Dn[3] + Dn[2];
	RzzD = Rzz * Dn[3] + Dn[2];
	twoD2 = 2.0 * Dn[2];
	
	/* calculate Cijk */
	D[XXX] = R[X] * (twoD2 + RxxD);
	D[XXY] = R[Y] * RxxD;
	D[XXZ] = R[Z] * RxxD;
	D[XYY] = R[X] * RyyD;
	D[XYZ] = R[X] * Ryz * Dn[3];
	D[XZZ] = R[X] * RzzD;
	D[YYY] = R[Y] * (twoD2 + RyyD);
	D[YYZ] = R[Z] * RyyD;
	D[YZZ] = R[Y] * RzzD;
	D[ZZZ] = R[Z] * (twoD2 + RzzD);

	/* calculate Ai for worldA */
	D1Q[X] = R[X] * Dn[1] + 0.5 * (QUAD(nodeB)[XX] * D[XXX] + QUAD(nodeB)[YY] * D[XYY] +
	         QUAD(nodeB)[ZZ] * D[XZZ]) + QUAD(nodeB)[XY] * D[XXY] + 
	         QUAD(nodeB)[XZ] * D[XXZ] + QUAD(nodeB)[YZ] * D[XYZ];
	D1Q[Y] = R[Y] * Dn[1] + 0.5 * (QUAD(nodeB)[XX] * D[YXX] + QUAD(nodeB)[YY] * D[YYY] +
		 QUAD(nodeB)[ZZ] * D[YZZ]) + QUAD(nodeB)[XY] * D[YXY] + 
		 QUAD(nodeB)[XZ] * D[YXZ] + QUAD(nodeB)[YZ] * D[YYZ];
	D1Q[Z] = R[Z] * Dn[1] + 0.5 * (QUAD(nodeB)[XX] * D[ZXX] + QUAD(nodeB)[YY] * D[ZYY] +
		 QUAD(nodeB)[ZZ] * D[ZZZ]) + QUAD(nodeB)[XY] * D[ZXY] + 
		 QUAD(nodeB)[XZ] * D[ZXZ] + QUAD(nodeB)[YZ] * D[ZYZ];

	VECTOR_MUL_ADD(D1(nodeA), MASS(nodeB), D1Q);
	VECTOR6_MUL_ADD(D2(nodeA), -MASS(nodeB), Dtwo);
	VECTOR10_MUL_ADD(D3(nodeA), MASS(nodeB), D);

	/* calculate Ai for worldB */
	D1Q[X] = R[X] * Dn[1] + 0.5 * (QUAD(nodeA)[XX] * D[XXX] + QUAD(nodeA)[YY] * D[XYY] +
		 QUAD(nodeA)[ZZ] * D[XZZ]) + QUAD(nodeA)[XY] * D[XXY] + 
		 QUAD(nodeA)[XZ] * D[XXZ] + QUAD(nodeA)[YZ] * D[XYZ];
	D1Q[Y] = R[Y] * Dn[1] + 0.5 * (QUAD(nodeA)[XX] * D[YXX] + QUAD(nodeA)[YY] * D[YYY] +
		 QUAD(nodeA)[ZZ] * D[YZZ]) + QUAD(nodeA)[XY] * D[YXY] + 
		 QUAD(nodeA)[XZ] * D[YXZ] + QUAD(nodeA)[YZ] * D[YYZ];
	D1Q[Z] = R[Z] * Dn[1] + 0.5 * (QUAD(nodeA)[XX] * D[ZXX] + QUAD(nodeA)[YY] * D[ZYY] +
	         QUAD(nodeA)[ZZ] * D[ZZZ]) + QUAD(nodeA)[XY] * D[ZXY] + 
	         QUAD(nodeA)[XZ] * D[ZXZ] + QUAD(nodeA)[YZ] * D[ZYZ];

	VECTOR_MUL_ADD(D1(nodeB), -MASS(nodeA), D1Q);
	VECTOR6_MUL_ADD(D2(nodeB), -MASS(nodeA), Dtwo);
	VECTOR10_MUL_ADD(D3(nodeB), -MASS(nodeA), D);
}

//------------------------------------------------------------------------------------
// this combines the coefficients from the current node, down into the children nodes
//------------------------------------------------------------------------------------
void merge_down(Node *n)
{
	float delta[3];
	Node *sn;

	/* propagate the taylor coefficients of the current node
	   down into the subnodes */
	sn = CHILD(n);
	while (sn != NULL) {
		VECTOR_SUB(delta, POS(sn), POS(n));

		/* both cells and particles need D1 */
		D1(sn)[X] +=
		    D1(n)[X] + delta[X] * D2(n)[XX] + delta[Y] * D2(n)[XY] +
		    delta[Z] * D2(n)[XZ]
		    + 0.5 * (delta[X] * delta[X] * D3(n)[XXX] +
			     delta[Y] * delta[Y] * D3(n)[XYY] +
			     delta[Z] * delta[Z] * D3(n)[XZZ]) +
		    delta[X] * delta[Y] * D3(n)[XXY] +
		    delta[X] * delta[Z] * D3(n)[XXZ] +
		    delta[Y] * delta[Z] * D3(n)[XYZ];

		D1(sn)[Y] +=
		    D1(n)[Y] + delta[X] * D2(n)[YX] + delta[Y] * D2(n)[YY] +
		    delta[Z] * D2(n)[YZ]
		    + 0.5 * (delta[X] * delta[X] * D3(n)[YXX] +
			     delta[Y] * delta[Y] * D3(n)[YYY] +
			     delta[Z] * delta[Z] * D3(n)[YZZ]) +
		    delta[X] * delta[Y] * D3(n)[YXY] +
		    delta[X] * delta[Z] * D3(n)[YXZ] +
		    delta[Y] * delta[Z] * D3(n)[YYZ];

		D1(sn)[Z] +=
		    D1(n)[Z] + delta[X] * D2(n)[ZX] + delta[Y] * D2(n)[ZY] +
		    delta[Z] * D2(n)[ZZ]
		    + 0.5 * (delta[X] * delta[X] * D3(n)[ZXX] +
			     delta[Y] * delta[Y] * D3(n)[ZYY] +
			     delta[Z] * delta[Z] * D3(n)[ZZZ]) +
		    delta[X] * delta[Y] * D3(n)[ZXY] +
		    delta[X] * delta[Z] * D3(n)[ZXZ] +
		    delta[Y] * delta[Z] * D3(n)[ZYZ];

		/* if subnode is a cell, it needs D2 and D3 as well */
		if (IS_CELL(sn)) {
			/* D2 */ 
			D2(sn)[XX] +=
			    D2(n)[XX] + delta[X] * D3(n)[XXX] +
			    delta[Y] * D3(n)[XXY] + delta[Z] * D3(n)[XXZ];
			D2(sn)[XY] +=
			    D2(n)[XY] + delta[X] * D3(n)[XYX] +
			    delta[Y] * D3(n)[XYY] + delta[Z] * D3(n)[XYZ];
			D2(sn)[XZ] +=
			    D2(n)[XZ] + delta[X] * D3(n)[XZX] +
			    delta[Y] * D3(n)[XZY] + delta[Z] * D3(n)[XZZ];
			D2(sn)[YY] +=
			    D2(n)[YY] + delta[X] * D3(n)[YYX] +
			    delta[Y] * D3(n)[YYY] + delta[Z] * D3(n)[YYZ];
			D2(sn)[YZ] +=
			    D2(n)[YZ] + delta[X] * D3(n)[YZX] +
			    delta[Y] * D3(n)[YZY] + delta[Z] * D3(n)[YZZ];
			D2(sn)[ZZ] +=
			    D2(n)[ZZ] + delta[X] * D3(n)[ZZX] +
			    delta[Y] * D3(n)[ZZY] + delta[Z] * D3(n)[ZZZ];

			/* D3 */
			VECTOR10_INCREMENT(D3(sn), D3(n));
			/* since this is a cell, propagate down as well */
			merge_down(sn);
		}

		sn = NEXT(sn);
	}
}

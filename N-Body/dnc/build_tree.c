#include "dnc.h"

/* local prototypes */
const float calc_split_value(Box *bounds, Node *particles, int N, int dim);
void split_particles(Node *p, int N, float *mid, int dim, int *nb, int *na);
Node *get_new_node(int N);
void set_rmax(Node *n, Box *b);
float calc_theta(float m);
int split_node(Node *parent, Box *bounds, Node *particles, int N);

/* global variables */
extern float theta;
extern int maxdepth;

/* local variables */
Node **NodeBuf = NULL;
Cell **CellBuf = NULL;
int nBuf;
int curBuf, nUsed;
const int bufSize = 1000;
int depth = 0;
int bucketSize;
float mTotal;
int lastCellCount = 0;

/*-----------------------------------------------------------------
 * build_tree: make the tree
 *----------------------------------------------------------------*/
int buildTree(Node *parent, Box *bounds, Node *particles, int N, int bSize)
{
	int nC;
	int i;

	/* initialize the rootnode */
	TYPE(parent) = CELL;
	NEXT(parent) = NULL;
	VECTOR_CLEAR(D1(parent));
	VECTOR6_CLEAR(QUAD(parent));
	VECTOR6_CLEAR(D2(parent));
	VECTOR10_CLEAR(D3(parent));
	
	/* get the total mass of the particles */
	mTotal = 0;
	for (i = 0; i < N; i++) {
		mTotal += MASS(particles + i);
	}
	
	bucketSize = bSize;
	nC = split_node(parent, bounds, particles, N);
	lastCellCount = nC;
	return nC;
}

/*-----------------------------------------------------------------
 * linkNode: link a node to it's parent
 *----------------------------------------------------------------*/
Node *linkNode(Node *curnode, Node *lastnode, Node *parent)
{
	/* if this is the first node, it must be the child of the parent */
	if (lastnode == NULL) {
		CHILD(parent) = curnode;
		NEXT(curnode) = NULL;
		return curnode;
	} else {
		NEXT(lastnode) = curnode;
		NEXT(curnode) = NULL;
		return curnode;
	}
}

/*-----------------------------------------------------------------
 * split_node: split the parent node into max of 8 child nodes
 * NOTE: split_node is called recursively so there is a minimum
 *       of variable declarations
 *----------------------------------------------------------------*/
int split_node(Node *parent, Box *bounds, Node *particles, int N)
{
	int i, j, d;
	int offset, octcount[8];
	int nCell, nReal;
	float mid[3];
	Box newBounds;
	Node *node, *lastNode;

	depth++;
	if (depth > maxdepth)
		maxdepth = depth;
		
	/* calculate the split values */
	for (i = 0; i < 3; i++) {
		mid[i] = calc_split_value(bounds, particles, N, i);
	}
	
	octcount[0] = N;
	for (i = 0; i < 3; i++) {
		d = 1 << (3 - i);
		offset = 0;
		for (j = 0; j < 8; j += d) {
			split_particles(particles + offset, octcount[j], mid, i,
					octcount + j, octcount + j + (d >> 1));
			offset += octcount[j] + octcount[j + (d >> 1)];
		}

	}

	offset = 0;
	nCell = 0;
	lastNode = NULL;
	for (i = 0; i < 8; i++) {
		if (octcount[i] == 0)
			continue;

		/* see how many particles are in this node */
		nReal = octcount[i];

		/* link in particle to the tree */
		if (nReal <= bucketSize) {
			node = particles + offset;

			/* link all of these particles into the tree */
			for (j = 0; j < octcount[i]; j++) {
				lastNode = linkNode(node, lastNode, parent);
				node++;
			}
		} else {
			/* get bounds for new node */
			newBounds = *bounds;
			for (j = 0; j < 3; j++) {
				if (i & (1 << (2 - j)))
					newBounds.min[j] = mid[j];
				else
					newBounds.max[j] = mid[j];
			}

			/* get a new node */
			node = get_new_node(1);

			nCell++;
			nCell += split_node(node, &newBounds, particles + offset, octcount[i]);

			/* link this cell into the tree */
			lastNode = linkNode(node, lastNode, parent);
		}

		offset += octcount[i];
	}

	/* next, calculate the COM  */
	node = CHILD(parent);
	VECTOR_CLEAR(POS(parent));
	MASS(parent) = 0;
	while (node != NULL) {
		VECTOR_MUL_ADD(POS(parent), MASS(node), POS(node));
		MASS(parent) += MASS(node);
		node = NEXT(node);
	}
	VECTOR_DIV(POS(parent), MASS(parent));

	/* and the quad. moment for this node */
	VECTOR6_CLEAR(QUAD(parent));
	node = CHILD(parent);
	while (node != NULL) {
		VECTOR_SUB(mid, POS(parent), POS(node));
		if (IS_CELL(node)) {
			QUAD(parent)[XX] += MASS(node) * (QUAD(node)[XX] + mid[X] * mid[X]);
			QUAD(parent)[XY] += MASS(node) * (QUAD(node)[XY] + mid[X] * mid[Y]);
			QUAD(parent)[XZ] += MASS(node) * (QUAD(node)[XZ] + mid[X] * mid[Z]);
			QUAD(parent)[YY] += MASS(node) * (QUAD(node)[YY] + mid[Y] * mid[Y]);
			QUAD(parent)[YZ] += MASS(node) * (QUAD(node)[YZ] + mid[Y] * mid[Z]);
			QUAD(parent)[ZZ] += MASS(node) * (QUAD(node)[ZZ] + mid[Z] * mid[Z]);
		} else {
			QUAD(parent)[XX] += MASS(node) * mid[X] * mid[X];
			QUAD(parent)[XY] += MASS(node) * mid[X] * mid[Y];
			QUAD(parent)[XZ] += MASS(node) * mid[X] * mid[Z];
			QUAD(parent)[YY] += MASS(node) * mid[Y] * mid[Y];
			QUAD(parent)[YZ] += MASS(node) * mid[Y] * mid[Z];
			QUAD(parent)[ZZ] += MASS(node) * mid[Z] * mid[Z];
		}
		node = NEXT(node);
	}
	VECTOR6_DIV(QUAD(parent), MASS(parent));

	/* next, set the update flag for this node */
	node = CHILD(parent);
	while (node != NULL) {
		if (UPDATE(node)) {
			UPDATE_FORCE(parent);
			break;
		}
		node = NEXT(node);
	}

	/* finally, calculate the opening criteria for this node */
	THETA(parent) = calc_theta(MASS(parent));
	set_rmax(parent, bounds);

	depth--;

	return nCell;
}

/*---------------------------------
 * calc_theta: determine theta
 *-------------------------------*/
float calc_theta(float m)
{
	float C;
	float t, lastt;
	float TOL = 1e-4;
	int i = 0;
	
	if (theta < 0)
		return -theta;
		
	C = theta * pow(1.0 - theta, -0.4) * pow(m / mTotal, -1.0 / 15);

	t = theta;
	do {
		lastt = t;
		/* simple Newton's method to find root */
		t -= (t - pow(1.0 - t, 0.4) * C) / (1.0 + 0.4 * C * pow(1.0 - t, -0.6));
		i++;
	} while (fabs(t - lastt) > TOL);

	return t;
}

/*----------------------------------------------------
 * calc_split_value: determine the value to split the
 *                   node at
 *----------------------------------------------------*/
/* geometrical node splitting
   just divide each dimension in half */
const float calc_split_value(Box *bounds, Node *particles, int N, int dim)
{
	return 0.5 * (bounds->min[dim] + bounds->max[dim]);
}

/*---------------------------------------------------------------
 * split_particles: split the particles above and below mid[dim]
 *---------------------------------------------------------------*/
void split_particles(Node *p, int N, float *mid, int dim, int *nb, int *na)
{
	int i, j;
	Node tmp;
	
	i = 0;
	j = N - 1;
	while (i <= j) {
		while (POS(p + i)[dim] < mid[dim]) {
			i++;
			if (i == N)
				break;
		}
		while (POS(p + j)[dim] >= mid[dim]) {
			j--;
			if (j < 0)
				break;
		}
		if (i < j) {
			tmp = p[i];
			p[i] = p[j];
			p[j] = tmp;
		}
	}
	*nb = i;
	*na = N - i;
	
	return;
}

/*------------------------------------------------------
 * get_new_node(): return a pointer to a new node
 * This allocates memory in blocks and then doles
 * it out as required
 * (should guess at initial amount needed and allocate that much)
 *----------------------------------------------------*/
Node *get_new_node(int N)
{
	Node *b;
	
	/* initialize the buffer space */
	if (NodeBuf == NULL) {
		if (lastCellCount == 0) {
			nBuf = 100;
		} else {
			nBuf = ceil(1.0 * lastCellCount / bufSize);
		}
		NodeBuf = (Node **) calloc(nBuf, sizeof(Node *));
		NodeBuf[0] = (Node *) calloc(bufSize, sizeof(Node));
		CellBuf = (Cell **) calloc(nBuf, sizeof(Cell *));
		CellBuf[0] = (Cell *) calloc(bufSize, sizeof(Cell));
		nUsed = 0;
		curBuf = 0;
	}
	/* not enough room in current buffer */
	if ((nUsed + N) > bufSize) {
		/*need to alloc next buffer */
		curBuf++;
		
		/* do we need more buffers */
		if (curBuf >= nBuf) {
			/* the first tree construction may be slow since many
			   expansions are needed, but afterwords, the
			   size isn't going to change much */
			nBuf++;
			NodeBuf = (Node **) realloc(NodeBuf, sizeof(Node *) * nBuf);
			CellBuf = (Cell **) realloc(CellBuf, sizeof(Cell *) * nBuf);
		}
		if ((NodeBuf == NULL) || (CellBuf == NULL)) {
			printf("Unable to allocate enough memory for buffer index. Aborting...\n");
			exit(9);
		}

		NodeBuf[curBuf] = (Node *) calloc(bufSize, sizeof(Node));
		CellBuf[curBuf] = (Cell *) calloc(bufSize, sizeof(Cell));
		if ((NodeBuf[curBuf] == NULL) || (CellBuf[curBuf] == NULL)) {
			printf("Unable to allocate enough memory new buffer. Aborting...\n");
			exit(9);
		}
		b = NodeBuf[curBuf];
		b->d.c = CellBuf[curBuf];
		nUsed = N;
	}
	/* otherwise, there is enough space */
	else {
		b = NodeBuf[curBuf] + nUsed;
		b->d.c = CellBuf[curBuf] + nUsed;
		nUsed += N;
	}

	TYPE(b) = CELL;
	VECTOR_CLEAR(D1(b));
	VECTOR6_CLEAR(QUAD(b));
	VECTOR6_CLEAR(D2(b));
	VECTOR10_CLEAR(D3(b));

	return b;
}

/*------------------------------------------------------
 * set_rmax: calculates opening radius criteria for 
 * specified node
 *----------------------------------------------------*/
void set_rmax(Node *n, Box *b)
{
	int i;
	float d, d1;
	float bmax;
	float mid;
	float delta[3];
	Node *node;

	/* calculate max distance from com to a node corner */
	d = 0;
	for (i = 0; i < 3; i++) {
		mid = 0.5 * (b->max[i] + b->min[i]);
		d1 = (POS(n)[i] < mid ? b->max[i] - POS(n)[i] : POS(n)[i] - b->min[i]);
		d += d1 * d1;
	}
	bmax = sqrt(d);

	/* next, calculate it the second way */
	d = 0;
	node = CHILD(n);
	while (node != NULL) {
		VECTOR_SUB(delta, POS(node), POS(n));
		VECTOR_MAG2(d1, delta);

		d1 = sqrt(d1);
		if (IS_CELL(node)) {
			d1 += RMAX(node) * THETA(node);
		}
		if (d1 > d)
			d = d1;
		node = NEXT(node);
	}

	if (d < bmax) {
		RMAX(n) = d / THETA(n);
	} else {
		RMAX(n) = bmax / THETA(n);
	}
}

/* free all the buffers used */
void freeTree(Node * rootnode)
{
	int i;

	for (i = 0; i <= curBuf; i++) {
		free(NodeBuf[i]);
		free(CellBuf[i]);
	}
	free(NodeBuf);
	free(CellBuf);
	nBuf = 0;
	NodeBuf = NULL;
	CellBuf = NULL;
	nUsed = 0;
	curBuf = 0;

	CHILD(rootnode) = NULL;
	NEXT(rootnode) = NULL;
}

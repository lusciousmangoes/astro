#include "dnc.h"
#include <sys/stat.h>

double StartTime = 0.0;

Node *read_particles(char *fname, int *Np)
{
	struct stat finfo;
	Particle *p;
	Node *n;
	FILE *fp;
	int i, N, fsize, dumpint;

	fp = fopen(fname, "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: Unable to open input file %s\n", fname);
		exit(2);
	}
	
	/* Get file size */
	stat(fname, &finfo);
	fsize = finfo.st_size;

	N = (fsize - sizeof(int) - sizeof(double)) / (7 * sizeof(float));
	n = (Node *) calloc(N, sizeof(Node));
	p = (Particle *) calloc(N, sizeof(Particle));

	fread(&dumpint, sizeof(int), 1, fp);
	fread(&StartTime, sizeof(double), 1, fp);
	printf("Starting at time %f\n", StartTime);

	for (i = 0; i < N; i++) {
		/* link the particle to the node */
		n[i].d.p = p + i;
		
		/* read in the details */
		fread(&(MASS(n + i)), sizeof(float), 1, fp);
		fread(POS(n + i), sizeof(float), 3, fp);
		fread(VEL(n + i), sizeof(float), 3, fp);
		
		TYPE(n + i) = PARTICLE;
		ID(n + i) = i;

		/* assume all forces are to be calculated */
		UPDATE_FORCE(n + i);
	}
	
	fclose(fp);
	*Np = N;
	return n;
}

int cmp_id(const void *_n1, const void *_n2)
{
	Node *n1 = (Node *) _n1;
	Node *n2 = (Node *) _n2;
	if (ID(n1) > ID(n2))
		return 1;
	else
		return -1;

}

void save_particles(Node * p, int N, int step, char *base, float dt)
{
	double temp_t = StartTime + dt * step;
	char fname[1000];
	FILE *fp;
	int i;
	
	sprintf(fname, "%s.%04d", base, step);
	printf("Saving to %s\n", fname);
	fflush(stdout);
	
	fp = fopen(fname, "w");
	
	/* make sure original particle order is kept */
	qsort(p, N, sizeof(Node), cmp_id);

	fwrite(&N, sizeof(int), 1, fp);
	fwrite(&temp_t, sizeof(double), 1, fp);
	for (i = 0; i < N; i++) {
		fwrite(&(MASS(p + i)), sizeof(float), 1, fp);
		fwrite(POS(p + i), sizeof(float), 3, fp);
		fwrite(VEL(p + i), sizeof(float), 3, fp);
	}
	
	fclose(fp);

}

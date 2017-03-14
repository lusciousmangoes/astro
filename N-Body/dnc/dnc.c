#include "dnc.h"
#include "calc_forces.h"
#include "particle_io.h"
#include "move.h"

#include <time.h>
#include <sys/timeb.h>
#include <signal.h>
#include <getopt.h>
#include <string.h>

/* global variables */
float theta, reciptheta;
int maxdepth, bucketSize;
int nThread;
int dump_and_exit;


void print_help()
{
    	fprintf(stderr, "Usage: dnc [options]\n");
    	fprintf(stderr, "Options:\n");
    	fprintf(stderr, "   -i file     Initial data file (binary) (init.dat)\n");
    	fprintf(stderr, "   -o file     Output file pattern (output)\n");
    	fprintf(stderr, "   -t #        Length of timestep (0.1)\n");
    	fprintf(stderr, "   -s #        Timestep to stop simulation at (100)\n");
    	fprintf(stderr, "   -a #        Opening angle (0.5)\n");
    	fprintf(stderr, "   -e #        Softening length (0.1)\n");
    	fprintf(stderr, "   -b #        Size of tree bucket (1)\n");
    	fprintf(stderr, "   -d #        Output timestep interval (100)\n");
}

void set_dump_and_exit(int sig)
{
	fprintf(stderr, "Finishing current step..... ");
	fprintf(stderr, "(Once more to stop now)\n");
	dump_and_exit = 1;
	signal(sig, SIG_DFL);
	return;
}

/*-----------------------------------------------
 * M A I N : 
 *----------------------------------------------*/
int main(int argc, char *argv[])
{
	int N, i;
	Node *particles;
	float dt, eps;
	int step, startStep, stopStep;
	int outInterval;
	char outputfile[200], inputfile[200];

	struct timeb tstart2f, tstop2f;
	float tinterval;
	
	int option;

	/* default input parameters */
	sprintf(inputfile, "init.dat");
	sprintf(outputfile, "output");
	dt = 0.1;
	stopStep = 100;
	theta = 0.5;
	eps = 0.1;
	bucketSize = 1;
	outInterval = 100;

	printf("Welcome to Dave's N-Body Code\n");
	printf("-----------------------------\n");
	printf("This is a C implementation of Dehnen's (2002) ");
	printf("multipole expansion algorithm.\n\n");

	if (argc < 2) {
		print_help();
		exit(9);
	} else {
		/* process flags */
		while ((option = getopt(argc, argv, "hi:o:t:s:a:e:b:d:r")) != EOF){
			switch (option){
				case 'i':
					strncpy(inputfile, optarg, 200);
					break;				
				case 'o':
					strncpy(outputfile, optarg, 200);
					break;				
				case 't':
					dt = atof(optarg);
					break;					
				case 's':
					stopStep = atoi(optarg);
					break;					
				case 'a':
					theta = atof(optarg);
					break;					
				case 'e':
					eps = atof(optarg);
					break;					
				case 'b':
					bucketSize = atoi(optarg);
					break;					
				case 'd':
					outInterval = atoi(optarg);
					break;		
				case 'h':
				case '?':
					print_help();
					exit(9);					
			}
		}
		
		particles = read_particles(inputfile, &N);
	}

	fprintf(stderr, "Using single precision\n");

	printf("Theta:\t%f\n", theta);
	printf("eps: \t%f\n", eps);
	printf("dt: \t%f\n", dt);
	printf("bSize:\t%d\n", bucketSize);

	dump_and_exit = 0;
	signal(SIGINT, set_dump_and_exit);
	signal(SIGTERM, set_dump_and_exit);

	/* ----------------------------------------------------------
	   M A I N L O O P S T A R T S H E R E
	   ---------------------------------------------------------- */	   
	startStep = 0;
	for (step = startStep; step <= stopStep; step++) {

		ftime(&tstart2f);

		if (step > startStep)
			AdvancePositions(particles, N, dt);

		for (i = 0; i < N; i++) {
			/* if adding external potential, set it here */
			VECTOR_CLEAR(ACC(particles + i));
		}
		
		calc_forces(particles, N, eps);
		
		if (step > startStep)
			AdvanceVelocities(particles, N, 0.5 * dt);

		ftime(&tstop2f);
		tinterval = 1.0 * (tstop2f.time - tstart2f.time) + 
		            0.001 * (tstop2f.millitm - tstart2f.millitm);
		printf("Step %d: \t%.3f s, %.1f particles/s\n",
		       step, tinterval, (tinterval > 0 ? 1.0 * N / tinterval : 0));

		if (step == startStep) {
			float amax, a;

			amax = 0;
			for (i = 0; i < N; i++) {
				VECTOR_MAG2(a, ACC(particles + i));
				if (a > amax)
					amax = a;
			}
			amax = sqrt(amax);
			printf("Maximum acceleration in %e\n", amax);
			printf("Time step should be about %e\n", 0.1 * sqrt(eps / amax));

		} else if ((!(step % outInterval)) || (step == stopStep)
			   || (dump_and_exit)) {
			save_particles(particles, N, step, outputfile, dt);
			if (dump_and_exit)
				exit(1);
		}

		AdvanceVelocities(particles, N, 0.5 * dt);
	}

	return 0;
}


#include "dnc.h"

void AdvancePositions(Node * p, int N, float dt)
{
	int i;
	
	for (i = 0; i < N; i++) {
		VECTOR_MUL_ADD(POS(p + i), dt, VEL(p + i));
	}
}

void AdvanceVelocities(Node * p, int N, float dt)
{
	int i;

	for (i = 0; i < N; i++) {
		VECTOR_MUL_ADD(VEL(p + i), dt, ACC(p + i));
	}
}

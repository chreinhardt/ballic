/*
 ** This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*granite->nTableMax) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy. We lookup u of v from the table
	** and print it to a file for debugging splint(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	// Try 100 only
	int nTableMax = 100;
	double v, u;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableMax: %i\n", granite->nTableMax);

	v = 0.0;
	u = 0.0;

	i = (granite->vmax/granite->delta)*0.5;
	i = 10;
//	for (i=0;i<granite->nTableMax;i++)
//	{
		for (j=0;j<granite->nTableMax;j++)
		{
			//v = j*granite->delta
			// Choose a point in the middle of the interval
			v = (j + 0.5)*granite->delta;
			u = tillSplineIntv(granite, v, i);

			fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);

			printf("%g %g %g\n",i*granite->delta, v, u);
		}
		printf("\n");
//	}

	tillFinalizeMaterial(granite);
}

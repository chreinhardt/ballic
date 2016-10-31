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
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double v, u;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableRho: %i, nTableV: %i\n", granite->nTableRho,granite->nTableV);

	v = 0.0;
	u = 0.0;

	i = (granite->vmax/granite->dv)*0.5;
	i = 10;
//	for (i=0;i<granite->nTableRho;i++)
//	{
		for (j=0;j<granite->nTableV;j++)
		{
			//v = j*granite->dv
			// Choose a point in the middle of the interval
			v = (j + 0.5)*granite->dv;
			u = tillSplineIntv(granite, v, i);

			fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);

			printf("%g %g %g\n",i*granite->drho, v, u);
		}
		printf("\n");
//	}

	tillFinalizeMaterial(granite);
}

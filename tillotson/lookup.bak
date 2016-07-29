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

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableMax = 1000;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableMax: %i\n", granite->nTableMax);

	/* Print the lookup table to a file. */
	for (i=0;i<granite->nTableMax;i++)
	{
		printf("%.30f",j*granite->delta);
	
		for (j=0;j<granite->nTableMax;j++)
		{
			printf(" %.30f",j*granite->delta,granite->Lookup[j*granite->nTableMax+i]);
		}
		printf("\n");
	}

	/* Solve for the cold curve.
	isentrope = tillSolveIsentrope(granite,0);
	for (j=0;j<granite->nTable;j++)
	{
		printf("%.30f %.30f %.30f\n",j*granite->delta,isentrope[j].u,isentrope[j].rho);
	}
	*/

	/* Debug the function tillColdULookup().
	u = tillColdULookup(granite,rho);
	printf("rho: %.30f u: %.30f u: %g\n", rho, tillColdULookup(granite, rho),u);
	*/
	tillFinalizeMaterial(granite);
}

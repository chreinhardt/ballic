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

#define INDEX(i, j) (((i)*granite->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy. We store the second derivative
	** of u with respect to rho in udv2 just for debugging
	** splint(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
<<<<<<< HEAD
	double rhomax = 100.0;
=======
	double rhomax =25.0;
>>>>>>> master
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;

	double rho, u;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	/* Create an output file for the look up table */
	FILE *fp = NULL;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	tillInitSplineRho(granite);

	fprintf(stderr, "Done.\n");

	rho = 0.0;
	u = 0.0;

	/*
	** Print the look up table to a file first.
	*/	
//#if 0

	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		fprintf(fp,"%g",rho);
		for (j=0;j<granite->nTableV;j+=1)
		{
			// v = j*granite->dv
			u = granite->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
//#endif

	j = 51;
	 	
<<<<<<< HEAD
	for (j=0;j<granite->nTableMax;j++)
	{
		for (i=0;i<granite->nTableMax;i++)
=======
//	for (j=0;j<granite->nTableV;j++)
//	{
		for (i=0;i<granite->nTableRho-1;i++)
>>>>>>> master
		{
			//rho = granite->Lookup[INDEX(i,j)].rho;
			// Choose a point in the middle of the interval
			rho = granite->Lookup[INDEX(i,j)].rho + 0.5*granite->drho;
//			fprintf(stderr,"%g  ", rho);
//			rho = (i+0.5)*granite->drho;
//			fprintf(stderr,"%g\n", rho);
			u = tillSplineIntrho(granite, rho, j);

			printf("%g %g %g\n", rho, u, tillCubicIntRho(granite, rho, j));
		}
<<<<<<< HEAD
		printf("\n");
	}
=======
//		printf("\n");
//	}
>>>>>>> master

	tillFinalizeMaterial(granite);
}

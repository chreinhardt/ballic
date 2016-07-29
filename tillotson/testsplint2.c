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
	** the internal energy. We now evolve a given point
	** (rho1, u1) along an isentrope to (rho2,u2) to debug
	** tillCubicInt(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// Need about 1000x1000 to get good results
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u, v;
	double rho1, u1, rho2, u2;

	int i = 0;
	int j = 0;
	double k = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableRho: %i nTableV: %g\n", granite->nTableRho, granite->nTableV);

	v = 0.0;
	u1 = 0.0;
	u2 = 0.0;

	/* If we set rho = rho0 we start on an isentrope (for u1=v). */
	rho1 = granite->rho0;
	rho2 = 1e-2;

	
	/* Create an output file for the look up table */
	FILE *fp = NULL;

	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);
	
	/* Print the look up table first */	
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

#if 0	
	/* Interpolate values between the isentropes */
	for (i=0;i<granite->nTableMax-1;i+=1)
	{
		// Middle of the interval (i,i+1)
		rho = (i + 0.5)*granite->drho;
		printf("%g", rho); 
		for (j=0;j<granite->nTableMax-1;j+=1)
		{
			// Middle of the interval (i,i+1)
			// v = (j + 0.5)*granite->dv;
			k = 0.1;
			while (k < 0.9)
			{
				v = (j + k)*granite->dv;
				u = tillCubicInt(granite, rho, v);

				//fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);
				printf("  %g", u);
				k+=0.1;
			}
		}
		printf("\n");
	}
#endif
	
	/* Do the isentropic evolution. */
	u1 = 0.0;
	for (j=0;j<granite->nTableV-1;j+=10)
	{
		/* Set u1 to v_i */
		u1 = j*granite->dv;
		fprintf(stderr,"Step %i: rho1=%g u1=%g rho2=%g\n",j,rho1,u1,rho2);

		/* From (rho1,u1)=(rho0,v) to (rho2,u2). */
		u2 = tillLookupU(granite, rho1, u1, rho2, 0);
//		printf("%g  %g  %g  %g\n", rho1, u1, rho2, u2);
		printf("%g  %g  %g  %g ", rho1, u1, rho2, u2);
		/* And back. */
		u1 = 0.0;
		u1 = tillLookupU(granite, rho2, u2, rho1, 0);
		printf("%g  %g  %g  %g\n", rho2, u2, rho1, u1);
	}

	tillFinalizeMaterial(granite);
}

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
	** the internal energy. We lookup u of rho and v from
	** the table and print it to a file to debug tillCubicInt(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	// vmax=25.0, rhomax=100.0, nTableV=10, nTableRho=4000
	//int nTableMax = 1000;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, v, u;

	int i = 0;
	int j = 0;
	int n = 1;
	double k = 0.0;
	double l = 0.0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, n);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	
	rho = 0.0;
	v = 0.0;
	u = 0.0;

	/* Create an output file for the look up table */
	FILE *fp = NULL;

	
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
	
	/* Interpolate values between the isentropes */
	for (i=0;i<granite->nTableRho-1;i+=1)
	{
		// Middle of the interval (i,i+1)
		//rho = (i + 0.5)*granite->drho;
		l = 0.0;
		while (l < 0.9)
		{
			// Try
			rho = (i + l)*granite->drho;
			//rho = granite->Lookup[INDEX(i, 0)].rho;
			//rho += l*fabs((granite->Lookup[INDEX(i, 0)].rho-granite->Lookup[INDEX(i+1, 0)].rho));

			printf("%g", rho);
			for (j=0;j<granite->nTableV-1;j+=1)
			{
				// Middle of the interval (i,i+1)
				// v = (j + 0.5)*granite->dv;
				k = 0.5;
				while (k < 0.9)
				{
					// This does not work for non uniform steps in v
					//v = (j + k)*granite->dv;
					v = granite->vmax/pow(granite->nTableV-1,n)*pow(j+k,n);

					u = tillCubicInt(granite, rho, v);

					//fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);
					printf("  %.8g", u);
					k+=0.5;
				}
			}
		printf("\n");
		l+=0.5;
		}
	}
	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

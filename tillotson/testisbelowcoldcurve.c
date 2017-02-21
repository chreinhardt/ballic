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
	** Debug tillIsBelowColdCurve(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 800.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	// vmax=25.0, rhomax=100.0, nTableV=10, nTableRho=4000
	//int nTableMax = 1000;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, v, u;
	double umax;

	int i = 0;
	int j = 0;
	int n = 1;

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
	u = 0.0;

	/* Create an output file for the look up table */
	FILE *fp = NULL;
	
	/*
	** Print the look up table to a file first.
	*/	
#if 0
	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		fprintf(fp,"%g",rho);
		for (j=0;j<granite->nTableV;j+=1)
		{
			u = granite->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
#endif

	/* Test if different points (rho,u) are below the cold curve */
	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		for (j=0;j<granite->nTableV;j+=1)
		{
			u = granite->Lookup[INDEX(i, j)].u;
			if (!tillIsInTable(granite, rho, u))
			{
				printf("tillIsInTable: value (%g,%g) below not in the look up table (iMat=%i)!\n",rho,u,granite->iMaterial);
				break;
			}

			fprintf(stderr,"i=%i j=%i rho=%g u=%g",i,j,rho,u);
			fprintf(stderr," uc=%g",tillColdULookup(granite, rho));
			fprintf(stderr," bBelowColdCurve=%i\n",tillIsBelowColdCurve(granite, rho,u));
			printf("  %i", tillIsBelowColdCurve(granite,rho,u));
		}
		printf("\n");
	}

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

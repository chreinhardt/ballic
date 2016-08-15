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

#define INDEX(i, j) ((i*mat->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 26.32;
//	double vmax = 28.0;
	double rho, u;
	int nTableRho = 1000;
	int nTableV = 1000;

	int i = 0;
	int j = 0;

	TILLMATERIAL *mat;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	mat = tillInitMaterial(BASALT, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(mat);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", mat->rhomax, mat->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", mat->nTableRho, mat->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", mat->drho, mat->dv);
	fprintf(stderr,"\n");

//#if 0
	/* Create an output file for the look up table */
	FILE *fp = NULL;

	/*
	** Print the look up table to a file first.
	*/	
	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup_cold.txt","w");
	assert(fp != NULL);

	/* Print the lookup table to a file. */
	for (i=0;i<mat->nTableRho;i+=1)
	{
		/* j=0 corresponds to the cold curve. */
		rho = i*mat->drho;
		// v = j*mat->dv
		u = mat->Lookup[INDEX(i, 0)].u;
		fprintf(fp,"%g  %g\n",rho,u);
	}
	fclose(fp);
//#endif
	for (i=0;i<mat->nTableRho;i++)
	{
		// Lookup(i, j) = Lookup(rho, v)
		printf("%.30f",mat->Lookup[INDEX(i, 0)].rho);
	
		for (j=0;j<mat->nTableV;j++)
		{
			printf(" %.30f",mat->Lookup[INDEX(i,j)].u);
		}
		printf("\n");
	}

	/* Solve for the cold curve.
	isentrope = tillSolveIsentrope(mat,0);
	for (j=0;j<mat->nTable;j++)
	{
		printf("%.30f %.30f %.30f\n",j*mat->delta,isentrope[j].u,isentrope[j].rho);
	}
	*/

	/* Debug the function tillColdULookup().
	u = tillColdULookup(mat,rho);
	printf("rho: %.30f u: %.30f u: %g\n", rho, tillColdULookup(mat, rho),u);
	*/
	tillFinalizeMaterial(mat);
}

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

#define INDEX(i, j) ((i*granite->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	** Calculate the isentrope for a given (rho,u) using tillLookupU().
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
//	double vmax = 28.0;
	double rho, u, v;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rhos, us, drho;
	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	/*
	** For now we just integrate from a given rho to zero.
	*/
	if (argc != 3) {
		fprintf(stderr,"Usage: calcisentrope <rho> <u>\n");
		exit(1);
	}
	rhos = atof(argv[1]);
	us = atof(argv[2]);

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"\n");

	fprintf(stderr,"Starting from rho=%g u=%g\n", rhos, us);
//#if 0
	/* Create an output file for the look up table */
	FILE *fp = NULL;

	/*
	** Print the look up table to a file first.
	*/	
	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	/* Print the lookup table to a file. */
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

#if 0
	for (i=0;i<granite->nTableRho;i++)
	{
		// Lookup(i, j) = Lookup(rho, v)
		printf("%.30f",granite->Lookup[INDEX(i, 0)].rho);
	
		for (j=0;j<granite->nTableV;j++)
		{
			printf(" %.30f",granite->Lookup[INDEX(i,j)].u);
		}
		printf("\n");
	}
#endif

	fprintf(stderr,"Solving for isentrope\n");
	drho = (rhos-0.0)/100.0;

	rho = rhos;
	u = us;
#if 0
	while (rho >= drho)
	{
		if (rho < drho) break;
		u = tillLookupU(granite,rho,u,rho-drho,0);
		rho-=drho;
		printf("%g  %g\n",rho, u);
	}
#endif

	/* Find which isentrope (rhos,us) is on. */
	v = tillFindEntropyCurve(granite,rhos,us,0);

	rho = rhos;
	u = us;
	
	while (rho >= 0.0)
	{
		if (rho < 0.0) break;
		u = tillFindUonIsentrope(granite,v,rho);
		printf("%g  %g\n",rho, u);
		rho-=drho;
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

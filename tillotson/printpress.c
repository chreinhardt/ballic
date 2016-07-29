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
	** Print the derivatives of u (in rho and v) used for the cubic
	** spline interpolator.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	int nTableRho = 1000;
	int nTableV = 100;
	double rho, v, u;

	double k = 0.0;
	double l = 0.0;

	int i = 0;
	int j = 0;
	int n = 2;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"n: %i\n", granite->n);
	
	/* Create an output file */
	FILE *fp = NULL;
	
	/*
	** Print the look up table to a file first.
	*/	
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

	/*
	** Now P(rho,u).
	*/	
	fp = fopen("printpress.txt","w");
	assert(fp != NULL);
#if 0
	for (i=0;i<2;i+=1)
	{
		for (j=0;j<granite->nTableV;j+=1)
		{
			rho = i*granite->drho;
			v = granite->vmax/pow(granite->nTableV-1,n)*pow(j,n);
			fprintf(fp,"rho=%g, v=%g, u=%g, ", rho, v, granite->Lookup[INDEX(i, j)].u);
			fprintf(fp,"u1=%g, ", granite->Lookup[INDEX(i, j)].u1);
			fprintf(fp,"udv2=%g, ", granite->Lookup[INDEX(i, j)].udv2);
			fprintf(fp,"u1dv2=%g, ", granite->Lookup[INDEX(i, j)].u1dv2);
			fprintf(fp,"P=%g\n", tillPressure(granite, rho, granite->Lookup[INDEX(i, j)].u));
		}
		fprintf(fp,"\n");
	}
#endif
	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		fprintf(fp,"%g",rho);
		for (j=0;j<granite->nTableV;j+=1)
		{
			fprintf(fp,"  %g", tillPressure(granite, rho, granite->Lookup[INDEX(i, j)].u));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	/*
	** Now print the pressure at u(rho=0).
	*/
	fp = fopen("printpress.i0.txt","w");
	assert(fp != NULL);

	for (i=0;i<1;i+=1)
	{
		fprintf(stderr,"i=%i\n",i);
		for (j=0;j<granite->nTableV-1;j+=1)
		{
			rho = granite->Lookup[INDEX(i, j)].rho;
			rho = 0.0; 
			u = granite->Lookup[INDEX(i, j)].u; 
			fprintf(fp,"%g %g  %g\n", rho, u, tillPressure(granite, rho, u));
		}
	}
	fclose(fp);

	/*
	** Now interpolate values on the grid.
	*/

	fp = fopen("testsplint.txt","w");
	assert(fp != NULL);

	for (i=0;i<granite->nTableRho-1;i+=1)
	{
		// Middle of the interval (i,i+1)
		rho = (i + 0.5)*granite->drho;
		fprintf(fp,"%g", rho);
		for (j=0;j<granite->nTableV-1;j+=1)
		{
			// Middle of the interval (i,i+1)
			//v = (j + 0.5)*granite->dv;
			v = granite->vmax/pow(granite->nTableV-1,n)*pow(j+0.5,n);
				
	//		fprintf(stderr,"j=%i,v=%g,n=%i,vmax=%g\n",j,v,n,granite->vmax);
			u = tillCubicInt(granite, rho, v);
		
			//fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fprintf(stderr,"Done.\n");
//	tillFinalizeMaterial(granite);
}

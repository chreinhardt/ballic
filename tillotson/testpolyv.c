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
	** Try how uneven spacing in dv works.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 40.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	int nTableRho = 100;
	int nTableV = 100;
	double rho, v, u;

	double k = 0.0;
	double l = 0.0;

	int i = 0;
	int j = 0;
	int n = 5;

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
#if 0
	double *v_array;

	v_array=malloc(nTableV*sizeof(double));

	rho = 0.0;
	v = 0.0;
	u = 0.0;

	for (i=0; i < nTableV; i++)
	{
		v_array[i] = vmax/pow(nTableV-1,n)*pow(i,n);
		printf("%i   %g\n", i, v_array[i]);
	}

	// Now find i so that v(i) <= v <= v(i+1)
	v = vmax*0.33;
	
	i = floor(pow(v/vmax,1.0/n)*(nTableV-1));

	fprintf(stderr,"v=%g, i=%i, v(i)=%g, v(i+1)=%g\n", v, i, v_array[i], v_array[i+1]);

	exit(0);
#endif
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
		rho = (i + 0.5)*granite->drho;
		printf("%g", rho);
		for (j=0;j<granite->nTableV-1;j+=1)
		{
			// Middle of the interval (i,i+1)
			// v = (j + 0.5)*granite->dv;
			v = granite->vmax/pow(granite->nTableV-1,n)*pow(j+0.5,n);
				
	//		fprintf(stderr,"j=%i,v=%g,n=%i,vmax=%g\n",j,v,n,granite->vmax);
			u = tillCubicInt(granite, rho, v);
		
			//fprintf(stderr,"i: %i, j: %i, v: %g, u: %g\n",i,j,v,u);
			printf("  %g", u);
		}
		printf("\n");
	}

	fprintf(stderr,"Done.\n");
//	tillFinalizeMaterial(granite);
}

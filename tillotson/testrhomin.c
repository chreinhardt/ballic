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
	** Check if rhomin is set properly and the integration boundaries work. 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 1200.0;
	// For vmax=25, rhomax=25,  nTableV=100,  nTableRho=1000 we get excellent results.
	// Or  vmax=800,rhomax=100, nTableV=1000, nTableRho=1000
	int nTableRho = 1000;
	int nTableV = 100;
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
	fprintf(stderr,"rhomin: %g, rhomax: %g, vmax: %g \n", granite->rhomin, granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"\n");

	/* Make sure that rho0 is indeed on our grid */	
	fprintf(stderr,"rho0: %g, n: %i, n*drho: %g, rho0-i*drho: %g \n",granite->rho0, granite->n, granite->n*granite->drho, granite->rho0-granite->n*granite->drho);
	fprintf(stderr,"\n");

	fprintf(stderr,"Look up table:\n");

	for (j=0;j<granite->nTableV;j+=1)
	{
//		printf("j=%i :",j);
//		printf("n=%i rho(%i,%i)=%g rho0=%g\n",granite->n, granite->n, j, granite->Lookup[INDEX(granite->n,j)].rho, granite->rho0);
	}

	for (i=0;i<granite->nTableRho;i+=1)
	{
//		printf("j=%i :",j);
		printf("i=%i rho(%i,0)=%g rho0=%g, u1=%g\n", i, i, granite->Lookup[INDEX(i,99)].rho, granite->rho0,granite->Lookup[INDEX(i,99)].u1);
	}
#if 0
	for (i=0;i<granite->nTableRho;i+=1)
	{
		for (j=0;j<granite->nTableV;j+=1)
		{
//			printf("j=%i :",j);
			printf("%g ",j, granite->Lookup[INDEX(i,j)].rho);
		}
		printf("\n");
	}
#endif
	/*
	** Print the look up table to a file first.
	*/	
	exit(1);
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
					printf("  %g", u);
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

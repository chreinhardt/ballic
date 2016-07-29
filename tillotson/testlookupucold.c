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
	** Debug tillColdULookup(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 25.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u;

	int i = 0;
	int j = 0;
	int n = 1;
	double l = 0.0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(IRON, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, n);
	
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

// (CR) Added this to check if u = uc +cv*T works
	rho = 10.0;
	u = 30.0;
	double P = tillPressure(granite, rho, u);
	double T = tillTempRhoU(granite, rho, u);

	fprintf(stderr,"Starting values: rho=%g u=%g P=%g T=%g\n",rho,u,P,T);

//	rho = 0.0;
	u = 0.0;
	P = 0.0;
//	T = 0.0;

	double rho2 = rho;
	double u2 = tillColdULookup(granite,rho2) + granite->cv*T;
	double P2 = tillPressure(granite, rho2, u2);
	double T2 = tillTempRhoU(granite, rho2, u2);

	fprintf(stderr,"Results: rho2=%g u2=%g P2=%g T2=%g\n",rho2,u2,P2,T2);
	
	exit(1);
// Just delete this code later

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
	/* Interpolate values along the cold curve (v=0) */
	for (i=0;i<granite->nTableRho-1;i+=1)
	{
		// Middle of the interval (i,i+1)
		//rho = (i + 0.5)*granite->drho;
		l = 0.0;
		while (l < 0.9)
		{
			rho = (i + l)*granite->drho;
			
			u = tillColdULookup(granite, rho);
			printf("%g  %g\n",rho, u);
			l+=0.2;
		}
	}
	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

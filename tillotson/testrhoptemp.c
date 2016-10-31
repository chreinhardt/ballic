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
	** Check, if the function tillRhoPTemp() produces reasonable
	** results.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	int nTableRho = 1000;
	int nTableV = 1000;
	double rho, u, P, T;
	double rhoPT;

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
	u = 0.0;

	rho = granite->rhomin*1.05;
	u = 10.0;

	P = tillPressure(granite, rho, u);
	T = tillTempRhoU(granite, rho, u);

	fprintf(stderr,"rho=%g u=%g P=%g T=%g\n",rho,u,P,T);
	rhoPT = tillRhoPTemp(granite, P, T);

	/* Create an output file for the look up table */
	FILE *fp = NULL;

	/*
	** Print the look up table to a file first.
	*/	
//#if 0

	//sprintf(achFile,"%s.log",msrOutName(msr));
	//fp = fopen("testrhoptemp.txt","w");
	//assert(fp != NULL);

	rho = granite->rhomin*1.05;

	while (rho < granite->rhomax*0.99)
	{
			u=6.0;
			while (u<30.0)
			{
				P = tillPressure(granite, rho, u);
				T = tillTempRhoU(granite, rho, u);
				fprintf(stderr,"rho=%g u=%g P=%g T=%g\n",rho,u,P,T);
				rhoPT = tillRhoPTemp(granite, P, T);
//				fprintf(fp,"%.8g %.8g %.8g %.8g\n",rho,u,P,T);
				if (fabs(rhoPT-rho)>1e-3) printf("rho=%g u=%g P=%g T=%g rhoPT=%g\n",rho,u,P,T,rhoPT);
				u +=1.0;
			}
			rho += granite->drho;
	}

#if 0
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
#endif
	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

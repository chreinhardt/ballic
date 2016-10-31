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
	** Find the region where the pressure becomes negative in the expanded states.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 1200.0;
	// For vmax=rhomax=25 and nTableV=100, nTableRho=1000 we get excellent results.
	int nTableRho = 1000;
	int nTableV = 100;
	double rho, a, b, c, Pa, Pb, Pc;

	int i = 0;
	int j = 0;
	int n = 2;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(BASALT, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, n);
	
	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"n: %i\n", granite->n);
	
	/*
	** Find where P<0 for 0 < rho < rho0.
	*/
	rho=TILL_RHO_MIN;

#ifdef TILL_PRESS_MELOSH
	/* In this case the pressure is set to zero for eta<0.8. */
	fprintf(stderr,"TILL_PRESS_MELOSH defined!\n");
	exit(1);
#endif
	
	while (rho <= granite->rho0)
	{
		/* Do bisection to find where P<0. */
		a = 1e-10;
		b = granite->us2;
		c = 0.0;

		Pa = tillPressureNP(granite, rho, a);
		Pb = tillPressureNP(granite, rho, b);
		Pc = 0.0;

		while (b-a > 1e-14)
		{
//			printf("rho=%g a=%g b=%g c=%g Pa=%g Pb=%g Pc=%g\n",rho,a,b,c,Pa,Pb,Pc);
			c = 0.5*(a + b);
			Pc = tillPressureNP(granite, rho, c);

			if (Pc < 0)
			{
				// Set a = c
				a = c;
				Pa = Pc;
			} else {
				// Set b = c
				b = c;
				Pb = Pc;
			}
		}

		printf("%g  %g\n", rho, c);
		rho += 0.01;
	}

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

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

#define INDEX(i, j) (((i)*granite->nTableMax) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy. We store the second derivative
	** of u with respect to rho in udv2 just for debugging
	** spline(). 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableMax = 1000;
	double rho, u, P, c2, udrho2, d2udrho2_n, d2udrho2_a;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableMax: %i\n", granite->nTableMax);
	
	for (j=0;j<granite->nTableMax;j++)
	{
		for (i=0;i<granite->nTableMax;i++)
		{
			// For the analyic expression
			rho = granite->Lookup[INDEX(i,j)].rho;
			u = granite->Lookup[INDEX(i,j)].u;
			P = tillPressureSound(granite, rho, u, &c2);

			// Calculate d2udrho2 numerically
			// d2y/dx2 = (u[j+1] - 2*u[j] + u[j-1])/(dx*dx)

			if (i > 0 && i < granite->nTableMax-1)
			{
				d2udrho2_n = (granite->Lookup[INDEX(i+1,j)].u - 2*granite->Lookup[INDEX(i,j)].u +granite->Lookup[INDEX(i-1,j)].u)/(granite->delta*granite->delta);
			} else {
				d2udrho2_n = 0;
			}
			
			// Calculate d2udrho2 analytically
			// d2u/drho2 = 1/rho^2 * (c^2 - 2*P/rho)
			d2udrho2_a = 1.0/(rho*rho)*(c2 - 2.0*P/rho);

			// We temporarily store d2u/drho2 in udv2
			udrho2 = granite->Lookup[INDEX(i,j)].udv2;

			// For small rho the analytic expression deviates a lot from the numerical solution!
			printf("%i %i %g %g %g %g %g\n", i, j, rho, u, udrho2, d2udrho2_n, d2udrho2_a);
		}
		printf("\n");
	}
	
	/* Print the lookup table to a file. */
//	for (i=0;i<granite->nTableMax;i++)
//	{
//		// Lookup(i, j) = Lookup(rho, v)
//		printf("%.30f",granite->Lookup[INDEX(i, 0)].rho);
//	
//		for (j=0;j<granite->nTableMax;j++)
//		{
//			printf(" %.30f",granite->Lookup[INDEX(i,j)].u);
//		}
//		printf("\n");
//	}

	tillFinalizeMaterial(granite);
}

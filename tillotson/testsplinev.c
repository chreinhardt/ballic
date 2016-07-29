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
	** the internal energy. We compare the second derivative
	** of u with respect to v from tillInitSplinev with 
	** the result of numerical differentiation. 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableMax = 1000;
	double rho, u, P, c2, d2udv2;

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
	
	for (i=0;i<granite->nTableMax;i++)
	{
		for (j=0;j<granite->nTableMax;j++)
		{
			// For the analyic expression
			rho = granite->Lookup[INDEX(i,j)].rho;
			u = granite->Lookup[INDEX(i,j)].u;
			P = tillPressureSound(granite, rho, u, &c2);

			// Calculate d2udrho2 numerically
			// d2y/dx2 = (u[j+1] - 2*u[j] + u[j-1])/(dx*dx)

			if (j > 0 && j < granite->nTableMax-1)
			{
				d2udv2 = (granite->Lookup[INDEX(i,j+1)].u - 2*granite->Lookup[INDEX(i,j)].u +granite->Lookup[INDEX(i,j-1)].u)/(granite->delta*granite->delta);
			} else {
				d2udv2 = 0;
			}
	
			printf("%i %i %g %g %g %g\n", i, j, rho, u, granite->Lookup[INDEX(i,j)].udv2, d2udv2);
		}
		printf("\n");
	}
	
	tillFinalizeMaterial(granite);
}

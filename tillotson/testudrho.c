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
	** Test to check of the first derivatives derived using
	** the cubic spline work. We calculate dudrho and compare
	** it to the analytic solution. 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableMax = 1000;
	double rho, drho, u, P, dudrho, dudrho_a;
	double A, B;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);

	/* Also solve the splines in rho */
	tillInitSplineRho(granite);

	fprintf(stderr, "Done.\n");

	fprintf(stderr,"nTableMax: %i\n", granite->nTableMax);
	
	for (j=0;j<granite->nTableMax;j++)
	{
		for (i=0;i<granite->nTableMax;i++)
		{
			// For the analyic expression
			rho = granite->Lookup[INDEX(i,j)].rho;
			u = granite->Lookup[INDEX(i,j)].u;
			P = tillPressure(granite, rho, u);

			// Analyic solution
			dudrho_a = P/(rho*rho);

			// Calculate the first derivative using splines

			/*
			** dx = x[j+1] - x[j]
			** A = (x[j+1]-x)/(x[j+1]-x[j])
			** B = (x-x[j])/(x[j+1-x[xj])
			** dy/dx = (y[j+1]-y[j])/(x[j+1]-x[j]) - (3.0*A*A-1.0)/6.0*(x[j+1]-x[j])*y2[j] + (3.0*B*B-1.0)/6.0*(x[j+1]-x[j])*y2[j+1]
			** 
			*/
			
			if (i > 0 && i < granite->nTableMax-1)
			{
				drho = granite->delta;

				// dudrho is at (i,j):   dudrho = u[j+1]-u[j])/dv-1.0/6.0*dv*(2.0*udv2[j]+udv2[j+1])
				dudrho = (granite->Lookup[INDEX(i+1, j)].u-granite->Lookup[INDEX(i, j)].u)/drho-1.0/6.0*drho*(2.0*granite->Lookup[INDEX(i, j)].udrho2+granite->Lookup[INDEX(i+1, j)].udrho2);
			} else {
				dudrho = 0;
			}

			// For small rho the analytic expression deviates a lot from the numerical solution!
			printf("%i %i %g %g %g %g\n", i, j, rho, u, dudrho, dudrho_a);
#if 0
			if (fabs(dudrho-dudrho_a)<0.01) 
			{
			} else {
				printf("%i %i %g\n", i, j, dudrho-dudrho_a);
			}
#endif
		}
		printf("\n");
	}

	tillFinalizeMaterial(granite);
}

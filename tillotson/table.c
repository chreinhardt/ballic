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

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 25.0;
	double vmax = 25.0;
	int nTableMax = 1000;

	double rho1 = 10.0;
	double u1 = 15.0;
	double rho2 = 10.0;

	int i = 0;
	int j = 0;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	if (argc != 4) {
		fprintf(stderr,"Usage: table <rho1> <u1> <rho2>\n");
		exit(1);
	}

    rho1 = atof(argv[1]);
    u1 = atof(argv[2]);
    rho2 = atof(argv[3]);

	assert(rho1 < rhomax && rho2 < rhomax);

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	/* Check if we get the same result for a very small difference. */
/*
	rho1 = granite->rho0;
	u1 = granite->delta*5;
	rho2 = granite->rho0+5*granite->delta;
	fprintf(stderr,"nTableMax: %i\n",granite->nTableMax);
*/

	fprintf(stderr,"rho1: %g u1: %g rho2: %g u2: %g (interpol) %g (integrated)\n",rho1,u1,rho2,tillLookupU(granite,rho1,u1,rho2,0),tillCalcU(granite,rho1,u1,rho2));

	/* Debug the function tillColdULookup().
	u = tillColdULookup(granite,rho);
	printf("rho: %.30f u: %.30f u: %g\n", rho, tillColdULookup(granite, rho),u);
	*/
	tillFinalizeMaterial(granite);
}

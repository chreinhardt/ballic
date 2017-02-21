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

#define INDEX(i, j) ((i*tillmat->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	// These values are arbitrary
	double rhomax = 100.0;
	double vmax = 1200.0;
	double rho, u, P;
	int nTableRho = 1000;
	int nTableV = 1000;
	int iMat = 0;

	int i = 0;
	int j = 0;

	TILLMATERIAL *tillmat;
	struct lookup *isentrope;

	if (argc != 4)
	{
			fprintf(stderr, "calcpressure <iMat> <rho> <u>\n");
			exit(1);
	}

	iMat = atoi(argv[1]);
	assert(iMat >= 0);

	rho = atof(argv[2]);
	assert(rho >= 0.0);

	u = atof(argv[3]);
	assert(u >= 0.0);

	tillmat = tillInitMaterial(iMat, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax, 1);
	
	P = tillPressure(tillmat, rho, u);
	printf("%i %.8g %.8g %.8g\n",iMat,rho,u,P);

	tillFinalizeMaterial(tillmat);
}

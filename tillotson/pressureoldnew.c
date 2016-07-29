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
	** Compare the two functions that calculate the pressure.
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 1000.0;
	double vmax = 1000.0;
	int nTableMax = 1000;

	double rho, u, Pold, Pnew;

	double umax = 1000.0;
	double delta = 0.01;

	TILLMATERIAL *granite;

/*
	if (argc != 5) {
	fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <Tcore> <gamma>  >myball.std\n");
	exit(1);
	}
    nDesired = atoi(argv[1]);
    mTot = atof(argv[2]);
    Tcore = atof(argv[3]);
    gamma = atof(argv[4]);
*/
    //gamma = 2.0; /* Choose a gamma between 0.2 and 0.6 */

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableMax, rhomax, vmax);

	rho = 0.0;
	u = 0.0;
	Pold = 0.0;
	Pnew = 0.0;
	
	while (u < umax)
	{
		rho = 0.0;
		while (rho < rhomax)
		{
			Pold = tillPressureSoundold(granite, rho, u, NULL);
			Pnew = tillPressureSound(granite, rho, u, NULL);

			//printf("%g %g %.20f %.20f %g\n", rho, u, Pold, Pnew, Pold-Pnew);

			//if (fabs(Pold-Pnew) > 10e-12) fprintf(stderr,"Pold not equal to Pnew: %g %g %f %f %g\n",rho,u,Pold,Pnew,Pold-Pnew);	
			if (fabs(Pold-Pnew) > 10e-12) printf("Pold-Pnew > 10e-12: %g %g %f %f %g\n",rho,u,Pold,Pnew,Pold-Pnew);	
	
			rho += delta;
		}
		u += delta;
	}

	tillFinalizeMaterial(granite);
}

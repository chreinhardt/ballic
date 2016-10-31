/*
 ** Copyright (c) 2014-2016 Christian Reinhardt and Joachim Stadel.
 **
 ** This file provides all the functions to do the density correction
 ** at a material interface proposed in Woolfson 2007
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"
#include "tillwoolfson.h"

/* Basic functions:
 *
 * tillInitInterfLookup: initialise the lookup table with the material coefficients f_ij.
 * tillFinalizeInterfLookup: free memory.
 */

//TILLINTERFACE *tillInitInterface(TILLMATERIAL **pMat, int nMat, int nTableP, int nTableT, double Pmin, double Pmax, double Tmin, double Tmax)
TILLINTERFACE *tillInitInterface(TILLMATERIAL **pMat, int nMat, int nTableP, int nTableT, double Pmin, double Pmax, double Tmin, double Tmax)
{
	/*
	 * Initialise the material interface data structure
	 *
	 * We do:
	 * Initialize variables and allocate memory
	 */
    TILLINTERFACE *interface;
	int i;
	 
    interface = malloc(sizeof(TILLINTERFACE));
    assert(interface != NULL);

	interface->nMat = nMat;
	
	/* Number of grid points for the look up table */
	interface->nTableP = nTableP;
	interface->nTableT = nTableT;
	
	/* Min and max values for P and T. */
	interface->Pmin = Pmin;
	interface->Pmax = Pmax;

	interface->Tmin = Tmin;
	interface->Tmax = Tmax;

	/* The memory for the lookup table is allocated when tillInitLookup is called. */
    return(interface);
}

void tillFinalizeInterface(TILLINTERFACE *interface)
{
	int i;
	/* Free the memory */
	for (i=0; i < interface->nMat; i++)
	{
		if (interface->InterfLookup[i] != NULL) free(interface->InterfLookup[i]);
	}

	free(interface);
}

void tillInitInterfLookup(TILLINTERFACE *interface)
{
	/*
	** Make a look up table for the correction factors f_ij for each combination of materials.
	*/
	double P, T;
	int i;

	/* Make sure that the materials are properly initialized. */
	for (i=0; i<interface->nMat; i++)
	{
		assert(interface->tillMat[i] != NULL);
	}

	/* Make a table for each material i and j. */
	for (i=0; i<interface->nMat; i++)
	{

		assert(interface->tillMat[i] != NULL);
	}
}

void tillSolveInterfLookup(TILLINTERFACE *interface, int Mati, int Matj)
{
	/*
	** Generate a lookup table for the correction factor f_ij for material i and j.
	*/
	double P, T, dP, dT;
	double rhoi,rhoj;

	/* Make sure that the materials are properly initialized. */
	assert(interface->tillMat[Mati] != NULL);
	assert(interface->tillMat[Matj] != NULL);

	dP = (interface->Pmax-interface->Pmin)/interface->nTableP;
	dT = (interface->Tmax-interface->Tmin)/interface->nTableT;

	P = interface->Pmin;
	T = interface->Tmin;
	
	while (T<interface->Tmax)
	{
		/* Reset P to Pmin. */
		P=interface->Pmin;
		while (P<interface->Pmax)
		{
			rhoi = tillRhoPTemp(interface->tillMat[Mati], P, T);
			rhoj = tillRhoPTemp(interface->tillMat[Matj], P, T);

		}
	}
#if 0
	Pc = 0.0;

	/* Calculate P and T in material 1. */
	P = tillPressure(mat1, rho1, u1);
	T = tillTempRhoU(mat1, rho1, u1);

	/*
	** We use rho1 as an upper limit for rho2 assuming that the denser component is in the inner shell.
	*/
	a = rho1;
	ua = tillURhoTemp(mat2, a, T);
	Pa = tillPressure(mat2, a, ua);

	b = 0.0;
	ub = tillURhoTemp(mat2, b, T);
	Pb = tillPressure(mat2, b, ub);
	
	assert (Pa > P && Pb < P);	
	fprintf(stderr,"modelSolveBC: starting with a=%g ua=%g Pa=%g b=%g ub=%g Pb=%g\n",a,ua,Pa,b,ub,Pb);

    /*
    ** Root bracketed by (a,b).
    */
    while (Pa-Pb > 1e-10) {
		c = 0.5*(a + b);
		uc = tillURhoTemp(mat2,c,T);
		Pc = tillPressure(mat2,c, uc);
		
		if (Pc < P) {
			b = c;
			Pb = Pc;
		}
		else {
			a = c;
			Pa = Pc;
		}
//		fprintf(stderr,"c:%.10g Pc:%.10g\n",c,Pc);
	}

//	fprintf(stderr,"modelSolveBC: rho1: %g, u1: %g, rho2:%g, u2:%g\n",rho1,u1,c,uc);
//	fprintf(stderr,"modelSolveBC: P1: %g, T1: %g, P2:%g, T2:%g\n",P,T,tillPressure(mat2,c,uc),tillTempRhoU(mat2,c,uc));
	/*
	** Return values.
	*/
	*prho2 = c;
	*pu2 = uc; 
#endif
}


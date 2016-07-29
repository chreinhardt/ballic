/*
 ** Copyright (c) 2014-2015 Christian Reinhardt and Joachim Stadel.
 **
 ** This file provides all the functions for the Tillotson EOS library.
 ** The Tillotson EOS (e.g. Benz 1986) is a relatively simple but reliable
 ** and convenient to use equation of state that can describe matter over
 ** a large range of pressures, densities and internal energies.
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "tillotson.h"

#include "interpol/coeff.h"
#include "interpol/interpol.h"

/* This will cut the pressure in the cold expanded states for rho/rho0 < 0.8 as suggested in Melosh1989. */
//#define TILL_PRESS_MELOSH

/* Basic functions:
 *
 * tillInit: initialise the library
 *
 * tillInitLookup: make the look up table for a) the cold curve and b) the isentropes (allocate memory!)
 *
 * tillPress: calculate the pressure for a given rho and u from the Tillotson EOS
 *  (also needed for the look up table)
 */

TILLMATERIAL *tillInitMaterial(int iMaterial, double dKpcUnit, double dMsolUnit, int nTableRho, int nTableV, double rhomax, double vmax, int iExpV)
{
	/*
	 * Initialise a material from the Tillotson library
	 *
	 * We do:
	 * Initialize variables
	 * Convert quantities to code units
	 * The memory for the look up table is allocated in tillInitLookup()
	 */

    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */

    TILLMATERIAL *material;
	int i;
	 
    material = malloc(sizeof(TILLMATERIAL));
    assert(material != NULL);

	material->iMaterial = iMaterial;
/*
    material->dKpcUnit = 2.06701e-13;
    material->dMsolUnit = 4.80438e-08;
*/
	/* This two parameters define the unit system we use */
    material->dKpcUnit = dKpcUnit;
    material->dMsolUnit = dMsolUnit;
	material->rhomax = rhomax;
	material->vmax = vmax;

	if (material->vmax == 0)
	{
		/* Just as a first step we have equal steps in rho and v. */
		material->vmax = material->rhomax;
	}

	/* Needs about 800M memory. */
//	material->nTableMax = 10000;
	/* Number of grid points for the look up table */
	material->nTableRho = nTableRho;
	material->nTableV = nTableV;
    /*
    ** Convert kboltz/mhydrogen to system units, assuming that
    ** G == 1.
    */
    material->dGasConst = material->dKpcUnit*KPCCM*KBOLTZ
	/MHYDR/GCGS/material->dMsolUnit/MSOLG;
    /* code energy per unit mass --> erg per g */
    material->dErgPerGmUnit = GCGS*material->dMsolUnit*MSOLG/(material->dKpcUnit*KPCCM);
    /* code density --> g per cc */
    material->dGmPerCcUnit = (material->dMsolUnit*MSOLG)/pow(material->dKpcUnit*KPCCM,3.0);
    /* code time --> seconds */
    material->dSecUnit = sqrt(1/(material->dGmPerCcUnit*GCGS));


	/* The memory for the lookup table is allocated when tillInitLookup is called. */
		
	/*
	** Set the Tillotson parameters for the material.
	*/
	switch(iMaterial)
	{
		case GRANITE:
			/*
			** Just set granite at the moment.
			*/
			material->a = 0.5;
			material->b = 1.3;
			material->u0 = 1.6e11; /* in ergs/g */
			material->rho0 = 2.7; /* g/cc */
			material->A = 1.8e11; /* ergs/cc */
			material->B = 1.8e11; /* ergs/cc */
			material->us = 3.5e10; /* ergs/g */
			material->us2 = 1.8e11; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.79e7; /* ergs/g K (or 790 J/kg K) */ 
			break;
		case IRON:
			/*
			** Material parameters from Benz 1987.
			*/
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 9.5e10; /* in ergs/g */
			material->rho0 = 7.86; /* g/cc */
			material->A = 1.28e12; /* ergs/cc */
			material->B = 1.05e12; /* ergs/cc */
			material->us = 1.425e10; /* ergs/g */
			material->us2 = 8.45e10; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.449e7; /* ergs/g K */ 
			break;
		default:
			/* Unknown material */
			assert(0);
	}

    /*
    ** Convert energies and densities to code units!
    */
    material->u0 /= material->dErgPerGmUnit;
    material->us /= material->dErgPerGmUnit;
    material->us2 /= material->dErgPerGmUnit;
    material->rho0 /= material->dGmPerCcUnit;
    material->A /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
    material->B /= (material->dGmPerCcUnit*material->dErgPerGmUnit);
 
 	material->cv /= material->dErgPerGmUnit;

	/* Set drho so that rho0 lies on the grid. */
	material->n = floor(material->rho0/material->rhomax*material->nTableRho);
 	material->drho =  material->rho0/material->n;
	
	/* Set the actual rhomax. */ 
  	material->rhomax = material->drho*(material->nTableRho-1);
	material->dv = material->vmax/(material->nTableV-1);
 	// (CR) set vmax to rhomax just to check if dv = drho = delta works. */
	// material->vmax = material->rhomax;

	// (CR) 15.11.15: try non uniform steps in v
	// But careful: material->n is used to integrate the isentropes
	//	material->n = 5;
	material->iExpV = iExpV;
	// (CR) 15.11.15: Change this back later
    return(material);
}

void tillFinalizeMaterial(TILLMATERIAL *material)
{
	/* Free the memory */
	if (material->Lookup != NULL) free(material->Lookup);
	if (material->cold != NULL) free(material->cold);

	free(material);
}

double tilldPdrho_s(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dP/drho at constant entropy.
	*/

	return (1.0/(rho*rho)*(tillSoundSpeed2old(material,rho, u)-2.0*tillPressure(material,rho,u)/rho));
}

double tillSoundSpeed2old(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate the sound speed for the Tillotson EOS like we did in
	** Gasoline. In the intermediate states its better however to do a
	** linear interpolation in the sound speed because we have no
	** discontinuity when we change from expanded hot to condensed states.
	*/

	double eta, mu;
	double Pc, Pe;
	double c2c, c2e;
	double Gammac, Gammae, w0, y, z;

	eta = rho/material->rho0;
	mu = eta - 1.0;
	z = (1.0 - eta)/eta;
	w0 = u/(material->u0*eta*eta)+1.0;
	
	/*
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;

		return ((Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho));
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		
		return ((Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho));
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
		Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

		return ((Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho)));
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
		Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));
	
		return ((Gammac*(1.0-y)+Gammae*y+1.0)*(Pc*(1.0-y)+Pe*y)/rho+((material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho))*(1.0-y)+(material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho)))*y);
	}
}

double tillPressureSoundold(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
	/*
	** Calculate the pressure and sound speed from the Tillotson EOS for a material.
	** Set pcSound = NULL to only calculate the pressure forces. Here we used the old
	** function that was originally implemented in Gasoline.
	*/

	double eta, mu;
	double P,c2;
	double Gamma, epsilon0, y, x;

	eta = rho/material->rho0;
	mu = eta - 1.0;
	x  = 1.0 - eta;
	
	/*
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		y = 0.0;
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		y = 0.0;
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		y = 1.0;
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);
	}

	Gamma = material->a + material->b/(u/(material->u0*eta*eta)+1.0)*((1.0-y) + y*exp(-material->beta*x*x/eta/eta));
	epsilon0 = 1.0/(Gamma*rho)*((material->A*mu + material->B*mu*mu)*(1.0-y)+y*material->A*mu*exp(-(material->alpha*x/eta + material->beta*x*x/eta/eta)));
	P = Gamma*rho*(u + epsilon0);

	if (pcSound != NULL)
	{
		/* calculate the sound speed */
		c2 = (Gamma+1.0)*P/rho + 1.0/rho*((material->A+material->B*(eta*eta-1.0))*(1.0-y)+y*material->A*exp(-(material->alpha*x/eta+material->beta*x*x/eta/eta))*(1+mu/eta*(material->alpha+2*material->beta*x/eta)));
		/* make sure that c^2 > 0 for rho < rho0 */
		if (c2 < material->A/material->rho0){
			c2 = material->A/material->rho0;
		}
		
		*pcSound = sqrt(c2);
	}

	return (P);
}

double tillPressureSound(TILLMATERIAL *material, double rho, double u, double *pcSound)
{
	/*
	** Calculate the pressure and sound speed from the Tillotson EOS for a material.
	** Set pcSound = NULL to only calculate the pressure forces.
	*/

	double eta, mu;
	double Pc, Pe;
	double c2c, c2e;
	double Gammac, Gammae, w0, y, z;

	eta = rho/material->rho0;
	mu = eta - 1.0;
	z = (1.0 - eta)/eta;
	w0 = u/(material->u0*eta*eta)+1.0;
	
	/*
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		
		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			//*pcSound = sqrt(c2c);
			*pcSound = c2c;
		}
		return (Pc);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		
		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			//*pcSound = sqrt(c2c);
			*pcSound = c2c;
		}
#ifdef TILL_PRESS_MELOSH
		/* Melosh 1989 suggests to cut the pressure in the expanded cold states for rho/rho < 0.8 */
		if (eta < 0.8)
		{
			fprintf(stderr,"Setting pressure to zero for eta=%g\n",eta);
			Pc = 0.0;
		}
#endif
		return (Pc);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
		Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));
			//*pcSound = sqrt(c2e);
			*pcSound = c2e;
		}
		
		return (Pe);
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		Gammac = material->a + material->b/w0;
		Pc = Gammac*u*rho + material->A*mu + material->B*mu*mu;
		Gammae = material->a + material->b/w0*exp(-material->beta*z*z);
		Pe = Gammae*u*rho + material->A*mu*exp(-(material->alpha*z+material->beta*z*z));

		if (pcSound != NULL)
		{
			/* calculate the sound speed */
			c2c = (Gammac+1.0)*Pc/rho + (material->A+material->B*(eta*eta-1.0))/rho + material->b/(w0*w0)*(w0-1.0)*(2*u-Pc/rho);
			c2e = (Gammae+1.0)*Pe/rho + material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta)) + material->b*rho*u/(w0*w0*eta*eta)*exp(-material->beta*z*z)*(2.0*material->beta*z*w0/material->rho0+1.0/(material->u0*rho)*(2.0*u-Pe/rho));

			//*pcSound = sqrt(c2c*(1.0-y)+c2e*y);
			*pcSound = c2c*(1.0-y)+c2e*y;
		}
	
#ifdef TILL_PRESS_MELOSH
		/* Melosh 1989 suggests to cut the pressure in the expanded cold states for rho/rho < 0.8 */
		if (eta < 0.8)
		{
			fprintf(stderr,"Setting pressure to zero for eta=%g\n",eta);
			Pc = 0.0;
		}
#endif
		return (Pc*(1.0-y)+Pe*y);
	}
}

double tillPressure(TILLMATERIAL *material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */
	double P = tillPressureSound(material, rho, u, NULL);
	if (P < 0.0 ) P = 0.0;
	return (P);
}

double tilldPdrho(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dP/drho at u=const.
	*/

	double eta, mu;
	double dPcdrho, dPedrho;
	double w0, y, z;

	eta = rho/material->rho0;
	mu = eta - 1.0;
	z = (1.0 - eta)/eta;
	w0 = u/(material->u0*eta*eta)+1.0;
	
	/*
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		dPcdrho = (material->a+material->b/w0*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;
		
		return (dPcdrho);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		dPcdrho = (material->a+material->b/w0*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;

		return (dPcdrho);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		dPedrho = (material->a + material->b/w0*exp(-material->beta*z*z)*(2.0*material->beta*z/eta+3.0-2.0/w0))*u+material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta));		
		return (dPedrho);
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		dPcdrho = (material->a+material->b/(w0*w0)*(3.0-2.0/w0))*u + (material->A+2.0*material->B*mu)/material->rho0;
		dPedrho = (material->a + material->b/w0*exp(-material->beta*z*z)*(2.0*material->beta*z/eta+3.0-2.0/w0))*u+material->A/material->rho0*exp(-(material->alpha*z+material->beta*z*z))*(1.0+mu/(eta*eta)*(material->alpha+2.0*material->beta*z-eta));		
		return (dPcdrho*(1.0-y)+dPedrho*y);
	}
}

double tilldPdu(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dP/du at rho=const.
	*/

	double eta, mu;
	double dPcdu, dPedu;
	double w0, y, z;

	eta = rho/material->rho0;
	mu = eta - 1.0;
	z = (1.0 - eta)/eta;
	w0 = u/(material->u0*eta*eta)+1.0;
	
	/*
	**  Here we evaluate, which part of the equation of state we need.
	*/
	if (rho >= material->rho0) {
		/*
		**  condensed states (rho > rho0)
		*/
		dPcdu = (material->a+material->b/(w0*w0))*rho;
		
		return (dPcdu);
	} else if (u < material->us) {
		/* 
		** cold expanded states (rho < rho0 and u < us)
		** P is like for the condensed states
		*/
		dPcdu = (material->a+material->b/(w0*w0))*rho;
		return (dPcdu);
	} else if (u > material->us2) {
		/*
		** expanded hot states (rho < rho0 and u > us2)
		*/
		dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		
		return (dPedu);
	} else {
		/*
		**  intermediate states (rho < rho0 and us < u < us2)
		*/
		y = (u - material->us)/(material->us2 - material->us);

		dPcdu = (material->a+material->b/(w0*w0))*rho;
		dPedu = (material->a + material->b/(w0*w0)*exp(-material->beta*z*z))*rho;		

		return (dPcdu*(1.0-y)+dPedu*y);
	}
}

double tilldTdrho(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dT/drho at u=const.
	*/
	assert(material->cv > 0.0);
	return (-1.0/material->cv*tillPressure(material,rho,tillColdULookup(material,rho))*(rho*rho));
}

double tilldTdu(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate dT/du at rho=const.
	*/

	assert(material->cv > 0.0);
	return (-1.0/material->cv);
}

double tillPressureRhoU(TILLMATERIAL material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */

}

double tillTempRhoU(TILLMATERIAL *material, double rho, double u)
{
	/*
	** Calculate T(rho,u) for a material. As an approximation
	** we use u(rho,T) = uc(rho) + cv*T.
	*/
	assert(material->cv > 0.0);
	return ((u-tillColdULookup(material,rho))/material->cv);
}

double tillTempRhoP(TILLMATERIAL *material, double rho, double P)
{
	/* Calculate T(rho,P) for a material */
}

double tillSoundSpeed(TILLMATERIAL *material, double rho, double u)
{
	/* Calculate sound speed for a material */
	double c;

	tillPressureSound(material, rho, u, &c);
	return (c);
}

double tillDensRatio(TILLMATERIAL material1, TILLMATERIAL material2, double P, double T)
{
	/* From Woolfson 2007 */
}

double tilldudrho(TILLMATERIAL *material, double rho, double u)
{
	return(tillPressure(material,rho,u)/(rho*rho));
}

int comparerho(const void* a, const void* b)
/*
** This function compares two entries in the look up table
** and returns -1 if a1.rho < a2.rho, 1 if a1.rho > a2.rho or 0 if
** they are equal (needed to sort the particles with qsort).
*/
{
	TILL_LOOKUP_ENTRY a1 = *(const TILL_LOOKUP_ENTRY*)(a);
    TILL_LOOKUP_ENTRY a2 = *(const TILL_LOOKUP_ENTRY*)(b);
    
    if (a1.rho < a2.rho) return -1;
    if (a1.rho > a2.rho) return 1;
    return 0;
}

void tillInitColdCurve(TILLMATERIAL *material)
{
	/* Generate the look up table for the cold curve */
	TILL_LOOKUP_ENTRY *isentrope;

	/* Lookup table */
    material->cold = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));
    assert(material->cold != NULL);

	/* v = u(rho=rho0) */
	isentrope = tillSolveIsentrope(material,0.0);

	material->cold = isentrope;
	/* Now sort the look up table. */
	// qsort(material->cold,material->nTable,sizeof(TILL_LOOKUP_ENTRY),comparerho);
}

void tillInitLookup(TILLMATERIAL *material)
{
	/*
	** Generate the look up table for the isentropic evolution of
	** the internal energy.
	*/
#define TILL_USE_RK4
	TILL_LOOKUP_ENTRY *isentrope;
	double v, dv;
    int i,j;

	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
//    material->Lookup = malloc(material->nTableMax*material->nTableMax*sizeof(double));
    material->Lookup = malloc(material->nTableRho*material->nTableV*sizeof(TILL_LOOKUP_ENTRY));
    assert(material->Lookup != NULL);

	v = 0.0;
	//fprintf(stderr, "Starting integration...\n");
	// (CR) This doesnt work, we need dv=vmax/(nTableMax-1)
	//dv = material->vmax/material->nTableMax;
	//dv = material->vmax/(material->nTableMax-1);
	dv = material->dv;

	/*
	** Integrate the isentropes for different v.
	*/
	for (j=0; j<material->nTableV; j++)
	{
//		printf("%15.7E%15.7E%15.7E%15.7E%15.7E%15.7E\n",v,j*material->delta,v-material->vmax, dv, material->delta, dv-material->delta);
//		printf("%15.7E%15.7E%15.7E%15.7E%15.7E%15.7E\n",v,j*material->dv,v-material->vmax, dv, material->dv, dv-material->dv);

//		(CR) 15.11.15: Try non uniform spacing in v
//		v = material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
//		(CR) 15.11.15: Until here	
		isentrope = tillSolveIsentrope(material,v);
		
		/* Copy one row to the look up table. This is of course not efficient at all. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[INDEX(i,j)] = isentrope[i];		
		}
		
		free(isentrope);
		v += dv;
	}

    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */
	tillInitSplines(material);

	/* Initialize the coefficients for the interpolation. */
//	SamplesToCoefficients(material->Lookup,material->nTableMax,material->nTableMax, TILL_SPLINE_DEGREE);
}

TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v)
{
	/*
	** Integrate one isentrope for the look up table. The parameter v corresponds to
	** u(rho=rho0) so v=0 gives the cold curve.
	*/
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	TILL_LOOKUP_ENTRY *isentrope;
    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));

	rho = material->rho0;
	u = v;
	h = material->drho;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].u = u;
	isentrope[i].u1 = tilldudrho(material, rho, u); // du/drho
	
	/* Output some information. */
#ifdef TILL_USE_RK4
	fprintf(stderr,"Using RK4.\n");
#else
	fprintf(stderr,"Using RK2.\n");
#endif

	/*
	** Integrate the condensed and expanded states separately.
	*/
#ifdef TILL_USE_RK4
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*tilldudrho(material,rho,u);
			k2u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}
#else
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = hs*tilldudrho(material,rho,u);
			k2u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k1u);

			u += k2u;
			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}
#endif
	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = v;

#ifdef TILL_USE_RK4
	for (i=material->n-1;i>=0;i--)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*-tilldudrho(material,rho,u);
			k2u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*-tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
#if 0
		if (i == 0)
		{
//			fprintf(stderr,"i=%i,rho(0)=%g,u(0)=%g,u1(0)=%g",i,isentrope[i].rho,isentrope[i].u,isentrope[i].u1);
//			fprintf(stderr,"  tilldudrho=%g P=%g\n",tilldudrho(material, rho, u),tillPressure(material, rho, u));
//			isentrope[i].u1 = 0.0;
//			isentrope[i].u1 = isentrope[i+1].u1; // Set u1(0)=u1(drho)
		}
#endif
	}
#else
	for (i=material->n-1;i>=0;i--)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = hs*-tilldudrho(material,rho,u);
			k2u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k1u);

			u += k2u;
			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}
#endif
	return isentrope;
}

TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v)
{
	/*
	** Integrate one isentrope for the look up table. The parameter v corresponds to
	** u(rho=rho0) so v=0 gives the cold curve. We use the BS algorithm from
	** the Numerical Recipes.
	*/
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	TILL_LOOKUP_ENTRY *isentrope;
    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));

	rho = material->rho0;
	u = v;
	h = material->drho;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].u = u;
	isentrope[i].u1 = tilldudrho(material, rho, u); // du/drho
	
	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*tilldudrho(material,rho,u);
			k2u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho += hs;
		}

	    isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}

	/*
	** Now the expanded states. Careful about the negative sign.
	*/
	rho = material->rho0;
	u = v;
	fprintf(stderr,"Expanded states\n");

	for (i=material->n-1;i>=0;i--)
	{
		float hs = h/100.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<100;s++)
		{
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = hs*-tilldudrho(material,rho,u);
			k2u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k1u);
			k3u = hs*-tilldudrho(material,rho+0.5*hs,u+0.5*k2u);
			k4u = hs*-tilldudrho(material,rho+hs,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
			rho -= hs;
		}
	
		isentrope[i].u = u;
	    isentrope[i].rho = rho;
		isentrope[i].u1 = tilldudrho(material, rho, u);
	}
	return isentrope;
}
void tillInitSplines(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u and u1 in v.
	** For this we use routines from Numerical Recipes.
	*/
//	tillInitSplineRho(material);
	tillInitSplineU1(material);
	tillInitSplineU(material);
}

void tillInitSplineRho(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u and u1 in rho.
	** For this we use routines from Numerical Recipes.
	*/
	float *x;
	float *y;
	int i,j,n;
	float yp1, ypn;
	float *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableRho;
	x = malloc(material->nTableRho*sizeof(float));
	y = malloc(material->nTableRho*sizeof(float));
	y2 = malloc(material->nTableRho*sizeof(float));

	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	/*
	** Just to debug the code we calculate d2u/drho2 and store it in udrho2
	** because for this we have an analytic expression.
	*/
	for (j=0; j<material->nTableV; j++)
	{
/*
		// Try this: udrho2 = 1/rho^2*(c2 - 2*P/rho)
		double c2, u, rho, P;
		rho = material->Lookup[INDEX(0,j)].rho;
		u = material->Lookup[INDEX(0,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		yp1 = 1.0/(rho*rho)*(c2 - 2.0*P/rho);


		rho = material->Lookup[INDEX(material->nTableMax-1,j)].rho;
		u = material->Lookup[INDEX(material->nTableMax-1,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		ypn = 1.0/(rho*rho)*(c2 - 2.0*P/rho);
*/
		/* Copy one row of the look up table to the temporary array. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			x[i] = 	material->Lookup[INDEX(i,j)].rho;
			y[i] = 	material->Lookup[INDEX(i,j)].u;
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[INDEX(i,j)].udrho2 = y2[i];
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}

void tillInitSplinev(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u in v.
	** For this we use routines from Numerical Recipes.
	*/
	float *x;
	float *y;
	int i,j,n;
	float yp1, ypn;
	float *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableV;
	x = malloc(material->nTableV*sizeof(float));
	y = malloc(material->nTableV*sizeof(float));
	y2 = malloc(material->nTableV*sizeof(float));

	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/* Copy one row of the look up table to the temporary array. */
		for (j=0; j<material->nTableV; j++)
		{
			/*
			** Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			** v = u(rho0)
			*/
			x[j] = 	j*material->dv;
			// (CR) 15.11.15: Try non uniform steps in v
			// x[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
			// (CR) 15.11.15: Done
			y[j] = 	material->Lookup[INDEX(i,j)].u;
			//printf("j: %i % g %g\n",j, x[j], y[j]);
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (j=0; j<material->nTableV; j++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[INDEX(i,j)].udv2 = y2[j];
			//printf("i: %i j: %i %g\n",i,j, y2[j]);
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}

void tillInitSplineU(TILLMATERIAL *material)
{
	tillInitSplinev(material);
}

void tillInitSplineU1(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u1 in v.
	** For this we use routines from Numerical Recipes.
	*/
	float *x;
	float *y;
	int i,j,n;
	float yp1, ypn;
	float *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableV;
	x = malloc(material->nTableV*sizeof(float));
	y = malloc(material->nTableV*sizeof(float));
	y2 = malloc(material->nTableV*sizeof(float));

	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/* Copy one row of the look up table to the temporary array. */
		for (j=0; j<material->nTableV; j++)
		{
			/*
			** Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			** v = u(rho0)
			*/
			x[j] = 	j*material->dv;
			// (CR) 15.11.15: Try non uniform steps in v
			// x[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
			// (CR) 15.11.15: Done
			y[j] = 	material->Lookup[INDEX(i,j)].u1;
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (j=0; j<material->nTableV; j++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[INDEX(i,j)].u1dv2 = y2[j];
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}

/*
** Just to debug the code we do an interpolation in u of rho.
*/
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv)
{
	/*
	** Do a cubic spline interpolation.
	*/
	float *xa;
	float *ya;
	float *y2a;
	float x,y;
	int i,n;
	
	/* Temporarily store data in an array to use splint() */
	n = material->nTableRho;
	xa = malloc(material->nTableRho*sizeof(float));
	ya = malloc(material->nTableRho*sizeof(float));
	y2a = malloc(material->nTableRho*sizeof(float));

	//void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)

	x = rho;

	/* Copy one row of the look up table to the temporary array. */
	for (i=0; i<material->nTableRho; i++)
	{
		/* Careful with the indices!
		** Lookup[i][j] = Lookup(rho,v)
		*/
		xa[i] =	material->Lookup[INDEX(i,iv)].rho;
		ya[i] =	material->Lookup[INDEX(i,iv)].u;
		y2a[i]=	material->Lookup[INDEX(i,iv)].udrho2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}

/*
** dx = x[j+1] - x[j]
** A = (x[j+1]-x)/(x[j+1]-x[j])
** B = (x-x[j])/(x[j+1-x[xj])
** dy/dx = (y[j+1]-y[j])/(x[j+1]-x[j]) - (3.0*A*A-1.0)/6.0*(x[j+1]-x[j])*y2[j] + (3.0*B*B-1.0)/6.0*(x[j+1]-x[j])*y2[j+1]
** 
*/
/*
** Just to debug the code we do an interpolation in u of v.
*/
double tillSplineIntv(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation.
	*/
	float *xa;
	float *ya;
	float *y2a;
	float x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(float));
	ya = malloc(material->nTableV*sizeof(float));
	y2a = malloc(material->nTableV*sizeof(float));

	//void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)

	x = v;

	/* Copy one row of the look up table to the temporary array. */
	for (j=0; j<material->nTableV; j++)
	{
		/* Careful with the indices!
		** Lookup[i][j] = Lookup(rho,v)
		*/
		xa[j] =	j*material->dv;
		// (CR) 15.11.15: Try non uniform steps in v
		// xa[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
		// (CR) 15.11.15: Done
		ya[j] =	material->Lookup[INDEX(irho,j)].u;
		y2a[j]=	material->Lookup[INDEX(irho,j)].udv2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}

double tillSplineIntU(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation of u in v.
	*/
	float *xa;
	float *ya;
	float *y2a;
	float x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(float));
	ya = malloc(material->nTableV*sizeof(float));
	y2a = malloc(material->nTableV*sizeof(float));

	//void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)

	x = v;

	/* Copy one row of the look up table to the temporary array. */
	for (j=0; j<material->nTableV; j++)
	{
		/* Careful with the indices!
		** Lookup[i][j] = Lookup(rho,v)
		*/
		xa[j] =	j*material->dv;
		// (CR) 15.11.15: Try non uniform steps in v
		// xa[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
		// (CR) 15.11.15: Done
		ya[j] =	material->Lookup[INDEX(irho,j)].u;
		y2a[j]=	material->Lookup[INDEX(irho,j)].udv2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}

double tillSplineIntU1(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation of u1 in v.
	*/
	float *xa;
	float *ya;
	float *y2a;
	float x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(float));
	ya = malloc(material->nTableV*sizeof(float));
	y2a = malloc(material->nTableV*sizeof(float));

	//void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)

	x = v;

	/* Copy one row of the look up table to the temporary array. */
	for (j=0; j<material->nTableV; j++)
	{
		/* Careful with the indices!
		** Lookup[i][j] = Lookup(rho,v)
		*/
		xa[j] =	j*material->dv;
		// (CR) 15.11.15: Try non uniform steps in v
		//xa[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
		// (CR) 15.11.15: Done
		ya[j] =	material->Lookup[INDEX(irho,j)].u1;
		y2a[j]=	material->Lookup[INDEX(irho,j)].u1dv2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}

double tillCubicInt(TILLMATERIAL *material, double rhoint, double vint) {
	/*
	** Use cubicint to interpolate u for a given rho and v.
	*/
	double dv, A, B;
	int i, j;
	double *u, *dudrho, *dudv, *dudvdrho, *rho, *intvalues;

	/* Allocate memory */
	u = malloc(2*sizeof(double));
	dudrho = malloc(2*sizeof(double));
	dudv = malloc(2*sizeof(double));
	dudvdrho = malloc(2*sizeof(double));
	rho = malloc(2*sizeof(double));

	intvalues = malloc(4*sizeof(double));

	// rhoint is between rho[i] and rho[i+1]
	i = floor(rhoint/material->drho);
	assert(i < material->nTableRho-1);

	// vint is between v[j] and v[j+1]
	j = floor(vint/material->dv);
	// (CR) Debug info
	// (CR) 15.11.15: Try non uniform spacing in v
	// j = floor(pow(vint/material->vmax,1.0/material->n)*(material->nTableV-1));
	//fprintf(stderr,"j: %i\n",j);
	// (CR) 15.11.15: Until here
	assert(j < material->nTableV-1);

	// dv = v[j+1]-v[j]
	// Uniform steps in v
	dv = material->dv;
	// (CR) 15.11.15: Try non uniform steps in v. Needs work!!
	//dv = (material->vmax/pow(material->nTableV-1,material->n))*(pow(j+1,material->n)-pow(j,material->n));
	// (CR) 15.11.15: Done
	
	rho[0] = i*material->drho;
	rho[1] = (i+1)*material->drho;

	/* Do the spline look up in v. */
	u[0] = tillSplineIntU(material, vint, i);
	u[1] = tillSplineIntU(material, vint, i+1);

	dudrho[0] = tillSplineIntU1(material, vint, i);
	dudrho[1] = tillSplineIntU1(material, vint, i+1);

	/*
	** Calculate dudv from udv2
	** dudv = (u[j+1]-u[j])/dv - (3.0*A*A-1.0)/6.0*dv*udv2[j] + (3.0*B*B-1.0)/6.0*dv*udv2[j+1]
	**
	** Calculate dudvdrho from u1dv2
	** dudvdrho = (u1[j+1]-u1[j])/dv - (3.0*A*A-1.0)/6.0*dv*u1dv2[j] + (3.0*B*B-1.0)/6.0*dv*u1dv2[j+1]
	*/

	/* Calculate dudv from udv2 */
	A = ((j+1)*material->dv-vint)/dv;
	B = (vint-material->dv)/dv;
	
	dudv[0] = (material->Lookup[INDEX(i, j+1)].u-material->Lookup[INDEX(i, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[INDEX(i, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[INDEX(i, j+1)].udv2;
	dudv[1] = (material->Lookup[INDEX(i+1, j+1)].u-material->Lookup[INDEX(i+1, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[INDEX(i+1, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[INDEX(i+1, j+1)].udv2;

	/* Calculate dudvdrho from u1dv2 */
	dudvdrho[0] = (material->Lookup[INDEX(i, j+1)].u1-material->Lookup[INDEX(i, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[INDEX(i, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[INDEX(i, j+1)].u1dv2;
	dudvdrho[1] = (material->Lookup[INDEX(i+1, j+1)].u1-material->Lookup[INDEX(i+1, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[INDEX(i+1, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[INDEX(i+1, j+1)].u1dv2;

	/*
	** Do the interpolation for u(i,v), u(i+1,v), udrho(i,v), udrho(i+1,v)
	*/
	cubicint(u, dudrho, dudv, dudvdrho, rho, rhoint, intvalues);

	return(intvalues[0]);

	/* Free memory */
	free(u);
	free(dudrho);
	free(dudv);
	free(dudvdrho);
	free(rho);
	free(intvalues);
}

/*
** u[0] = u(v,0)
** u[1] = u(v,1)
** dudrho[0] = dudrho(v,0)
** dudrho[1] = dudrho(v,1)
** dudv[0] = dudv(v,0)
** dudv[1] = dudv(v,1)
** dudvdrho[0] = dudvdrho(v,0)
** dudvdrho[1] = dudvdrho(v,1)
** Where v is interpolated between vj and vj+1 using a cubic spline.
*/
void cubicint(double u[2],double dudrho[2], double dudv[2], double dudvdrho[2], double rho[2], double rhoint, double *intvalues) {
	double dx, e, e1;
	double *ce;
	
	assert(intvalues != NULL);

	/* Allocate memory */
	ce = malloc(4*sizeof(double));
	assert(ce != NULL);

	dx = rho[1] - rho[0];
	e = (rhoint - rho[0])/dx;
	e1 = e - 1;

	// these are the 4 Hermite functions
	ce[0] = (2*e + 1)*e1*e1;
	ce[1] = e*e1*e1;
	ce[2] = e*e*(3 - 2*e);
	ce[3] = e*e*e1;

	/*
	**    = ce[0]*u(v,0) + ce[1]*dudrho(v,0)*dx + ce[2]*u(v,1) + ce[3]*dudrho(v,1);
	** the above is written as 4 independent spline lookups in the table v lies between some j and j+1
	*/
	intvalues[0] = (ce[0]*u[0] + ce[1]*dudrho[0]*dx + ce[2]*u[1] + ce[3]*dudrho[1]*dx);
	intvalues[1] = (ce[0]*dudv[0] + ce[1]*dudvdrho[0]*dx + ce[2]*dudv[1] + ce[3]*dudvdrho[1]*dx);
	
	// free memory
	free(ce);

//	return(ce[0]*u[0] + ce[1]*dudrho[0]*dx + ce[2]*u[1] + ce[3]*dudrho[1]*dx);
//	return(ce[0]*dudv[0] + ce[1]*dudvdrho[0]*dx + ce[2]*dudv[1] + ce[3]*dudvdrho[1]*dx);
}

#if 0
double u(rho,v) {
}
double tilldudv(rho,v) {
}
void uandudv(rho,v,double *u,double *dudv);
#endif

/* Prasenjits root finder */
float brent(float (*func)(TILLMATERIAL *,float,float,float),TILLMATERIAL *material,float a,float b,float rho,float u,float tol,int iOrder);

float tillFindUonIsentrope(TILLMATERIAL *material,float v,float rho)
{
	float iv,irho,u;
	/* Needed for the interpolation function. */
	//iv = (material->nTableMax-1)*v/material->vmax;
	//irho = (material->nTableMax-1)*rho/material->rhomax;

//	return InterpolatedValue(material->Lookup,material->nTableMax,material->nTableMax,iv,irho,TILL_SPLINE_DEGREE);
	return (tillCubicInt(material, rho, v));
}

float denergy(TILLMATERIAL *material,float v,float rho,float u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
float tillFindEntropyCurve(TILLMATERIAL *material,float rho,float u,int iOrder)
{
	float tol=1e-6;

	return brent(denergy,material,0.0,material->vmax-material->dv,rho,u,tol,iOrder);
//	return brent(denergy,material,0.0,material->vmax,rho,u,tol,iOrder);
}

double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder)
{
	/* Calculates u2 for a given rho1,u2,rho2. */
	double v;

	v = tillFindEntropyCurve(material,rho1,u1,iOrder);

	return tillFindUonIsentrope(material,v,rho2);
}

#ifdef TILL_USE_OLD_BCINT
float tillFindUonIsentrope(TILLMATERIAL *material,float v,float rho)
{
	float iv,irho,u;
	/* Needed for the interpolation function. */
	//iv = (material->nTableMax-1)*v/material->vmax;
	//irho = (material->nTableMax-1)*rho/material->rhomax;

//	return InterpolatedValue(material->Lookup,material->nTableMax,material->nTableMax,iv,irho,TILL_SPLINE_DEGREE);
}

float denergy(TILLMATERIAL *material,float v,float rho,float u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
float tillFindEntropyCurve(TILLMATERIAL *material,float rho,float u,int iOrder)
{
	float tol=1e-6;
	return brent(denergy,material,0,material->vmax,rho,u,tol,iOrder);
}

double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder)
{
	/* Calculates u2 for a given rho1,u2,rho2. */
	double v;

	v = tillFindEntropyCurve(material,rho1,u1,iOrder);

	return tillFindUonIsentrope(material,v,rho2);
}
#endif

double tillColdULookup(TILLMATERIAL *material,double rho)
{
    double x,xi;
	double drho;
    int i;

    i = material->nTableRho-1;
	
	/* What do we do if rho > rhomax */
	/* if (r >= material->r[i]) return(material->rho[i]*exp(-(r-material->r[i]))); */
	
	x = rho/material->drho;
	xi = floor(x);
	assert(xi >= 0.0);
	x -= xi;

	i = (int)xi;
	if (i < 0)
	{
		fprintf(stderr,"ERROR rho:%.14g x:%.14g xi:%.14g i:%d\n",rho,x,xi,i);
	}
    assert(i >= 0);

	if (i > material->nTableRho-2) fprintf(stderr,"ERROR: out of bounds rho:%.14g rhomax:%.14g i:%i nTableRho: %i\n",rho,material->rhomax,i,material->nTableRho);
	assert(i < material->nTableRho-1);

	if (i <= material->nTableRho-2)
	{
		/* linear interpolation for now. */
		return(material->cold[i].u*(1.0-x) + material->cold[i+1].u*x);
	}
	
	/* This would only be needed if we cut off the model between the last two steps.	
	if (i == material->nTable-2)
	{
		dr = material->r[i+1] - material->r[i];
		x = r/dr;
		xi = floor(x);
		x -= xi;
		return(material->rho[i]*(1.0-x) + material->rho[i+1]*x);
	*/
	/* What do we do if i >= nTable-1
	} else {
		i = material->nTable - 1;
		return(material->rho[i]*exp(-(r-material->r[i])));
	}*/
}
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2)
{
	/* Calculate u2 by solving the ODE */
    double rho;
    double u;
    double k1u,k2u;
	double h;

	rho = rho1;
	u = u1;
	/* Make smaller steps than we used for look up table. */
	h = material->drho/100.0;

	if (rho1 < rho2)
	{
		while (rho < rho2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*tilldudrho(material,rho,u);
			k2u = h*tilldudrho(material,rho+0.5*h,u+0.5*k1u);
	
			u += k2u;
			rho += h;
		}
	} else if (rho1 > rho2) {
		while (rho > rho2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*-tilldudrho(material,rho,u);
			k2u = h*-tilldudrho(material,rho+0.5*h,u+0.5*k1u);

			u += k2u;
			rho -= h;
		}
	}
	return u;
}



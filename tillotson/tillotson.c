/*
 ** Copyright (c) 2014-2016 Christian Reinhardt and Joachim Stadel.
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
//#include "tillinitlookup.h"
//#include "tillsplint.h"

/* This will cut the pressure in the cold expanded states for rho/rho0 < 0.8 as suggested in Melosh1989. */
//#define TILL_PRESS_MELOSH

/* Basic functions:
 *
 * tillInit: initialise the library (contains all materials and converts to a given unit system).
 *
 * tillPressureSound: calculate the pressure and soundspeed for a given density and internal enegry using the Tillotson EOS.
 *
 * tillPressure: calculate the pressure for a given rho and u from the Tillotson EOS (uses tillPressureSound).
 *
 * tillSoundSpeed: calculate the sound speed for a given rho and u from the Tillotson EOS (uses tillPressureSound).
 *
 * tillFinalize: free memory.
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
//			printf("Using larger reference density for granite!\n");
//			material->rho0 = 3.32; /* g/cc */
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
			material->A = 1.279e12; /* ergs/cc */
			material->B = 1.05e12; /* ergs/cc */
			material->us = 1.42e10; /* ergs/g */
			material->us2 = 8.45e10; /* ergs/g */
//			printf("Using much larger values for us and us2\n");
//			material->us = 5e10; /* ergs/g */
//			material->us2 = 8.45e11; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.449e7; /* ergs/g K */ 
			break;
		case BASALT:
			/*
			** Material parameters from Benz 1999.
			*/
			material->a = 0.5;
			material->b = 1.5;
			material->u0 = 4.87e12; /* in ergs/g */
			material->rho0 = 2.7; /* g/cc */
//			Using a much larger rho0, to compare to dunite in ANEOS.
//			material->rho0 = 3.32; /* g/cc */
			material->A = 2.67e11; /* ergs/cc */
			material->B = 2.67e11; /* ergs/cc */
			material->us = 4.72e10; /* ergs/g */
			material->us2 = 1.82e11; /* ergs/g */
			material->alpha = 5.0;
			material->beta = 5.0;
			material->cv = 0.84e7; /* ergs/g K */ 
			break;
		case ICE:
			/*
			** Material parameters from Benz 1999.
			*/
			material->a = 0.3;
			material->b = 0.1;
			material->u0 = 1.0e11; /* in ergs/g */
			material->rho0 = 0.917; /* g/cc */
//			Using a much larger rho0, to compare to dunite in ANEOS.
//			material->rho0 = 3.32; /* g/cc */
			material->A = 9.47e10; /* ergs/cc */
			material->B = 9.47e10; /* ergs/cc */
			material->us = 7.73e9; /* ergs/g */
			material->us2 = 3.04e10; /* ergs/g */
			material->alpha = 10.0;
			material->beta = 5.0;
			// Have to look up in more details?
			material->cv = 1.0e7; /* ergs/g K */ 
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

	/* Set rhomin */
	material->rhomin = TILL_RHO_MIN;
	
	/* Set drho so that rho0 lies on the grid. */
	//material->n = floor(material->rho0/material->rhomax*material->nTableRho);
 	//material->drho =  material->rho0/material->n;
	material->n = floor((material->rho0-material->rhomin)/(material->rhomax-material->rhomin)*material->nTableRho);
 	material->drho =  (material->rho0-material->rhomin)/material->n;

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
//	if (material->cold != NULL) free(material->cold);

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

double tillPressureNP(TILLMATERIAL *material, double rho, double u)
{
	/* Calculate the pressure from the Tillotson EOS for a material */
	double P = tillPressureSound(material, rho, u, NULL);

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
	assert(0);
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
	assert(0);
}

double tillURhoTemp(TILLMATERIAL *material, double rho, double T)
{
	/* Calculate u(rho,T) for a material */
	return(tillColdULookup(material,rho) + material->cv*T);
}

double tillRhoPTemp(TILLMATERIAL *material, double P, double T)
{
	/* Calculate rho(P,T) for a material */
	double a, ua, Pa, b, ub, Pb, c, uc, Pc;

	Pc = 0.0;

	/*
	** We use rhoa=0 and rhob=rhomax as limits.
	*/
	a = material->rhomin;
	ua = tillURhoTemp(material, a, T);
	Pa = tillPressure(material, a, ua);

	b = material->rhomax*0.99;
	ub = tillURhoTemp(material, b, T);
	Pb = tillPressure(material, b, ub);
	
	/* What do we do for P=0 in the expanded cold states?*/
	fprintf(stderr,"tillRhoPTemp: starting with a=%g ua=%g Pa=%g b=%g ub=%g Pb=%g\n",a,ua,Pa,b,ub,Pb);
	assert (Pa < P && P < Pb);	
	//fprintf(stderr,"tillRhoPTemp: starting with a=%g ua=%g Pa=%g b=%g ub=%g Pb=%g\n",a,ua,Pa,b,ub,Pb);

    /*
    ** Root bracketed by (a,b).
    */
    while (Pb-Pa > 1e-10) {
		c = 0.5*(a + b);
		uc = tillURhoTemp(material,c,T);
		Pc = tillPressure(material,c, uc);
		
		if (Pc < P) {
			a = c;
			Pa = Pc;
		}
		else {
			b = c;
			Pb = Pc;
		}
//		fprintf(stderr,"c:%.10g Pc:%.10g\n",c,Pc);
	}

	fprintf(stderr,"tillRhoPTemp: rhoc=%g uc=%g Pc=%g P=%g T=%g\n",c,uc,Pc,P,T);
	/*
	** Return values.
	*/
	return(c);
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

void tillSolveBC(TILLMATERIAL *mat1, TILLMATERIAL *mat2, double rho1, double u1, double *prho2, double *pu2)
{
	/*
	** Given rho1, u1 (material 1) solve for rho2, u2 (material 2) at the interface between two material.
	** The b.c. are P1(rho1,u1)=P2(rho2,u2) and T1(rho1,u1)=T2(rho2,u2). We solve for P2(rho2)-P1=0.
	*/
	double P, T;
	double a, ua, Pa, b, ub, Pb, c, uc, Pc;

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

//	float brent(float ax, float bx, float cx, float (*f)(float), float tol,
//	float *xmin)



}


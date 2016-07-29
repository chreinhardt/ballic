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
//#include "tillinitlookup.h"
//#include "tillsplint.h"

/* Use Runge-Kutta 4th order to solve the ODEs */
#define TILL_USE_RK4

/* Basic functions:
 *
 * tillInitLookup: generate the look up table for the isentropes (allocate memory!)
 *
 * tillSolveIsentrope: solve the ODE du/drho=P/(rho*rho) for a given initial value v=u(rho0).
 */

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
	TILL_LOOKUP_ENTRY *isentrope;
	double v, dv;
    int i,j;

	/* We arrange the look up table as a 1D array with Lookup[i][j] = Lookup[i*Ntable+j] */
//    material->Lookup = malloc(material->nTableMax*material->nTableMax*sizeof(double));
    material->Lookup = malloc(material->nTableRho*material->nTableV*sizeof(TILL_LOOKUP_ENTRY));
    assert(material->Lookup != NULL);

	v = 0.0;
	// (CR) There was a bug before. We need dv=vmax/(nTableMax-1).
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
			material->Lookup[TILL_INDEX(i,j)] = isentrope[i];		
		}
		
		free(isentrope);
		v += dv;
	}

	fprintf(stderr,"Init splines\n",j);	
    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */
	tillInitSplines(material);
	fprintf(stderr,"Splines initialised\n",j);	
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
//	fprintf(stderr,"Using RK4.\n");
#else
	fprintf(stderr,"Using RK2.\n");
#endif
	/*
	** Integrate the condensed and expanded states separately.
	*/
#ifdef TILL_USE_RK4
	for (i=material->n+1;i<material->nTableRho;i++)
	{
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
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
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
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
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
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
//#if 0
		if (i == 0)
		{
//			fprintf(stderr,"i=%i,rho(0)=%g,u(0)=%g,u1(0)=%g",i,isentrope[i].rho,isentrope[i].u,isentrope[i].u1);
//			fprintf(stderr,"  tilldudrho=%g P=%g\n",tilldudrho(material, rho, u),tillPressure(material, rho, u));
//			isentrope[i].u1 = 0.0;
			isentrope[i].u1 = isentrope[i+1].u1; // Set u1(0)=u1(drho)
		}
//#endif
	}
#else
	for (i=material->n-1;i>=0;i--)
	{
		double hs = h/10.0;
		
		/* We do substeps that saved to increase the accuracy. */
		for (s=0;s<10;s++)
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

void tillBSderivs(TILLMATERIAL *material, float x, float y[], float dydx[])
{
	/*
	** A function that provides the first derivative of u for bsstep.
	*/
}
#if 0
TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v)
{
	/*
	** Integrate one isentrope for the look up table. The parameter v corresponds to
	** u(rho=rho0) so v=0 gives the cold curve. We use the Bulirsch-Stoer method from
	** the Numerical Recipes.
	*/
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;
    int i,s;

	/* Use this as a temporary data structure because it is easy to sort with qsort. */
	TILL_LOOKUP_ENTRY *isentrope;

	/* Variables for bsstep. */
	float *y;
	float *xx;
	float *dydx;

	float *yscal;
	
	int nv = 1;

	float eps;

	float htry;
	float hdid;
	float hnext;

	/* Allocate memory */
	y = malloc(sizeof(float));
	xx = malloc(sizeof(float));
	dydx = malloc(sizeof(float));

	eps = 1e-8; 
	
	// Carefull: nr expects unit offset arrays!
    isentrope = malloc(material->nTableRho*sizeof(TILL_LOOKUP_ENTRY));

	rho = material->rho0;
	u = v;
	htry = material->drho;

	i = material->n;

	isentrope[i].rho = rho;
	isentrope[i].u = u;
	isentrope[i].u1 = tilldudrho(material, rho, u); // du/drho
	
	/*
	** Integrate the condensed and expanded states separately.
	*/
	for (i=material->n+1;i<material->nTableRho;i++)
	{

		
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
		double hs = h/100.0;
		
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

	/* Free memory */
	free(y);
	free(xx);
	free(dydx);
	return isentrope;
}
#endif 
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2)
{
	/* Calculate u2 by solving the ODE */
    double rho;
    double u;
    double k1u,k2u,k3u,k4u;
	double h;

	rho = rho1;
	u = u1;
	/* Make smaller steps than we used for look up table. */
	h = material->drho/100.0;

	if (rho1 < rho2)
	{
		while (rho < rho2) {
#ifdef TILL_USE_RK2
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*tilldudrho(material,rho,u);
			k2u = h*tilldudrho(material,rho+0.5*h,u+0.5*k1u);
	
			u += k2u;
#else
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = h*tilldudrho(material,rho,u);
			k2u = h*tilldudrho(material,rho+0.5*h,u+0.5*k1u);
			k3u = h*tilldudrho(material,rho+0.5*h,u+0.5*k2u);
			k4u = h*tilldudrho(material,rho+h,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
#endif
			rho += h;
		}
	} else if (rho1 > rho2) {
		while (rho > rho2) {
#ifdef TILL_USE_RK2
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1u = h*-tilldudrho(material,rho,u);
			k2u = h*-tilldudrho(material,rho+0.5*h,u+0.5*k1u);

			u += k2u;
#else
			/*
			** Midpoint Runga-Kutta (4nd order).
			*/
			k1u = h*tilldudrho(material,rho,u);
			k2u = h*tilldudrho(material,rho+0.5*h,u+0.5*k1u);
			k3u = h*tilldudrho(material,rho+0.5*h,u+0.5*k2u);
			k4u = h*tilldudrho(material,rho+h,u+k3u);

			u += k1u/6.0+k2u/3.0+k3u/3.0+k4u/6.0;
#endif
			rho -= h;
		}
	}
	return u;
}

int tillIsInTable(TILLMATERIAL *material,double rho,double u)
{
	/*
	** This function checks if a given (rho,u) is in our look up
	** table or not.
	**
	** Returns 0 if (rho,u) is in the table and 1 if not.
	*/
	int iRet = 1;

	/* Check if rho < rhomin or rho > rhomax */
	if (rho < material->rhomin || rho > material->rhomax)
	{
		return(iRet);
	}
	
	/* Check if u > u(rho,vmax) */
	if (u > tillCubicIntRho(material, rho, material->nTableV-1))
	{
		return(iRet);
	}

	/* Check if u < u(rho,0) where v=0 is the cold curve */
	if (u < tillCubicIntRho(material, rho, 0))
	{
		/* We are in the unphysical region below the cold curve */
//		fprintf(stderr,"tillIsInTable: value (%g,%g) below the cold curve!\n",rho,u);
		printf("tillIsInTable: value (%g,%g) below the cold curve (iMat=%i)!\n",rho,u,material->iMaterial);
		assert(0);
		return(iRet);
	}
	iRet = 0;
	return(iRet);
}


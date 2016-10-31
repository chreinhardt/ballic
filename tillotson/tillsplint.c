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
//#include "tillsplint.h"

/* Header file for the Numerical Recipes routines */
#include "nr/nrcubicspline.h"

#if 0
/* Code for the old interpolator */
#include "interpol/coeff.h"
#include "interpol/interpol.h"
#endif

/* Basic functions:
 *
 * tillInitSplines: initialise the cubic splines in u and u1.
 *
 * tillSplineIntU: do a cubic spline interpolation of u in v.
 *
 * tillSplineIntU1: do a cubic spline interpolation of u1 in v.
 *
 * tillCubicInt: interpolate a value (i,i+1) in rho and (j,j+1) in v (uses Joachim's 2D interpolator).
 */

void tillInitSplines(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u and u1 in v.
	** For this we use routines from Numerical Recipes.
	*/
	tillInitSplineU1(material);
	tillInitSplineU(material);
}

#if 0
/*
** Old function that uses standard Numerical Recipes functions.
*/
#ifdef TILL_DEBUG_SPLINT
void tillInitSplineRho(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u and u1 in rho.
	** For this we use routines from Numerical Recipes.
	*/
	double *x;
	double *y;
	int i,j,n;
	double yp1, ypn;
	double *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableRho;
	x = malloc(material->nTableRho*sizeof(double));
	y = malloc(material->nTableRho*sizeof(double));
	y2 = malloc(material->nTableRho*sizeof(double));

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
		rho = material->Lookup[TILL_INDEX(0,j)].rho;
		u = material->Lookup[TILL_INDEX(0,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		yp1 = 1.0/(rho*rho)*(c2 - 2.0*P/rho);


		rho = material->Lookup[TILL_INDEX(material->nTableRho-1,j)].rho;
		u = material->Lookup[TILL_INDEX(material->nTableRho-1,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		ypn = 1.0/(rho*rho)*(c2 - 2.0*P/rho);
*/
		/* Copy one row of the look up table to the temporary array. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			x[i] = 	material->Lookup[TILL_INDEX(i,j)].rho;
			y[i] = 	material->Lookup[TILL_INDEX(i,j)].u;
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[TILL_INDEX(i,j)].udrho2 = y2[i];
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}
#endif
void tillInitSplinev(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u in v.
	** For this we use routines from Numerical Recipes.
	*/
	double *x;
	double *y;
	int i,j,n;
	double yp1, ypn;
	double *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableV;
	x = malloc(material->nTableV*sizeof(double));
	y = malloc(material->nTableV*sizeof(double));
	y2 = malloc(material->nTableV*sizeof(double));

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
			y[j] = 	material->Lookup[TILL_INDEX(i,j)].u;
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
			material->Lookup[TILL_INDEX(i,j)].udv2 = y2[j];
			//printf("i: %i j: %i %g\n",i,j, y2[j]);
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}

void tillInitSplineU1(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u1 in v.
	** For this we use routines from Numerical Recipes.
	*/
	double *x;
	double *y;
	int i,j,n;
	double yp1, ypn;
	double *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableV;
	x = malloc(material->nTableV*sizeof(double));
	y = malloc(material->nTableV*sizeof(double));
	y2 = malloc(material->nTableV*sizeof(double));

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
			y[j] = 	material->Lookup[TILL_INDEX(i,j)].u1;
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (j=0; j<material->nTableV; j++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[TILL_INDEX(i,j)].u1dv2 = y2[j];
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}
#endif

/*
** (CR) 10.08.16: Adapting the lookup algorithm to our data structure.
*/
#ifdef TILL_DEBUG_SPLINT
void tillInitSplineRho(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u in rho.
	*/
	double *x;
	double *y;
	int i,j,n;
	double yp1, ypn;
	double *y2;
	
	/* Temporarily store data in an array to use spline() */
	n = material->nTableRho;
	x = malloc(material->nTableRho*sizeof(double));
	y = malloc(material->nTableRho*sizeof(double));
	y2 = malloc(material->nTableRho*sizeof(double));

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
		rho = material->Lookup[TILL_INDEX(0,j)].rho;
		u = material->Lookup[TILL_INDEX(0,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		yp1 = 1.0/(rho*rho)*(c2 - 2.0*P/rho);


		rho = material->Lookup[TILL_INDEX(material->nTableMax-1,j)].rho;
		u = material->Lookup[TILL_INDEX(material->nTableMax-1,j)].u;
		P = tillPressureSound(material, rho, u, &c2);
		ypn = 1.0/(rho*rho)*(c2 - 2.0*P/rho);
*/
		/* Copy one row of the look up table to the temporary array. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			x[i] = 	material->Lookup[TILL_INDEX(i,j)].rho;
			y[i] = 	material->Lookup[TILL_INDEX(i,j)].u;
		}
		
		// Careful nr expects unit offset
		spline(x-1, y-1, n, yp1, ypn, y2-1);
	
		/* Copy the second dervative back to the look up table. */
		for (i=0; i<material->nTableRho; i++)
		{
			/* Careful with the indices!
			** Lookup[i][j] = Lookup(rho,v)
			*/
			material->Lookup[TILL_INDEX(i,j)].udrho2 = y2[i];
		}
	}

	/* Free memory */
	free(x);
	free(y);
	free(y2);
}
#endif

void tillInitSplinev(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u in v.
	*/
	int i,j,k,n;
	double yp1, ypn;	// dudv(v=0) and dudv(v=n-1)
	double p,qn,sig,un,*u;

	n = material->nTableV;
	/* Allocate memory for temporary array */
	u = malloc(n*sizeof(double));
		
	// (CR) 15.11.15: Try non uniform steps in v
	// x[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
	// (CR) 15.11.15: Done
	
	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/*
		** Set up splines in u(v) for a given rho.
		*/

		/* Set b.c */	
		if (yp1 > 0.99e30)
			material->Lookup[TILL_INDEX(i,0)].udv2=u[0]=0.0;
		else {
			material->Lookup[TILL_INDEX(i,0)].udv2 = -0.5;
//			u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
			/* (CR) This part has to be changed if dv is not constant! */
			u[0]=(3.0/(material->dv))*((material->Lookup[TILL_INDEX(i,1)].u-material->Lookup[TILL_INDEX(i,0)].u)/(material->dv)-yp1);
		}

		/* Calculate second derivatives of u in v */
		for (j=1;j<=n-2;j++) {
//			sig=(x[j]-x[j-1])/(x[j+1]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			sig=material->dv/(2.0*material->dv);
			p=sig*material->Lookup[TILL_INDEX(i,j-1)].udv2+2.0;
			material->Lookup[TILL_INDEX(i,j)].udv2=(sig-1.0)/p;
//			u[j]=(y[j+1]-y[j])/(x[j+1]-x[j]) - (y[j]-y[j-1])/(x[j]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(material->Lookup[TILL_INDEX(i,j+1)].u-material->Lookup[TILL_INDEX(i,j)].u)/(material->dv) - (material->Lookup[TILL_INDEX(i,j)].u-material->Lookup[TILL_INDEX(i,j-1)].u)/(material->dv);
//			u[j]=(6.0*u[j]/(x[j+1]-x[j-1])-sig*u[j-1])/p;
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(6.0*u[j]/(2.0*material->dv)-sig*u[j-1])/p;
		}

		/* Set b.c. */
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
//			Careful with index!! This array ends at n-1, so this code is already using the right indices.			
//			un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
			/* (CR) This part has to be changed if dv is not constant! */
			un=(3.0/(material->dv))*(ypn-(material->Lookup[TILL_INDEX(i,n-1)].u-material->Lookup[TILL_INDEX(i,n-2)].u)/(material->dv));
		}
		material->Lookup[TILL_INDEX(i,n-1)].udv2=(un-qn*u[n-2])/(qn*material->Lookup[TILL_INDEX(i,n-2)].udv2+1.0);
		for (k=n-2;k>=0;k--)
			material->Lookup[TILL_INDEX(i,k)].udv2=material->Lookup[TILL_INDEX(i,k)].udv2*material->Lookup[TILL_INDEX(i,k+1)].udv2+u[k];
	}
}

void tillInitSplineU1(TILLMATERIAL *material)
{
	/*
	** Calculate the second derivatives for u1 in v.
	** For this we use routines from Numerical Recipes.
	*/
	int i,j,k,n;
	double yp1, ypn;	// du1dv(v=0) and d1udv(v=n-1)
	double p,qn,sig,un,*u;

	n = material->nTableV;
	/* Allocate memory for temporary array */
	u = malloc(n*sizeof(double));
		
	// (CR) 15.11.15: Try non uniform steps in v
	// x[j] =  material->vmax/pow(material->nTableV-1,material->iExpV)*pow(j,material->iExpV);
	// (CR) 15.11.15: Done
	
	/* Set b.c. for natural cubic spline */
	yp1 = 1e30;
	ypn = 1e30;

	for (i=0; i<material->nTableRho; i++)
	{
		/*
		** Set up splines in u(v) for a given rho.
		*/

		/* Set b.c */	
		if (yp1 > 0.99e30)
			material->Lookup[TILL_INDEX(i,0)].u1dv2=u[0]=0.0;
		else {
			material->Lookup[TILL_INDEX(i,0)].u1dv2 = -0.5;
//			u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
			/* (CR) This part has to be changed if dv is not constant! */
			u[0]=(3.0/(material->dv))*((material->Lookup[TILL_INDEX(i,1)].u-material->Lookup[TILL_INDEX(i,0)].u)/(material->dv)-yp1);
		}

		/* Calculate second derivatives of u in v */
		for (j=1;j<=n-2;j++) {
//			sig=(x[j]-x[j-1])/(x[j+1]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			sig=material->dv/(2.0*material->dv);
			p=sig*material->Lookup[TILL_INDEX(i,j-1)].u1dv2+2.0;
			material->Lookup[TILL_INDEX(i,j)].u1dv2=(sig-1.0)/p;
//			u[j]=(y[j+1]-y[j])/(x[j+1]-x[j]) - (y[j]-y[j-1])/(x[j]-x[j-1]);
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(material->Lookup[TILL_INDEX(i,j+1)].u-material->Lookup[TILL_INDEX(i,j)].u)/(material->dv) - (material->Lookup[TILL_INDEX(i,j)].u-material->Lookup[TILL_INDEX(i,j-1)].u)/(material->dv);
//			u[j]=(6.0*u[j]/(x[j+1]-x[j-1])-sig*u[j-1])/p;
			/* (CR) This part has to be changed if dv is not constant! */
			u[j]=(6.0*u[j]/(2.0*material->dv)-sig*u[j-1])/p;
		}

		/* Set b.c. */
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
//			Careful with index!! This array ends at n-1, so this code is already using the right indices.			
//			un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
			/* (CR) This part has to be changed if dv is not constant! */
			un=(3.0/(material->dv))*(ypn-(material->Lookup[TILL_INDEX(i,n-1)].u-material->Lookup[TILL_INDEX(i,n-2)].u)/(material->dv));
		}
		material->Lookup[TILL_INDEX(i,n-1)].u1dv2=(un-qn*u[n-2])/(qn*material->Lookup[TILL_INDEX(i,n-2)].u1dv2+1.0);
		for (k=n-2;k>=0;k--)
			material->Lookup[TILL_INDEX(i,k)].u1dv2=material->Lookup[TILL_INDEX(i,k)].u1dv2*material->Lookup[TILL_INDEX(i,k+1)].u1dv2+u[k];
	}
}

void tillInitSplineU(TILLMATERIAL *material)
{
	tillInitSplinev(material);
}

#if 0
/*
** This is an old version that has to copy all data to an extra array for each lookup
** in order to be compatible with the standard Numerical Recipes routines.
*/

/*
** Just to debug the code we do an interpolation in u of rho.
*/
#ifdef TILL_DEBUG_SPLINT
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv)
{
	/*
	** Do a cubic spline interpolation.
	*/
	double *xa;
	double *ya;
	double *y2a;
	double x,y;
	int i,n;
	
	/* Temporarily store data in an array to use splint() */
	n = material->nTableRho;
	xa = malloc(material->nTableRho*sizeof(double));
	ya = malloc(material->nTableRho*sizeof(double));
	y2a = malloc(material->nTableRho*sizeof(double));

	//void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)

	x = rho;

	/* Copy one row of the look up table to the temporary array. */
	for (i=0; i<material->nTableRho; i++)
	{
		/* Careful with the indices!
		** Lookup[i][j] = Lookup(rho,v)
		*/
		xa[i] =	material->Lookup[TILL_INDEX(i,iv)].rho;
		ya[i] =	material->Lookup[TILL_INDEX(i,iv)].u;
		y2a[i]=	material->Lookup[TILL_INDEX(i,iv)].udrho2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}
#endif
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
	double *xa;
	double *ya;
	double *y2a;
	double x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(double));
	ya = malloc(material->nTableV*sizeof(double));
	y2a = malloc(material->nTableV*sizeof(double));

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
		ya[j] =	material->Lookup[TILL_INDEX(irho,j)].u;
		y2a[j]=	material->Lookup[TILL_INDEX(irho,j)].udv2;
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
	double *xa;
	double *ya;
	double *y2a;
	double x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(double));
	ya = malloc(material->nTableV*sizeof(double));
	y2a = malloc(material->nTableV*sizeof(double));

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
		ya[j] =	material->Lookup[TILL_INDEX(irho,j)].u;
		y2a[j]=	material->Lookup[TILL_INDEX(irho,j)].udv2;
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
	double *xa;
	double *ya;
	double *y2a;
	double x,y;
	int j,n;

	/* Temporarily store data in an array to use splint() */
	n = material->nTableV;
	xa = malloc(material->nTableV*sizeof(double));
	ya = malloc(material->nTableV*sizeof(double));
	y2a = malloc(material->nTableV*sizeof(double));

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
		ya[j] =	material->Lookup[TILL_INDEX(irho,j)].u1;
		y2a[j]=	material->Lookup[TILL_INDEX(irho,j)].u1dv2;
	}
	
	// Careful nr expects unit offset
	splint(xa-1, ya-1, y2a-1, n, x, &y);

	/* Free memory */
	free(xa);
	free(ya);
	free(y2a);

	return y;
}
#endif

/*
** (CR) 10.08.2016: Integrated splint() into the function, allowing to use the
** TILL_LOOKUP_ENTRY structure for the interpolation.
*/
double tillSplineIntv(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation of u in v.
	*/
	int klo,khi,k;
	double h,b,a;
	double u;
/*
	klo=1;
	khi=material->nTableV;
	
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (material->Lookup[TILL_INDEX(irho,k)].v > v) khi=k;
		else klo=k;
	}
	
	h=material->Lookup[TILL_INDEX(irho,khi)].v-material->Lookup[TILL_INDEX(irho,klo)].v;
//	if (h == 0.0) nrerror("Bad xa input to routine splint");
	assert(h != 0.0);
	a=(material->Lookup[TILL_INDEX(irho,khi)].v-v)/h;
	b=(v-material->Lookup[TILL_INDEX(irho,klo)])/h;
*/
	/* Use v=k*material->dv */
	
	klo = floor(v/material->dv);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;

	u=a*material->Lookup[TILL_INDEX(irho,klo)].u+b*material->Lookup[TILL_INDEX(irho,khi)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(irho,klo)].udv2+(b*b*b-b)*material->Lookup[TILL_INDEX(irho,khi)].udv2)*(h*h)/6.0;

	return u;
}

double tillSplineIntU(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation of u in v.
	*/
	int klo,khi,k;
	double h,b,a;
	double u;
/*
	klo=1;
	khi=material->nTableV;
	
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (material->Lookup[TILL_INDEX(irho,k)].v > v) khi=k;
		else klo=k;
	}
	
	h=material->Lookup[TILL_INDEX(irho,khi)].v-material->Lookup[TILL_INDEX(irho,klo)].v;
//	if (h == 0.0) nrerror("Bad xa input to routine splint");
	assert(h != 0.0);
	a=(material->Lookup[TILL_INDEX(irho,khi)].v-v)/h;
	b=(v-material->Lookup[TILL_INDEX(irho,klo)])/h;
*/
	/* Use v=k*material->dv */
	
	klo = floor(v/material->dv);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;

	u=a*material->Lookup[TILL_INDEX(irho,klo)].u+b*material->Lookup[TILL_INDEX(irho,khi)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(irho,klo)].udv2+(b*b*b-b)*material->Lookup[TILL_INDEX(irho,khi)].udv2)*(h*h)/6.0;

	return u;
}

double tillSplineIntU1(TILLMATERIAL *material, double v, int irho)
{
	/*
	** Do a cubic spline interpolation of u1 in v.
	*/
	int klo,khi,k;
	double h,b,a;
	double u1;
	
	klo = floor(v/material->dv);
	khi = klo+1;

	h = (khi-klo)*material->dv;
	assert(h != 0.0);
	assert(h == material->dv);

	a=(khi*material->dv-v)/h;
	b=(v-klo*material->dv)/h;
	
	u1=a*material->Lookup[TILL_INDEX(irho,klo)].u1+b*material->Lookup[TILL_INDEX(irho,khi)].u1+((a*a*a-a)*material->Lookup[TILL_INDEX(irho,klo)].u1dv2+(b*b*b-b)*material->Lookup[TILL_INDEX(irho,khi)].u1dv2)*(h*h)/6.0;

	return u1;
}

#ifdef TILL_DEBUG_SPLINT
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv)
{
	/*
	** Do a cubic spline interpolation in rho.
	*/
	int klo,khi,k;
	double h,b,a;
	double u;
	
	klo = floor(rho/material->drho);
	khi = klo+1;

	h = (khi-klo)*material->drho;
	assert(h != 0.0);
	assert(h == material->drho);

	a=(khi*material->drho-rho)/h;
	b=(rho-klo*material->drho)/h;
	
	return a*material->Lookup[TILL_INDEX(klo,iv)].u+b*material->Lookup[TILL_INDEX(khi,iv)].u+((a*a*a-a)*material->Lookup[TILL_INDEX(klo,iv)].udrho2+(b*b*b-b)*material->Lookup[TILL_INDEX(khi,iv)].udrho2)*(h*h)/6.0;
}
#endif

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

double tillCubicInt(TILLMATERIAL *material, double rhoint, double vint) {
	/*
	** Use cubicint to interpolate u for a given rho and v.
	*/
	double dv, A, B;
	int i, j;
	double *u, *dudrho, *dudv, *dudvdrho, *rho, *intvalues;
	double uint;

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
	
	dudv[0] = (material->Lookup[TILL_INDEX(i, j+1)].u-material->Lookup[TILL_INDEX(i, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j+1)].udv2;
	dudv[1] = (material->Lookup[TILL_INDEX(i+1, j+1)].u-material->Lookup[TILL_INDEX(i+1, j)].u)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j)].udv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j+1)].udv2;

	/* Calculate dudvdrho from u1dv2 */
	dudvdrho[0] = (material->Lookup[TILL_INDEX(i, j+1)].u1-material->Lookup[TILL_INDEX(i, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i, j+1)].u1dv2;
	dudvdrho[1] = (material->Lookup[TILL_INDEX(i+1, j+1)].u1-material->Lookup[TILL_INDEX(i+1, j)].u1)/dv-(3.0*A*A-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j)].u1dv2+(3.0*B*B-1.0)/6.0*dv*material->Lookup[TILL_INDEX(i+1, j+1)].u1dv2;

	/*
	** Do the interpolation for u(i,v), u(i+1,v), udrho(i,v), udrho(i+1,v)
	*/
	cubicint(u, dudrho, dudv, dudvdrho, rho, rhoint, intvalues);
	
	uint = intvalues[0];

	/* Free memory */
	free(u);
	free(dudrho);
	free(dudv);
	free(dudvdrho);
	free(rho);
	free(intvalues);

	return(uint);
}

/* Prasenjits root finder */
double brent(double (*func)(TILLMATERIAL *,double,double,double),TILLMATERIAL *material,double a,double b,double rho,double u,double tol,int iOrder);

double tillFindUonIsentrope(TILLMATERIAL *material,double v,double rho)
{
	double iv,irho,u;
	/* Needed for the interpolation function. */
	//iv = (material->nTableMax-1)*v/material->vmax;
	//irho = (material->nTableMax-1)*rho/material->rhomax;

//	return InterpolatedValue(material->Lookup,material->nTableMax,material->nTableMax,iv,irho,TILL_SPLINE_DEGREE);
	return (tillCubicInt(material, rho, v));
}

double denergy(TILLMATERIAL *material,double v,double rho,double u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder)
{
	double tol=1e-6;

	return brent(denergy,material,0.0,material->vmax-material->dv,rho,u,tol,iOrder);
//	return brent(denergy,material,0.0,material->vmax,rho,u,tol,iOrder);
}

double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder)
{
	/* Calculates u2 for a given rho1,u2,rho2. */
	double v, u;

	/* Check if the starting and end point are inside of the look up table */
	if (tillIsInTable(material, rho1, u1) != 0 || rho2 < material->rhomin || rho2 > material->rhomax)
	{
//		fprintf(stderr,"tillLookupU: values outside of look up table, doing direct integration.\n");
#ifdef TILL_OUTPUT_ALL_WARNINGS
		printf("tillLookupU: values outside of look up table, doing direct integration.\n");
#endif
		/* In this case we either about with an error or do a direct integration using tillCalcU(). */
		u = tillCalcU(material, rho1, u1, rho2);
	} else {
		/* Interpolate using the look up table */
		v = tillFindEntropyCurve(material,rho1,u1,iOrder);
		u = tillFindUonIsentrope(material,v,rho2);
	}
	return u;
//	return tillFindUonIsentrope(material,v,rho2);
}

#ifdef TILL_USE_OLD_BCINT
/*
** This code uses the old interpolator we used in the first version of the code.
*/
double tillFindUonIsentrope(TILLMATERIAL *material,double v,double rho)
{
	double iv,irho,u;
	/* Needed for the interpolation function. */
	//iv = (material->nTableMax-1)*v/material->vmax;
	//irho = (material->nTableMax-1)*rho/material->rhomax;

	return InterpolatedValue(material->Lookup,material->nTableMax,material->nTableMax,iv,irho,TILL_SPLINE_DEGREE);
}

double denergy(TILLMATERIAL *material,double v,double rho,double u)
{
	return (tillFindUonIsentrope(material,v,rho)-u);
}

/* Find isentrope for a given rho and u */
double tillFindEntropyCurve(TILLMATERIAL *material,double rho,double u,int iOrder)
{
	double tol=1e-6;
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

double tillCubicIntRho(TILLMATERIAL *material, double rhoint, int iv) {
	/* Do an interpolation of u in rho for a given isentrope v. */
	double dx, e, e1;
	double *ce;
	double *u, *dudrho;
	double uint;
	int i;

	// rhoint is between rho[i] and rho[i+1]
	i = floor(rhoint/material->drho);
	assert(i < material->nTableRho-1);

	/* Allocate memory */
	ce = malloc(4*sizeof(double));
	assert(ce != NULL);

	u = malloc(2*sizeof(double));
	assert(u != NULL);

	dudrho = malloc(2*sizeof(double));
	assert(dudrho != NULL);

	/* Get u and dudrho from the lookup table */
	u[0] = material->Lookup[TILL_INDEX(i, iv)].u;
	u[1] = material->Lookup[TILL_INDEX(i+1, iv)].u;
	dudrho[0] = material->Lookup[TILL_INDEX(i, iv)].u1;
	dudrho[1] = material->Lookup[TILL_INDEX(i+1, iv)].u1;
	
	//dx = rho[1] - rho[0];
	//e = (rhoint - rho[0])/dx;
	//e1 = e - 1;
	
	// Only works for uniform steps in rho
	dx = material->drho;
	e = (rhoint - i*material->drho)/dx;
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
	uint = (ce[0]*u[0] + ce[1]*dudrho[0]*dx + ce[2]*u[1] + ce[3]*dudrho[1]*dx);
	
	// free memory
	free(ce);
	free(u);
	free(dudrho);
	return(uint);
}

double tillColdULookup(TILLMATERIAL *material,double rho)
{
	// (CR) 6.1.2016: Updated the code and tested it!
	// However we have to still do something in case that we are not inside of the look up table
	// Make sure that the look up table is initialized.
	assert(material->Lookup != NULL);

	// iv = 0 corresponds to the cold curve
	return(tillCubicIntRho(material, rho, 0));
}


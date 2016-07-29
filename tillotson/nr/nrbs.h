/*
 ** The header file for the Bulirsch-Stoer integrator from the Numerical Recipes book.
 */
#ifndef TILLOTSON_NRBS_HINCLUDED
#define TILLOTSON_NRBS_HINCLUDED

#include "nrutil.h"

void mmid(float y[], float dydx[], int nvar, float xs, float htot, int nstep,
	float yout[], void (*derivs)(float, float[], float[]));

void pzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);

/*
** Bulirsch-Stoer method to integrate ODEs.
**
** y[]:			dependent variable
** dydx[]:		its first derivative at the starting value x
** nv:			number of variables y1,...,yn
** xx:			
** htry:		step size (the algorithm can use a smaller value if needed)
** eps:			required accuracy
** yscal[]:		vector to scale the error
** hdid:		return actual step size
** hnext:		return estimated next step size
** derivs:		function to calculate the right hand side derivatives
*/
void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));

void odeint(float ystart[], int nvar, float x1, float x2, float eps, float h1,
	float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*rkqs)(float [], float [], int, float *, float, float, float [],
	float *, float *, void (*)(float, float [], float [])));

#endif


/*
 ** The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_INITLOOKUP_HINCLUDED
#define TILLOTSON_INITLOOKUP_HINCLUDED

#include "tillotson.h"
#include "tillsplint.h"

#include "nr/nrcubicspline.h"

void tillInitColdCurve(TILLMATERIAL *material);
void tillInitLookup(TILLMATERIAL *material);
TILL_LOOKUP_ENTRY *tillSolveIsentrope(TILLMATERIAL *material, double v);
/* Use bsstep.c from the Numerical Recipes */
TILL_LOOKUP_ENTRY *tillSolveIsentropeBS(TILLMATERIAL *material, double v);
double tillCalcU(TILLMATERIAL *material,double rho1,double u1,double rho2);

/* Defines for the Numerical Recipes routines */

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

#endif


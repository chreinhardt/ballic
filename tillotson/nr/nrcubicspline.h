/*
 ** The header file for cubic spline interpolator from the Numerical Recipes book.
 */
#ifndef TILLOTSON_NRCUBICSPLINE_HINCLUDED
#define TILLOTSON_NRCUBICSPLINE_HINCLUDED

#include "nrutil.h"
#include "../tillotson.h"

/* Defines for the Numerical Recipes routines */

/*
** Solver for a tridiagonal set of equations (needed to do calculate the second derivatives
** for the cubic spline).
**
** a[], b[], c[]:	input vectors containing the non zero elements of the matrix A
** r[]:				right side of the equation
** n:				number of elements (1 to n)
** u[]:				return vector that contains the solution
*/
void tridag(double a[], double b[], double c[], double r[], double u[], unsigned long n);

/*
** Calculate the second derivatives of an interpolated function y(x) for splint().
**
** x[], y[]:	tabulated function y(x)
** n:			number of elements (1 to n)
** yp1[]:		first derivative at x[1] (define boundary conditions for the algorithm)
** ypn[]:		first derivative at x[n] (use natural cubic spline is > 1e30)
** y2[]:		return second derivatives of y(x) 
*/
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);

/*
** Do a cubic spline interpolation.
**
** xa[], ya[]:	tabulated function y(x)
** y2a[]:		second derivatives of y(x) (caluclated with spline())
** n:			number of elements (1 to n)
** x, y:		return interpolated value at x
*/
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
//void splint(TILLMATERIAL *material, int v, double rho, double *u);

#endif


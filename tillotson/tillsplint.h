/*
 ** The header file for the Tillotson EOS library.
 */
#ifndef TILLOTSON_SPLINT_HINCLUDED
#define TILLOTSON_SPLINT_HINCLUDED

#include "tillotson.h"
/* Header file for the Numerical Recipes routines */
#include "nr/nrcubicspline.h"

#include "interpol/coeff.h"
#include "interpol/interpol.h"

/* Stuff for the cubic spline interpolator */
void tillInitSplines(TILLMATERIAL *material);
void tillInitSplineRho(TILLMATERIAL *material);
void tillInitSplinev(TILLMATERIAL *material);
void tillInitSplineU(TILLMATERIAL *material);
void tillInitSplineU1(TILLMATERIAL *material);
// Just for debugging
double tillSplineIntrho(TILLMATERIAL *material, double rho, int iv);
double tillSplineIntv(TILLMATERIAL *material, double v, int irho);
double tillSplineIntU(TILLMATERIAL *material, double v, int irho);
double tillSplineIntU1(TILLMATERIAL *material, double v, int irho);

void cubicint(double u[2],double dudrho[2], double dudv[2], double dudvdrho[2], double rho[2], double rhoint, double *intvalues);
double tillCubicInt(TILLMATERIAL *material, double rho, double v);

float tillFindUonIsentrope(TILLMATERIAL *material,float v,float rho);
/* Used for the root finder */
float denergy(TILLMATERIAL *material,float v,float rho,float u);
float tillFindEntropyCurve(TILLMATERIAL *material,float rho,float u,int iOrder);
double tillLookupU(TILLMATERIAL *material,double rho1,double u1,double rho2,int iOrder);
double tillColdULookup(TILLMATERIAL *material,double rho);

#endif


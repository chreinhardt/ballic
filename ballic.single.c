/*
 ** Copyright (c) 2014-2016 Joachim Stadel and Christian Reinhardt.
 **
 ** ballic provides a low noise particle representation of equilibrium
 ** models.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsy.h"
#include "tillotson/tillotson.h"
//#include "ballic.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

// Isentropic thermal profile
#define BALLIC_U_ISENTROPIC

typedef struct icosa_struct {
    float R[180];
    float v[36];
    } ICOSA;

ICOSA *icosaInit(void) {
    ICOSA *ctx;
    ctx = malloc(sizeof(ICOSA));
    assert(ctx != NULL);
    compute_matrices_(&ctx->R);
    compute_corners_(&ctx->v);
    return(ctx);
    }

void icosaPix2Vec(ICOSA *ctx,int i,int resolution,double *vec) {
    float v[3];
    pixel2vector_(&i,&resolution,&ctx->R,&ctx->v,v);
    vec[0] = v[0];
    vec[1] = v[1];
    vec[2] = v[2];
    }


/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
void pix2vec_ring( long nside, long ipix, double *vec) {
  /*
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING) 
    c     for a parameter nside
    c=======================================================================
  */
  
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double  fact1, fact2, fodd, hip, fihip;
  double PI=M_PI;
  double z, sz, phi;
  //      PARAMETER (pi     = 3.1415926535897932384626434d0)
  //      parameter (ns_max = 8192) ! 2^13 : largest nside available
  
  int ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12*nside*nside;      // ! total number of points
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
  }
  
  ipix1 = ipix + 1; // in {1, npix}
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if( ipix1 <= ncap ) {  //! North Polar cap -------------
    
    hip   = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1);
    
    z = 1. - iring*iring / fact2 ;
    phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
  else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
    
    ip    = ipix1 - ncap - 1;
    iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
    iphi  = (int)fmod(ip,nl4) + 1;
    
    fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
    z = (nl2 - iring) / fact1;
    phi   = (1.*iphi - fodd) * PI /(2.*nside);
  }
  else {//! South Polar cap -----------------------------------
    
    ip    = npix - ipix1 + 1;
    hip   = ip/2.;
/* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
    
    z = -1. + iring*iring / fact2 ;
    phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }

  sz = sqrt( 1.0 - z*z );
  vec[0] = sz * cos(phi);
  vec[1] = sz * sin(phi);
  vec[2] = z;
  
}

double Packed49[49][3] = {
    {1.0820379E-008,-3.2300459E-006,    0.7790824},
    { 0.0851419,      -0.0604194,       0.3497294},
    {-0.0914288,       0.4137149,       0.6537967},
    {-0.4233704,      -0.0167043,       0.6537949},
    { 0.4161716,       0.0795128,       0.6537953},
    { 0.2606476,      -0.3340458,       0.6537932},
    {-0.1783143,      -0.3843534,       0.6537934},
    { 0.3214161,       0.4881783,       0.5151146},
    {-0.4918671,       0.3793286,       0.4702615},
    {-0.0929470,       0.3254455,       0.2208710},
    {-0.5385916,      -0.3977282,       0.3983726},
    { 0.6544342,      -0.1767253,       0.3839966},
    {-0.0464433,      -0.6884808,       0.3616719},
    { 0.6543453 ,      0.2637903 ,      0.3304789},
    { 0.3839234,      -0.5971428,       0.3209247},
    {-0.7108776,       0.0012718,       0.3187802},
    {-0.1876019,      -0.3352532,       0.1368439},
    { 0.0704146,       0.7351211,       0.2481284},
    {-0.3638534,       0.6695231,       0.1622308},
    { 0.3033485,       0.2022102,       0.0692749},
    { 0.4742969,       0.6067250,       0.1178835},
    {-0.6795831,       0.3744220,       0.0703152},
    {-0.3794767,       0.0482692,       0.0303653},
    { 0.6538247,      -0.4232979,       0.0173632},
    { 0.2332925,      -0.2885923,       0.0015460},
    { 0.7790813, -5.7994478E-013,      -0.0012951},
    {-0.5429419,      -0.5585797,      -0.0131201},
    {-0.1452212,      -0.7628716,      -0.0625065},
    {-0.7541588,      -0.1768841,      -0.0832219},
    { 0.2928920,      -0.7156702,      -0.0948672},
    {-0.0815266,       0.7567586,      -0.1662504},
    { 0.6534047,       0.3539813,      -0.2339385},
    {-0.4662671,       0.5642231,      -0.2668646},
    { 0.3250845,       0.6429856,      -0.2964101},
    {-0.6822678,       0.1837015,      -0.3282282},
    { 0.4930282,      -0.4578927,      -0.3927173},
    {-0.3428409,      -0.5690840,      -0.4069064},
    { 0.6530941,      -0.0471244,      -0.4221572},
    {-0.1036092,       0.2867325,      -0.2191358},
    { 0.0858176,      -0.5795982,      -0.5134887},
    {-0.5513357,      -0.1947856,      -0.5148367},
    { 0.0032634,       0.5253046,      -0.5753380},
    {-0.1660916,      -0.1851457,      -0.2781839},
    { 0.3945018,       0.3205518,      -0.5904102},
    {-0.3693146,       0.2921726,      -0.6206539},
    { 0.2268184,       0.0121013,      -0.3221586},
    { 0.2750635,      -0.2203113,      -0.6948182},
    {-0.1633231,      -0.2729136,      -0.7112054},
    { 0.0133638,       0.1279994,      -0.7683794}
    };


typedef struct model_ctx {
	/* Material coefficients from the Tillotson EOS. */
	TILLMATERIAL *tillMat;
	/*
	** Some unit conversion factors.
	*/
	double dKpcUnit;
	double dMsolUnit;
	/*
	** The lookup table for the equilibrium model.
	*/
	int nTableMax;
	int nTable;
	double uc; /* u at r = 0 */
	double *M;
	double *rho;
	double *u;
	double *r;
	double dr;
	double R;
	} MODEL;


MODEL *modelInit(double ucore, int iMat) {
    /* Initialize the model */
	MODEL *model;
    
    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
	model->tillMat = malloc(sizeof(TILLMATERIAL));

	/*
	** Initialize one material.
	** i=0: Granite
	** i=1: Iron
	** i=2: Basalt
	** i=3: Ice
	*/
	model->tillMat = tillInitMaterial(iMat, model->dKpcUnit, model->dMsolUnit, 100, 100, 50.0, 50.0, 1);

	// Debug information

	fprintf(stderr,"\n");	
	fprintf(stderr,"Material: %i\n",iMat);	
	fprintf(stderr,"a: %g\n", model->tillMat->a);
	fprintf(stderr,"b: %g\n", model->tillMat->b);
	fprintf(stderr,"A: %g\n", model->tillMat->A);
	fprintf(stderr,"B: %g\n", model->tillMat->B);
	fprintf(stderr,"rho0: %g\n", model->tillMat->rho0);
	fprintf(stderr,"u0: %g\n", model->tillMat->u0);
	fprintf(stderr,"us: %g\n", model->tillMat->us);
	fprintf(stderr,"us2: %g\n", model->tillMat->us2);
	fprintf(stderr,"alpha: %g\n", model->tillMat->alpha);
	fprintf(stderr,"beta: %g\n", model->tillMat->beta);
	fprintf(stderr,"cv: %g\n", model->tillMat->cv);
	fprintf(stderr,"\n");

    /* model->uFixed = uFixed/model->dErgPerGmUnit; */
    model->uc = ucore;

    model->nTableMax = 10000; 
    model->M = malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);
    model->rho = malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);
    model->u = malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);
    model->r = malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);
    model->dr =  0.0;
    model->nTable = 0;
    
    return(model);
    }

#if 0
double Gamma(MODEL *model,double rho,double u) {
    double eta = rho/model->par.rho0;
    double w0 = u/(model->par.u0*eta*eta) + 1.0;

    return(model->par.a + model->par.b/w0);
    }

/*
** Currently for condensed states only!
*/
double Pressure(MODEL *model,double rho,double u) {
    double mu = rho/model->par.rho0 - 1.0;
    return(Gamma(model,rho,u)*rho*u + model->par.A*mu + model->par.B*mu*mu);
    }

double drhodr(MODEL *model,double r,double rho,double M,double u);
double dudrho(MODEL *model,double rho,double u);
double dudr(MODEL *model,double r,double rho,double M,double u);

/*
** Currently for condensed states only! Changed this 8.2.13 to include a polytropic temperature profile.
*/
double drhodr(MODEL *model,double r,double rho,double M,double u) {
    double eta = rho/model->par.rho0;
    double w0 = u/(model->par.u0*eta*eta) + 1.0;
    double dPdrho,dPdu;

//    fprintf(stderr,"u in system units: %.14g\n",u);
	dPdrho=tilldPdrho; // dP/drho at u=const.
    dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*u + 
	(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

    dPdu = (model->par.a + model->par.b/(w0*w0))*rho;
	dPdu = tilldPdu; // dP/du at rho=const.

    /*
    ** Here is the old version of dPdrho which assumes Gamma is constant for the derivative and
    ** as such is not quite self consistent. This makes about a 3% difference in the final 
    ** radius of the planet. Also u is assumed to be constant here!
    **
    dPdrho = Pressure(model,rho,model->uFixed)/rho + (model->par.A + model->par.B*(eta*eta - 1.0))/rho;
    */
    assert(r >= 0.0);
    if (r > 0.0) {
      return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
	}
    else {
	return(0.0);
	}
    }

/*
** We assume an isentropic internal energy profile!
*/
double dudrho(MODEL *model,double rho,double u) {
  return(Pressure(model,rho,u)/(rho*rho));
}

double dudr(MODEL *model,double r,double rho,double M,double u) {
  return(dudrho(model,rho,u)*drhodr(model,r,rho,M,u));
}
#endif

double drhodr(MODEL *model,double r,double rho,double M,double u);
double dudr(MODEL *model,double r,double rho,double M,double u);

/*
** dudrho depends on the internal energy profile that we choose!
*/
double dudrho(MODEL *model,double rho,double u) {
#ifdef BALLIC_U_POLYTROPIC
	/*
	** We assume a polytropic internal energy profile!
	*/
	// Not implemented yet
	assert(0);
	return(Pressure(model,rho,u)/(rho*rho));
#else
#ifdef BALLIC_U_ISENTROPIC
	/*
	** We assume an isentropic internal energy profile!
	*/
	return(tillPressure(model->tillMat,rho,u)/(rho*rho));
#else
  	fprintf(stderr,"No thermal profile defined when compiled!\n");
	assert(0);
#endif
#endif
}

/*
** Calculate dudrho to solve for the equilibrium model.
*/
double drhodr(MODEL *model,double r,double rho,double M,double u) {
    double dPdrho,dPdu;

	dPdrho=tilldPdrho(model->tillMat, rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat, rho, u);; // dP/du at rho=const.

	/*
	** drho/dr = -G*M*rho/(dPdrho+dPdu*dudrho)
	*/
	assert(r >= 0.0);
	if (r > 0.0) {
	// We assume G=1
		return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,rho,u))));
	}
	else {
		return(0.0);
	}
}


double dudr(MODEL *model,double r,double rho,double M,double u) {
	return(dudrho(model,rho,u)*drhodr(model,r,rho,M,u));
//	return(0.0);
}

/*
** This derivative is independent of the model and only involves geometry.
*/
double dMdr(double r,double rho) {
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

const double fact = 1.0;

/*
** This function solves the model as an initial value problem with rho_initial = rho and 
** M_initial = 0 at r = 0. This function returns the mass when rho == model->tillMat[i]->rho0.
*/
double midPtRK(MODEL *model,int bSetModel,double rho,double h,double *pR) {
    FILE *fp;
    double M = 0.0;
    double r = 0.0;
	double u = model->uc;
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
    int i;

    if (bSetModel) {
		i = 0;
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		fp = fopen("ballic.model","w");
		assert(fp != NULL);
		/* Output in temperature! */
		fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
		++i;
	}

    while (rho > fact*model->tillMat->rho0) {
	/*
	** Midpoint Runga-Kutta (2nd order).
	*/
	k1rho = h*drhodr(model,r,rho,M,u);
	k1M = h*dMdr(r,rho);
	k1u = h*dudr(model,r,rho,M,u);

	k2rho = h*drhodr(model,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
	k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
	k2u = h*dudr(model,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

	rho += k2rho;
	M += k2M;
	u += k2u;
	r += h;

	if (bSetModel) {
	    model->rho[i] = rho;
	    model->M[i] = M;
	    model->u[i] = u;
	    model->r[i] = r;
	    fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
	    ++i;
	    }
	}
    /*
    ** Now do a linear interpolation to rho == fact*rho0.
    */
    x = (fact*model->tillMat->rho0 - rho)/k2rho;
    assert(x <= 0.0);
    r += h*x;
    M += k2M*x;
    rho += k2rho*x;
	u += k2u*x;

    if (bSetModel) {
	--i;
	model->M[i] = M;
	model->r[i] = r;
	model->rho[i] = rho;
	model->u[i] = u;

	fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
    fclose(fp);
	++i;
	model->nTable = i;
	model->dr = h;
	}
    *pR = r;
    return(M);
    }


double modelSolve(MODEL *model,double M) {
    const int nStepsMax = 10000;
    int bSetModel;
    double rmax;
    double dr,R;
    double a,Ma,b,Mb,c,Mc;

    /*
    ** First estimate the maximum possible radius.
    */
    R = cbrt(3.0*M/(4.0*M_PI*model->tillMat->rho0));
    dr = R/nStepsMax;
    a = 1.01*model->tillMat->rho0; /* starts with 1% larger central density */
    Ma = midPtRK(model,bSetModel=0,a,dr,&R);
    fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
    b = a;
    Mb = 0.5*M;
    while (Ma > M) {
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat->rho0 + a);
		Ma = midPtRK(model,bSetModel=0,a,dr,&R);
	}
    while (Mb < M) {
		b = 2.0*b;
	   	Mb = midPtRK(model,bSetModel=0,b,dr,&R);	
		fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}

	// (CR) Debug
	fprintf(stderr,"Root bracketed.\n");

    /*
    ** Root bracketed by (a,b).
    */
    while (Mb-Ma > 1e-10*Mc) {
	c = 0.5*(a + b);
        Mc = midPtRK(model,bSetModel=0,c,dr,&R);	
	if (Mc < M) {
	    a = c;
	    Ma = Mc;
	    }
	else {
	    b = c;
	    Mb = Mc;
	    }
//	fprintf(stderr,"c:%.10g Mc:%.10g R:%.10g\n",c/model->tillMat[0]->rho0,Mc,R);
	}
    /*
    ** Solve it once more setting up the lookup table.
    */
    fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->tillMat->cv,model->uc);
    Mc = midPtRK(model,bSetModel=1,c,dr,&R);
    model->R = R;
    return c;
    }


double MLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->M[i]*(1.0 + log(r-model->r[i]+1)));
    x = r/model->dr;
    xi = floor(x);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
	fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
	return(model->M[i]*(1.0-x) + model->M[i+1]*x);
	}
    if (i == model->nTable-2) {
	dr = model->r[i+1] - model->r[i];
	x = r/dr;
	xi = floor(x);
	x -= xi;
	return(model->M[i]*(1.0-x) + model->M[i+1]*x);
	}
    else {
	i = model->nTable - 1;
	return(model->M[i]*(1.0 + log(r-model->r[i]+1)));
	}
    }

double rhoLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->rho[i]*exp(-(r-model->r[i])));
    x = r/model->dr;
    xi = floor(x);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
	fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
	return(model->rho[i]*(1.0-x) + model->rho[i+1]*x);
	}
    if (i == model->nTable-2) {
	dr = model->r[i+1] - model->r[i];
	x = r/dr;
	xi = floor(x);
	x -= xi;
	return(model->rho[i]*(1.0-x) + model->rho[i+1]*x);
	}
    else {
	i = model->nTable - 1;
	return(model->rho[i]*exp(-(r-model->r[i])));
	}
    }

double uLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->u[i]*exp(-(r-model->r[i])));
    x = r/model->dr;
    xi = floor(x);
    assert(xi >= 0.0);
    x -= xi;
    i = (int)xi;
    if (i < 0) {
	fprintf(stderr,"ERROR r:%.14g x:%.14g xi:%.14g i:%d\n",r,x,xi,i);
	}
    assert(i >= 0);
    if (i < model->nTable-2) {
	return(model->u[i]*(1.0-x) + model->u[i+1]*x);
	}
    if (i == model->nTable-2) {
	dr = model->r[i+1] - model->r[i];
	x = r/dr;
	xi = floor(x);
	x -= xi;
	return(model->u[i]*(1.0-x) + model->u[i+1]*x);
	}
    else {
	i = model->nTable - 1;
	return(model->u[i]*exp(-(r-model->r[i])));
	}
    }

double Fzero(MODEL *model,int bIcosa,double r,double ri,double m,int ns) {
    long npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
    return(MLookup(model,r)-MLookup(model,ri)-npix*m);
    }


double rShell(MODEL *model,double m,double ri) {
    double a = ri;
    double b = 1.0;
    double c,Mc,Ma;

    Ma = MLookup(model,a);
    Mc = 0.0;
    while (m > (MLookup(model,b)-Ma)) b *= 2.0;
    while (fabs(m-(Mc-Ma)) > 1e-10*m) {
	c = 0.5*(a + b);
	Mc = MLookup(model,c);
//	fprintf(stderr,"c:%.7g Mc:%.7g\n",c,Mc);
	if (m > (Mc-Ma)) a = c;
	else b = c;
	}
    return c;
    }


double rShell2(MODEL *model,int bIcosa,double m,double ri,int ns) {
    double a = ri;
    double b = 1.0;
    double c;
    double z;

    z = 1.0;
    while (Fzero(model,bIcosa,b,ri,m,ns) < 0) b *= 2.0;
    while ((b-a)/c > 1e-10) {
	c = 0.5*(a + b);
	z = Fzero(model,bIcosa,c,ri,m,ns);
	if (z < 0) a = c;
	else b = c;
//	printf("c:%.14g M(c):%.14g\n",c,MLookup(model,c));
	}
    return c;
    }


void main(int argc, char **argv) {
    const int bCentral = 1;
    int bIcosa = 0;
    const int bRandomRotate = 1;
    TCTX out;
    long ns,npix,ipix,na,nb;
    struct gas_particle gp;
    double r[3];
    double rhoCenter;
    double ri,ro,rs,roa,rob,ros,rta,rtb,rts;
    double m,mTot,l1,l2,nsf;
    double x,y,ang1,ang2,ang3;
    int j,iShell,nDesired,nReached,nLast;
    double theta,phi;
    float rr[3];
    ICOSA *ctx;
    int nShell,nMaxShell;
    double *rsShell;
    long *nsShell;
    long *isShell;
    double *xyz;
	// Model
    MODEL *model;
    double ucore;
	int iMat;
    int iter;
    int nSmooth;
    double d,rho,u,eta,w0,dPdrho;
    FILE *fpi,*fpo;
    int iRet,i;
	/*
	** These variables are used to find the optimal softening.
	*/
	double l1max = 0.0;
	double l2max = 0.0;

    if (argc != 5) {
	fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <ucore> <iMat> >myball.std\n");
	exit(1);
	}
    nDesired = atoi(argv[1]);
    mTot = atof(argv[2]);
    ucore = atof(argv[3]);
    iMat = atoi(argv[4]);

    model = modelInit(ucore, iMat);
    rhoCenter = modelSolve(model,mTot);

    m = mTot/nDesired;   /* a first guess at the particle mass */
    /*
    ** Initialize icosahedral parameters.
    */
    ctx = icosaInit();

#if (0)
    /*
    ** Test the icosahedral grid.
    */
    ns = 4;
    npix = 40*ns*(ns-1) + 12;
    for (ipix=0;ipix<nipix;++ipix) {
	icosaPix2Vec(ctx,ipix,ns,r);
	printf("%d %g %g %g\n",i,r[0],r[1],r[2]);
	}
#endif
    /*
    ** Initialize the array of shell radii and resolutions.
    */
    nMaxShell = 1000;
    nShell = 0;
    rsShell = malloc(nMaxShell*sizeof(double));
    assert(rsShell != NULL);
    nsShell = malloc(nMaxShell*sizeof(long));
    assert(nsShell != NULL);
    isShell = malloc(nMaxShell*sizeof(long));
    assert(nsShell != NULL);
    /*
    ** Phase 1: We first determine the number of particles per shell such that
    ** the ratio of radial to tangential lengths of their volumes is as close
    ** to 1 as possible. This fixes the total number of particles in the model.
    */
    fprintf(stderr,"PHASE 1: ODE approach\n");
    for (iter=0;iter<2;++iter) {
	nReached = 0;
 	if (bCentral) {
	    ro = rShell(model,m,0.0);
	    }
	else {
	    ro = 0.0;
	    }
	iShell = 0;
	for (;;) {
	    ri = ro;

	    na = 1;
	    roa = rShell2(model,bIcosa,m,ri,na);
	    npix = (bIcosa)?(40*na*(na-1)+12):(12*na*na);
	    l1 = roa-ri;
	    l2 = sqrt(M_PI/npix)*(roa+ri);
	    rta = l1/l2;
	    nb = 16;
	    do {
		nb *= 2;
		rob = rShell2(model,bIcosa,m,ri,nb);
		npix = (bIcosa)?(40*nb*(nb-1)+12):(12*nb*nb);
		l1 = rob-ri;
		l2 = sqrt(M_PI/npix)*(rob+ri);
		rtb = l1/l2;
		} while (rtb < 1.0);
	    while (nb - na > 1) {
		ns = (na+nb)/2;
		ros = rShell2(model,bIcosa,m,ri,ns);
		npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
		l1 = ros-ri;
		l2 = sqrt(M_PI/npix)*(ros+ri);
		rts = l1/l2;
/*		fprintf(stderr,"ns:%d rts:%g\n",ns,rts); */
		if (rts < 1.0) {
		    na = ns;
		    roa = ros;
		    rta = rts;
		    }
		else {
		    nb = ns;
		    rob = ros;
		    rtb = rts;
		    }
		}
/*
	    if (iShell >= 5) {
		fprintf(stderr,"na:%d rta:%g\n",na,rta);
		fprintf(stderr,"nb:%d rtb:%g\n",nb,rtb);
		}
*/
	    /*
	    ** if the two possible ratios differ by less that 1% then we favour
	    ** the higher resolution spherical grid (nb).
	    */
	    if (1/rta+0.01 < rtb) {
		ro = roa;
		ns = na;
		rts = rta;
		}
	    else {
		ro = rob;
		ns = nb;
		rts = rtb;
		}
	    npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	    if (iShell == nMaxShell) {
		nMaxShell *= 2;
		rsShell = realloc(rsShell,nMaxShell*sizeof(double));
		assert(rsShell != NULL);
		nsShell = realloc(nsShell,nMaxShell*sizeof(long));
		assert(nsShell != NULL);
		}
	    nsShell[iShell] = ns;
/*	fprintf(stderr,"nReached:%d npix:%d\n",nReached,npix);*/
	    if ((nReached + npix) < nDesired) {
		nReached += npix;
		fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g\n",iShell,ns,rts);
		++iShell;
		}
	    else {
		nShell = iShell;
		break;
		}
	    }  /* end of iShell loop */
	ns = nsShell[iShell-1];
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);	
	if (nDesired - nReached > npix/2) {
	    nReached += npix;
	    nsShell[nShell] = ns;
	    fprintf(stderr,"iShell:%d ns:%d radial/tangential:??? (added)\n",nShell,ns);
	    nShell++;
	    }
	fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
	    nReached,m,mTot/nReached);
	m = mTot/nReached;
	nDesired = nReached+1;
	}
    /*
    ** Phase 2: With the numbers of particles in each of the shells now fixed
    ** we continue by recalculating the radii of the shells based on the updated
    ** particle mass.
    */
    fprintf(stderr,"PHASE 2\n");
    if (bCentral) {
	ro = rShell(model,m,0.0);
	}
    else {
	ro = 0.0;
	}
    for (iShell=0;iShell<nShell;++iShell) {

	ri = ro;
	ns = nsShell[iShell];
	ro = rShell2(model,bIcosa,m,ri,ns);
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	l1 = ro-ri;
	l2 = sqrt(M_PI/npix)*(ro+ri);
	rts = l1/l2;
	/*
	** Here we calculate the optimal softening for gravity using the
	** largest bin size. l1: radial size, l2: tangental size
	*/
	if (l1 > l1max) {
			l1max = l1;
	}

	if (l2 > l2max) {
			l2max = l2;
	}

	if (bIcosa) {
	    d = 1.0;
	    if (iShell == 0) d = 2.5;
	    if (iShell == 1) d = 1.0;
	    if (iShell == 2) d = 3.0;
	    }
	else {
	    d = 1.9;
	    if (iShell == 0) d = 2.0;
	    if (iShell == 1) d = 1.5;
	    }

	rs = rsShell[iShell] = pow(0.5*(pow(ri,d) + pow(ro,d)),1/d);

	rho = rhoLookup(model,rs);

    u = uLookup(model,rs); /* We also have to look up u from a table */

	eta = rho/model->tillMat->rho0;
	    /* This was the old code using a constant internal energy uFixed.
	    w0 = model->uFixed/(model->par.u0*eta*eta) + 1.0;
		dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*model->uFixed + 
		(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

		fprintf(stderr,"iShell:%d r:%g M:%g rho:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,model->uFixed));
        */	
    	w0 = u/(model->tillMat->u0*eta*eta) + 1.0;
        dPdrho = (model->tillMat->a + (model->tillMat->b/w0)*(3 - 2/w0))*u + 
        (model->tillMat->A + 2*model->tillMat->B*(eta - 1))/model->tillMat->rho0;

//        fprintf(stderr,"iShell:%d r:%g M:%g rho:%g u:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,u,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,u));
        }
    /*
    ** Now generate the coordinates of all the particles in each shell as they are on the unit
    ** sphere. This simplifies the later adjusting of the radii of the shells (unit sphere 
    ** coordinates are independent of this.
    */

    /*
    ** Phase 3: With the masses and numbers of the particles fixed and the 
    ** as well as their angular positions, we now attempt to use an SPH density
    ** estimate to calculate the mean radial pressure force on the shell 
    ** and converge on the radius of the shell such that it is balanced by
    ** gravity on the shell (analytic).
    */
    if (bCentral) {
	/*
	** First adjust the density of the central particle.
	*/
	}

    /*
    ** Now output the particles.
    */
    TipsyInitialize(&out,0,NULL);

    for (j=0;j<3;++j) gp.vel[j] = 0.0;
    gp.phi = 0.0;
	/*
	** As a first guess we use max(l1max,l2max) for the softening.
	*/
	fprintf(stderr,"l1max: %g l2max: %g epsilon: %g\n",l1max,l2max,MAX(l1max,l2max));
    gp.hsmooth = MAX(l1max,l2max);  /* is actually eps for gasoline */
    gp.mass = m;

	fprintf(stderr,"hsmooth=%g\n",gp.hsmooth);

	gp.hsmooth=0.001;
    
	//gp.temp = model->uFixed;   /* Christian's version of gasoline uses thermal energy instead of temperature as input! */
    nLast = nReached;
    nReached = 0;
    if (bCentral) {
		for (j=0;j<3;++j) gp.pos[j] = 0.0;
		gp.temp = uLookup(model, 0);
		// Dont forget to set the material for the central particle
		gp.metals = iMat;
		TipsyAddGas(out,&gp);
	}
    for (iShell=0;iShell<nShell;++iShell) {
	rs = rsShell[iShell];
	ns = nsShell[iShell];
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	nReached += npix;
	ang1 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
	ang2 = 2.0*M_PI*rand()/(RAND_MAX+1.0);
	ang3 = 2.0*M_PI*rand()/(RAND_MAX+1.0);

	/* Experiment with grav. softening. */
	if (iShell < nShell-1)
	{
		// Inner shell
		gp.hsmooth = 0.001;
	} else {
		// Most outer shell
		gp.hsmooth = 0.01;
	}

	for (ipix = 0;ipix < npix;++ipix) {
	    if (bIcosa) icosaPix2Vec(ctx,ipix,ns,r);
	    else pix2vec_ring(ns,ipix,r);
	    if (bRandomRotate) {
		y = r[1]*cos(ang1) - r[2]*sin(ang1);
		r[2] = r[1]*sin(ang1) + r[2]*cos(ang1);
		
		x = r[0]*cos(ang2) - r[2]*sin(ang2);
		r[2] = r[0]*sin(ang2) + r[2]*cos(ang2);
		r[0] = x;

		r[1] = y*cos(ang3) - r[2]*sin(ang3);
		r[2] = y*sin(ang3) + r[2]*cos(ang3);
		}

	    for (j=0;j<3;++j) gp.pos[j] = rs*r[j];
	    
//	    rho = rhoLookup(model,rs);
	    gp.temp = uLookup(model,rs);
		// Save Material
		gp.metals = iMat;
	    TipsyAddGas(out,&gp);		
	    }
	}
    fprintf(stderr,"Writing %d particles. Model R:%g Last Shell r:%g\n",nReached,model->R,rsShell[nShell-1]);
    TipsyWriteAll(out,0.0,"ballic.std");
    TipsyFinish(out);    
#if 0
    /* Smooth with gather doesn work with my version! */
    system("sleep 1; smooth -s 128 density <ballic.std; sleep 1");

    fpi = fopen("smooth.den","r");
    assert(fpi != NULL);
    fpo = fopen("ballic.den","w");
    assert(fpo != NULL);
    iRet = fscanf(fpi,"%d",&nReached);
    assert(iRet == 1);
    if (bCentral) {
	iRet = fscanf(fpi,"%lf",&rho);
	assert(iRet == 1);
	fprintf(fpo,"0.0 %g\n",rho);
	nReached -= 1;
	}
    for (iShell=0;iShell<nShell;++iShell) {
	rs = rsShell[iShell];
	ns = nsShell[iShell];
	npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
	nReached -= npix;
	for (ipix = 0;ipix < npix;++ipix) {
	    iRet = fscanf(fpi,"%lf",&rho);
	    assert(iRet == 1);
	    fprintf(fpo,"%g %g\n",rs,rho);
	    }
	}
    fclose(fpi);
    fclose(fpo);
    assert(nReached == 0);
#endif
    }

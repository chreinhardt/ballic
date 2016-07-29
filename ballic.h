#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tipsy.h"

#ifndef BALLIC
#define BALLIC

typedef struct icosa_struct {
    float R[180];
    float v[36];
    } ICOSA;

ICOSA *icosaInit(void);

void icosaPix2Vec(ICOSA *ctx,int i,int resolution,double *vec);

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
void pix2vec_ring( long nside, long ipix, double *vec);

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
    struct params {
	double a;
	double b;
	double A;
	double B;
	double rho0;
	double u0;
	double us;
	double us2;
	double alpha;
	double beta;
    double cv; /* specific heat capacity at constant volume (assumed to be constant) */ 
	} par;
    /*
    ** Some unit conversion factors.
    */
    double dKpcUnit;
    double dMsolUnit;
    double dGasConst;
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;
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


MODEL *modelInit(double ucore) {
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */

    MODEL *model;
    
    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
    /*
    ** Convert kboltz/mhydrogen to system units, assuming that
    ** G == 1.
    */
    model->dGasConst = model->dKpcUnit*KPCCM*KBOLTZ
	/MHYDR/GCGS/model->dMsolUnit/MSOLG;
    /* code energy per unit mass --> erg per g */
    model->dErgPerGmUnit = GCGS*model->dMsolUnit*MSOLG/(model->dKpcUnit*KPCCM);
    /* code density --> g per cc */
    model->dGmPerCcUnit = (model->dMsolUnit*MSOLG)/pow(model->dKpcUnit*KPCCM,3.0);
    /* code time --> seconds */
    model->dSecUnit = sqrt(1/(model->dGmPerCcUnit*GCGS));
    /*
    ** Just set granite at the moment.
    */
    model->par.a = 0.5;
    model->par.b = 1.3;
    model->par.u0 = 1.6e11; /* in ergs/g */
    model->par.rho0 = 2.7; /* g/cc */
    model->par.A = 1.8e11; /* ergs/cc */
    model->par.B = 1.8e11; /* ergs/cc */
    model->par.us = 3.5e10; /* ergs/g */
    model->par.us2 = 1.8e11; /* ergs/g */
    model->par.alpha = 5.0;
    model->par.beta = 5.0;
    model->par.cv = 0.79e7; /* ergs/g K (or 790 J/kg K) */ 
    /*
    ** Convert energies and densities to code units!
    */
    model->par.u0 /= model->dErgPerGmUnit;
    model->par.us /= model->dErgPerGmUnit;
    model->par.us2 /= model->dErgPerGmUnit;
    model->par.rho0 /= model->dGmPerCcUnit;
    model->par.A /= (model->dGmPerCcUnit*model->dErgPerGmUnit);
    model->par.B /= (model->dGmPerCcUnit*model->dErgPerGmUnit);
    /* convert cv to code units here !! */
    model->par.cv /= model->dErgPerGmUnit;

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

/*
** Currently for condensed states only! Changed this 8.2.13 to include a polytropic temperature profile.
*/
double drhodr(MODEL *model,double r,double rho,double M,double u) {
    double eta = rho/model->par.rho0;
    double w0 = u/(model->par.u0*eta*eta) + 1.0;
    double dPdrho,dPdu;

//    fprintf(stderr,"u in system units: %.14g\n",u);

    dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*u + 
	(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

    dPdu = (model->par.a + model->par.b/(w0*w0))*rho;
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
** M_initial = 0 at r = 0. This function returns the mass when rho == model->par.rho0.
*/
double midPtRK(MODEL *model,int bSetModel,double rho,double h,double *pR) {
    FILE *fp;
    double M = 0.0;
    double r = 0.0;
	double u = model->par.uc;
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
    while (rho > fact*model->par.rho0) {
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
    x = (fact*model->par.rho0 - rho)/k2rho;
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
    R = cbrt(3.0*M/(4.0*M_PI*model->par.rho0));
    dr = R/nStepsMax;
    a = 1.01*model->par.rho0; /* starts with 1% larger central density */
    Ma = midPtRK(model,bSetModel=0,a,dr,&R);
    fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
    b = a;
    Mb = 0.5*M;
    while (Ma > M) {
	b = a;
	Mb = Ma;
	a = 0.5*(model->par.rho0 + a);
	Ma = midPtRK(model,bSetModel=0,a,dr,&R);
	}
    while (Mb < M) {
	b = 2.0*b;
   	Mb = midPtRK(model,bSetModel=0,b,dr,&R);	
	fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}
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
//	fprintf(stderr,"c:%.10g Mc:%.10g R:%.10g\n",c/model->par.rho0,Mc,R);
	}
    /*
    ** Solve it once more setting up the lookup table.
    */
    fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->par.cv,model->uc);
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
    MODEL *model;
    double ucore;
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

    if (argc != 4) {
	fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <ucore> >myball.std\n");
	exit(1);
	}
    nDesired = atoi(argv[1]);
    mTot = atof(argv[2]);
    ucore = atof(argv[3]);

    model = modelInit(ucore);
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

    u = uFunc(model,rho); /* We also have to look up u from a table */

	eta = rho/model->par.rho0;
	    /* This was the old code using a constant internal energy uFixed.
	    w0 = model->uFixed/(model->par.u0*eta*eta) + 1.0;
		dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*model->uFixed + 
		(model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

		fprintf(stderr,"iShell:%d r:%g M:%g rho:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,model->uFixed));
        */	
    	w0 = u/(model->par.u0*eta*eta) + 1.0;
        dPdrho = (model->par.a + (model->par.b/w0)*(3 - 2/w0))*u + 
        (model->par.A + 2*model->par.B*(eta - 1))/model->par.rho0;

        fprintf(stderr,"iShell:%d r:%g M:%g rho:%g u:%g ns:%d radial/tangential:%g dr:%g <? Jeans:%g Gamma:%g\n",iShell,rs,MLookup(model,rs),rho,u,ns,rts,ro-ri,sqrt(dPdrho/rho),Gamma(model,rho,u));
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
    //gp.temp = model->uFixed;   /* Christian's version of gasoline uses thermal energy instead of temperature as input! */
    nLast = nReached;
    nReached = 0;
    if (bCentral) {
	for (j=0;j<3;++j) gp.pos[j] = 0.0;
	gp.temp = uFunc(model, rhoLookup(model,0));   /* The new version uses a polytropic thermal profile! */
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
	    
	    rho = rhoLookup(model,rs);
	    gp.temp = uFunc(model,rho);
	    TipsyAddGas(out,&gp);		
	    }
	}
    fprintf(stderr,"Writing %d particles. Model R:%g Last Shell r:%g\n",nReached,model->R,rsShell[nShell-1]);
    TipsyWriteAll(out,0.0,"ballic.std");
    TipsyFinish(out);    
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
    }

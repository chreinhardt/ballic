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
#include "tipsydefs/tipsy.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

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
    int nMat;
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

    double *r;
    double *rho;
    double *M;
    double *u;
    double *P;
    double *T;

    int *mat;
    double dr;
    double R;
} MODEL;


MODEL *modelInit() {
    int i;
    /* Initialize the model */
    MODEL *model;

    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;

    /* Calculate the conversion factors. */
    model->nTableMax = 10000; 
    model->r = malloc(model->nTableMax*sizeof(double));
    assert(model->r != NULL);
    model->rho = malloc(model->nTableMax*sizeof(double));
    assert(model->rho != NULL);
    model->M = malloc(model->nTableMax*sizeof(double));
    assert(model->M != NULL);
    model->u = malloc(model->nTableMax*sizeof(double));
    assert(model->u != NULL);
    model->P = malloc(model->nTableMax*sizeof(double));
    assert(model->P != NULL);
    model->T = malloc(model->nTableMax*sizeof(double));
    assert(model->T != NULL);

    model->mat = malloc(model->nTableMax*sizeof(double));
    assert(model->mat != NULL);
    model->dr =  0.0;
    model->nTable = 0;

    return(model);
}

/*
 ** Read a text file containing an equilibrium model.
 */
double modelRead(MODEL *model, char *file) {
    FILE *fp;
    int iRet, i;

    assert(model != NULL);

    /* Open the file. */
    fp = fopen(file,"r");
    assert(fp != NULL);

    /* Read the model. */
    for (i=0; i<model->nTableMax;i++)
    {
        /* Currently only single material models implemented so iMat is read but ignored. */
        iRet = fscanf(fp, "%lf %lf %lf %lf %d %lf %lf", &model->r[i], &model->rho[i],
                &model->M[i], &model->u[i],  &model->mat[i], &model->P[i], &model->T[i]);

        if (iRet <= 0 && feof(fp)) break;

        assert(iRet > 0);

        assert(model->r[i] >= 0.0);
        assert(model->rho[i] >= 0.0);
        assert(model->M[i] >= 0.0);
        assert(model->u[i] >= 0.0);
        assert(model->mat[i] >= 0.0);

        model->nTable = i;
        //printf("%g  %g  %g  %g  %i\n",model->r[i],model->M[i],model->rho[i],model->u[i],model->mat[i]);
    }

    model->nTable += 1;

    fprintf(stderr,"nTable=%i\n",model->nTable);

}
#if 0
double modelRead(MODEL *model, char *file) {
    FILE *fp;
    int iRet, i;

    assert(model != NULL);

    /* Open the file. */
    fp = fopen(file,"r");
    assert(fp != NULL);

    /* Read the model. */
    for (i=0; i<model->nTableMax;i++)
    {

        //		fprintf(stderr,"i=%i\n",i);
        // R  M(R)  rho  u  iMat
        iRet = fscanf(fp, "%lf %lf %lf %lf %d", &model->r[i],&model->M[i],&model->rho[i],&model->u[i],&model->mat[i]);

        if (iRet <= 0 && feof(fp)) break;

        assert(iRet > 0);

        assert(model->r[i] >= 0.0);
        assert(model->M[i] >= 0.0);
        assert(model->rho[i] >= 0.0);
        assert(model->u[i] >= 0.0);
        assert(model->mat[i] >= 0.0);

        model->nTable = i;
        //printf("%g  %g  %g  %g  %i\n",model->r[i],model->M[i],model->rho[i],model->u[i],model->mat[i]);
    }

    model->nTable += 1;

    fprintf(stderr,"nTable=%i\n",model->nTable);

}
#endif
/*
 ** Do bisection to find the values r_i and r_i+1 that
 ** bracket r in the lookup table. The function returns
 ** r_i.
 */
int rLookup(MODEL *model,double r) {
    int iLower,iUpper,i;

    iLower = 0;
    iUpper = model->nTable-1;

    fprintf(stderr,"rLookup: r=%g\n", r);
    /* Make sure that r is in the lookup table. */
    assert(r >= model->r[iLower]);
    assert(r <= model->r[iUpper]);

    assert(model->r[iLower] < model->r[iUpper]);

    /* Do bisection. */
    while (iUpper-iLower > 1)
    {
        i = (iUpper + iLower) >> 1;

        if (r >= model->r[i])
        {
            iLower = i;
        } else {
            iUpper = i;
        }
    }

    if (r == model->r[0]) return 0;
    if (r == model->r[model->nTable-1]) return 0;

    return iLower;
}

/*
 ** Find the mass enclosed for a given r.
 */
double MLookup(MODEL *model,double r) {
    int i;
    double A;

    /* If r is outside of the lookup table, return something. */
    i = model->nTable-1;
    if (r >= model->r[i]) return(model->M[i]*(1.0 + log(r-model->r[i]+1)));

    i = rLookup(model,r);

    /* Do a linear interpolation. */
    A = (model->r[i+1]-r)/(model->r[i+1]-model->r[i]);
    return (A*model->M[i]+(1.0-A)*model->M[i+1]);
}

/*
 ** Find the density for a given r.
 */
double rhoLookup(MODEL *model,double r) {
    int i;
    double A;

    /* If r is outside of the lookup table, return something. */
    i = model->nTable-1;
    if (r >= model->r[i]) return(model->rho[i]*exp(-(r-model->r[i])));

    i = rLookup(model,r);

    /* Do a linear interpolation. */
    A = (model->r[i+1]-r)/(model->r[i+1]-model->r[i]);
    return (A*model->rho[i]+(1.0-A)*model->rho[i+1]);
}

/*
 ** Find the internal energy for a given r.
 */
double uLookup(MODEL *model,double r) {
    int i;
    double A;

    /* If r is outside of the lookup table, return something. */
    i = model->nTable-1;
    if (r >= model->r[i]) return(model->u[i]*exp(-(r-model->r[i])));

    i = rLookup(model,r);

    /* Do a linear interpolation. */
    A = (model->r[i+1]-r)/(model->r[i+1]-model->r[i]);
    return (A*model->u[i]+(1.0-A)*model->u[i+1]);
}

/*
 ** Find the temperature for a given r.
 */
double TLookup(MODEL *model,double r) {
    int i;
    double A;

    /* If r is outside of the lookup table, return something. */
    i = model->nTable-1;
    if (r >= model->r[i]) return(model->T[i]*exp(-(r-model->r[i])));

    i = rLookup(model,r);

    /* Do a linear interpolation. */
    A = (model->r[i+1]-r)/(model->r[i+1]-model->r[i]);
    return (A*model->T[i]+(1.0-A)*model->T[i+1]);
}

/* Needs work. */
int matLookup(MODEL *model,double r) {
    return 2;
#if 0
    double x,xi,dr;
    int i;

    assert(0);

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
#endif
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
        //		fprintf(stderr,"c:%.7g Mc:%.7g\n",c,Mc);
        if (m > (Mc-Ma)) a = c;
        else b = c;
    }
    printf("rShell: r=%g ri=%g m=%g\n", c,ri,m);
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
        //		printf("c:%.14g M(c):%.14g\n",c,MLookup(model,c));
    }
    printf("rShell2: r=%g ri=%g\n", c,ri);
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

    if (argc != 3) {
        fprintf(stderr,"Usage: ballic <nDesired> <model.txt> >myball.std\n");
        exit(1);
    }
    nDesired = atoi(argv[1]);

    fprintf(stderr,"Initializing model...\n");
    model = modelInit();

    fprintf(stderr,"Reading...\n");
    modelRead(model,argv[2]);

    fprintf(stderr,"Done.\n");
#if 0
    rs = 0.1;
    rs = 1.14236;
    rs = 1.14246;
    fprintf(stderr,"M(r=%g)=%g\n",rs, MLookup(model,rs));

    rs=0.0;
    while (rs < 2.0)
    {
        printf("%g  %g  %g  %g  %i\n",rs,MLookup(model,rs),rhoLookup(model,rs),uLookup(model,rs),iRet);
        rs += 0.01;
    }

    exit(1);
#endif
    /*
     ** This sould work, but in general a rho lookup might be problematic as
     ** the values can be discontinous.
     */
    rhoCenter = rhoLookup(model,0.0);

    mTot = model->M[model->nTable-1];
    m = mTot/nDesired;   /* a first guess at the particle mass */
    fprintf(stderr,"m=%g mTot=%g\n",m,mTot);
    assert(m>0.0);

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
        } else {
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
                /*				fprintf(stderr,"ns:%d rts:%g\n",ns,rts); */
                if (rts < 1.0) {
                    na = ns;
                    roa = ros;
                    rta = rts;
                } else {
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

            fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g ri:%g ro:%g rta:%g\n",iShell,ns,rts,ri,ro,rta);
            /*
             ** if the two possible ratios differ by less that 1% then we favour
             ** the higher resolution spherical grid (nb).
             */
            if (1/rta+0.01 < rtb) {
                ro = roa;
                ns = na;
                rts = rta;
            } else {
                ro = rob;
                ns = nb;
                rts = rtb;
            }

            fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g ri:%g ro:%g rta:%g\n",iShell,ns,rts,ri,ro,rta);
            npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
            if (iShell == nMaxShell) {
                nMaxShell *= 2;
                rsShell = realloc(rsShell,nMaxShell*sizeof(double));
                assert(rsShell != NULL);
                nsShell = realloc(nsShell,nMaxShell*sizeof(long));
                assert(nsShell != NULL);
            }
            nsShell[iShell] = ns;
            /*			fprintf(stderr,"nReached:%d npix:%d\n",nReached,npix);*/
            if ((nReached + npix) < nDesired) {
                nReached += npix;
                // (CR) 23.10.2016
                fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g ri:%g ro:%g rta:%g\n",iShell,ns,rts,ri,ro,rta);
                fprintf(stderr,"iShell:%d ns:%d radial/tangential:%g\n",iShell,ns,rts);
                ++iShell;
            } else {
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
    } else {
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

        if (bWriteTemp) {
            u = TLookup(model,rs); /* pkdgrav3 needs temperature */
        } else {
            u = uLookup(model,rs); /* We also have to look up u from a table */
        }

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
    //gp.temp = model->uFixed;   /* Christian's version of gasoline uses thermal energy instead of temperature as input! */
    nLast = nReached;
    nReached = 0;
    if (bCentral) {
        for (j=0;j<3;++j) gp.pos[j] = 0.0;
        gp.temp = uLookup(model, 0);
        // Dont forget to set the material for the central particle
        gp.metals = matLookup(model,0);
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

            //	    rho = rhoLookup(model,rs);
            gp.temp = uLookup(model,rs);
            // Save Material
            gp.metals = matLookup(model,rs);
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

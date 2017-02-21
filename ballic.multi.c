#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "ballic.h"
#include "tipsy.h"
#include "tillotson/tillotson.h"

/* Functions for Icosahedron. */
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

MODEL *modelInit(double M,double ucore) {
	int i;
    /* Initialize the model */
	MODEL *model;
    
    model = malloc(sizeof(MODEL));
    assert(model != NULL);

    model->dKpcUnit = 2.06701e-13;
    model->dMsolUnit = 4.80438e-08;
	
	/* Hard coded for the moment */
	model->nLayer = 2;
	assert(model->nLayer <= TILL_N_MATERIAL_MAX);

	model->tillMat = malloc(model->nLayer*sizeof(TILLMATERIAL *));
	assert(model->tillMat != NULL);
	
	model->iLayer = malloc(model->nLayer*sizeof(int));
	assert(model->iLayer != NULL);
	
	model->MLayer = malloc(model->nLayer*sizeof(double));
	assert(model->MLayer != NULL);

	/* Hard coded too */
	model->iLayer[0] = IRON;
//	model->iLayer[1] = GRANITE;
	model->iLayer[1] = BASALT;
	// It might be better so save M in model and use only mass fractions in MLayer
	model->MLayer[0] = 0.3*M;
	model->MLayer[1] = 0.7*M;

#if 0
	/* Debugging ice. */
	model->iLayer[0] = GRANITE;
	model->iLayer[1] = ICE;
	model->MLayer[0] = 0.25*M;
	model->MLayer[1] = 0.75*M;
#endif

#if 0
	/* Single component model. */
	model->iLayer[0] = GRANITE;
	model->MLayer[0] = 1.0*M;
#endif
	fprintf(stderr,"Initializing model:\n");
	fprintf(stderr,"Mtot=%g ucore=%g\n",M,ucore);

	fprintf(stderr,"iLayer[ ");
	for (i=0; i<model->nLayer; i++)
		fprintf(stderr,"%i ",model->iLayer[i]);
	fprintf(stderr,"]\n");

	fprintf(stderr,"MLayer[ ");
	for (i=0; i<model->nLayer; i++)
		fprintf(stderr,"%g ",model->MLayer[i]);
	fprintf(stderr,"]\n");

	for (i=0; i<model->nLayer; i++)
	{
		/*
		** Initialize one material.
		*/
		model->tillMat[i] = tillInitMaterial(model->iLayer[i], model->dKpcUnit, model->dMsolUnit, 100, 100, 50.0, 50.0, 1);

		// Debug information

		fprintf(stderr,"\n");	
		fprintf(stderr,"Material: %i\n",model->iLayer[i]);	
		fprintf(stderr,"a: %g\n", model->tillMat[i]->a);
		fprintf(stderr,"b: %g\n", model->tillMat[i]->b);
		fprintf(stderr,"A: %g\n", model->tillMat[i]->A);
		fprintf(stderr,"B: %g\n", model->tillMat[i]->B);
		fprintf(stderr,"rho0: %g\n", model->tillMat[i]->rho0);
		fprintf(stderr,"u0: %g\n", model->tillMat[i]->u0);
		fprintf(stderr,"us: %g\n", model->tillMat[i]->us);
		fprintf(stderr,"us2: %g\n", model->tillMat[i]->us2);
		fprintf(stderr,"alpha: %g\n", model->tillMat[i]->alpha);
		fprintf(stderr,"beta: %g\n", model->tillMat[i]->beta);
   		fprintf(stderr,"cv: %g\n", model->tillMat[i]->cv);
		fprintf(stderr,"\n");

		/* Generate the look up table needed for tillColdULookup(). */
		tillInitLookup(model->tillMat[i]);
	}


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
    model->mat = malloc(model->nTableMax*sizeof(int));
    assert(model->mat != NULL);
    model->dr =  0.0;
    model->nTable = 0;
    
	return(model);
    }

/*
** dudrho depends on the internal energy profile that we choose!
*/
double dudrho(MODEL *model,int iLayer,double rho,double u) {
	/*
	** We assume an isentropic internal energy profile!
	*/
	return(tillPressure(model->tillMat[iLayer],rho,u)/(rho*rho));
}

/*
** Calculate dudrho to solve for the equilibrium model.
*/
double drhodr(MODEL *model, int iLayer, double r,double rho,double M,double u) {
    double dPdrho,dPdu;

	dPdrho=tilldPdrho(model->tillMat[iLayer], rho, u); // dP/drho at u=const.
	dPdu = tilldPdu(model->tillMat[iLayer], rho, u);; // dP/du at rho=const.

	// (CR) Debug
//	fprintf(stderr,"dPdrho:%g, dPdu:%g,dudrho:%g\n",dPdrho,dPdu,dudrho(model,iLayer,rho,u));

	/*
	** drho/dr = -G*M*rho/(dPdrho+dPdu*dudrho)
	*/
	assert(r >= 0.0);
	if (r > 0.0) {
		// We assume G=1
		return(-M*rho/(r*r*(dPdrho + dPdu*dudrho(model,iLayer,rho,u))));
	}
	else {
		return(0.0);
	}
}

double dudr(MODEL *model,int iLayer,double r,double rho,double M,double u) {
  return(dudrho(model,iLayer,rho,u)*drhodr(model,iLayer,r,rho,M,u));
}

/*
** This derivative is independent of the model and only involves geometry.
*/
double dMdr(double r,double rho) {
	assert(r >= 0.0);
	return(4.0*M_PI*r*r*rho);
}

/*
** Solve for rho2 and u2 using the b.c. P1=P2 and T1=T2.
*/
void modelSolveBC(MODEL *model, double *prho, double *pu, int iLayer1, int iLayer2) {
	double P1, T1, P2, T2;
    double rho1,u1,rho2,u2;
    double rhoa,ua,rhob,ub;

	assert(prho != NULL);
	assert(pu != NULL);
	
	// Make sure that (rho,u) is not below the cold curve (not implemented yet!)
	
	rho1 = *prho;
	u1 = *pu;

	P1 = tillPressure(model->tillMat[model->iLayer[iLayer1]], rho1, u1);
	T1 = tillTempRhoU(model->tillMat[model->iLayer[iLayer1]], rho1, u1);

	/*
	** We use rho1 as an upper limit for rho2 assuming that the denser component is in the inner shell.
	*/
	rhoa = rho1;
	ua = tillURhoTemp(model->tillMat[model->iLayer[iLayer2]],rhoa,T1);

	rhob = 0.0;
	ub = tillURhoTemp(model->tillMat[model->iLayer[iLayer2]],rhob,T1);

	fprintf(stderr,"\n");
	fprintf(stderr,"*****************************************************************\n");
	fprintf(stderr,"modelSolveBC:\n");
	fprintf(stderr,"iLayer1=%i (iMat=%i) iLayer2=%i (iMat=%i)\n",iLayer1,model->iLayer[iLayer1],iLayer2,model->iLayer[iLayer2]);
	fprintf(stderr,"rho1=%g u1=%g P1=%g T1=%g\n",rho1,u1,P1,T1);
	fprintf(stderr,"rhoa=%g ua=%g Pa=%g Ta=%g\n",rhoa,ua,tillPressure(model->tillMat[model->iLayer[iLayer2]], rhoa, ua),tillTempRhoU(model->tillMat[model->iLayer[iLayer2]], rhoa, ua));
   	fprintf(stderr,"rhob=%g ub=%g Pb=%g Tb=%g\n",rhob,ub,tillPressure(model->tillMat[model->iLayer[iLayer2]], rhob, ub),tillTempRhoU(model->tillMat[model->iLayer[iLayer2]], rhob, ub));
	
	tillSolveBC(model->tillMat[model->iLayer[iLayer1]],model->tillMat[model->iLayer[iLayer2]],rho1,u1,&rho2,&u2);
	fprintf(stderr,"modelSolveBC: rho1: %g, u1: %g, rho2:%g, u2:%g\n",rho1,u1,rho2,u2);
	fprintf(stderr,"*****************************************************************\n");
	fprintf(stderr,"\n");

	/*
	** Return values.
	*/
	*prho = rho2;
	*pu = u2; 
}

/*
** This function integrates the ODEs with b.c. rho_initial=rho1, u_initial=u1
** and M_initial=M1 until M=M2. The final values for rho and u are returned
** (so the original value will be overwritten!!). The parameter h sets the
** stepsize for the RK2 algorithm and for bSetModel=1 the results are saved
** in the lookup table (starting from index pIndex). If bLastLayer=1 then
** the algorithm enforces rho(r=R)=rho0.
*/
void modelSolveComponent(MODEL *model,int iLayer,int bSetModel,int bLastLayer,int *pIndex,double h,double *prho1,double *pu1,double *pM1,double M2,double *pR)
{
    FILE *fp;
	// Set inital values rho1, u1, M1
	double rho=*prho1;
	double u = *pu1;
    double M = *pM1;
	double r = *pR;

    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveComponent (inital values):\n");
	fprintf(stderr,"iLayer: %i (iMat: %i), Index: %i, h: %g\n",iLayer,model->iLayer[iLayer],*pIndex,h);
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho,u,M,r);
	fprintf(stderr,"bSetModel: %i, bLastLayer: %i\n",bSetModel,bLastLayer);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

	if (bSetModel) {
		model->rho[*pIndex] = rho;
		model->M[*pIndex] = M;
		model->u[*pIndex] = u;
		model->r[*pIndex] = r;
		model->mat[*pIndex] = model->iLayer[iLayer];
		++*pIndex;
	}

	if (bLastLayer != 1)
	{
		/* Integrate from M1 to M2 for the inner layers. */
		while (M < M2) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1rho = h*drhodr(model,iLayer,r,rho,M,u);
			k1M = h*dMdr(r,rho);
			k1u = h*dudr(model,iLayer,r,rho,M,u);

			k2rho = h*drhodr(model,iLayer,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
			k2u = h*dudr(model,iLayer,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;
		
			if (bSetModel) {
				model->rho[*pIndex] = rho;
				model->M[*pIndex] = M;
				model->u[*pIndex] = u;
				model->r[*pIndex] = r;
				model->mat[*pIndex] = model->iLayer[iLayer];
	//			fprintf(fp,"%g %g %g %g\n",r,rho,M,u);
				++*pIndex;
		    }
		}

		/*
		** Now do a linear interpolation to M == M2.
		*/
		x = (M2 - M)/k2M;
		fprintf(stderr,"M2=%g, M=%g, x=%g\n",M2,M,x);
		assert(x <= 0.0);
		r += h*x;
		M += k2M*x;
		rho += k2rho*x;
		u += k2u*x;

		if (bSetModel) {
			--*pIndex;
			model->M[*pIndex] = M;
			model->r[*pIndex] = r;
			model->rho[*pIndex] = rho;
			model->u[*pIndex] = u;
			model->mat[*pIndex] = model->iLayer[iLayer];
			++*pIndex;
		}
		/* Make sure that the material is in the condensed state */
		assert(rho >= model->tillMat[model->iLayer[iLayer]]->rho0);
	} else {
		// For the last layer we integrate until rho == rho0. */
		while (rho > fact*model->tillMat[model->iLayer[iLayer]]->rho0) {
			/*
			** Midpoint Runga-Kutta (2nd order).
			*/
			k1rho = h*drhodr(model,iLayer,r,rho,M,u);
			k1M = h*dMdr(r,rho);
			k1u = h*dudr(model,iLayer,r,rho,M,u);

			k2rho = h*drhodr(model,iLayer,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
			k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
			k2u = h*dudr(model,iLayer,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

			rho += k2rho;
			M += k2M;
			u += k2u;
			r += h;
	
//			fprintf(stderr,"r=%g rho=%g M=%g u=%g k1rho=%g k1M=%g k1u=%g k2rho=%g k2M=%g k2u=%g\n",
//							r, rho, M, u, k1rho, k1M, k1u, k2rho, k2M, k2u);	
			if (bSetModel) {
				model->M[*pIndex] = M;
				model->r[*pIndex] = r;
				model->rho[*pIndex] = rho;
				model->u[*pIndex] = u;
				model->mat[*pIndex] = model->iLayer[iLayer];
				++*pIndex;
			}
		}
		
		/*
		** Now do a linear interpolation to rho == fact*rho0.
		*/
		x = (fact*model->tillMat[model->iLayer[iLayer]]->rho0 - rho)/k2rho;
		fprintf(stderr,"iLayer=%i (iMat=%i), rho0=%g, rho=%g, x=%g\n",iLayer,model->iLayer[iLayer],model->tillMat[model->iLayer[iLayer]]->rho0,rho,x);
		fprintf(stderr,"rho0-rho=%g, k2rho0=%g\n",model->tillMat[model->iLayer[iLayer]]->rho0-rho,k2rho);
		assert(x <= 0.0);
		r += h*x;
		M += k2M*x;
		rho += k2rho*x;
		u += k2u*x;
		
		if (bSetModel) {
			--*pIndex;
			model->M[*pIndex] = M;
			model->r[*pIndex] = r;
			model->rho[*pIndex] = rho;
			model->u[*pIndex] = u;
			model->mat[*pIndex] = model->iLayer[iLayer];
			++*pIndex;
		}	
	}

	// Return values
	*prho1 = rho;
	*pu1 = u;
	*pM1 = M;
    *pR = r;
}

/*
** This function integrates the ODEs for a two component model with b.c.
** rho_initial=rho, u_initial=u until rho(r=R)=rho0. It returns the total
** mass of the model. The parameter h sets the stepsize for the RK2 algorithm
** and for bSetModel=1 the results are saved.
*/
double modelSolveTwoComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR)
{
//    FILE *fp;
	// Set inital values for M1 and R
    double M = 0.0;
	double r = 0.0;
	double Mc = model->MLayer[0];
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
	int i = 0;

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveTwoComponent (inital values):\n");
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho,u,M,r);
	fprintf(stderr,"bSetModel: %i, h: %g\n",bSetModel,h);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

	/*
	** Start with the core
	*/
	if (bSetModel) {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = model->iLayer[0];
		++i;
	}

	/*
	** First integrate the core until M == Mcore.
	*/
	while (M < Mc) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,0,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,0,r,rho,M,u);

		k2rho = h*drhodr(model,0,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,0,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;
		
		if (bSetModel) {
			model->rho[i] = rho;
			model->M[i] = M;
			model->u[i] = u;
			model->r[i] = r;
			model->mat[i] = model->iLayer[0];
			++i;
		}
	}

	/*
	** Now do a linear interpolation to M == Mcore.
	*/
	x = (Mc - M)/k2M;
	fprintf(stderr,"M2=%g, M=%g, x=%g\n",Mc,M,x);
	assert(x <= 0.0);
	r += h*x;
	M += k2M*x;
	rho += k2rho*x;
	u += k2u*x;
	fprintf(stderr,"After correction: M2=%g, M=%g, x=%g\n",Mc,M,x);

	if (bSetModel) {
		--i;
		model->M[i] = M;
		model->r[i] = r;
		model->rho[i] = rho;
		model->u[i] = u;
		model->mat[i] = model->iLayer[0];
		++i;
	}
	
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveTwoComponent (core/mantle boundary):\n");
	fprintf(stderr,"Core: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",model->iLayer[0],rho,u,M,r,tillPressure(model->tillMat[0],rho,u),tillTempRhoU(model->tillMat[0],rho,u));

	/*
	** Now calculate rho and u for the mantle using b.c. T=const and P=const.
	*/
	tillSolveBC(model->tillMat[0],model->tillMat[1],rho,u,&rho,&u);

	fprintf(stderr,"Mantle: iMat=%i rho=%g u=%g M=%g r=%g P=%g T=%g\n",model->iLayer[1],rho,u,M,r,tillPressure(model->tillMat[1],rho,u),tillTempRhoU(model->tillMat[1],rho,u));
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");
	
	/*
	** The integrate the mantle until rho(r=R)=rho0.
	*/
	while (rho > fact*model->tillMat[1]->rho0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,1,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,1,r,rho,M,u);

		k2rho = h*drhodr(model,1,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,1,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;
	
//		fprintf(stderr,"r=%g rho=%g M=%g u=%g k1rho=%g k1M=%g k1u=%g k2rho=%g k2M=%g k2u=%g\n",
//							r, rho, M, u, k1rho, k1M, k1u, k2rho, k2M, k2u);	
		if (bSetModel) {
			model->M[i] = M;
			model->r[i] = r;
			model->rho[i] = rho;
			model->u[i] = u;
			model->mat[i] = model->iLayer[1];
			++i;
			model->nTable = i;
			model->dr = h;
		}
	}

	/*
	** Now do a linear interpolation to rho == fact*rho0.
	*/
	x = (fact*model->tillMat[1]->rho0 - rho)/k2rho;
	fprintf(stderr,"iMat=%i, rho0=%g, rho=%g, x=%g\n",model->iLayer[1],model->tillMat[1]->rho0,rho,x);
	fprintf(stderr,"rho0-rho=%g, k2rho0=%g\n",model->tillMat[1]->rho0-rho,k2rho);
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
		model->mat[i] = model->iLayer[1];
		++i;
	}	

	// Return values
    *pR = r;
	return(M);
}


/*
** This function integrates the ODEs for a single component model with b.c.
** rho_initial=rho, u_initial=u until rho(r=R)=rho0. It returns the total
** mass of the model. The parameter h sets the stepsize for the RK2 algorithm
** and for bSetModel=1 the results are saved.
*/
double modelSolveSingleComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR)
{
//    FILE *fp;
	// Set inital values for M1 and R
    double M = 0.0;
	double r = 0.0;
	double Mc = model->MLayer[0];
    double k1rho,k1M,k1u,k2rho,k2M,k2u,x;
	int i = 0;

	assert(model->nLayer == 1);

	// (CR) Debug information
	fprintf(stderr,"\n");
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"modelSolveSingleComponent (inital values):\n");
	fprintf(stderr,"rho: %g, u: %g, M:%g, r:%g\n",rho,u,M,r);
	fprintf(stderr,"bSetModel: %i, h: %g\n",bSetModel,h);
	fprintf(stderr,"******************************************************************\n");
	fprintf(stderr,"\n");

	/*
	** Save the values at the core.
	*/
	if (bSetModel) {
		model->rho[i] = rho;
		model->M[i] = M;
		model->u[i] = u;
		model->r[i] = r;
		model->mat[i] = model->iLayer[0];
		++i;
	}
	
	/*
	** The integrate the mantle until rho(r=R)=rho0.
	*/
	while (rho > fact*model->tillMat[0]->rho0) {
		/*
		** Midpoint Runga-Kutta (2nd order).
		*/
		k1rho = h*drhodr(model,0,r,rho,M,u);
		k1M = h*dMdr(r,rho);
		k1u = h*dudr(model,0,r,rho,M,u);

		k2rho = h*drhodr(model,0,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);
		k2M = h*dMdr(r+0.5*h,rho+0.5*k1rho);
		k2u = h*dudr(model,0,r+0.5*h,rho+0.5*k1rho,M+0.5*k1M,u+0.5*k1u);

		rho += k2rho;
		M += k2M;
		u += k2u;
		r += h;
	
		if (bSetModel) {
			model->M[i] = M;
			model->r[i] = r;
			model->rho[i] = rho;
			model->u[i] = u;
			model->mat[i] = model->iLayer[0];
			++i;
			model->nTable = i;
			model->dr = h;
		}
	}

	/*
	** Now do a linear interpolation to rho == fact*rho0.
	*/
	x = (fact*model->tillMat[0]->rho0 - rho)/k2rho;
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
		model->mat[i] = model->iLayer[0];
		++i;
	}	

	// Return values
    *pR = r;
	return(M);
}

/*
** Write the lookup table to ballic.model.
*/
void modelWriteToFile(MODEL *model)
{
	FILE *fp;
	int iLayer;
	int i;

	fprintf(stderr,"Writing model to file.\n");
	fp = fopen("ballic.model","w");
	assert(fp != NULL);

	// (CR) Some problems here with tillMat[model->mat[i]] as the materials are now stored in the order that they appear
	// in the model, e.g., tillMat[iMatLayer1, iMatLayer2,...].	
	fprintf(fp,"#R  rho  M  u  mat tillPressure  tillTempRhoU\n");
	for (i=0; i<model->nTable;i++)
	{
		/* Check which Layer it is from the material. This should be optimized and generalized for more layers!! */
		if (model->mat[i] == model->iLayer[0])
		{
			iLayer = 0;
		} else {
			iLayer = 1;
		}
/*		fprintf(fp,"%g %g %g %g %i %g %g\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i],
			tillPressure(model->tillMat[model->mat[i]], model->rho[i], model->u[i]),
			tillTempRhoU(model->tillMat[model->mat[i]], model->rho[i], model->u[i]));
*/
		fprintf(fp,"%g %g %g %g %i %g %g\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i],
			tillPressure(model->tillMat[iLayer], model->rho[i], model->u[i]),
			tillTempRhoU(model->tillMat[iLayer], model->rho[i], model->u[i]));
//		fprintf(fp,"%g %g %g %g %i\n",model->r[i],model->rho[i],model->M[i],model->u[i],model->mat[i]);

	
	}
	fclose(fp);
}

/*
** This function calls modelSolveComponent for every material for a given
** density rhoc and internal energy uc in the code and returns the
** total mass of the resulting planet. If bSetModel=1 the results are saved in
** the lookup table and both nTable and dr are set.
*/
double modelSolveAll(MODEL *model,int bSetModel,double rhoc,double uc,double h,double *pR)
{
    FILE *fp;
	// Set inital values rhoc, uc and M
	double rho=rhoc;
	double u = uc;
	double M = 0.0;
	// Radius of the model
	double R = 0.0;
	int Index = 0;

	int i;

	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"modelSolveAll:\n");
	fprintf(stderr,"modelSolveAll: bSetModel=%i rhoc=%g uc=%g h=%g\n",bSetModel,rhoc,uc,h);
	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"\n");

	M = modelSolveTwoComponent(model,bSetModel,rho,u,h,&R);
/*
	for (i=0; i<model->nLayer; i++)
	{
		// Solve for component i
		if (i == model->nMat-1)
		{
			fprintf(stderr,"modelSolveAll: Component:%i Index:%i h:%g rho:%g u:%g M1:%g M2:%g R:%g\n",model->iLayer[i], Index, h, rho, u, M,M+model->MLayer[i], R);
			// Enforce rho(r=R)=rho0 for the last layer
			modelSolveComponent(model,model->iLayer[i], bSetModel, 1, &Index, h, &rho, &u, &M,M+model->MLayer[i], &R);
		} else {
			fprintf(stderr,"modelSolveAll: Component:%i Index:%i h:%g rho:%g u:%g M1:%g M2:%g R:%g\n",model->iLayer[i], Index, h, rho, u, M,M+model->MLayer[i], R);
			modelSolveComponent(model,model->iLayer[i], bSetModel, 0, &Index, h, &rho, &u, &M,M+model->MLayer[i], &R);
			// Determine rho2 and u2 for the next material
			modelSolveBC(model, &rho, &u, model->iLayer[i], model->iLayer[i+1]);
			// tillSolveBC(model->tillMat[model->iLayer[i]],model->tillMat[i],rho,u,&rho,&u);

		}
	}
*/
	fprintf(stderr,"***********************************************************\n");
	fprintf(stderr,"modelSolveAll: Done. Index:%i h:%g rho:%g u:%g M:%g R:%g\n", Index, h, rho, u, M, R);
	fprintf(stderr,"***********************************************************\n");

	if (bSetModel) {
		// Set nTable and dr
		model->nTable = Index;
		model->dr = h;
		// Write the lookup table to the file ballic.model
		modelWriteToFile(model);
	}
	
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
    R = cbrt(3.0*M/(4.0*M_PI*model->tillMat[model->nLayer-1]->rho0)); // Use lower density material for max radius
    dr = R/nStepsMax;
    a = 2.0*model->tillMat[0]->rho0; /* starts with 100% larger central density */

	// (CR) Debug
	fprintf(stderr,"R: %g a: %g\n",R,a);

//	Ma = modelSolveSingleComponent(model,bSetModel=0,a,model->uc,dr,&R);
	Ma = modelSolveTwoComponent(model,bSetModel=0,a,model->uc,dr,&R);
//	Ma = modelSolveAll(model,bSetModel=0,a,model->uc,dr,&R);
    fprintf(stderr,"first Ma:%g R:%g\n",Ma,R);
    b = a;
    Mb = 0.5*M;
    while (Ma > M) {
		b = a;
		Mb = Ma;
		a = 0.5*(model->tillMat[model->iLayer[0]]->rho0 + a);
//		Ma = modelSolveSingleComponent(model,bSetModel=0,a,model->uc,dr,&R);
		Ma = modelSolveTwoComponent(model,bSetModel=0,a,model->uc,dr,&R);
//		Ma = modelSolveAll(model,bSetModel=0,a,model->uc,dr,&R);
	}
    while (Mb < M) {
		b = 2.0*b;
//	   	Mb = modelSolveSingleComponent(model,bSetModel=0,b,model->uc,dr,&R);	
		Mb = modelSolveTwoComponent(model,bSetModel=0,b,model->uc,dr,&R);	
//	   	Mb = modelSolveAll(model,bSetModel=0,b,model->uc,dr,&R);	
		fprintf(stderr,"first Mb:%g R:%g\n",Mb,R);
	}

	// (CR) Debug
	fprintf(stderr,"Root bracketed.\n");

    /*
    ** Root bracketed by (a,b).
    */
    while (Mb-Ma > 1e-10*Mc) {
		c = 0.5*(a + b);
//		Mc = modelSolveSingleComponent(model,bSetModel=0,c,model->uc,dr,&R);
		Mc = modelSolveTwoComponent(model,bSetModel=0,c,model->uc,dr,&R);
//		Mc = modelSolveAll(model,bSetModel=0,c,model->uc,dr,&R);	
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
    fprintf(stderr,"rho_core: %g cv: %g uc: %g (in system units)\n",c,model->tillMat[0]->cv,model->uc);
//    Mc = modelSolveSingleComponent(model,bSetModel=1,c,model->uc,dr,&R);
	Mc = modelSolveTwoComponent(model,bSetModel=1,c,model->uc,dr,&R);
	// This is only needed for modelSolveTwoComponent
	modelWriteToFile(model);
//    Mc = modelSolveAll(model,bSetModel=1,c,model->uc,dr,&R);
    model->R = R;
    return c;
    }
/*
** This code depends on the arrays being ordered in r but does not change
** if M, rho or u are not monotonic. However the discontinuities in rho
** u must be handled appropriately.
*/
double MLookup(MODEL *model,double r) {
    double x,xi,dr;
    int i;

    i = model->nTable-1;
    if (r >= model->r[i]) return(model->M[i]*(1.0 + log(r-model->r[i]+1)));
    x = r/model->dr;
    xi = floor(x);
	//(CR) Debugging
//	fprintf(stderr,"r=%g dr=%g x=%g xi=%g\n",r,model->dr,x,xi);
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
//		fprintf(stderr,"c:%.7g Mc:%.7g\n",c,Mc);
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
//		printf("c:%.14g M(c):%.14g\n",c,MLookup(model,c));
	}
	return c;
	}

/*
** Generate a HEALPix or Icosahedron grid and distribute the particles on it.
** mTot:		Desired mass total mass
** nDesired:	Desired number of particles (can vary due to constraints from the grid)
** ucore:		Value of the total energy in the center of the model (ucore = u(r=0))
** bCentral:	Do we want a central particle
** bIcosa:		Use Icosahedron grid (if 0 use HEALPIX) 
**
** This is a two component version of ballic that solves the problem of distributing the particles separately for the core and the mantle.
*/
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
	double mCore;
	/* A first guess for the number of particles of a layer. */
	int nDesiredLayer;
	int nReachedLastLayer = 0;
	int iShellLastLayer = 0;
	int iLayer;
	double rLayer; /* Final Radius of one layer. */

    if (argc != 4) {
		fprintf(stderr,"Usage: ballic <nDesired> <TotalMass> <ucore> >myball.std\n");
		exit(1);
	}
	
	/* Get command line parameters. */
	nDesired = atoi(argv[1]);
	mTot = atof(argv[2]);
	ucore = atof(argv[3]);

    model = modelInit(mTot,ucore);
	double R = 0.0;

	fprintf(stderr,"Model initialized.\n"); // CR
    rhoCenter = modelSolve(model,mTot);
	fprintf(stderr,"Model solved.\n");	// CR
#if 0
	/*
	** Desired number of particles in the core (Nc=fc*N).
	*/
	mCore = model->MLayer[0];
	nDesiredCore = mCore/mTot*nDesired;
#endif
	//m = mTot/nDesired;   /* a first guess at the particle mass */

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
	
	TipsyInitialize(&out,0,NULL);

	fprintf(stderr,"\n");
	fprintf(stderr,"*******************************************************\n");
	fprintf(stderr,"Distributing shells.\n");
	fprintf(stderr,"nLayer=%i nDesired=%i mTot=%g\n",model->nLayer,nDesired,mTot);
	fprintf(stderr,"*******************************************************\n");
	fprintf(stderr,"\n");
	for (iLayer=0;iLayer<model->nLayer;iLayer++)
	{
		fprintf(stderr,"Layer %i (Material %i):\n",iLayer,model->iLayer[iLayer]);
		fprintf(stderr,"tillMat=%i MLayer=%g nDesiredMat=%g m=%g (estimate)\n",
		model->iLayer[iLayer],
		model->MLayer[iLayer],
		model->MLayer[iLayer]/mTot*nDesired,
		model->MLayer[iLayer]/(model->MLayer[iLayer]/mTot*nDesired));
	}

	/*
	** Distribute the particles in shells for each material separately.
	*/
	for (iLayer=0;iLayer<model->nLayer;iLayer++)
	{
		/* Careful, iLayer is the number in which the materials are ordered in model->iLayer[]. */
		nDesiredLayer = model->MLayer[iLayer]/mTot*nDesired+nReachedLastLayer;
		m = model->MLayer[iLayer]/nDesiredLayer;   /* a first guess at the particle mass */
		
		fprintf(stderr,"\n");
		fprintf(stderr,"Layer: %i Material %i: nDesiredLayer=%i M=%g m=%g\n",iLayer,model->iLayer[iLayer], nDesiredLayer, model->MLayer[iLayer],m);

		/*
		** Phase 1: We first determine the number of particles per shell such that
		** the ratio of radial to tangential lengths of their volumes is as close
		** to 1 as possible. This fixes the total number of particles in the model.
		*/
		fprintf(stderr,"PHASE 1: ODE approach (iLayer=%i)\n",iLayer);
		for (iter=0;iter<2;++iter) {
			if (iLayer == 0)
			{
				/* Initialize nReached, ro and iShell for material 0. */
				nReached = 0;
				/* Set starting radius */
				if (bCentral) {
					ro = rShell(model,m,0.0);
				} else {
					ro = 0.0;
				}
				iShell = 0;
				rLayer = 0.0;
			} else {
				/* For the other layers we start from rLayer of the previous material. */
				ro = rLayer;
				nReached = nReachedLastLayer;
				iShell = iShellLastLayer;
			}

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

//					(CR) Debugging
//					fprintf(stderr,"npix=%i ros=%g ri=%g ns=%i\n",npix,ros,ri,ns);
					l1 = ros-ri;
					l2 = sqrt(M_PI/npix)*(ros+ri);
					rts = l1/l2;
/*					fprintf(stderr,"ns:%d rts:%g\n",ns,rts); */
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
/*				fprintf(stderr,"nReached:%d npix:%d\n",nReached,npix);*/
				if ((nReached + npix) < nDesiredLayer) {
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
			if (nDesiredLayer - nReached > npix/2) {
				nReached += npix;
				nsShell[nShell] = ns;
				fprintf(stderr,"iShell:%d ns:%d radial/tangential:??? (added)\n",nShell,ns);
				nShell++;
			}

			fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
				nReached,m,model->MLayer[iLayer]/nReached);
//			fprintf(stderr,"nReached:%d old mass:%.7g new mass:%.7g\n",
//				nReached,m,mTot/nReached);
			/* Update particles mass for the reached number of particles. */
			m = model->MLayer[iLayer]/(nReached-nReachedLastLayer);
			fprintf(stderr,"PHASE 1 update m: nReached=%i nReachedLastLayer=%i dReached=%i iShell=%i iShellLastLayer=%i nDesired=%i M=%g m=%g\n",nReached,nReachedLastLayer,nReached-nReachedLastLayer,iShell,iShellLastLayer,nDesiredLayer,model->MLayer[iLayer],m);
//			m = mTot/nReached;
			nDesiredLayer = nReached+1;
		}
//assert(0);
		fprintf(stderr,"PHASE 1 done: nReached=%i nReachedLastLayer=%i iShell=%i iShellLastLayer=%i nDesired=%i M=%g m=%g\n",nReached,nReachedLastLayer,iShell,iShellLastLayer,nDesiredLayer,model->MLayer[iLayer],m);
		/* Save how many particles we distributed for the last material. */
//		nReachedLastLayer = nReached;
//		iShellLastLayer = iShell;

		/*
		** Phase 2: With the numbers of particles in each of the shells now fixed
		** we continue by recalculating the radii of the shells based on the updated
		** particle mass.
		*/
		fprintf(stderr,"PHASE 2 (iLayer=%i)\n",iLayer);
		// (CR) For iLayer > 0 we have to set ro properly...
		if (iLayer == 0)
		{
			if (bCentral) {
				ro = rShell(model,m,0.0);
			} else {
				ro = 0.0;
			}
		} else {
			/* For the other components we start from rLayer of the previous material. */
			ro = rLayer;
		}
		
		for (iShell=iShellLastLayer;iShell<nShell;++iShell) {
			ri = ro;
			ns = nsShell[iShell];
			ro = rShell2(model,bIcosa,m,ri,ns);
			npix = (bIcosa)?(40*ns*(ns-1)+12):(12*ns*ns);
			l1 = ro-ri;
			l2 = sqrt(M_PI/npix)*(ro+ri);
			rts = l1/l2;

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

			// Set rLayer to ro
			rLayer = ro;
		}
		fprintf(stderr,"PHASE 2 done: nReached=%i nReachedLastLayer=%i iShell=%i iShellLastLayer=%i nDesired=%i M=%g m=%g rLayer=%g ro=%g\n",nReached,nReachedLastLayer,iShell,iShellLastLayer,nDesiredLayer,model->MLayer[iLayer],m,rLayer,ro);
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
		for (j=0;j<3;++j) gp.vel[j] = 0.0;
		gp.phi = 0.0;
		gp.mass = m;
		/* Save the material */
		gp.metals = model->iLayer[iLayer];
		nLast = nReached;
		/* Start writing from nReachedLastLayer otherwise we write the inner components more than once. */
		nReached = nReachedLastLayer;
//		nReached = 0;

//*****************************************************************************************************************************************
//		(CR) 19.1.16: Writing a central particle that was not originally included in the calculation could cause the mass to deviate from M
//*****************************************************************************************************************************************
		if (bCentral && iLayer==0) {
			/* Write the central particle (only for iLayer=0). */
			for (j=0;j<3;++j) gp.pos[j] = 0.0;
			gp.temp = uLookup(model, 0);
			// Dont forget to set the material for the central particle
			//gp.metals = model->iLayer[iLayer];
			TipsyAddGas(out,&gp);
		}

		/* Start writing from iShellLastLayer otherwise we write the inner components more than once. */
		for (iShell=iShellLastLayer;iShell<nShell;++iShell) {
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
	    
//				rho = rhoLookup(model,rs);
				gp.temp = uLookup(model,rs);
				TipsyAddGas(out,&gp);
			}
		}

		/* Save how many particles we distributed for the last material. */
		nReachedLastLayer = nReached;
		iShellLastLayer = iShell;

		fprintf(stderr,"iLayer %i: iShell=%i nReached=%i\n",iLayer, iShell, nReached);
	} /* iLayer */
 	fprintf(stderr,"\n");
	fprintf(stderr,"Writing %d particles. Model R:%g Last Shell r:%g\n",nReached,model->R,rsShell[nShell-1]);

	/* Write all particles to ballic.std */
	TipsyWriteAll(out,0.0,"ballic.std");
	TipsyFinish(out);

	/*
	** This code was added by Joachim to have a radial density profile (smoothed) for debugging.
	*/
#if 0	
	/* Smooth with gather doesn work with my version! */
	system("sleep 1; smooth -s 128 density <ballic.std; sleep 1");

	/* Write the smoothed quantities to a file. */
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


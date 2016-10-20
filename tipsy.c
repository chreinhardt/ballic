#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <malloc.h>
#include <assert.h>
#include "tipsy.h"

int xdr_header(XDR *xdrs, struct dump *header)
{
	int pad=0;
  
	if (xdr_double(xdrs,&header->time) != TRUE) return 0;
	if (xdr_u_int(xdrs,&header->nbodies) != TRUE) return 0;
	if (xdr_u_int(xdrs,&header->ndim) != TRUE) return 0;
	if (xdr_u_int(xdrs,&header->nsph) != TRUE) return 0;
	if (xdr_u_int(xdrs,&header->ndark) != TRUE) return 0;
	if (xdr_u_int(xdrs,&header->nstar) != TRUE) return 0;
	if (xdr_u_int(xdrs,&pad) != TRUE) return 0;
	return 1;
	}


int xdr_gas(XDR *xdrs,struct gas_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
#if defined(DOUBLE_POS) || defined(DOUBLE_POS_VEL)
	if (xdr_double(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
#endif
#ifdef DOUBLE_POS_VEL
	if (xdr_double(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
#endif
	if (xdr_float(xdrs,&p->rho) != TRUE) return 0;
	if (xdr_float(xdrs,&p->temp) != TRUE) return 0;
	if (xdr_float(xdrs,&p->hsmooth) != TRUE) return 0;
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_dark(XDR *xdrs,struct dark_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
#if defined(DOUBLE_POS) || defined(DOUBLE_POS_VEL)
	if (xdr_double(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
#endif
#ifdef DOUBLE_POS_VEL
	if (xdr_double(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
#endif
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_star(XDR *xdrs,struct star_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
#if defined(DOUBLE_POS) || defined(DOUBLE_POS_VEL)
	if (xdr_double(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->pos[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
#endif
#ifdef DOUBLE_POS_VEL
	if (xdr_double(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_double(xdrs,&p->vel[2]) != TRUE) return 0;
#else
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
#endif
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->tform) != TRUE) return 0;
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


#define TIPSY_EXTEND_ARRAY		10000


void TipsyInitialize(TCTX *pctx,int bNative,char *pchInFile)
{
	TCTX ctx;
	int ret;

	ctx = malloc(sizeof(struct TipsyContext));
	assert(ctx != NULL);
	ctx->bNative = bNative;
	if (pchInFile != NULL) {
		ctx->pchInFile = malloc(strlen(pchInFile)+1);
		assert(ctx->pchInFile != NULL);
		strcpy(ctx->pchInFile,pchInFile);
		if (!strcmp(pchInFile,"stdin")) {
			ctx->fp = stdin;
			}
		else {
			ctx->fp = fopen(ctx->pchInFile,"r");
			if (ctx->fp == NULL) {
				fprintf(stderr,"ERROR(TipsyInitialize): Could not open input file: %s\n",
						ctx->pchInFile);
				exit(1);
				}
			}
		if (bNative) {
			ret = fread(&ctx->head,sizeof(struct dump),1,ctx->fp);
			ctx->nRemain = ctx->head.nsph;
			ctx->iCurrType = TIPSY_TYPE_GAS;
			ctx->iBlock = 0;
			}
		else {
			xdrstdio_create(&ctx->xdr,ctx->fp,XDR_DECODE);
			ret = xdr_header(&ctx->xdr,&ctx->head);
			fprintf(stderr,"nbodies:%u\n",ctx->head.nbodies);
			}
		if (ret != 1) {
			fprintf(stderr,"ERROR(TipsyInitialize): Could not read header from: %s\n",
					ctx->pchInFile);
			exit(1);
			}
		}
	else {
		ctx->fp = NULL;
		ctx->pchInFile = NULL;
		ctx->head.time = 0;
		ctx->head.nbodies = 0;
		ctx->head.ndim = 0;
		ctx->head.nsph = 0;
		ctx->head.ndark = 0;
		ctx->head.nstar = 0;
		}
	ctx->nMaxGas = 0;
	ctx->nGas = 0;
	ctx->nMaxDark = 0;
	ctx->nDark = 0;
	ctx->nMaxStar = 0;
	ctx->nStar = 0;
	ctx->gp = NULL;
	ctx->dp = NULL;
	ctx->sp = NULL;
	ctx->iCurr = 0;
	*pctx = ctx;
	}


void TipsyFinish(TCTX ctx)
{
	if (ctx->fp != NULL) {
		if (strcmp(ctx->pchInFile,"stdin") != 0) { 
			fclose(ctx->fp);
			}
		if (!ctx->bNative) {
			xdr_destroy(&ctx->xdr);
			}
		}
	if (ctx->pchInFile) free(ctx->pchInFile);
	if (ctx->gp) free(ctx->gp);
	if (ctx->dp) free(ctx->dp);
	if (ctx->sp) free(ctx->sp);
	free(ctx);
	}


struct base_particle *pTipsyReadNative(TCTX ctx,int *piType,double *pdSoft)
{
	int i = ctx->iBlock;
	int ret;

	switch (ctx->iCurrType) {
	case TIPSY_TYPE_GAS:
		if (i == 0) {
			ctx->nBlock = TIPSY_BLOCK;
			if (ctx->nBlock > ctx->nRemain) ctx->nBlock = ctx->nRemain;
			if (ctx->nBlock > 0) {
				if (ctx->gp == NULL) {
					ctx->gp = malloc(TIPSY_BLOCK*sizeof(struct gas_particle));
					assert(ctx->gp != NULL);
					}
				ret = fread(ctx->gp,sizeof(struct gas_particle),ctx->nBlock,ctx->fp);
				if (ret != ctx->nBlock) {
					fprintf(stderr,"ERROR(pTispyReadNative): Could not read gas particles from: %s\n",
							ctx->pchInFile);
					exit(1);
					}
				ctx->nRemain -= ctx->nBlock;
				}
			else {
				ctx->nRemain = ctx->head.ndark;
				ctx->iCurrType = TIPSY_TYPE_DARK;
				if (ctx->gp != NULL) {
					free(ctx->gp);
					ctx->gp == NULL;
					};
				goto jump_dark;
				}
			}
		if (++ctx->iBlock == ctx->nBlock) ctx->iBlock = 0;
		*pdSoft = ctx->gp[i].hsmooth;
		*piType = ctx->iCurrType;
		return((struct base_particle *)&ctx->gp[i]);
	case TIPSY_TYPE_DARK:
		if (i == 0) {
		jump_dark:
			ctx->nBlock = TIPSY_BLOCK;
			if (ctx->nBlock > ctx->nRemain) ctx->nBlock = ctx->nRemain;
			if (ctx->nBlock > 0) {
				if (ctx->dp == NULL) {
					ctx->dp = malloc(TIPSY_BLOCK*sizeof(struct dark_particle));
					assert(ctx->dp != NULL);
					}
				ret = fread(ctx->dp,sizeof(struct dark_particle),ctx->nBlock,ctx->fp);
				if (ret != ctx->nBlock) {
					fprintf(stderr,"ERROR(pTispyReadNative): Could not read dark particles from: %s\n",
							ctx->pchInFile);
					exit(1);
					}
				ctx->nRemain -= ctx->nBlock;
				}
			else {
				ctx->nRemain = ctx->head.nstar;
				ctx->iCurrType = TIPSY_TYPE_STAR;
				if (ctx->dp != NULL) {
					free(ctx->dp);
					ctx->dp == NULL;
					};
				goto jump_star;
				}
			}
		if (++ctx->iBlock == ctx->nBlock) ctx->iBlock = 0;
		*pdSoft = ctx->dp[i].eps;
		*piType = ctx->iCurrType;
		return((struct base_particle *)&ctx->dp[i]);
	case TIPSY_TYPE_STAR:
		if (i == 0) {
		jump_star:
			ctx->nBlock = TIPSY_BLOCK;
			if (ctx->nBlock > ctx->nRemain) ctx->nBlock = ctx->nRemain;
			if (ctx->nBlock > 0) {
				if (ctx->sp == NULL) {
					ctx->sp = malloc(TIPSY_BLOCK*sizeof(struct star_particle));
					assert(ctx->sp != NULL);
					}
				ret = fread(ctx->sp,sizeof(struct star_particle),ctx->nBlock,ctx->fp);
				if (ret != ctx->nBlock) {
					fprintf(stderr,"ERROR(pTispyReadNative): Could not read star particles from: %s\n",
							ctx->pchInFile);
					exit(1);
					}
				ctx->nRemain -= ctx->nBlock;
				}
			else {
				if (ctx->sp != NULL) {
					free(ctx->sp);
					ctx->sp == NULL;
					};
				*piType = 0;
				return(NULL);
				}
			}
		if (++ctx->iBlock == ctx->nBlock) ctx->iBlock = 0;
		*pdSoft = ctx->sp[i].eps;
		*piType = ctx->iCurrType;
		return((struct base_particle *)&ctx->sp[i]);
		}
	}


struct base_particle *pTipsyRead(TCTX ctx,int *piType,double *pdSoft)
{
	int ret;

	if (ctx->iCurr < ctx->head.nsph) {
		if (ctx->bNative) {
			ret = fread(&ctx->gpCurr,sizeof(struct gas_particle),1,ctx->fp);
			}
		else {
			ret = xdr_gas(&ctx->xdr,&ctx->gpCurr);
			}
		if (ret != 1) {
			fprintf(stderr,"ERROR(pTispyRead): Could not read gas particle from: %s\n",
					ctx->pchInFile);
			exit(1);
			}
		++ctx->iCurr;
		*pdSoft = ctx->gpCurr.hsmooth;
		*piType = TIPSY_TYPE_GAS;
		return((struct base_particle *)&ctx->gpCurr);
		}
	else if (ctx->iCurr < ctx->head.nsph+ctx->head.ndark) {
		if (ctx->bNative) {
			ret = fread(&ctx->dpCurr,sizeof(struct dark_particle),1,ctx->fp);
			}
		else {
			ret = xdr_dark(&ctx->xdr,&ctx->dpCurr);
			}
		if (ret != 1) {
			fprintf(stderr,"ERROR(pTispyRead): Could not read dark particle from: %s\n",
					ctx->pchInFile);
			exit(1);
			}
		++ctx->iCurr;
		*pdSoft = ctx->dpCurr.eps;
		*piType = TIPSY_TYPE_DARK;
		return((struct base_particle *)&ctx->dpCurr);
		}
	else if (ctx->iCurr < ctx->head.nbodies) {
		if (ctx->bNative) {
			ret = fread(&ctx->spCurr,sizeof(struct star_particle),1,ctx->fp);
			}
		else {
			ret = xdr_star(&ctx->xdr,&ctx->spCurr);
			}
		if (ret != 1) {
			fprintf(stderr,"ERROR(pTispyRead): Could not read star particle from: %s\n",
					ctx->pchInFile);
			exit(1);
			}
		++ctx->iCurr;
		*pdSoft = ctx->spCurr.eps;
		*piType = TIPSY_TYPE_STAR;
		return((struct base_particle *)&ctx->spCurr);
		}
	else {
		*piType = 0;
		return(NULL);
		}
	}


struct base_particle *pTipsyParticle(TCTX ctx,unsigned iIndex,int *piType,double *pdSoft)
{
	*piType = 0;
	if (iIndex < 0) return(NULL);
	if (iIndex < ctx->nGas) {
		*pdSoft = ctx->gp[iIndex].hsmooth;
		*piType = TIPSY_TYPE_GAS;
		return((struct base_particle *)&ctx->gp[iIndex]);
		}
	iIndex -= ctx->nGas;
	if (iIndex < ctx->nDark) {
		*pdSoft = ctx->dp[iIndex].eps;
		*piType = TIPSY_TYPE_DARK;
		return((struct base_particle *)&ctx->dp[iIndex]);
		}
	iIndex -= ctx->nDark;
	if (iIndex < ctx->nStar) {
		*pdSoft = ctx->sp[iIndex].eps;
		*piType = TIPSY_TYPE_STAR;
		return((struct base_particle *)&ctx->sp[iIndex]);
		}
	return(NULL);
	}


void TipsyReadAll(TCTX ctx)
{
	int ret,i;

	/*
	 ** Make sure we are at the start of this file, otherwise we will not
	 ** be reading all of the particles which this function implies.
	 */
	assert(ctx->iCurr == 0);
	if (ctx->nGas != 0 || ctx->nDark != 0 || ctx->nStar != 0) {
	  fprintf(stderr,"WARNING(TispyReadAll): Reading file after adding particles, destroys added particles!\n");
	}
	if (ctx->head.nsph > ctx->nMaxGas) {
		ctx->nMaxGas = ctx->head.nsph;
		ctx->gp = realloc(ctx->gp,ctx->nMaxGas*sizeof(struct gas_particle));
		if (ctx->gp == NULL) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not extend gas particle array.\n");
			exit(1);
			}
		}
	if (ctx->head.ndark > ctx->nMaxDark) {
		ctx->nMaxDark = ctx->head.ndark;
		ctx->dp = realloc(ctx->dp,ctx->nMaxDark*sizeof(struct dark_particle));
		if (ctx->dp == NULL) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not extend dark particle array.\n");
			exit(1);
			}
		}
	if (ctx->head.nstar > ctx->nMaxStar) {
		ctx->nMaxStar = ctx->head.nstar;
		ctx->sp = realloc(ctx->sp,ctx->nMaxStar*sizeof(struct star_particle));
		if (ctx->sp == NULL) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not extend star particle array.\n");
			exit(1);
			}
		}
	if (ctx->bNative) {
		ret = fread(&ctx->gp[ctx->nGas],sizeof(struct gas_particle),ctx->head.nsph,ctx->fp);
		if (ret != ctx->head.nsph) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not read all gas particles from: %s\n",
					ctx->pchInFile);
			fprintf(stderr,"                     Gas paricles read: %d, expected: %d\n",
					ret,ctx->head.nsph);
			exit(1);
			}
		ctx->nGas = ctx->head.nsph;
		ret = fread(&ctx->dp[ctx->nDark],sizeof(struct dark_particle),ctx->head.ndark,ctx->fp);
		if (ret != ctx->head.ndark) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not read all dark particles from: %s\n",
					ctx->pchInFile);
			fprintf(stderr,"                     Dark paricles read: %d, expected: %d\n",
					ret,ctx->head.ndark);
			exit(1);
			}
		ctx->nDark = ctx->head.ndark;
		ret = fread(&ctx->sp[ctx->nStar],sizeof(struct star_particle),ctx->head.nstar,ctx->fp);
		if (ret != ctx->head.nstar) {
			fprintf(stderr,"ERROR(TispyReadAll): Could not read all star particles from: %s\n",
					ctx->pchInFile);
			fprintf(stderr,"                     Star paricles read: %d, expected: %d\n",
					ret,ctx->head.nstar);
			exit(1);
			}
		ctx->nStar = ctx->head.nstar;
		}
	else {
		for(i=0;i<ctx->head.nsph;++i) {
			ret = xdr_gas(&ctx->xdr,&ctx->gp[ctx->nGas]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TispyReadAll): Could not read all gas particles from: %s\n",
						ctx->pchInFile);
				fprintf(stderr,"                     Gas paricles read: %d, expected: %d\n",
						i,ctx->head.nsph);
				exit(1);
				}
			++ctx->nGas;
			}
		for(i=0;i<ctx->head.ndark;++i) {
			ret = xdr_dark(&ctx->xdr,&ctx->dp[ctx->nDark]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TispyReadAll): Could not read all dark particles from: %s\n",
						ctx->pchInFile);
				fprintf(stderr,"                     Dark paricles read: %d, expected: %d\n",
						i,ctx->head.ndark);
				exit(1);
				}
			++ctx->nDark;
			}
		for(i=0;i<ctx->head.nstar;++i) {
			ret = xdr_star(&ctx->xdr,&ctx->sp[ctx->nStar]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TispyReadAll): Could not read all star particles from: %s\n",
						ctx->pchInFile);
				fprintf(stderr,"                     Star paricles read: %d, expected: %d\n",
						i,ctx->head.nstar);
				exit(1);
				}
			++ctx->nStar;
			}
		}
	ctx->iCurr = ctx->head.nbodies;
	}


void TipsyAddGas(TCTX ctx,struct gas_particle *pGas)
{
	if (ctx->nGas == ctx->nMaxGas) {
		ctx->nMaxGas += TIPSY_EXTEND_ARRAY;
		ctx->gp = realloc(ctx->gp,ctx->nMaxGas*sizeof(struct gas_particle));
		if (ctx->gp == NULL) {
			fprintf(stderr,"ERROR(TispyAddGas): Could not extend gas particle array.\n");
			exit(1);
			}
		}
	ctx->gp[ctx->nGas++] = *pGas;
	}

void TipsyAddDark(TCTX ctx,struct dark_particle *pDark)
{
	if (ctx->nDark == ctx->nMaxDark) {
		ctx->nMaxDark += TIPSY_EXTEND_ARRAY;
		ctx->dp = realloc(ctx->dp,ctx->nMaxDark*sizeof(struct dark_particle));
		if (ctx->dp == NULL) {
			fprintf(stderr,"ERROR(TispyAddDark): Could not extend dark particle array.\n");
			exit(1);
			}
		}
	ctx->dp[ctx->nDark++] = *pDark;
	}

void TipsyAddStar(TCTX ctx,struct star_particle *pStar)
{
	if (ctx->nStar == ctx->nMaxStar) {
		ctx->nMaxStar += TIPSY_EXTEND_ARRAY;
		ctx->sp = realloc(ctx->sp,ctx->nMaxStar*sizeof(struct star_particle));
		if (ctx->sp == NULL) {
			fprintf(stderr,"ERROR(TispyAddStar): Could not extend star particle array.\n");
			exit(1);
			}
		}
	ctx->sp[ctx->nStar++] = *pStar;
	}

/*
 ** If an input file was specified, can we add the added particles?
 ** This isn't really needed though.
 */
void TipsyWriteAll(TCTX ctx,double dTime,char *pchFileName)
{
	FILE *fp;
	XDR xdr;
	struct dump h;
	int ret,i;

	if (pchFileName == NULL) {
		fp = stdout;
		}
	else {
		fp = fopen(pchFileName,"w");
		if (fp == NULL) {
			fprintf(stderr,"ERROR(TipsyWrite): Could not open output file: %s\n",
					pchFileName);
			exit(1);
			}
		}
	h.time = dTime;
	h.nbodies = ctx->nGas + ctx->nDark + ctx->nStar;
   	h.ndim = 3;
	h.nsph = ctx->nGas;
	h.ndark = ctx->nDark;
	h.nstar = ctx->nStar;
	if (ctx->bNative) {
		ret = fwrite(&h,sizeof(struct dump),1,fp);
		}
	else {
		xdrstdio_create(&xdr,fp,XDR_ENCODE);
		ret = xdr_header(&xdr,&h);
		}
	if (ret != 1) {
		fprintf(stderr,"ERROR(TipsyWrite): Could not write header to ");
		if (pchFileName == NULL) {
			fprintf(stderr,"stdout\n");
			}
		else {
			fprintf(stderr,"output file: %s\n",pchFileName);
			}
		exit(1);
		}
	if (ctx->bNative) {
		ret = fwrite(ctx->gp,sizeof(struct gas_particle),ctx->nGas,fp);
		if (ret != ctx->nGas) {
			fprintf(stderr,"ERROR(TipsyWrite): Could not write gas particles to ");
			if (pchFileName == NULL) {
				fprintf(stderr,"stdout\n");
				}
			else {
				fprintf(stderr,"output file: %s\n",pchFileName);
				}
			exit(1);
			}
		ret = fwrite(ctx->dp,sizeof(struct dark_particle),ctx->nDark,fp);
		if (ret != ctx->nDark) {
			fprintf(stderr,"ERROR(TipsyWrite): Could not write dark particles to ");
			if (pchFileName == NULL) {
				fprintf(stderr,"stdout\n");
				}
			else {
				fprintf(stderr,"output file: %s\n",pchFileName);
				}
			exit(1);
			}
		ret = fwrite(ctx->sp,sizeof(struct star_particle),ctx->nStar,fp);
		if (ret != ctx->nStar) {
			fprintf(stderr,"ERROR(TipsyWrite): Could not write star particles to ");
			if (pchFileName == NULL) {
				fprintf(stderr,"stdout\n");
				}
			else {
				fprintf(stderr,"output file: %s\n",pchFileName);
				}
			exit(1);
			}
		}
	else {
		for (i=0;i<ctx->nGas;++i) {
			ret = xdr_gas(&xdr,&ctx->gp[i]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TipsyWrite): Could not write gas particle to ");
				if (pchFileName == NULL) {
					fprintf(stderr,"stdout\n");
					}
				else {
					fprintf(stderr,"output file: %s\n",pchFileName);
					}
				exit(1);
				}
			}
		for (i=0;i<ctx->nDark;++i) {
			ret = xdr_dark(&xdr,&ctx->dp[i]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TipsyWrite): Could not write dark particle to ");
				if (pchFileName == NULL) {
					fprintf(stderr,"stdout\n");
					}
				else {
					fprintf(stderr,"output file: %s\n",pchFileName);
					}
				exit(1);
				}
			}
		for (i=0;i<ctx->nStar;++i) {
			ret = xdr_star(&xdr,&ctx->sp[i]);
			if (ret != 1) {
				fprintf(stderr,"ERROR(TipsyWrite): Could not write star particle to ");
				if (pchFileName == NULL) {
					fprintf(stderr,"stdout\n");
					}
				else {
					fprintf(stderr,"output file: %s\n",pchFileName);
					}
				exit(1);
				}
			}
		}
	if (!ctx->bNative) xdr_destroy(&xdr);
	if (pchFileName != NULL) fclose(fp);
	}


unsigned iTipsyNumParticles(TCTX ctx)
{
  if (ctx->nGas + ctx->nDark + ctx->nStar == 0) return(ctx->head.nbodies);
  else return(ctx->nGas + ctx->nDark + ctx->nStar);
}


double dTipsyTime(TCTX ctx)
{
	return(ctx->head.time);
	}

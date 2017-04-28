#ifndef FLOATTYPE_INCLUDED
#define FLOATTYPE_INCLUDED

#include <limits.h>
#include <float.h>
#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.402823466E+38F
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

#ifndef SINGLE
#define FLOAT double
#define FLOAT_MAXVAL DBL_MAX
/*
** Wrapper function to allow fully double precision tipsy standard files!
** NOTE: Gasoline is never compiled with SINGLE!
*/
#ifdef DOUBLE_TIPSY
static xdr_FLOAT(XDR *xdrs,FLOAT *pf) {
    return xdr_double(xdrs,pf);
    }
#else
static int xdr_FLOAT(XDR *xdrs,FLOAT *pf) {
    float fTmp = *pf;
    int iRet;
	/*
	 * (CR) and (MH): Changed the code, so that the value in fTmp is stored in pf before
	 * returning the error code.
	 */
	iRet = xdr_float(xdrs,&fTmp);
	*pf = fTmp;
	return iRet;
    }
#endif
#else
#define FLOAT float
#define FLOAT_MAXVAL FLT_MAX
/*
** NOTE: Gasoline is never compiled with SINGLE! BUT just in case it is,
** make sure it doesn't break.
*/
static int xdr_FLOAT(XDR *xdrs,FLOAT *pf) {
    return xdr_float(xdrs,pf);
    }

#endif

#endif

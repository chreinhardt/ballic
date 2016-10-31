/*
 ** The header file for the Tillotson EOS library.
 ** This file contains all the code needed to implement
 ** the correction at a material interface proposed by
 ** Woolfson2007.
 */
#ifndef TILLWOOLFSON_HINCLUDED
#define TILLWOOLFSON_HINCLUDED

#include "tillotson.h"

typedef struct till_interf_lookup_entry
{
//	double u;
//	double rho;
	double P;
	double T;
	double f;
} TILL_INTERF_LOOKUP_ENTRY;

typedef struct tillinterface
{
	TILLMATERIAL **tillMat;		/* Array that contains pointers to the materials */
	int nMat;					/* How many materials are there? */
	int nTableP;				/* Number of entries in the look up table in P */
	int nTableT;				/* and in T */
	double Pmin;				/* Min value for the lookup table */
	double Pmax;				/* Max value for the lookup table */
	double Tmin;
	double Tmax;

	double dP;
	double dT;
	
	/* An array of lookup tables for f_ij */
	TILL_INTERF_LOOKUP_ENTRY **InterfLookup;
} TILLINTERFACE;

//double tillInitInterfaceLookup(TILLMATERIAL *material1 TILLMATERIAL *material2, double rho, double u);


#endif


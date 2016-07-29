/*
 ** This is a simple program to test the Tillotson EOS library.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*granite->nTableV) + (j))

void main(int argc, char **argv) {
	/*
	** Debug the look up table for the isentropic evolution
	** the internal energy. We search for the smallest and
	** largest du(v) in the lookup table and print it to a
	** file. 
	*/
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	double rhomax = 100.0;
	double vmax = 40.0;
	// Try 100 only but we need 1000 x 1000 to get good results
	// vmax=25.0, rhomax=100.0, nTableV=10, nTableRho=4000
	//int nTableMax = 1000;
	int nTableRho = 1000;
	int nTableV = 100;
	double rho, v, u;

	/* Initialise du_min and du_max with crazy values */
	double du_min = 1e30;
	double du_max = -1.0; 

	int i = 0;
	int j = 0;
	int i_min = -1;
	int j_min = -1;
	int i_max = -1;
	int j_max = -1;

	TILLMATERIAL *granite;
	struct lookup *isentrope;

	fprintf(stderr, "Initializing material...\n");

	granite = tillInitMaterial(GRANITE, dKpcUnit, dMsolUnit, nTableRho, nTableV, rhomax, vmax);
	
	fprintf(stderr, "Initializing the look up table...\n");
	/* Solve ODE and splines */
	tillInitLookup(granite);
	fprintf(stderr, "Done.\n");

	fprintf(stderr,"\n");
	fprintf(stderr,"rhomax: %g, vmax: %g \n", granite->rhomax, granite->vmax);
	fprintf(stderr,"nTableRho: %i, nTableV: %i \n", granite->nTableRho, granite->nTableV);
	fprintf(stderr,"drho: %g, dv: %g \n", granite->drho, granite->dv);
	fprintf(stderr,"\n");
	
	rho = 0.0;
	v = 0.0;
	u = 0.0;

	/* Create an output file for the look up table */
	FILE *fp = NULL;

	/*
	** Print the look up table to a file first.
	*/	
#if 0
	//sprintf(achFile,"%s.log",msrOutName(msr));
	fp = fopen("lookup.txt","w");
	assert(fp != NULL);

	for (i=0;i<granite->nTableRho;i+=1)
	{
		rho = i*granite->drho;
		fprintf(fp,"%g",rho);
		for (j=0;j<granite->nTableV;j+=1)
		{
			// v = j*granite->dv
			u = granite->Lookup[INDEX(i, j)].u;
			fprintf(fp,"  %g", u);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
#endif

	/* Search the grid for du_min and du_max */
	for (i=0;i<granite->nTableRho;i+=1)
	{
		for (j=0;j<granite->nTableV-1;j+=1)
		{
			double du = fabs(granite->Lookup[INDEX(i, j)].u-granite->Lookup[INDEX(i, j+1)].u);
			
			if (du <= du_min)
			{
				du_min = du;
				i_min = i;
				j_min = j;
			}
			if (du >= du_max)
			{
				du_max = du;
				i_max = i;
				j_max = j;
			}
		}
	}

	fprintf(stderr,"du_min=%g, i_min=%i, j_min=%i\n",du_min,i_min,j_min);
	fprintf(stderr,"du_max=%g, i_max=%i, j_max=%i\n",du_max,i_max,j_max);
	fprintf(stderr,"\n");

	fprintf(stderr,"Done.\n");
	tillFinalizeMaterial(granite);
}

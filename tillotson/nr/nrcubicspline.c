#define NRANSI

#include "nrcubicspline.h"

/* (CR) Changed all functions from float to double on 25.11.15 */
void tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n)
{
	unsigned long j;
	double bet,*gam;

//	gam=vector(1,n);
	gam=dvector(1,n);
	if (b[1] == 0.0) nrerror("Error 1 in tridag");
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0)	nrerror("Error 2 in tridag");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
//	free_vector(gam,1,n);
	free_dvector(gam,1,n);
}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

//	u=vector(1,n-1);
	u=dvector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
//	free_vector(u,1,n-1);
	free_dvector(u,1,n-1);
}
#undef NRANSI


/*
 * This is a modified version so that we do not have to copy our data into
 * an array for every lookup. As a test we do a lookup in rho.
 */
/*
void splint(TILLMATERIAL *material, int v, double rho, double *u)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=material->nTableMax;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		// For a lookup along one isentrope k = i else k = j
		if (material->Lookup[INDEX(k, v)].rho > rho) khi=k;
		else klo=k;
	}
	// Modified this
	h=material->Lookup[INDEX(khi,v)].rho-material->Lookup[INDEX(klo,v)].rho;

	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(material->Lookup[INDEX(khi,v)].rho-rho)/h;
	b=(rho-material->Lookup[INDEX(klo,v)].rho)/h;
	// Careful in this version udv2 is d2u/drho2
	*u=a*material->Lookup[INDEX(klo,v)].u+b*material->Lookup[INDEX(khi,v)].u+((a*a*a-a)*material->Lookup[INDEX(klo,v)].udv2+(b*b*b-b)*material->Lookup[INDEX(khi,v)].rho)*(h*h)/6.0;
}
*/
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

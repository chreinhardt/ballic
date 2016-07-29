#include "../tillotson.h"

/*
 * This is a modified version so that we do not have to copy our data into
 * an array for every lookup. As a test we do a lookup in rho.
 */
/*
void splint(TILLMATERIAL *material, int v, float rho, float *u)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	float h,b,a;

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
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	float h,b,a;

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

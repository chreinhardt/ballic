#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../tillotson.h"

double brent(double (*func)(TILLMATERIAL *,double,double,double), TILLMATERIAL *material, double a, double b, double rho, double u, double tol,int iOrder)
{ int iter;
  double c,d,e,min1,min2,eps=1e-7;
  double fa,fb,fc,p,q,r,s,tol1,xm;
  fa = (*func)(material,a,rho,u);  fb = (*func)(material,b,rho,u);
  if (fa*fb > 0)
    {
		fprintf(stderr,"Particle: %i\n",iOrder);
		fprintf(stderr,"Root must be bracketed in zbrent\n");  exit(3);
    }
  c = b; fc = fb;
  for (iter=1; iter<=100; iter++)
    { if (fb*fc > 0)
        { c=a; fc=fa; e=d=b-a;
        }
      if (fabs(fc) < fabs(fb))
        { a=b; b=c; c=a;
          fa=fb; fb=fc; fc=fa;
        }
      tol1=2*eps*fabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0)
        return b;
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        { s=fb/fa;
          if (a == c) 
            { p=2.0*xm*s; q=1.0-s;
            } 
          else
            { q=fa/fc;  r=fb/fc;
              p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
              q=(q-1.0)*(r-1.0)*(s-1.0);
            }
          if (p > 0.0) q = -q;
          p=fabs(p);
          min1=3.0*xm*q-fabs(tol1*q);
          min2=fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2))
            { e=d; d=p/q;
            }
          else
            { d=xm; e=d;
            }
        }
      else
        { d=xm; e=d;
        }
      a=b; fa=fb;
      if (fabs(d) > tol1) b += d;
      else b += (xm >= 0.0 ? fabs(tol1) : -fabs(tol1));
      fb=(*func)(material,b,rho,u);
    }
  fprintf(stderr,"Maximum number of iterations exceeded in zbrent\n"); exit(0);
}



double cubicint(double u[2],double dudrho[2],double rho[2], rhoint) {
  dx = rho[1] - rho[0];
  e = (rhoint - rho[0])/dx;
  e1 = e - 1;

  // these are the 4 Hermite functions
  ce[0] = (2*e + 1)*e1*e1;
  ce[1] = e*e1*e1;
  ce[2] = e*e*(3 - 2*e);
  ce[3] = e*e*e1;

  /*
  **    = ce[0]*u(v,0) + ce[1]*dudrho(v,0)*dx + ce[2]*u(v,1) + ce[3]*dudrho(v,1);
  ** the above is written as 4 independent spline lookups in the table v lies between some j and j+1
  */
  return(ce[0]*u[0] + ce[1]*dudrho[0]*dx + ce[2]*u[1] + ce[3]*dudrho[1]*dx);
  return(ce[0]*dudv[0] + ce[1]*dudvdrho[0]*dx + ce[2]*dudv[1] + ce[3]*dudvdrho[1]*dx);
}


double u(rho,v) {
}

double dudv(rho,v) {
}

void uandudv(rho,v,double *u,double *dudv);

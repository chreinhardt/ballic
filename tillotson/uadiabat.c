/*
** Note: we could have done the spline fit along the rho axis instead of
** the v-axis and then done a piecewise cubic fit in the v axis. I don't
** know which is better, but generally we want the spline fit in the axis
** of less "polynomial" behaviour. As long as enough points are used in
** the ODE solution along rho, it should not matter. -Joachim Stadel
*/

typedef struct table_entry {
    double u;
    double u1;
    /*
    ** The following 2 variables are the second derivatives of the above
    ** variables (u and u1) with respect to v (which is the value of the
    ** constant entropy curve (adiabat) at rho_0). Both of these are obtained
    ** by fitting splines to u and u1 runs in the v axis.
    */
    double udv2;
    double u1dv2;
    } UADTE;


typedef struct ulookup_table {

    } * UAD;




/*
** This function sets up the table given the dudrho function.
** This is the right hand side of the differential equation that needs
** to be solved to find the values of u. It will also be directly
** tabulated. The step sizes hrho and hu0 are suggestions as the code
** will choose the next smallest value that will land exactly at rho=0
** and at u0=u0max respectively.
*/
void uadSetupTable(UAD *ctx,double dudrho(double u,double rho),
    double rho0,double hrho,double rhomax,double hu0,double u0max) {

    /* Solve ODE du/drho for each v */


    /* Solve splines for both u and u1 in v storing the 2nd derivatives wrt v */

    };

double u(double rhoStart,double uStart,double rhoEnd) {

    /* find i and i+1 splines which bound rhoStart */

    /* find v(uStart,rhoStart) using multiple interpolations using
       index i in a root finder */

    /* find k and k+1 splines which bound rhoEnd */
    /* find u using k and k+1 values at the determined v. */
    }


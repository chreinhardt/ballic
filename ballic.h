/*
 ** The header file for ballic.
 */
#ifndef BALLIC_HINCLUDED
#define BALLIC_HINCLUDED

#include "tipsy.h"
#include "tillotson/tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

// Isentropic thermal profile
#define BALLIC_U_ISENTROPIC

typedef struct icosa_struct {
    float R[180];
    float v[36];
    } ICOSA;

double Packed49[49][3] = {
    {1.0820379E-008,-3.2300459E-006,    0.7790824},
    { 0.0851419,      -0.0604194,       0.3497294},
    {-0.0914288,       0.4137149,       0.6537967},
    {-0.4233704,      -0.0167043,       0.6537949},
    { 0.4161716,       0.0795128,       0.6537953},
    { 0.2606476,      -0.3340458,       0.6537932},
    {-0.1783143,      -0.3843534,       0.6537934},
    { 0.3214161,       0.4881783,       0.5151146},
    {-0.4918671,       0.3793286,       0.4702615},
    {-0.0929470,       0.3254455,       0.2208710},
    {-0.5385916,      -0.3977282,       0.3983726},
    { 0.6544342,      -0.1767253,       0.3839966},
    {-0.0464433,      -0.6884808,       0.3616719},
    { 0.6543453 ,      0.2637903 ,      0.3304789},
    { 0.3839234,      -0.5971428,       0.3209247},
    {-0.7108776,       0.0012718,       0.3187802},
    {-0.1876019,      -0.3352532,       0.1368439},
    { 0.0704146,       0.7351211,       0.2481284},
    {-0.3638534,       0.6695231,       0.1622308},
    { 0.3033485,       0.2022102,       0.0692749},
    { 0.4742969,       0.6067250,       0.1178835},
    {-0.6795831,       0.3744220,       0.0703152},
    {-0.3794767,       0.0482692,       0.0303653},
    { 0.6538247,      -0.4232979,       0.0173632},
    { 0.2332925,      -0.2885923,       0.0015460},
    { 0.7790813, -5.7994478E-013,      -0.0012951},
    {-0.5429419,      -0.5585797,      -0.0131201},
    {-0.1452212,      -0.7628716,      -0.0625065},
    {-0.7541588,      -0.1768841,      -0.0832219},
    { 0.2928920,      -0.7156702,      -0.0948672},
    {-0.0815266,       0.7567586,      -0.1662504},
    { 0.6534047,       0.3539813,      -0.2339385},
    {-0.4662671,       0.5642231,      -0.2668646},
    { 0.3250845,       0.6429856,      -0.2964101},
    {-0.6822678,       0.1837015,      -0.3282282},
    { 0.4930282,      -0.4578927,      -0.3927173},
    {-0.3428409,      -0.5690840,      -0.4069064},
    { 0.6530941,      -0.0471244,      -0.4221572},
    {-0.1036092,       0.2867325,      -0.2191358},
    { 0.0858176,      -0.5795982,      -0.5134887},
    {-0.5513357,      -0.1947856,      -0.5148367},
    { 0.0032634,       0.5253046,      -0.5753380},
    {-0.1660916,      -0.1851457,      -0.2781839},
    { 0.3945018,       0.3205518,      -0.5904102},
    {-0.3693146,       0.2921726,      -0.6206539},
    { 0.2268184,       0.0121013,      -0.3221586},
    { 0.2750635,      -0.2203113,      -0.6948182},
    {-0.1633231,      -0.2729136,      -0.7112054},
    { 0.0133638,       0.1279994,      -0.7683794}
    };


typedef struct model_ctx {
	/* Material coefficients from the Tillotson EOS. */
	TILLMATERIAL **tillMat;
	int nLayer;
	/*
	** This array contains the material number for each layer (e.g., IRON, GRANITE).
	*/
	int *iLayer;	
	double *MLayer;		// Mass of each layer

	/*
	** Some unit conversion factors.
	*/
	double dKpcUnit;
	double dMsolUnit;

	/*
	** The lookup table for the equilibrium model.
	*/
	int nTableMax;
	int nTable;
	double uc; /* u at r = 0 */
	double *M;
	double *rho;
	double *u;
	double *r;
	int *mat;

	double dr;
	double R;
	} MODEL;

const double fact = 1.0;

/* Icosahedron. */
ICOSA *icosaInit(void);
void icosaPix2Vec(ICOSA *ctx,int i,int resolution,double *vec);
/* HEALPix. */
void pix2vec_ring( long nside, long ipix, double *vec);

MODEL *modelInit(double M,double ucore);
double drhodr(MODEL *model,int iLayer,double r,double rho,double M,double u);
double dudr(MODEL *model,int iLayer,double r,double rho,double M,double u);
double dudrho(MODEL *model,int iLayer,double rho,double u);
double dMdr(double r,double rho);

void modelSolveBC(MODEL *model, double *prho, double *pu, int iLayer1, int iLayer2);
void modelSolveComponent(MODEL *model,int iLayer,int bSetModel,int bLastLayer,int *pIndex,double h,double *prho1,double *pu1,double *pM1,double M2,double *pR);

double modelSolveTwoComponent(MODEL *model,int bSetModel,double rho,double u,double h,double *pR);
void modelWriteToFile(MODEL *model);

double modelSolveAll(MODEL *model,int bSetModel,double rhoc,double uc,double h,double *pR);

double modelSolve(MODEL *model,double M);

double MLookup(MODEL *model,double r);

double rhoLookup(MODEL *model,double r);

double uLookup(MODEL *model,double r);

/* Basic ballic routines. */
double Fzero(MODEL *model,int bIcosa,double r,double ri,double m,int ns);
double rShell(MODEL *model,double m,double ri);
double rShell2(MODEL *model,int bIcosa,double m,double ri,int ns);
#endif

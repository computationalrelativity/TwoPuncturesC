/* TwoPunctures.h 
   
   Files of "TwoPunctures":

   TwoPunctures.c
   TP_FuncAndJacobian.c
   TP_CoordTransf.c
   TP_Equations.c
   TP_Newton.c
   TP_utilities.c 
*/

#ifndef TP_HEADER_H
#define TP_HEADER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

/* Macros */
#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	/* Pi/2*/
#define Piq 0.78539816339744830961566084582	/* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

#define GFINDEX3D (i, j, k, n, nn) ( (i) + (n)*(j) + (nn)*(k) )
#define ERROR(s) {printf(s); exit(0);}
#define DBGSTOP(s) {printf("DEBUGSTOP: %s",s); exit(0);}

/* Various enum for options */

/* Interpolations options */
enum{
  taylor_expansion,// use a Taylor expansion about the nearest collocation point (fast, but might be inaccurate)
  evaluation, // evaluate using all spectral coefficients (slow)
  N_interp_opt,
};
static const char* str_interp_opt[N_interp_opt] = {
  "taylor", "spectral"
};

/* lapse options */
enum{
  antisymmetric, // antisymmetric lapse for two puncture black holes, -1 <= alpha <= +1
  averaged, // averaged lapse for two puncture black holes, 0 <= alpha <= +1"
  psin, // Based on the initial conformal factor
  brownsville, // See Phys. Rev. D 74, 041501 (2006)
  N_lapse_opt,
};
static const char* str_lapse_opt[N_lapse_opt] = {
  "antisymmetric", "averaged", "psin", "brownsville",
};

/* conformalfactor options */
enum{
  NONSTATIC, // conformal_state = 0 (?)
  FACTOR0, // CCTK_EQUALS(metric_type, "static conformal") && CCTK_EQUALS(conformal_storage, "factor")
  FACTOR1, // CCTK_EQUALS(metric_type, "static conformal") && CCTK_EQUALS(conformal_storage, "factor+derivs")
  FACTOR2, // CCTK_EQUALS(metric_type, "static conformal") && CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")
  N_cf_opt,
};
static const char* str_cf_opt[N_cf_opt] = {
  "nonstatic", "static0", "static01", "static012",
};

#define use_sources (0) //SB: do not add the source terms
#define rescale_sources (1) // If sources are used - rescale them after solving?

/* Structure with field and derivatives */
typedef struct DERIVS
{
  int size;
  double *d0; // field
  double *d1, *d2, *d3; // field 1st derivatives in directions 1,2,3
  double *d11, *d12, *d13, *d22, *d23, *d33; // field 2nd derivatives
} derivs;

/* params routines are in TP_Utilities.c */
#define NPARAMS (100)
#define STRLEN (256)
#define STREQL(s1,s2) ((strcmp((s1),(s2))==0))
enum{
  INTEGER,
  REAL,
  STRING,
  N_PAR_TYPES,
};
static const char* str_par_type[N_PAR_TYPES] = {
  "INTEGER", "REAL", "STRING", 
};

typedef struct {
  int type[NPARAMS];
  char key[NPARAMS][STRLEN];
  char val[NPARAMS][STRLEN];
  int n;
  //parameters *next;
} parameters;

/* Prototypes */

/* TP_CoordTransf.c */
void AB_To_XR (int nvar, double A, double B, double *X, double *R, derivs *U);
void C_To_c (int nvar, double X, double R, double *x, double *r, double par_b, derivs *U);
void rx3_To_xyz (int nvar, double x, double r, double phi, double *y, double *z, derivs *U);

/* TP_Equations.c */
double BY_KKofxyz (double x, double y, double z);
void BY_Aijofxyz (double x, double y, double z, double Aij[3][3]);
void NonLinEquations (double rho_adm,
		      double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs *U, double *values);
void LinEquations (double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs *dU, derivs *U, double *values);

/* TP_FuncAndJacobian.c */
int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs (derivs **v, int n);
void free_derivs (derivs *v);
void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs * v);
void F_of_v (int nvar, int n1, int n2, int n3, derivs * v,
	     double *F, derivs * u);
void J_times_dv (int nvar, int n1, int n2, int n3, derivs * dv,
		 double *Jdv, derivs * u);
void JFD_times_dv (int i, int j, int k, int nvar,
		   int n1, int n2, int n3,
		   derivs * dv, derivs * u, double *values);
void SetMatrix_JFD (int nvar, int n1, int n2, int n3,
		    derivs * u, int *ncols, int **cols, double **Matrix);
double PunctEvalAtArbitPosition (double *v, int ivar, double A, double B, double phi,
				 int nvar, int n1, int n2, int n3);
void calculate_derivs (int i, int j, int k, int ivar, int nvar,
		       int n1, int n2, int n3,
		       derivs * v, derivs * vv);
double interpol (double a, double b, double c, derivs * v);
double PunctTaylorExpandAtArbitPosition (int ivar, int nvar,
					 int n1, int n2, int n3,
					 derivs * v,
					 double x, double y, double z);
double PunctIntPolAtArbitPosition (int ivar, int nvar,
				   int n1, int n2, int n3,
				   derivs * v,
				   double x, double y, double z);
void SpecCoef(int n1, int n2, int n3, int ivar, double *v, double *cf);
double PunctEvalAtArbitPositionFast (double *v, int ivar, double A, double B, double phi,
				     int nvar, int n1, int n2, int n3);
double PunctIntPolAtArbitPositionFast (int ivar, int nvar, int n1,
				       int n2, int n3, derivs * v,
				       double x, double y, double z);
void set_initial_guess(derivs *v);

/* TP_Newton.c */
void Newton (int nvar, int n1, int n2, int n3, derivs *v,
	     double tol, int itmax);

/* TP_Utilies.c */
void nrerror (char error_text[]);
int *ivector (long nl, long nh);
double *dvector (long nl, long nh);
int **imatrix (long nrl, long nrh, long ncl, long nch);
double **dmatrix (long nrl, long nrh, long ncl, long nch);
double ***d3tensor (long nrl, long nrh, long ncl, long nch, long ndl,
		    long ndh);
void free_ivector (int *v, long nl, long nh);
void free_dvector (double *v, long nl, long nh);
void free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor (double ***t, long nrl, long nrh, long ncl, long nch,
		    long ndl, long ndh);

int minimum2 (int i, int j);
int minimum3 (int i, int j, int k);
int maximum2 (int i, int j);
int maximum3 (int i, int j, int k);
int pow_int (int mantisse, int exponent);

void chebft_Zeros (double u[], int n, int inv);
void chebft_Extremes (double u[], int n, int inv);
void chder (double *c, double *cder, int n);
double chebev (double a, double b, double c[], int m, double x);
void fourft (double *u, int N, int inv);
void fourder (double u[], double du[], int N);
void fourder2 (double u[], double d2u[], int N);
double fourev (double *u, int N, double x);

double norm1 (double *v, int n);
double norm2 (double *v, int n);
double scalarproduct (double *v, double *w, int n);

void params_alloc();
void params_free();
void params_read(char *fname);
void params_write(char * fname);
double params_get_real(char * key);
int params_get_int(char * key);
char *params_get_str(char * key);
void params_setadd(char * key, int type, char *val, int addpar);
void params_set_int(char * key, int val);
void params_set_bool(char * key, bool val);
void params_set_real(char * key, double val);
void params_set_str(char * key, char* val);
void params_add_int(char * key, int val);
void params_add_real(char * key, double val);
void params_add_str(char * key, char *val);

void make_output_dir();
void write_derivs(derivs *u, const int n1, const int n2, const int n3,
		  int include_derivatives_order, 
		  const char *fname);
void write_confact_atxyz(double x, double y, double z, double u,  
			 const char *fname);
void write_bam_inifile(derivs *u, const int n1, const int n2, const int n3,
		       const char *fname);

/* TwoPunctures.c */

// Expose for C++ intefacing
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    double *F;
    derivs *u;
    derivs *v;
    derivs *cf_v;
    int ntotal;
  } ini_data;

  // set default parameters
  void TwoPunctures_params_set_default();

  // set based on input file
  void TwoPunctures_params_set_inputfile(char *inputfile);
  void TwoPunctures_params_set_Real(char *key, double value);
  void TwoPunctures_params_set_Int(char *key, int value);
  void TwoPunctures_params_set_Boolean(char *key, bool value);

  ini_data * TwoPunctures_make_initial_data();

  void TwoPunctures_Cartesian_interpolation
  (ini_data *data,     // struct containing the previously calculated solution
   int *imin,         // min, max idxs of Cartesian Grid in the three directions
   int *imax,         // in the three dirs
   int *nxyz,
   double *x,         // Cartesian coordinates
   double *y,
   double *z,
   double *alp,       // lapse
   double *psi,       // conformal factor and derivatives
   double *psix,
   double *psiy,
   double *psiz,
   double *psixx,
   double *psixy,
   double *psixz,
   double *psiyy,
   double *psiyz,
   double *psizz,
   double *gxx,       // metric components
   double *gxy,
   double *gxz,
   double *gyy,
   double *gyz,
   double *gzz,
   double *kxx,       // extrinsic curvature components
   double *kxy,
   double *kxz,
   double *kyy,
   double *kyz,
   double *kzz);

  // for cleanup purposes
  void TwoPunctures_finalise(ini_data *data);

#ifdef __cplusplus
}
#endif

#endif /* TP_HEADER_H */

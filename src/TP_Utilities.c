/* TP_utilities.c */

/* Various shit from numerical recipes 
   Routines for managing params 
   Routines for output */

#include "TwoPunctures.h"

// static const char* str_interp_opt[N_interp_opt] = {
//   "taylor", "spectral"
// };

// static const char* str_lapse_opt[N_lapse_opt] = {
//   "antisymmetric", "averaged", "psin", "brownsville",
// };

// static const char* str_cf_opt[N_cf_opt] = {
//   "nonstatic", "static0", "static01", "static012",
// };

static const char* str_par_type[N_PAR_TYPES] = {
  "INTEGER", "REAL", "STRING", 
};

/* -------------------------------------------------------------------------*/

int *
TP_ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *retval;

  retval = malloc(sizeof(int)*(nh-nl+1));
  if(retval == NULL)
    ERROR ("allocation failure in TP_ivector()");
  return retval - nl;
}

double *
dvector (long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *retval;

  retval = malloc(sizeof(double)*(nh-nl+1));
  if(retval == NULL)
    ERROR ("allocation failure in dvector()");
  return retval - nl;
}

int **
imatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int **retval;

  retval = malloc(sizeof(int *)*(nrh-nrl+1));
  if(retval == NULL)
    ERROR ("allocation failure (1) in imatrix()");
  
  /* get all memory for the matrix in on chunk */
  retval[0] = malloc(sizeof(int)*(nrh-nrl+1)*(nch-ncl+1));
  if(retval[0] == NULL)
    ERROR ("allocation failure (2) in imatrix()");

  /* apply column and row offsets */
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows */
  long width = (nch-ncl+1);
  for(long i = nrl+1 ; i <= nrh ; i++)
    retval[i] = retval[i-1] + width;
  assert(retval[nrh]-retval[nrl] == (nrh-nrl)*width);

  return retval;
}

double **
dmatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  double **retval;

  retval = malloc(sizeof(double *)*(nrh-nrl+1));
  if(retval == NULL)
    ERROR ("allocation failure (1) in dmatrix()");

  /* get all memory for the matrix in on chunk */
  retval[0] = malloc(sizeof(double)*(nrh-nrl+1)*(nch-ncl+1));
  if(retval[0] == NULL)
    ERROR ("allocation failure (2) in dmatrix()");

  /* apply column and row offsets */
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows */
  long width = (nch-ncl+1);
  for(long i = nrl+1 ; i <= nrh ; i++)
    retval[i] = retval[i-1] + width;
  assert(retval[nrh]-retval[nrl] == (nrh-nrl)*width);

  return retval;
}

double ***
d3tensor (long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  double ***retval;

  /* get memory for index structures */
  retval = malloc(sizeof(double **)*(nrh-nrl+1));
  if(retval == NULL)
    ERROR ("allocation failure (1) in d3tensor()");

  retval[0] = malloc(sizeof(double *)*(nrh-nrl+1)*(nch-ncl+1));
  if(retval[0] == NULL)
    ERROR ("allocation failure (2) in d3tensor()");

  /* get all memory for the tensor in on chunk */
  retval[0][0] = malloc(sizeof(double)*(nrh-nrl+1)*(nch-ncl+1)*(ndh-ndl+1));
  if(retval[0][0] == NULL)
    ERROR ("allocation failure (3) in d3tensor()");

  /* apply all offsets */
  retval[0][0] -= ndl;
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows and columns */
  long width = (nch-ncl+1);
  long depth = (ndh-ndl+1);
  for(long j = ncl+1 ; j <= nch ; j++) { /* first row of columns */
    retval[nrl][j] = retval[nrl][j-1] + depth;
  }
  assert(retval[nrl][nch]-retval[nrl][ncl] == (nch-ncl)*depth);
  for(long i = nrl+1 ; i <= nrh ; i++) {
    retval[i] = retval[i-1] + width;
    retval[i][ncl] = retval[i-1][ncl] + width*depth; /* first cell in column */
    for(long j = ncl+1 ; j <= nch ; j++) {
      retval[i][j] = retval[i][j-1] + depth;
    }
    assert(retval[i][nch]-retval[i][ncl] == (nch-ncl)*depth);
  }
  assert(retval[nrh]-retval[nrl] == (nrh-nrl)*width);
  assert(&retval[nrh][nch][ndh]-&retval[nrl][ncl][ndl] == (nrh-nrl+1)*(nch-ncl+1)*(ndh-ndl+1)-1);

  return retval;
}

void
TP_free_ivector (int *v, long nl, long nh)
/* free an int vector allocated with TP_ivector() */
{
  free(v+nl);
}

void
free_dvector (double *v, long nl, long nh)
/* free an double vector allocated with dvector() */
{
  free(v+nl);
}

void
free_imatrix (int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free(m[nrl]+ncl);
  free(m+nrl);
}

void
free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free(m[nrl]+ncl);
  free(m+nrl);
}

void
free_d3tensor (double ***t, long nrl, long nrh, long ncl, long nch,
	       long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
  free(t[nrl][ncl]+ndl);
  free(t[nrl]+ncl);
  free(t+nrl);
}

int
minimum2 (int i, int j)
{
  int result = i;
  if (j < result)
    result = j;
  return result;
}

int
minimum3 (int i, int j, int k)
{
  int result = i;
  if (j < result)
    result = j;
  if (k < result)
    result = k;
  return result;
}

int
maximum2 (int i, int j)
{
  int result = i;
  if (j > result)
    result = j;
  return result;
}

int
maximum3 (int i, int j, int k)
{
  int result = i;
  if (j > result)
    result = j;
  if (k > result)
    result = k;
  return result;
}

int
pow_int (int mantisse, int exponent)
{
  int i, result = 1;
  
  for (i = 1; i <= exponent; i++)
    result *= mantisse;
  
  return result;
}

void
chebft_Zeros (double u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.4) of 2nd edition C++ NR */
{
  int k, j, isignum;
  double fac, sum, Pion, *c;
  
  c = dvector (0, n);
  Pion = Pi / n;
  if (inv == 0)
    {
      fac = 2.0 / n;
      isignum = 1;
      for (j = 0; j < n; j++)
	{
	  sum = 0.0;
	  for (k = 0; k < n; k++)
	    sum += u[k] * cos (Pion * j * (k + 0.5));
	  c[j] = fac * sum * isignum;
	  isignum = -isignum;
	}
    }
  else
    {
      for (j = 0; j < n; j++)
	{
	  sum = -0.5 * u[0];
	  isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * cos (Pion * (j + 0.5) * k) * isignum;
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
#if 0
    if (fabs(c[j]) < 5.e-16)
      u[j] = 0.0;
    else
#endif
      u[j] = c[j];
  free_dvector (c, 0, n);
}

void
chebft_Extremes (double u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.5) of 2nd edition C++ NR */
{
  int k, j, isignum, N = n - 1;
  double fac, sum, PioN, *c;

  c = dvector (0, N);
  PioN = Pi / N;
  if (inv == 0)
  {
    fac = 2.0 / N;
    isignum = 1;
    for (j = 0; j < n; j++)
    {
      sum = 0.5 * (u[0] + u[N] * isignum);
      for (k = 1; k < N; k++)
	sum += u[k] * cos (PioN * j * k);
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
    c[N] = 0.5 * c[N];
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * cos (PioN * j * k) * isignum;
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
    u[j] = c[j];
  free_dvector (c, 0, N);
}

void
chder (double *c, double *cder, int n)
{
  int j;

  cder[n] = 0.0;
  cder[n - 1] = 0.0;
  for (j = n - 2; j >= 0; j--)
    cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
}

double
chebev (double a, double b, double c[], int m, double x)
/* eq. 5.8.11 of C++ NR (2nd ed) */
{
  int j;
  double djp2, djp1, dj; /* d_{j+2}, d_{j+1} and d_j */
  double y;

  /* rescale input to lie within [-1,1] */
  y = 2*(x - 0.5*(b+a))/(b-a);

  dj = djp1 = 0;
  for(j = m-1 ; j >= 1; j--)
  {
    /* advance the coefficients */
    djp2 = djp1;
    djp1 = dj;
    dj   = 2*y*djp1 - djp2 + c[j];
  }

  return y*dj - djp1 + 0.5*c[0];
}

void
fourft (double *u, int N, int inv)
/* a (slow) Fourier transform, seems to be just eq. 12.1.6 and 12.1.9 of C++ NR (2nd ed) */
{
  int l, k, iy, M;
  double x, x1, fac, Pi_fac, *a, *b;
  
  // this check resolves warning that ‘a[0]’ and ‘a[M]’ may be used uninitialized in the inv==0 code block
  if (N < 0) ERROR("fourft: number of points must be >= 0");

  M = N / 2;
  a = dvector (0, M);
  b = dvector (1, M);		/* Actually: b=vector(1,M-1) but this is problematic if M=1*/
  fac = 1. / M;
  Pi_fac = Pi * fac;
  if (inv == 0)
  {
    for (l = 0; l <= M; l++)
    {
      a[l] = 0;
      if (l > 0 && l < M)
	b[l] = 0;
      x1 = Pi_fac * l;
      for (k = 0; k < N; k++)
      {
	x = x1 * k;
	a[l] += fac * u[k] * cos (x);
	if (l > 0 && l < M)
	  b[l] += fac * u[k] * sin (x);
      }
    }
    u[0] = a[0];
    u[M] = a[M];
    for (l = 1; l < M; l++)
    {
      u[l] = a[l];
      u[l + M] = b[l];
    }
  }
  else
  {
    a[0] = u[0];
    a[M] = u[M];
    for (l = 1; l < M; l++)
    {
      a[l] = u[l];
      b[l] = u[M + l];
    }
    iy = 1;
    for (k = 0; k < N; k++)
    {
      u[k] = 0.5 * (a[0] + a[M] * iy);
      x1 = Pi_fac * k;
      for (l = 1; l < M; l++)
      {
	x = x1 * l;
	u[k] += a[l] * cos (x) + b[l] * sin (x);
      }
      iy = -iy;
    }
  }
  free_dvector (a, 0, M);
  free_dvector (b, 1, M);
}

void
fourder (double u[], double du[], int N)
{
  int l, M, lpM;

  M = N / 2;
  du[0] = 0.;
  du[M] = 0.;
  for (l = 1; l < M; l++)
  {
    lpM = l + M;
    du[l] = u[lpM] * l;
    du[lpM] = -u[l] * l;
  }
}

void
fourder2 (double u[], double d2u[], int N)
{
  int l, l2, M, lpM;

  d2u[0] = 0.;
  M = N / 2;
  for (l = 1; l <= M; l++)
  {
    l2 = l * l;
    lpM = l + M;
    d2u[l] = -u[l] * l2;
    if (l < M)
      d2u[lpM] = -u[lpM] * l2;
  }
}

double
fourev (double *u, int N, double x)
{
  int l, M = N / 2;
  double xl, result;

  result = 0.5 * (u[0] + u[M] * cos (x * M));
  for (l = 1; l < M; l++)
  {
    xl = x * l;
    result += u[l] * cos (xl) + u[M + l] * sin (xl);
  }
  return result;
}

double
norm1 (double *v, int n)
{
  int i;
  double result = -1;

  for (i = 0; i < n; i++)
    if (fabs (v[i]) > result)
      result = fabs (v[i]);

  return result;
}

double
norm2 (double *v, int n)
{
  int i;
  double result = 0;
  
  for (i = 0; i < n; i++)
    result += v[i] * v[i];
  
  return sqrt (result);
}

double
scalarproduct (double *v, double *w, int n)
{
  int i;
  double result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * w[i];

  return result;
}

/* -------------------------------------------------------------------------*/

/* Dummy/minimal stuff for managing parameters
   - Params are stored as double and can be retrived only as int or real
   - Input file parser assumes input file is correctly written with lines
   'key=val\n'
   or
   '# comment'
*/

parameters *params;

/* alloc/free  */
void params_alloc()
{
  params = calloc(sizeof(parameters), 1);
  params->n=0;
}
void params_free()
{
  free(params);
}

/* parse parameter file */
void params_read(char *fname) {

#define verbose (0)

  FILE *fp;
  char line[STRLEN], key[STRLEN], val[STRLEN];
  char *s;

  fp = fopen(fname, "r");
  if (fp == NULL)
    ERROR("[params_read] Failed to open input file");
  while ( fgets(line,sizeof line,fp) != NULL ) {
    /* TODO: better parsing */
    if (line[0]=='#') continue;
    s = strtok(line, "=");
    strcpy(key,s);
    /* s = strtok(NULL, s); */
    s = strtok(0, "=");
    strcpy(val,s);
    strtok(val, "\n");
    if (verbose && 0) printf("|%s|%s|\n",key,val);
    // Do not add parameters from input file, make sure they exist already
    for (int i =0; i < params->n; i++) {
      if (STREQL(key,params->key[i])) {
	strcpy(params->val[i],val);
	if (verbose) printf("[params_read] read: %s = %s\n",key,val);
	break; // goto next line
      }
    }

  }
  fclose(fp);
}

/* dump parameters to file */
void params_write(char * fname)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (fp == NULL)
    ERROR("[params_write] Failed to open file\n");
  for (int i =0; i < params->n; i++) {
    fprintf(fp,"%s=%s\n",params->key[i],params->val[i]);
  }
  fclose(fp);
}

/* get double */
double params_get_real(char * key)
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==REAL) {
	return atof(params->val[i]);
      }  else
	ERROR("[params_get] parameter is not of type REAL\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* get int */
int params_get_int(char * key) 
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==INTEGER) 
	return atoi(params->val[i]);
      else
	ERROR("[params_get] parameter is not of type INTEGER\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* get string */
char * params_get_str(char * key)
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==STRING) 
	return params->val[i];
      else
	ERROR("[params_get] parameter is not of type REAL\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* set parameter */
void params_setadd(char * key, int type, char* val, int addpar)
{
#define check_unique (0) // addpar=1 assumes parameter is not there,
			 // activate this if you need to enforce this.
  
  if (addpar) {
#if (check_unique)
    //Do nothing if present
    for (int i =0; i < params->n; i++)
      if (strcmp(params->key[i],key) == 0) return;
#endif
    int i = params->n;
    if (i>=NPARAMS) ERROR("Increase NPARAMS.");
    if (type>=N_PAR_TYPES) ERROR("Unknown parameter type.");
    strcpy(params->key[i],key);
    strcpy(params->val[i],val);    
    params->type[i] = type;
    params->n++;
    if (0) printf ("Add parameter: %s TYPE=%s VALUE=%s (%d)\n",
		   key,str_par_type[type],val,params->n);
    return;
  }

  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      strcpy(params->val[i],val);
      return;
    }
  }

  ERROR("[params_set] parameter not found\n");
}

/* 'set' wrappers 
   Replace values, do not append new (assume exists)
   Note: addpar=0 does not set the type
*/

void params_set_int(char * key, int val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val);
  params_setadd(key, INTEGER, cval, 0); 
}
 
void params_set_bool(char * key, bool val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val?1:0);
  params_setadd(key, INTEGER, cval, 0); 
}

void params_set_real(char * key, double val) {
  char cval[STRLEN];
  sprintf(cval, "%.16e", val);
  params_setadd(key, REAL, cval, 0); 
}

void params_set_str(char * key, char *val) {
  params_setadd(key, STRING, val, 0);
}

/* 'add' wrappers 
   New parameter, append
 */

void params_add_int(char * key, int val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val);
  params_setadd(key, INTEGER, cval, 1);
}

void params_add_bool(char * key, bool val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val?1:0);
  params_setadd(key, INTEGER, cval, 1); 
}

void params_add_real(char * key, double val) {
  char cval[STRLEN];
  sprintf(cval, "%.16e", val);
  params_setadd(key, REAL, cval, 1);
}

void params_add_str(char * key, char *val) {
  params_setadd(key, STRING, val, 1);
}


/* -------------------------------------------------------------------------*/

void make_output_dir()
/* note this is unsafe */
{
  char s[STRLEN*4];
  sprintf(s,"mkdir -p %s",params_get_str("outputdir"));
  if (system(s)==-1) ERROR("System call failed to mkdir");
  strcpy(s,params_get_str("outputdir"));
  strcat (s,"/params.par");
  params_write(s);
}

void write_derivs(derivs *u, const int n1, const int n2, const int n3,
		  int include_derivatives_order, 
		  const char *fname)
/* output routine for derivs */
{

  char f[STRLEN*4];
  strcpy(f,params_get_str("outputdir"));
  strcat(f,"/");
  strcat(f,fname);

  FILE *fp;
  fp = fopen(f, "w");
  assert(fp);
  const int size = n1*n2*n3;
  fprintf(fp, "# %d %d %d\n",n1,n2,n3);
  if (include_derivatives_order==0)
    for (int i=0; i<size; i++)
      fprintf(fp, "%.16e\n",u->d0[i]);
  else if (include_derivatives_order==1)
    for (int i=0; i<size; i++)
      fprintf(fp, "%.16e %.16e %.16e %.16e\n",
	      u->d0[i],
	      u->d1[i],u->d2[i],u->d3[i]);
  else if (include_derivatives_order==2)
    for (int i=0; i<size; i++)
      fprintf(fp, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      u->d0[i],
	      u->d1[i],u->d2[i],u->d3[i], 
	      u->d11[i],u->d12[i],u->d13[i], u->d22[i],u->d23[i],u->d33[i]);
  else  
    ERROR ("include_derivatives_order = 0,1,2\n");
  fclose(fp);
}

void write_confact_atxyz(double x, double y, double z, double u, 
			 const char *fname)
{
  double par_b = params_get_real("par_b");
  double par_m_plus = params_get_real("par_m_plus");
  double par_m_minus = params_get_real("par_m_minus");
  
  double r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  double r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  double psi = 1.+ 0.5 * par_m_plus  / r_plus + 0.5 * par_m_minus / r_minus + u;

  char f[STRLEN*4];
  strcpy(f,params_get_str("outputdir"));
  strcat(f,"/");
  strcat(f,fname);
  
  FILE *fp;
  fp = fopen(f, "w");
  assert(fp);
  fprintf(fp, "%.16e %.16e %.16e %.16e \n", x, y, z, psi);
  fclose(fp);
}

void write_bam_inifile(derivs *u, const int n1, const int n2, const int n3,
		       const char *fname)
/* output routine for derivs in BAM format */
{

  char f[STRLEN*4];
  strcpy(f,params_get_str("outputdir"));
  strcat(f,"/");
  strcat(f,fname);

  FILE *fp;
  fp = fopen(f, "w");
  assert(fp);
  const int size = n1*n2*n3;
  fprintf(fp, "# TwoPuncture.c\n");
  fprintf(fp, "# Black hole puncture data produced for bam's punctures_ps\n");
  fprintf(fp, "# Newton_tol = %e\n",params_get_real("Newton_tol"));
  fprintf(fp, "# nx = %d\n",n1);
  fprintf(fp, "# ny = %d\n",n2);
  fprintf(fp, "# nz = %d\n",n3);
  fprintf(fp, "bhmass1 = %e\n",params_get_real("par_m_plus"));
  fprintf(fp, "bhx1 = %e\n",params_get_real("par_b")); // FIXME
  fprintf(fp, "bhy1 = %e\n",0.); //
  fprintf(fp, "bhz1 = %e\n",0.); //
  fprintf(fp, "bhpx1 = %e\n",params_get_real("par_P_plus1")); // FIXME: might need rotation
  fprintf(fp, "bhpy1 = %e\n",params_get_real("par_P_plus2"));
  fprintf(fp, "bhpz1 = %e\n",params_get_real("par_P_plus3"));
  fprintf(fp, "bhsx1 = %e\n",params_get_real("par_S_plus1")); // FIXME: might need rotation
  fprintf(fp, "bhsy1 = %e\n",params_get_real("par_S_plus2"));
  fprintf(fp, "bhsz1 = %e\n",params_get_real("par_S_plus3"));
  fprintf(fp, "bhmass2 = %e\n",params_get_real("par_m_minus"));
  fprintf(fp, "bhx2 = -%e\n",params_get_real("par_b")); // FIXME
  fprintf(fp, "bhy2 = %e\n",0.); //
  fprintf(fp, "bhz2 = %e\n",0.); //
  fprintf(fp, "bhpx2 = %e\n",params_get_real("par_P_minus1")); // FIXME: might need rotation
  fprintf(fp, "bhpy2 = %e\n",params_get_real("par_P_minus2"));
  fprintf(fp, "bhpz2 = %e\n",params_get_real("par_P_minus3"));
  fprintf(fp, "bhsx2 = %e\n",params_get_real("par_S_minus1")); // FIXME: might need rotation
  fprintf(fp, "bhsy2 = %e\n",params_get_real("par_S_minus2"));
  fprintf(fp, "bhsz2 = %e\n",params_get_real("par_S_minus3"));
  fprintf(fp,"#\n"); // FIXME: Are the following data use dby BAM?
		     //  some must be calculated.
  fprintf(fp,"# d = %e\n",0.); // FIXME
  fprintf(fp,"# m1 = %e\n",params_get_real("par_m_plus"));
  fprintf(fp,"# m2 = %e\n",params_get_real("par_m_minus"));
  fprintf(fp,"# P1 = %e\n",0.);
  fprintf(fp,"# P2 = %e\n",0.);
  fprintf(fp,"# S1 = %e\n",0.);
  fprintf(fp,"# S2 = %e\n",0.);
  fprintf(fp,"#\n");
  fprintf(fp,"# M1 = %e\n",0.); // FIXME
  fprintf(fp,"# M2 = %e\n",0.);
  fprintf(fp,"# Madm = %e\n",0.);
  fprintf(fp,"# deltaMadm = %e\n",0.);
  fprintf(fp,"# Ebind = %e\n",0.); // FIXME
  fprintf(fp,"#\n");
  fprintf(fp,"# d/Madm = %e\n",0.);
  fprintf(fp,"# M1/Madm = %e\n",0.);
  fprintf(fp,"# M2/Madm = %e\n",0.);
  fprintf(fp,"# P1/Madm = %e\n",0.);
  fprintf(fp,"# P2/Madm = %e\n",0.);
  fprintf(fp,"# S1/Madm = %e\n",0.);
  fprintf(fp,"# S2/Madm = %e\n",0.);
  fprintf(fp,"# S1/M1^2 = %e\n",0.);
  fprintf(fp,"# S2/M2^2 = %e\n",0.);
  fprintf(fp,"#\n");
  fprintf(fp,"data %d %d %d\n",n1,n2,n3);
  for (int i=0; i<size; i++)
    fprintf(fp, "%.16e\n",u->d0[i]);
  fclose(fp);
}

/* output routine for ini_data */
//TODO

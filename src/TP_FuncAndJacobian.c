/* TP_FuncAndJacobian.c */

#include "TwoPunctures.h"

#define FAC (sin(al)*sin(be)*sin(al)*sin(be)*sin(al)*sin(be))

static inline double min (double const x, double const y)
{
  return x<y ? x : y;
}

int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3)
{
  int i1 = i, j1 = j, k1 = k;

  if (i1 < 0)
    i1 = -(i1 + 1);
  if (i1 >= n1)
    i1 = 2 * n1 - (i1 + 1);

  if (j1 < 0)
    j1 = -(j1 + 1);
  if (j1 >= n2)
    j1 = 2 * n2 - (j1 + 1);

  if (k1 < 0)
    k1 = k1 + n3;
  if (k1 >= n3)
    k1 = k1 - n3;

  return ivar + nvar * (i1 + n1 * (j1 + n2 * k1));
}

void allocate_derivs (derivs **v, const int n)
{
  *v = (derivs *) calloc(1, sizeof(derivs)); 
  if (v == NULL) ERROR("Out of memory");
  const int m = n - 1;
  (*v)->size = n;
  (*v)->d0 = dvector (0, m);
  (*v)->d1 = dvector (0, m);
  (*v)->d2 = dvector (0, m);
  (*v)->d3 = dvector (0, m);
  (*v)->d11 = dvector (0, m);
  (*v)->d12 = dvector (0, m);
  (*v)->d13 = dvector (0, m);
  (*v)->d22 = dvector (0, m);
  (*v)->d23 = dvector (0, m);
  (*v)->d33 = dvector (0, m);
}

void free_derivs (derivs * v)
{
  const int m = v->size;
  free_dvector (v->d0, 0, m);
  free_dvector (v->d1, 0, m);
  free_dvector (v->d2, 0, m);
  free_dvector (v->d3, 0, m);
  free_dvector (v->d11, 0, m);
  free_dvector (v->d12, 0, m);
  free_dvector (v->d13, 0, m);
  free_dvector (v->d22, 0, m);
  free_dvector (v->d23, 0, m);
  free_dvector (v->d33, 0, m);
  free(v);
}

void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs *v)
{
  int i, j, k, ivar, N, *indx;
  double *p, *dp, *d2p, *q, *dq, *r, *dr;

  N = maximum3 (n1, n2, n3);
  p = dvector (0, N);
  dp = dvector (0, N);
  d2p = dvector (0, N);
  q = dvector (0, N);
  dq = dvector (0, N);
  r = dvector (0, N);
  dr = dvector (0, N);
  indx = ivector (0, N);

  for (ivar = 0; ivar < nvar; ivar++)
  {
    for (k = 0; k < n3; k++)
    {				/* Calculation of Derivatives w.r.t. A-Dir. */
      for (j = 0; j < n2; j++)
      {				/* (Chebyshev_Zeros)*/
	for (i = 0; i < n1; i++)
	  {
	  indx[i] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[i] = v->d0[indx[i]];
	  }
	chebft_Zeros (p, n1, 0);
	chder (p, dp, n1);
	chder (dp, d2p, n1);
	chebft_Zeros (dp, n1, 1);
	chebft_Zeros (d2p, n1, 1);
	for (i = 0; i < n1; i++)
	{
	  v->d1[indx[i]] = dp[i];
	  v->d11[indx[i]] = d2p[i];
	}
      }
    }
    for (k = 0; k < n3; k++)
    {				/* Calculation of Derivatives w.r.t. B-Dir. */
      for (i = 0; i < n1; i++)
      {				/* (Chebyshev_Zeros)*/
	for (j = 0; j < n2; j++)
	{
	  indx[j] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[j] = v->d0[indx[j]];
	  q[j] = v->d1[indx[j]];
	}
	chebft_Zeros (p, n2, 0);
	chebft_Zeros (q, n2, 0);
	chder (p, dp, n2);
	chder (dp, d2p, n2);
	chder (q, dq, n2);
	chebft_Zeros (dp, n2, 1);
	chebft_Zeros (d2p, n2, 1);
	chebft_Zeros (dq, n2, 1);
	for (j = 0; j < n2; j++)
	{
	  v->d2[indx[j]] = dp[j];
	  v->d22[indx[j]] = d2p[j];
	  v->d12[indx[j]] = dq[j];
	}
      }
    }
    for (i = 0; i < n1; i++)
    {				/* Calculation of Derivatives w.r.t. phi-Dir. (Fourier)*/
      for (j = 0; j < n2; j++)
      {
	for (k = 0; k < n3; k++)
	{
	  indx[k] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[k] = v->d0[indx[k]];
	  q[k] = v->d1[indx[k]];
	  r[k] = v->d2[indx[k]];
	}
	fourft (p, n3, 0);
	fourder (p, dp, n3);
	fourder2 (p, d2p, n3);
	fourft (dp, n3, 1);
	fourft (d2p, n3, 1);
	fourft (q, n3, 0);
	fourder (q, dq, n3);
	fourft (dq, n3, 1);
	fourft (r, n3, 0);
	fourder (r, dr, n3);
	fourft (dr, n3, 1);
	for (k = 0; k < n3; k++)
	{
	  v->d3[indx[k]] = dp[k];
	  v->d33[indx[k]] = d2p[k];
	  v->d13[indx[k]] = dq[k];
	  v->d23[indx[k]] = dr[k];
	}
      }
    }
  }
  free_dvector (p, 0, N);
  free_dvector (dp, 0, N);
  free_dvector (d2p, 0, N);
  free_dvector (q, 0, N);
  free_dvector (dq, 0, N);
  free_dvector (r, 0, N);
  free_dvector (dr, 0, N);
  free_ivector (indx, 0, N);
}

void F_of_v (int nvar, int n1, int n2, int n3, derivs *v, double *F,
	     derivs *u)
{
  /* Calculates the left hand sides of the non-linear equations F_m(v_n)=0
     and the function u (u->d0[]) as well as its derivatives
     (u->d1[], u->d2[], u->d3[], u->d11[], u->d12[], u->d13[], u->d22[], u->d23[], u->d33[])
     at interior points and at the boundaries "+/-" 
  */

  double par_b = params_get_real("par_b");
  double par_m_plus = params_get_real("par_m_plus");
  double par_m_minus = params_get_real("par_m_minus");
  
  int i, j, k, ivar, indx;
  double al, be, A, B, X, R, x, r, phi, y, z, Am1, *values;
  derivs *U;

  values = dvector (0, nvar - 1);
  allocate_derivs (&U, nvar);

  double *sources;  
  sources = calloc(n1*n2*n3, sizeof(double));
  
#if (0) //SB: Keep sources = 0. , switch off code below
  if (use_sources)
  {
    double *s_x, *s_y, *s_z;
    int i3D;
    s_x    =calloc(n1*n2*n3, sizeof(double));
    s_y    =calloc(n1*n2*n3, sizeof(double));
    s_z    =calloc(n1*n2*n3, sizeof(double));
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          i3D = Index(0,i,j,k,1,n1,n2,n3);

          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;

          Am1 = A - 1;
          for (ivar = 0; ivar < nvar; ivar++)
          {
            indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
            U->d0[ivar] = Am1 * v->d0[indx];        /* U*/
            U->d1[ivar] = v->d0[indx] + Am1 * v->d1[indx];        /* U_A*/
            U->d2[ivar] = Am1 * v->d2[indx];        /* U_B*/
            U->d3[ivar] = Am1 * v->d3[indx];        /* U_3*/
            U->d11[ivar] = 2 * v->d1[indx] + Am1 * v->d11[indx];        /* U_AA*/
            U->d12[ivar] = v->d2[indx] + Am1 * v->d12[indx];        /* U_AB*/
            U->d13[ivar] = v->d3[indx] + Am1 * v->d13[indx];        /* U_AB*/
            U->d22[ivar] = Am1 * v->d22[indx];        /* U_BB*/
            U->d23[ivar] = Am1 * v->d23[indx];        /* U_B3*/
            U->d33[ivar] = Am1 * v->d33[indx];        /* U_33*/
          }
          /* Calculation of (X,R) and*/
          /* (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33)*/
          AB_To_XR (nvar, A, B, &X, &R, U);
          /* Calculation of (x,r) and*/
          /* (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33)*/
          C_To_c (nvar, X, R, &(s_x[i3D]), &r, par_b, U);
          /* Calculation of (y,z) and*/
          /* (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz)*/
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
        }
    Set_Rho_ADM(n1*n2*n3, sources, s_x, s_y, s_z);
    free(s_z);
    free(s_y);
    free(s_x);
  }
  else
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
          sources[Index(0,i,j,k,1,n1,n2,n3)]=0.0;
#endif //SB: Set sources = 0.
    
  Derivatives_AB3 (nvar, n1, n2, n3, v);

  double psi, psi2, psi4, psi7, r_plus, r_minus;

  char fname[STRLEN];
  FILE *debugfile = NULL;
  if (params_get_int("do_residuum_debug_output")) {
    strcpy(fname,params_get_str("outputdir"));
    strcat (fname,"/res.dat");
    debugfile = fopen(fname, "w");
    assert(debugfile);
  }
  
  for (i = 0; i < n1; i++)
    {
      for (j = 0; j < n2; j++)
	{
	  for (k = 0; k < n3; k++)
	    {
	      
	      al = Pih * (2 * i + 1) / n1;
	      A = -cos (al);
	      be = Pih * (2 * j + 1) / n2;
	      B = -cos (be);
	      phi = 2. * Pi * k / n3;
	      
	      Am1 = A - 1;
	      for (ivar = 0; ivar < nvar; ivar++)
		{
		  indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
		  U->d0[ivar] = Am1 * v->d0[indx];        /* U*/
		  U->d1[ivar] = v->d0[indx] + Am1 * v->d1[indx];        /* U_A*/
		  U->d2[ivar] = Am1 * v->d2[indx];        /* U_B*/
		  U->d3[ivar] = Am1 * v->d3[indx];        /* U_3*/
		  U->d11[ivar] = 2 * v->d1[indx] + Am1 * v->d11[indx];        /* U_AA*/
		  U->d12[ivar] = v->d2[indx] + Am1 * v->d12[indx];        /* U_AB*/
		  U->d13[ivar] = v->d3[indx] + Am1 * v->d13[indx];        /* U_AB*/
		  U->d22[ivar] = Am1 * v->d22[indx];        /* U_BB*/
		  U->d23[ivar] = Am1 * v->d23[indx];        /* U_B3*/
		  U->d33[ivar] = Am1 * v->d33[indx];        /* U_33*/
		}
	      
	      /* Calculation of (X,R) and*/
	      /* (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33)*/
	      AB_To_XR (nvar, A, B, &X, &R, U);

	      /* Calculation of (x,r) and*/
	      /* (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33)*/
	      C_To_c (nvar, X, R, &x, &r, par_b, U);

	      /* Calculation of (y,z) and*/
	      /* (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz)*/
	      rx3_To_xyz (nvar, x, r, phi, &y, &z, U);

	      NonLinEquations (sources[Index(0,i,j,k,1,n1,n2,n3)],
			       A, B, X, R, x, r, phi, y, z, U, values);
	      
	      for (ivar = 0; ivar < nvar; ivar++)
		{
		  indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
		  F[indx] = values[ivar] * FAC;
		  u->d0[indx] = U->d0[ivar];      
		  u->d1[indx] = U->d1[ivar];      
		  u->d2[indx] = U->d2[ivar];      
		  u->d3[indx] = U->d3[ivar];      
		  u->d11[indx] = U->d11[ivar];    
		  u->d12[indx] = U->d12[ivar];    
		  u->d13[indx] = U->d13[ivar];    
		  u->d22[indx] = U->d22[ivar];    
		  u->d23[indx] = U->d23[ivar];    
		  u->d33[indx] = U->d33[ivar];    
		}

	      
	      if (debugfile && (k==0)) {
		r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
		r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
		psi = 1.+
		  0.5 * par_m_plus  / r_plus +
		  0.5 * par_m_minus / r_minus +
		  U->d0[0];
		psi2 = psi * psi;
		psi4 = psi2 * psi2;
		psi7 = psi * psi2 * psi4;
		fprintf(debugfile,
			"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
			(double)x, (double)y, (double)A, (double)B,
			(double)(U->d11[0] +
				 U->d22[0] +
				 U->d33[0] +
				   /*                      0.125 * BY_KKofxyz (x, y, z) / psi7 +*/
				 (2.0 * Pi / psi2/psi * sources[indx]) * FAC),
			(double)((U->d11[0] +
				  U->d22[0] +
				  U->d33[0])*FAC),
			(double)(-(2.0 * Pi / psi2/psi * sources[indx]) * FAC),
			(double)sources[indx]
			/*(double)F[indx]*/
			);
	      }// debug file
	      
	    } // for k
	} // for j
    } // for i
    
  if (debugfile)
    fclose(debugfile);
  
  free(sources);
  free_dvector (values, 0, nvar - 1);
  free_derivs (U);
}

void J_times_dv (int nvar, int n1, int n2, int n3, derivs *dv,
		 double *Jdv, derivs *u)
{
  /* Calculates the left hand sides of the non-linear equations F_m(v_n)=0
     and the function u (u->d0[]) as well as its derivatives
     (u->d1[], u->d2[], u->d3[], u->d11[], u->d12[], u->d13[], u->d22[], u->d23[], u->d33[])
     at interior points and at the boundaries "+/-" 
  */

  double par_b = params_get_real("par_b");
  
  int i, j, k, ivar, indx;
  double al, be, A, B, X, R, x, r, phi, y, z, Am1, *values;
  derivs *dU, *U;
  
  Derivatives_AB3 (nvar, n1, n2, n3, dv);
  
#ifdef TP_OMP
#pragma omp parallel for private (values,dU,U,i,j,k,al,A,be,B,phi,X,R,x,r,y,z,Am1,ivar,indx,par_b) schedule(dynamic)
#endif
  for (i = 0; i < n1; i++)
    {
      values = dvector (0, nvar - 1);
      allocate_derivs (&dU, nvar);
      allocate_derivs (&U, nvar);
      for (j = 0; j < n2; j++)
	{
	  for (k = 0; k < n3; k++)
	    {
	      
	      al = Pih * (2 * i + 1) / n1;
	      A = -cos (al);
	      be = Pih * (2 * j + 1) / n2;
	      B = -cos (be);
	      phi = 2. * Pi * k / n3;
	      
	      Am1 = A - 1;
	      for (ivar = 0; ivar < nvar; ivar++)
		{
		  indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
		  dU->d0[ivar] = Am1 * dv->d0[indx];	/* dU*/
		  dU->d1[ivar] = dv->d0[indx] + Am1 * dv->d1[indx];	/* dU_A*/
		  dU->d2[ivar] = Am1 * dv->d2[indx];	/* dU_B*/
		  dU->d3[ivar] = Am1 * dv->d3[indx];	/* dU_3*/
		  dU->d11[ivar] = 2 * dv->d1[indx] + Am1 * dv->d11[indx];	/* dU_AA*/
		  dU->d12[ivar] = dv->d2[indx] + Am1 * dv->d12[indx];	/* dU_AB*/
		  dU->d13[ivar] = dv->d3[indx] + Am1 * dv->d13[indx];	/* dU_AB*/
		  dU->d22[ivar] = Am1 * dv->d22[indx];	/* dU_BB*/
		  dU->d23[ivar] = Am1 * dv->d23[indx];	/* dU_B3*/
		  dU->d33[ivar] = Am1 * dv->d33[indx];	/* dU_33*/
		  U->d0[ivar] = u->d0[indx];	/* U   */
		  U->d1[ivar] = u->d1[indx];	/* U_x*/
		  U->d2[ivar] = u->d2[indx];	/* U_y*/
		  U->d3[ivar] = u->d3[indx];	/* U_z*/
		  U->d11[ivar] = u->d11[indx];	/* U_xx*/
		  U->d12[ivar] = u->d12[indx];	/* U_xy*/
		  U->d13[ivar] = u->d13[indx];	/* U_xz*/
		  U->d22[ivar] = u->d22[indx];	/* U_yy*/
		  U->d23[ivar] = u->d23[indx];	/* U_yz*/
		  U->d33[ivar] = u->d33[indx];	/* U_zz*/
		}
	      
	      /* Calculation of (X,R) and*/
	      /* (dU_X, dU_R, dU_3, dU_XX, dU_XR, dU_X3, dU_RR, dU_R3, dU_33)*/
	      AB_To_XR (nvar, A, B, &X, &R, dU);

	      /* Calculation of (x,r) and*/
	      /* (dU, dU_x, dU_r, dU_3, dU_xx, dU_xr, dU_x3, dU_rr, dU_r3, dU_33)*/
	      C_To_c (nvar, X, R, &x, &r, par_b, dU);
	      
	      /* Calculation of (y,z) and*/
	      /* (dU, dU_x, dU_y, dU_z, dU_xx, dU_xy, dU_xz, dU_yy, dU_yz, dU_zz)*/
	      rx3_To_xyz (nvar, x, r, phi, &y, &z, dU);
	      LinEquations (A, B, X, R, x, r, phi, y, z, dU, U, values);
	      for (ivar = 0; ivar < nvar; ivar++)
		{
		  indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
		  Jdv[indx] = values[ivar] * FAC;
		}
	    }
	}
      
      free_dvector (values, 0, nvar - 1);
      free_derivs (dU);
      free_derivs (U);
    }
}

void JFD_times_dv (int i, int j, int k, int nvar, int n1, int n2,
		   int n3, derivs *dv, derivs *u, double *values)
{
  /* Calculates rows of the vector 'J(FD)*dv'.
     First row to be calculated: row = Index(0,      i, j, k; nvar, n1, n2, n3)
     Last  row to be calculated: row = Index(nvar-1, i, j, k; nvar, n1, n2, n3)
     These rows are stored in the vector JFDdv[0] ... JFDdv[nvar-1]. */

  double par_b = params_get_real("par_b");
  
  int ivar, indx;
  double al, be, A, B, X, R, x, r, phi, y, z, Am1;
  double sin_al, sin_al_i1, sin_al_i2, sin_al_i3, cos_al;
  double sin_be, sin_be_i1, sin_be_i2, sin_be_i3, cos_be;
  double dV0, dV1, dV2, dV3, dV11, dV12, dV13, dV22, dV23, dV33,
    ha, ga, ga2, hb, gb, gb2, hp, gp, gp2, gagb, gagp, gbgp;
  derivs *dU, *U;

  allocate_derivs (&dU, nvar);
  allocate_derivs (&U, nvar);
  
  if (k < 0)
    k = k + n3;
  if (k >= n3)
    k = k - n3;
  
  ha = Pi / n1;			/* ha: Stepsize with respect to (al)*/
  al = ha * (i + 0.5);
  A = -cos (al);
  ga = 1 / ha;
  ga2 = ga * ga;
  
  hb = Pi / n2;			/* hb: Stepsize with respect to (be)*/
  be = hb * (j + 0.5);
  B = -cos (be);
  gb = 1 / hb;
  gb2 = gb * gb;
  gagb = ga * gb;

  hp = 2 * Pi / n3;		/* hp: Stepsize with respect to (phi)*/
  phi = hp * k;
  gp = 1 / hp;
  gp2 = gp * gp;
  gagp = ga * gp;
  gbgp = gb * gp;


  sin_al = sin (al);
  sin_be = sin (be);
  sin_al_i1 = 1 / sin_al;
  sin_be_i1 = 1 / sin_be;
  sin_al_i2 = sin_al_i1 * sin_al_i1;
  sin_be_i2 = sin_be_i1 * sin_be_i1;
  sin_al_i3 = sin_al_i1 * sin_al_i2;
  sin_be_i3 = sin_be_i1 * sin_be_i2;
  cos_al = -A;
  cos_be = -B;

  Am1 = A - 1;
  for (ivar = 0; ivar < nvar; ivar++)
    {
      int iccc = Index (ivar, i, j, k, nvar, n1, n2, n3),
	ipcc = Index (ivar, i + 1, j, k, nvar, n1, n2, n3),
	imcc = Index (ivar, i - 1, j, k, nvar, n1, n2, n3),
	icpc = Index (ivar, i, j + 1, k, nvar, n1, n2, n3),
	icmc = Index (ivar, i, j - 1, k, nvar, n1, n2, n3),
	iccp = Index (ivar, i, j, k + 1, nvar, n1, n2, n3),
	iccm = Index (ivar, i, j, k - 1, nvar, n1, n2, n3),
	icpp = Index (ivar, i, j + 1, k + 1, nvar, n1, n2, n3),
	icmp = Index (ivar, i, j - 1, k + 1, nvar, n1, n2, n3),
	icpm = Index (ivar, i, j + 1, k - 1, nvar, n1, n2, n3),
	icmm = Index (ivar, i, j - 1, k - 1, nvar, n1, n2, n3),
	ipcp = Index (ivar, i + 1, j, k + 1, nvar, n1, n2, n3),
	imcp = Index (ivar, i - 1, j, k + 1, nvar, n1, n2, n3),
	ipcm = Index (ivar, i + 1, j, k - 1, nvar, n1, n2, n3),
	imcm = Index (ivar, i - 1, j, k - 1, nvar, n1, n2, n3),
	ippc = Index (ivar, i + 1, j + 1, k, nvar, n1, n2, n3),
	impc = Index (ivar, i - 1, j + 1, k, nvar, n1, n2, n3),
	ipmc = Index (ivar, i + 1, j - 1, k, nvar, n1, n2, n3),
	immc = Index (ivar, i - 1, j - 1, k, nvar, n1, n2, n3);
      /* Derivatives of (dv) w.r.t. (al,be,phi):*/
      dV0 = dv->d0[iccc];
      dV1 = 0.5 * ga * (dv->d0[ipcc] - dv->d0[imcc]);
      dV2 = 0.5 * gb * (dv->d0[icpc] - dv->d0[icmc]);
      dV3 = 0.5 * gp * (dv->d0[iccp] - dv->d0[iccm]);
      dV11 = ga2 * (dv->d0[ipcc] + dv->d0[imcc] - 2 * dv->d0[iccc]);
      dV22 = gb2 * (dv->d0[icpc] + dv->d0[icmc] - 2 * dv->d0[iccc]);
      dV33 = gp2 * (dv->d0[iccp] + dv->d0[iccm] - 2 * dv->d0[iccc]);
      dV12 =
	0.25 * gagb * (dv->d0[ippc] - dv->d0[ipmc] + dv->d0[immc] - dv->d0[impc]);
      dV13 =
	0.25 * gagp * (dv->d0[ipcp] - dv->d0[imcp] + dv->d0[imcm] - dv->d0[ipcm]);
      dV23 =
	0.25 * gbgp * (dv->d0[icpp] - dv->d0[icpm] + dv->d0[icmm] - dv->d0[icmp]);
      /* Derivatives of (dv) w.r.t. (A,B,phi):*/
      dV11 = sin_al_i3 * (sin_al * dV11 - cos_al * dV1);
      dV12 = sin_al_i1 * sin_be_i1 * dV12;
      dV13 = sin_al_i1 * dV13;
      dV22 = sin_be_i3 * (sin_be * dV22 - cos_be * dV2);
      dV23 = sin_be_i1 * dV23;
      dV1 = sin_al_i1 * dV1;
      dV2 = sin_be_i1 * dV2;
      /* Derivatives of (dU) w.r.t. (A,B,phi):*/
      dU->d0[ivar] = Am1 * dV0;
      dU->d1[ivar] = dV0 + Am1 * dV1;
      dU->d2[ivar] = Am1 * dV2;
      dU->d3[ivar] = Am1 * dV3;
      dU->d11[ivar] = 2 * dV1 + Am1 * dV11;
      dU->d12[ivar] = dV2 + Am1 * dV12;
      dU->d13[ivar] = dV3 + Am1 * dV13;
      dU->d22[ivar] = Am1 * dV22;
      dU->d23[ivar] = Am1 * dV23;
      dU->d33[ivar] = Am1 * dV33;
      
      indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
      U->d0[ivar] = u->d0[indx];	/* U   */
      U->d1[ivar] = u->d1[indx];	/* U_x*/
      U->d2[ivar] = u->d2[indx];	/* U_y*/
      U->d3[ivar] = u->d3[indx];	/* U_z*/
      U->d11[ivar] = u->d11[indx];	/* U_xx*/
      U->d12[ivar] = u->d12[indx];	/* U_xy*/
      U->d13[ivar] = u->d13[indx];	/* U_xz*/
      U->d22[ivar] = u->d22[indx];	/* U_yy*/
      U->d23[ivar] = u->d23[indx];	/* U_yz*/
      U->d33[ivar] = u->d33[indx];	/* U_zz*/
    }
  
  /* Calculation of (X,R) and*/
  /* (dU_X, dU_R, dU_3, dU_XX, dU_XR, dU_X3, dU_RR, dU_R3, dU_33)*/
  AB_To_XR (nvar, A, B, &X, &R, dU);

  /* Calculation of (x,r) and*/
  /* (dU, dU_x, dU_r, dU_3, dU_xx, dU_xr, dU_x3, dU_rr, dU_r3, dU_33)*/
  C_To_c (nvar, X, R, &x, &r, par_b, dU);

  /* Calculation of (y,z) and*/
  /* (dU, dU_x, dU_y, dU_z, dU_xx, dU_xy, dU_xz, dU_yy, dU_yz, dU_zz)*/
  rx3_To_xyz (nvar, x, r, phi, &y, &z, dU);

  LinEquations (A, B, X, R, x, r, phi, y, z, dU, U, values);

  for (ivar = 0; ivar < nvar; ivar++)
    values[ivar] *= FAC;
  
  free_derivs (dU);
  free_derivs (U);
}

void SetMatrix_JFD (int nvar, int n1, int n2, int n3, derivs *u,
		    int *ncols, int **cols, double **Matrix)
{
  int column, row, mcol;
  int i, i1, i_0, i_1, j, j1, j_0, j_1, k, k1, k_0, k_1, N1, N2, N3,
    ivar, ivar1, ntotal = nvar * n1 * n2 * n3;
  double *values;
  derivs *dv;

  values = dvector (0, nvar - 1);
  allocate_derivs (&dv, ntotal);

  N1 = n1 - 1;
  N2 = n2 - 1;
  N3 = n3 - 1;

#ifdef TP_OMP
#pragma omp parallel for private (i,j,k,ivar,row) schedule(dynamic)
#endif
  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < n2; j++)
      {
	for (k = 0; k < n3; k++)
	  {
	    for (ivar = 0; ivar < nvar; ivar++)
	      {
		row = Index (ivar, i, j, k, nvar, n1, n2, n3);
		ncols[row] = 0;
		dv->d0[row] = 0;
	      }
	  }
      }
  }
  for (i = 0; i < n1; i++)
    {
      for (j = 0; j < n2; j++)
	{
	  for (k = 0; k < n3; k++)
	    {
	      for (ivar = 0; ivar < nvar; ivar++)
		{
		  column = Index (ivar, i, j, k, nvar, n1, n2, n3);
		  dv->d0[column] = 1;
		  
		  i_0 = maximum2 (0, i - 1);
		  i_1 = minimum2 (N1, i + 1);
		  j_0 = maximum2 (0, j - 1);
		  j_1 = minimum2 (N2, j + 1);
		  k_0 = k - 1;
		  k_1 = k + 1;
		  /*
		    i_0 = 0;
		    i_1 = N1;
		    j_0 = 0;
		    j_1 = N2;
		    k_0 = 0;
		    k_1 = N3;
		  */
		  
		  for (i1 = i_0; i1 <= i_1; i1++)
		    {
		      for (j1 = j_0; j1 <= j_1; j1++)
			{
			  for (k1 = k_0; k1 <= k_1; k1++)
			    {
			      JFD_times_dv (i1, j1, k1, nvar, n1, n2, n3,
					    dv, u, values);
			      for (ivar1 = 0; ivar1 < nvar; ivar1++)
				{
				  if (values[ivar1] != 0)
				    {
				      row = Index (ivar1, i1, j1, k1, nvar, n1, n2, n3);
				      mcol = ncols[row];
				      cols[row][mcol] = column;
				      Matrix[row][mcol] = values[ivar1];
				      ncols[row] += 1;
				    }
				}
			    }
			}
		    }
		  
		  dv->d0[column] = 0;
		}
	    }
	}
    }
  free_derivs (dv);
  free_dvector (values, 0, nvar - 1);
}

double PunctEvalAtArbitPosition (double *v, int ivar, double A, double B, double phi,
				 int nvar, int n1, int n2, int n3)
{

  /* Calculates the value of v at an arbitrary position (A,B,phi)*/
  
  int i, j, k, N;
  double *p, *values1, **values2, result;
  
  N = maximum3 (n1, n2, n3);
  p = dvector (0, N);
  values1 = dvector (0, N);
  values2 = dmatrix (0, N, 0, N);
  
  for (k = 0; k < n3; k++)
    {
      for (j = 0; j < n2; j++)
	{
	  for (i = 0; i < n1; i++)
	    /* p[i] = v[ivar + nvar * (i + n1 * (j + n2 * k))]; */
	    p[i] = v[Index (ivar, i, j, k, nvar, n1, n2, n3)];
	  chebft_Zeros (p, n1, 0);
	  values2[j][k] = chebev (-1, 1, p, n1, A);
	}
    }
  
  for (k = 0; k < n3; k++)
    {
      for (j = 0; j < n2; j++)
	p[j] = values2[j][k];
      chebft_Zeros (p, n2, 0);
      values1[k] = chebev (-1, 1, p, n2, B);
    }
  
  fourft (values1, n3, 0);
  result = fourev (values1, n3, phi);
  
  free_dvector (p, 0, N);
  free_dvector (values1, 0, N);
  free_dmatrix (values2, 0, N, 0, N);
  
  return result;
}

void calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1, int n2,
		       int n3, derivs *v, derivs *vv)
{
  double al = Pih * (2 * i + 1) / n1, be = Pih * (2 * j + 1) / n2,
    sin_al = sin (al), sin2_al = sin_al * sin_al, cos_al = cos (al),
    sin_be = sin (be), sin2_be = sin_be * sin_be, cos_be = cos (be);
  
  vv->d0[0] = v->d0[Index (ivar, i, j, k, nvar, n1, n2, n3)];
  vv->d1[0] = v->d1[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al;
  vv->d2[0] = v->d2[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_be;
  vv->d3[0] = v->d3[Index (ivar, i, j, k, nvar, n1, n2, n3)];
  vv->d11[0] = v->d11[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin2_al
    + v->d1[Index (ivar, i, j, k, nvar, n1, n2, n3)] * cos_al;
  vv->d12[0] = v->d12[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al * sin_be;
  vv->d13[0] = v->d13[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al;
  vv->d22[0] = v->d22[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin2_be
    + v->d2[Index (ivar, i, j, k, nvar, n1, n2, n3)] * cos_be;
  vv->d23[0] = v->d23[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_be;
  vv->d33[0] = v->d33[Index (ivar, i, j, k, nvar, n1, n2, n3)];
}

double interpol (double a, double b, double c, derivs *v)
{
  return v->d0[0]
    + a * v->d1[0] + b * v->d2[0] + c * v->d3[0]
    + 0.5 * a * a * v->d11[0] + a * b * v->d12[0] + a * c * v->d13[0]
    + 0.5 * b * b * v->d22[0] + b * c * v->d23[0] + 0.5 * c * c * v->d33[0];
}

static double clamp_pm_one (double val)
{
  return val < -1 ? -1 : val > 1 ? 1 : val;
}

double PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1, int n2, int n3,
					 derivs *v,
					 double x, double y, double z)
{
  /* Calculates the value of v at an arbitrary position (x,y,z)*/
  
  double par_b = params_get_real("par_b");
  
  double xs, ys, zs, rs2, phi, X, R, A, B, al, be, aux1, aux2, a, b, c,
    result, Ui;
  int i, j, k;

  derivs *vv;
  allocate_derivs (&vv, 1);

  xs = x / par_b;
  ys = y / par_b;
  zs = z / par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;
  
  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  
  /* Note: Range of R = asin(Q) is [0,pi] for Q in [0,1] */
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;
  
  A = clamp_pm_one( 2 * tanh (0.5 * X) - 1 );

  /* Note: Range of R/2 - pi/4 is [ -pi/4, pi/4 ] and so range of tan
   * is [-1,1], for R in [0,pi]. */
  B = clamp_pm_one( tan (0.5 * R - Piq) );
  al = Pi - acos (A);
  be = Pi - acos (B);

  i = rint (al * n1 / Pi - 0.5);
  j = rint (be * n2 / Pi - 0.5);
  k = rint (0.5 * phi * n3 / Pi);

  a = al - Pi * (i + 0.5) / n1;
  b = be - Pi * (j + 0.5) / n2;
  c = phi - 2 * Pi * k / n3;

  calculate_derivs (i, j, k, ivar, nvar, n1, n2, n3, v, vv);
  result = interpol (a, b, c, vv);
  
  free_derivs (vv);
  
  Ui = (A - 1) * result;
  
  assert( isfinite( Ui ) );
  
  return Ui;
}

double PunctIntPolAtArbitPosition (int ivar, int nvar,
				   int n1, int n2,int n3,
				   derivs *v,
				   double x, double y, double z)
{
  /* Calculates the value of v at an arbitrary position (x,y,z)*/
  
  double par_b = params_get_real("par_b");

  double xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;

  xs = x / par_b;
  ys = y / par_b;
  zs = z / par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);
  
  result = PunctEvalAtArbitPosition (v->d0, ivar, A, B, phi, nvar, n1, n2, n3);

  Ui = (A - 1) * result;

  assert( isfinite( Ui ) );
  
  return Ui;
}


//////////////////////////////////////////////////////
/// Fast Spectral Interpolation Routine Stuff
//////////////////////////////////////////////////////

double  PunctEvalAtArbitPositionFast (double *v, int ivar, double A, double B, double phi, int nvar, int n1, int n2, int n3)
{
  /* Calculates the value of v at an arbitrary position (A,B,phi)* using the fast routine */

  int i, j, k, N;
  double *p, *values1, **values2, result;
  // VASILIS: Nothing should be changed in this routine. This is used by PunctIntPolAtArbitPositionFast

  N = maximum3 (n1, n2, n3);

  p = dvector (0, N);
  values1 = dvector (0, N);
  values2 = dmatrix (0, N, 0, N);

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++)
    {
      for (i = 0; i < n1; i++)
	/* p[i] = v[ivar + nvar * (i + n1 * (j + n2 * k))]; */
	p[i] = v[Index (ivar, i, j, k, nvar, n1, n2, n3)];
      //      chebft_Zeros (p, n1, 0);
      values2[j][k] = chebev (-1, 1, p, n1, A);
    }
  }

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++) p[j] = values2[j][k];
    //    chebft_Zeros (p, n2, 0);
    values1[k] = chebev (-1, 1, p, n2, B);
  }

  //  fourft (values1, n3, 0);
  result = fourev (values1, n3, phi);

  free_dvector (p, 0, N);
  free_dvector (values1, 0, N);
  free_dmatrix (values2, 0, N, 0, N);

  return result;
  //  */
  //  return 0.;
}

double PunctIntPolAtArbitPositionFast (int ivar, int nvar,
				       int n1, int n2, int n3,
				       derivs *v,
				       double x, double y, double z)
{
  
  /* Calculates the value of v at an arbitrary position (x,y,z) if the spectral coefficients are known */
  double par_b = params_get_real("par_b");
  
  double xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;
  // VASILIS: Here the struct derivs v refers to the spectral coeffiecients of variable v not the variable v itself

  xs = x / par_b;
  ys = y / par_b;
  zs = z / par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);

  result = PunctEvalAtArbitPositionFast (v->d0, ivar, A, B, phi, nvar, n1, n2, n3);

  Ui = (A - 1) * result;

  return Ui;
}

void SpecCoef(int n1, int n2, int n3, int nvar, double *v, double *cf)
{
  /* Evaluates the spectral expansion coefficients of v   */

  // VASILIS: Here v is a pointer to the values of the variable v at the collocation points and cf_v a pointer to the spectral coefficients that this routine calculates

  int i, j, k, N, n, l, m;
  double *p, ***values3, ***values4;

  N=maximum3(n1,n2,n3);
  p=dvector(0,N);
  values3=d3tensor(0,n1,0,n2,0,n3);
  values4=d3tensor(0,n1,0,n2,0,n3);

for (int ivar=0; ivar<nvar; ivar++){
  
  // Caclulate values3[n,j,k] = a_n^{j,k} = (sum_i^(n1-1) f(A_i,B_j,phi_k) Tn(-A_i))/k_n , k_n = N/2 or N
  for(k=0;k<n3;k++) {
    for(j=0;j<n2;j++) {

      for(i=0;i<n1;i++)
	/* p[i]=v[ivar + (i + n1 * (j + n2 * k))]; */
	p[i]=v[Index (ivar, i, j, k, nvar, n1, n2, n3)];
      
      chebft_Zeros(p,n1,0);
      for (n=0;n<n1;n++)	{
	values3[n][j][k] = p[n];
      }
    }
  }

  // Caclulate values4[n,l,k] = a_{n,l}^{k} = (sum_j^(n2-1) a_n^{j,k} Tn(B_j))/k_l , k_l = N/2 or N

  for (n = 0; n < n1; n++){
    for(k=0;k<n3;k++) {
      for(j=0;j<n2;j++) p[j]=values3[n][j][k];
      chebft_Zeros(p,n2,0);
      for (l = 0; l < n2; l++){
	values4[n][l][k] = p[l];
      }
    }
  }

  // Caclulate coefficients  a_{n,l,m} = (sum_k^(n3-1) a_{n,m}^{k} fourier(phi_k))/k_m , k_m = N/2 or N
  for (i = 0; i < n1; i++){
    for (j = 0; j < n2; j++){
      for(k=0;k<n3;k++) p[k]=values4[i][j][k];
      fourft(p,n3,0);
      for (k = 0; k<n3; k++){
	/* cf[ivar + (i + n1 * (j + n2 * k))] = p[k]; */
	cf[Index (ivar, i, j, k, nvar, n1, n2, n3)] = p[k];
      }
    }
  }

}

  free_dvector(p,0,N);
  free_d3tensor(values3,0,n1,0,n2,0,n3);
  free_d3tensor(values4,0,n1,0,n2,0,n3);

}

/* Initial guess */
void set_initial_guess(derivs *v)
{
  int nvar = 1,
    n1 = params_get_int("npoints_A"),
    n2 = params_get_int("npoints_B"),
    n3 = params_get_int("npoints_phi");

  double par_b = params_get_real("par_b");

  double *s_x, *s_y, *s_z;
  double al, A, Am1, be, B, phi, R, r, X;
  int ivar, i, j, k, i3D, indx;
  derivs *U;

  char fname[STRLEN];
  FILE *debug_file =NULL;
  
  if (params_get_int("solve_momentum_constraint"))
    nvar = 4;
  
  s_x = calloc(n1*n2*n3, sizeof(double));
  s_y = calloc(n1*n2*n3, sizeof(double));
  s_z = calloc(n1*n2*n3, sizeof(double));

  allocate_derivs (&U, nvar);
  
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++) {
	  
          i3D = Index(ivar,i,j,k,1,n1,n2,n3);
	  
          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;
	  
          /* Calculation of (X,R)*/
          AB_To_XR (nvar, A, B, &X, &R, U);
	  
	  /* Calculation of (x,r)*/
          C_To_c (nvar, X, R, &(s_x[i3D]), &r, par_b, U);
	  
	  /* Calculation of (y,z)*/
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
	  
	}
  
  //SB See: TwoPunctures/interface.ccl
  //SB /* Set_Initial_Guess_for_u(n1*n2*n3, v->d0, s_x, s_y, s_z); */ 
  //SB: In the BAM code this guess looks like  v = 0 , dv = 0
  
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++) {
	  indx = Index(ivar,i,j,k,1,n1,n2,n3);
          v->d0[indx]/=(-cos(Pih * (2 * i + 1) / n1)-1.0);
	}
  
  Derivatives_AB3 (nvar, n1, n2, n3, v);
  
  if (params_get_int("do_initial_debug_output")) {

    strcpy(fname,params_get_str("outputdir"));
    strcat (fname,"/initial.dat");
    debug_file=fopen(fname, "w");
    assert(debug_file);

    for (ivar = 0; ivar < nvar; ivar++)
      for (i = 0; i < n1; i++)
	for (j = 0; j < n2; j++) {
	  
	  al = Pih * (2 * i + 1) / n1;
	  A = -cos (al);
	  Am1 = A -1.0;
	  be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 0.0;
          indx = Index(ivar,i,j,0,1,n1,n2,n3);
          U->d0[0] = Am1 * v->d0[indx];        /* U*/
          U->d1[0] = v->d0[indx] + Am1 * v->d1[indx];        /* U_A*/
          U->d2[0] = Am1 * v->d2[indx];        /* U_B*/
          U->d3[0] = Am1 * v->d3[indx];        /* U_3*/
          U->d11[0] = 2 * v->d1[indx] + Am1 * v->d11[indx];        /* U_AA*/
          U->d12[0] = v->d2[indx] + Am1 * v->d12[indx];        /* U_AB*/
          U->d13[0] = v->d3[indx] + Am1 * v->d13[indx];        /* U_AB*/
          U->d22[0] = Am1 * v->d22[indx];        /* U_BB*/
          U->d23[0] = Am1 * v->d23[indx];        /* U_B3*/
          U->d33[0] = Am1 * v->d33[indx];        /* U_33*/

	  /* Calculation of (X,R)*/
	  AB_To_XR (nvar, A, B, &X, &R, U);

	  /* Calculation of (x,r)*/
	  C_To_c (nvar, X, R, &(s_x[indx]), &r, par_b, U);

	  /* Calculation of (y,z)*/
	  rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
	  fprintf(debug_file,
		  "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
		  "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
		  (double)s_x[indx], (double)s_y[indx],
		  (double)A,(double)B,
		  (double)U->d0[0],
		  (double)(-cos(Pih * (2 * i + 1) / n1)-1.0),
		  (double)U->d1[0],
		  (double)U->d2[0],
		  (double)U->d3[0],
		  (double)U->d11[0],
		  (double)U->d22[0],
		  (double)U->d33[0],
		  (double)v->d0[indx],
		  (double)v->d1[indx],
		  (double)v->d2[indx],
		  (double)v->d3[indx],
		  (double)v->d11[indx],
		  (double)v->d22[indx],
		  (double)v->d33[indx]
		  );
        }
    
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2; i++) {
      double d;
      indx = Index(0,0,i,0,1,n1,n2,n3);
      d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
				     s_x[indx], 0.0, 0.0);
      fprintf(debug_file, "%.16g %.16g\n",
	      (double)s_x[indx], (double)d);
    }
    
    fprintf(debug_file, "\n\n");

    for (i=n2-10; i<n2-1; i++) {
      double d;
      int ip= Index(0,0,i+1,0,1,n1,n2,n3);
      indx = Index(0,0,i,0,1,n1,n2,n3);
      for (j=-10; j<10; j++) {
        d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
				       s_x[indx]+(s_x[ip]-s_x[indx])*j/10,
				       0.0, 0.0);
        fprintf(debug_file, "%.16g %.16g\n",
                (double)(s_x[indx]+(s_x[ip]-s_x[indx])*j/10), (double)d);
      }
    }
    
    fprintf(debug_file, "\n\n");
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++) {
        X = 2*(2.0*i/n1-1.0);
        R = 2*(1.0*j/n2);
        if (X*X+R*R > 1.0) {
	  C_To_c (nvar, X, R, &(s_x[indx]), &r, par_b, U);
	  rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);

	  //SB: this is wrong:
	  *U->d0  = s_x[indx]*s_x[indx];
	  *U->d1  = 2*s_x[indx];
	  *U->d2  = 0.0;
	  *U->d3 = 0.0;
	  *U->d11 = 2.0;
	  *U->d22 = 0.0;
	  *U->d33 = *U->d12 = *U->d23 = *U->d13 = 0.0;

	  C_To_c (nvar, X, R, &(s_x[indx]), &r, par_b, U);
	  fprintf(debug_file,
		    "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                  (double)s_x[indx], (double)r, (double)X, (double)R, (double)U->d0[0],
		  (double)U->d1[0],
		  (double)U->d2[0],
		    (double)U->d3[0],
		  (double)U->d11[0],
		  (double)U->d22[0],
		  (double)U->d33[0]);
        }
      }
    fclose(debug_file);
  }
  
  free(s_z);
  free(s_y);
  free(s_x);
  free_derivs (U);
}

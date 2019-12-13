/* TP_Equations.c */

#include "TwoPunctures.h"

double BY_KKofxyz (double x, double y, double z)
{

  double par_b = params_getd("par_b");

  double par_P_plus[3], par_P_minus[3];
  double par_S_plus[3], par_S_minus[3];

  par_P_plus[0] = params_getd("par_P_plus1");
  par_P_plus[1] = params_getd("par_P_plus2");
  par_P_plus[2] = params_getd("par_P_plus3");
  par_P_minus[0] = params_getd("par_P_minus1");
  par_P_minus[1] = params_getd("par_P_minus2");
  par_P_minus[2] = params_getd("par_P_minus3");
  
  par_S_plus[0] = params_getd("par_S_plus1");
  par_S_plus[1] = params_getd("par_S_plus2");
  par_S_plus[2] = params_getd("par_S_plus3");
  par_S_minus[0] = params_getd("par_S_minus1");
  par_S_minus[1] = params_getd("par_S_minus2");
  par_S_minus[2] = params_getd("par_S_minus3");
  
  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];
  
  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;
  
  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;
  
  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++) {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }

  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];

  AijAij = 0;
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Bowen-York-Curvature :*/
      Aij =
	+ 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
                 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
      AijAij += Aij * Aij;
    }
  }
  
  return AijAij;
}

void BY_Aijofxyz (double x, double y, double z, double Aij[3][3])
{
  
  double par_b = params_getd("par_b");
  double TP_epsilon = params_getd("TP_epsilon");
  double TP_Tiny = params_getd("TP_Tiny");

  double par_P_plus[3], par_P_minus[3];
  double par_S_plus[3], par_S_minus[3];

  par_P_plus[0] = params_getd("par_P_plus1");
  par_P_plus[1] = params_getd("par_P_plus2");
  par_P_plus[2] = params_getd("par_P_plus3");
  par_P_minus[0] = params_getd("par_P_minus1");
  par_P_minus[1] = params_getd("par_P_minus2");
  par_P_minus[2] = params_getd("par_P_minus3");
  
  par_S_plus[0] = params_getd("par_S_plus1");
  par_S_plus[1] = params_getd("par_S_plus2");
  par_S_plus[2] = params_getd("par_S_plus3");
  par_S_minus[0] = params_getd("par_S_minus1");
  par_S_minus[1] = params_getd("par_S_minus2");
  par_S_minus[2] = params_getd("par_S_minus3");
  
  int i, j;
  double r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_plus = sqrt (pow (r2_plus, 2) + pow (TP_epsilon, 4));
  r2_minus = sqrt (pow (r2_minus, 2) + pow (TP_epsilon, 4));
  if (r2_plus < pow(TP_Tiny,2))
    r2_plus = pow(TP_Tiny,2);
  if (r2_minus < pow(TP_Tiny,2))
    r2_minus = pow(TP_Tiny,2);
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++) {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }
  
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
		 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }
}

void NonLinEquations (double rho_adm,
		      double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs *U, double *values)
{
  double par_b = params_getd("par_b");
  double par_m_plus = params_getd("par_m_plus");
  double par_m_minus = params_getd("par_m_minus");

  double r_plus, r_minus, psi, psi2, psi4, psi7;
  double mu;

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);

  psi = 1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U->d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;
  
  values[0] = U->d11[0] + U->d22[0] + U->d33[0] + 0.125 * BY_KKofxyz (x, y, z) / psi7 
    + 2.0 * Pi / psi2/psi * rho_adm;
}

void LinEquations (double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs *dU, derivs *U, double *values)
{
  double par_b = params_getd("par_b");
  double par_m_plus = params_getd("par_m_plus");
  double par_m_minus = params_getd("par_m_minus");
  
  double r_plus, r_minus, psi, psi2, psi4, psi8;
  
  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  
  psi = 1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U->d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;
  
  values[0] = dU->d11[0] + dU->d22[0] + dU->d33[0] - 0.875 * BY_KKofxyz (x, y, z) / psi8 * dU->d0[0];
}


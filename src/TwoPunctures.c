/* TwoPunctures.c */

#include "TwoPunctures.h"
#include "stdbool.h"

/* -------------------------------------------------------------------------*/

/* Manage memory internally */
static int allocated_params = 0;

static void _alloc_params_mem_if_req(){
  if(!allocated_params){
    params_alloc();
    allocated_params = 1;
  }
}

static void _dealloc_params_mem_if_req(){
  if(allocated_params){
    params_free();
    allocated_params = 0;
  }
}

/* Control interaction with internal parameters in the following */ //

void TwoPunctures_params_set_Real(char *key, double value){
  /*
    Set parameters according to input.
  */
  params_set_real(key, value);       // replace, don't append new
}

void TwoPunctures_params_set_Int(char *key, int value){
  /*
    Set parameters according to input.
  */
  params_set_int(key, value);   // replace, don't append new
}

void TwoPunctures_params_set_Boolean(char *key, bool value){
  /*
    Set parameters according to input.
  */
  params_set_bool(key, value);   // replace, don't append new
}

void TwoPunctures_params_set_String(char *key, char * value){
  /*
    Set parameters according to input.
  */
  params_set_str(key, value);   // replace, don't append new
}

void TwoPunctures_params_set_inputfile(char* inputfile){
  /*
    Set parameters from input file
  */

  // allocate and seed with defaults such that we have a guaranteee of params
  _alloc_params_mem_if_req();
  TwoPunctures_params_set_default();

  if(inputfile!=NULL){
    params_read(inputfile);

    char od[STRLEN];//SB: this can be improved.
    int n = strlen(inputfile);
    if (STREQL(inputfile+(n-4),".par")) {
      strncpy (od,inputfile,n-4);
      od[n-4] = '\0';
      //printf("%s\n",s);
      params_set_str("outputdir",od); 
    }
  }
}

void TwoPunctures_params_set_default(){
  /*
    Set default parameters.
  */
  _alloc_params_mem_if_req();

  /* Add parameters here */
  params_add_real("par_b",1.0); // x coordinate of the m+ puncture; 1.0
  params_add_real("par_m_plus",1.0); // mass of the m+ puncture; 1.0
  params_add_real("par_m_minus",0.0); // mass of the m- puncture; 0.0
  params_add_real("target_M_plus",0.5); // target ADM mass for m+; 0.5
  params_add_real("target_M_minus",0.5); // target ADM mass for m-; 0.5
  params_add_real("par_P_plus1",0.); // momentum of the m+ puncture
  params_add_real("par_P_plus2",0.);
  params_add_real("par_P_plus3",0.);
  params_add_real("par_P_minus1",0.); // momentum of the m- puncture
  params_add_real("par_P_minus2",0.);
  params_add_real("par_P_minus3",0.);
  params_add_real("par_S_plus1",0.); // spin of the m+ puncture
  params_add_real("par_S_plus2",0.);
  params_add_real("par_S_plus3",0.);
  params_add_real("par_S_minus1",0.); // spin of the m- puncture
  params_add_real("par_S_minus2",0.);
  params_add_real("par_S_minus3",0.);
  params_add_real("center_offset1",0.0); // offset b=0 to position (x,y,z) ; 0
  params_add_real("center_offset2",0.0);  // ; 0
  params_add_real("center_offset3",0.0); // ; 0

  params_add_int("give_bare_mass",1); // User provides bare masses rather than target ADM masses ; 1
  params_add_int("npoints_A",30); // Number of coefficients in the compactified radial direction
  params_add_int("npoints_B",30); // Number of coefficients in the angular direction
  params_add_int("npoints_phi",16); // Number of coefficients in the phi direction

  params_add_real("Newton_tol",1e-10); // Tolerance for Newton solver
  params_add_int("Newton_maxit",5); // Maximum number of Newton iterations

  params_add_real("TP_epsilon",0.);  // A small number to smooth out singularities at the puncture locations
  params_add_real("TP_Tiny",0.);  // Tiny number to avoid nans near or at the pucture locations
  params_add_real("TP_Extend_Radius",0); // Radius of an extended spacetime instead of the puncture ; 0
  params_add_real("adm_tol",1e-10); // Tolerance of ADM masses when give_bare_mass=no
  params_add_int("solve_momentum_constraint",0); // Solve for momentum constraint?
  params_add_int("use_external_initial_guess",0); // Set initial guess by external function?

  // output
  params_add_str("outputdir","./"); // Output debug information about the residuum
  params_add_int("do_residuum_debug_output",0); // Output debug information about the residuum
  params_add_int("do_initial_debug_output",0); // Output debug information about initial guess
  params_add_int("do_solution_file_output",0); // output .data files
  params_add_int("do_bam_file_output",0); // output input files for bam's puncture_ps
  
  // Interpolation
  params_add_int("grid_setup_method",taylor_expansion); // How to fill the 3D grid from the spectral grid ?
  params_add_int("initial_lapse",psin); // How to set the lapse ?
  params_add_real("initial_lapse_psi_exponent",-2.0); // Exponent n for psi^-n initial lapse profile (<0)
  params_add_int("conformal_state",1); // Ways to set the conformal factor
  params_add_int("swap_xz",0); // Swap x and z coordinates when interpolating, // so that the black holes are separated in the z direction; 0
  params_add_int("multiply_old_lapse",0); // Multiply the old lapse with the new one

  params_add_int("verbose",0);

}

/* -------------------------------------------------------------------------*/

/* Swap two variables */
static inline
void swap (double * restrict const a, double * restrict const b)
{
  double const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))

ini_data* TwoPunctures_make_initial_data() {
  
  /* Prepare initial data based on internal settings 
     This alloc the mem */

  const int verbose = params_get_int("verbose");

  char outdir[STRLEN];
  if (params_get_int("do_residuum_debug_output")+
      params_get_int("do_initial_debug_output")+
      params_get_int("do_solution_file_output")+
      params_get_int("do_bam_file_output")) {
    make_output_dir();
  }
  
  double par_P_plus[3], par_P_minus[3];
  par_P_plus[0] = params_get_real("par_P_plus1");
  par_P_plus[1] = params_get_real("par_P_plus2");
  par_P_plus[2] = params_get_real("par_P_plus3");
  par_P_minus[0] = params_get_real("par_P_minus1");
  par_P_minus[1] = params_get_real("par_P_minus2");
  par_P_minus[2] = params_get_real("par_P_minus3");

  double par_S_plus[3], par_S_minus[3];
  par_S_plus[0] = params_get_real("par_S_plus1");
  par_S_plus[1] = params_get_real("par_S_plus2");
  par_S_plus[2] = params_get_real("par_S_plus3");
  par_S_minus[0] = params_get_real("par_S_minus1");
  par_S_minus[1] = params_get_real("par_S_minus2");
  par_S_minus[2] = params_get_real("par_S_minus3");

  double par_b = params_get_real("par_b");

  double center_offset[3];
  center_offset[0] = params_get_real("center_offset1");
  center_offset[1] = params_get_real("center_offset2");
  center_offset[2] = params_get_real("center_offset3");

  // ---------------------------------------------------

  double E; // ADM energy of the Bowen-York spacetime"
  double J1, J2, J3; // Angular momentum of the Bowen-York spacetime"
  double mp = params_get_real("par_m_plus"), mm = params_get_real("par_m_minus"); // Bare masses of the punctures
  double mp_adm, mm_adm; // ADM masses of the punctures (measured at the other spatial infinities
  double admMass;

  // ---------------------------------------------------

  int const nvar = 1,
    n1 = params_get_int("npoints_A"),
    n2 = params_get_int("npoints_B"),
    n3 = params_get_int("npoints_phi");

  int const ntotal = n1 * n2 * n3 * nvar;
#if (0)
  int percent10 = 0;
#endif
  static double *F = NULL;
  static derivs *u, *v, *cf_v;

  if (! F) {
    double up, um;

    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);
    allocate_derivs (&cf_v, ntotal);

    //if (use_sources) {
    //  printf ("Solving puncture equation for BH-NS/NS-NS system");
    //} else {
    if (verbose) printf ("Solving puncture equation for BH-BH system\n");
    //}

    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++) {
      cf_v->d0[j] = 0.0;
      cf_v->d1[j] = 0.0;
      cf_v->d2[j] = 0.0;
      cf_v->d3[j] = 0.0;
      cf_v->d11[j] = 0.0;
      cf_v->d12[j] = 0.0;
      cf_v->d13[j] = 0.0;
      cf_v->d22[j] = 0.0;
      cf_v->d23[j] = 0.0;
      cf_v->d33[j] = 0.0;
      v->d0[j] = 0.0;
      v->d1[j] = 0.0;
      v->d2[j] = 0.0;
      v->d3[j] = 0.0;
      v->d11[j] = 0.0;
      v->d12[j] = 0.0;
      v->d13[j] = 0.0;
      v->d22[j] = 0.0;
      v->d23[j] = 0.0;
      v->d33[j] = 0.0;
    }

    /* call for external initial guess */
    if (params_get_int("use_external_initial_guess")) {
      set_initial_guess(v);
    }

    /* If bare masses are not given, iteratively solve for them given the
       target ADM masses target_M_plus and target_M_minus and with initial
       guesses given by par_m_plus and par_m_minus. */
    if(!(params_get_int("give_bare_mass"))) {

      double tmp, mp_adm_err, mm_adm_err;
      //char valbuf[100];

      double M_p = params_get_real("target_M_plus");
      double M_m = params_get_real("target_M_minus");

      if (verbose) {
	printf ("Attempting to find bare masses.\n");
	printf ("Target ADM masses: M_p=%g and M_m=%g\n", M_p, M_m);
	printf ("ADM mass tolerance: %g\n", params_get_real("adm_tol"));
      }

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        if (verbose) printf ("Bare masses: mp=%.15g, mm=%.15g\n", mp, mm);
        Newton (nvar, n1, n2, n3, v, params_get_real("Newton_tol"), 1);

        F_of_v (nvar, n1, n2, n3, v, F, u);

        up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        mp_adm = (1 + up) * mp + mp * mm / (4. * par_b);
        mm_adm = (1 + um) * mm + mp * mm / (4. * par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p-mp_adm);
        mm_adm_err = fabs(M_m-mm_adm);
	if (verbose) printf ("ADM mass error: M_p_err=%.15g, M_m_err=%.15g\n",
			     mp_adm_err, mm_adm_err);

        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*par_b*( 1 + um + up + um*up ) +
	  sqrt(16*par_b*M_m*(1 + um)*(1 + up) +
	       pow(-M_m + M_p + 4*par_b*(1 + um)*(1 + up),2));
        mp = (tmp + M_p - M_m)/(2.*(1 + up));
        mm = (tmp - M_p + M_m)/(2.*(1 + um));

      } while ( (mp_adm_err > params_get_real("adm_tol")) ||
                (mm_adm_err > params_get_real("adm_tol")) );

      if (verbose) printf ("Found bare masses.\n");
    }

    Newton (nvar, n1, n2, n3, v, params_get_real("Newton_tol"), params_get_int("Newton_maxit"));

    F_of_v (nvar, n1, n2, n3, v, F, u);

    /* SpecCoef(n1, n2, n3, 0, v->d0, cf_v->d0); */
    SpecCoef(n1, n2, n3, nvar, v->d0, cf_v->d0);
 
    if (verbose) printf ("The two puncture masses are mp=%.17g and mm=%.17g\n", mp, mm);

    up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
    um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    mp_adm = (1 + up) * mp + mp * mm / (4. * par_b);
    mm_adm = (1 + um) * mm + mp * mm / (4. * par_b);

    if (verbose) printf ("Puncture 1 ADM mass is %g\n", mp_adm);
    if (verbose) printf ("Puncture 2 ADM mass is %g\n", mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (mp + mm - 4*par_b*PunctEvalAtArbitPosition(v->d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    if (verbose) printf ("The total ADM mass is %g\n", admMass);

    E = admMass;

    J1 = -(center_offset[2]*par_P_minus[1]) + center_offset[1]*par_P_minus[2] - center_offset[2]*par_P_plus[1] + center_offset[1]*par_P_plus[2] + par_S_minus[0] + par_S_plus[0];
    J2 = center_offset[2]*par_P_minus[0] - center_offset[0]*par_P_minus[2] + par_b*par_P_minus[2] + center_offset[2]*par_P_plus[0] - center_offset[0]*par_P_plus[2] - par_b*par_P_plus[2] + par_S_minus[1] + par_S_plus[1];
    J3 = -(center_offset[1]*par_P_minus[0]) + center_offset[0]*par_P_minus[1] - par_b*par_P_minus[1] - center_offset[1]*par_P_plus[0] + center_offset[0]*par_P_plus[1] + par_b*par_P_plus[1] + par_S_minus[2] + par_S_plus[2];

  }

  if (params_get_int("do_solution_file_output")) {
    /* Output the solution */    
    write_derivs( u, n1,n2,n3, 
		  0, // =0,1,2 derivatives to output
		  "u.data");
    write_derivs( v, n1,n2,n3, 
		  0, // =0,1,2 derivatives to output
		  "v.data");
    write_derivs( cf_v, n1,n2,n3, 
		  0, // =0,1,2 derivatives to output
		  "cf_v.data");
  }

  if (params_get_int("do_bam_file_output")) {
    /* Output the solution for bam */
    write_bam_inifile( u, n1,n2,n3,"u.bamdata");
    write_bam_inifile( v, n1,n2,n3,"v.bamdata");
    write_bam_inifile( cf_v, n1,n2,n3,"cf_v.bamdata");
  }
  
  /*
    free_dvector (F, 0, ntotal - 1);
    free_derivs (u);
    free_derivs (v);
    free_derivs (cf_v);
  */  
  /*
  ini_data ret;  // ret_ptr = &ret;
  ret.F = F;     // (&ret)->F = F;
  ret.u = &u;
  ret.v = &v;
  ret.cf_v = &cf_v;
  ret.ntotal = ntotal;
  */
  //SB: shouldn't be something like:
  ini_data *data;
  data = (ini_data *) calloc(1, sizeof(ini_data)); 
  if (data == NULL) ERROR("Out of memory");  
  data->F = F; //SB: Or should we alloc data->arrays and copy?
  data->u = u;
  data->v = v;
  data->cf_v = cf_v;
  data->ntotal = ntotal;
  
  return data;
}

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
 double *kzz){

  // TODO:
  // add safety check for null ini_data to make it
  if (data == NULL) ERROR("No data to interpolate");
  
  /* // unpack initial data structure */
  /* double *F = data.F; */
  /* // derivs *u = data.u; */
  /* derivs v = *(data.v); */
  /* derivs cf_v = *(data.cf_v); */
  /* int ntotal = data.ntotal; */
  //SB:
  double *F = data->F; 
  derivs *u = data->u; 
  derivs *v = data->v;
  derivs *cf_v = data->cf_v; 
  int ntotal = data->ntotal; 
  
  // ---- prepare required parameters
  double par_b = params_get_real("par_b");
  double mp = params_get_real("par_m_plus"),
    mm = params_get_real("par_m_minus"); // Bare masses of the punctures

  double center_offset[3];
  center_offset[0] = params_get_real("center_offset1");
  center_offset[1] = params_get_real("center_offset2");
  center_offset[2] = params_get_real("center_offset3");

  int const nvar = 1,
    n1 = params_get_int("npoints_A"),
    n2 = params_get_int("npoints_B"),
    n3 = params_get_int("npoints_phi");

  const int verbose = params_get_int("verbose");

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  switch (params_get_int("grid_setup_method")) {
  case taylor_expansion:
    gsm = GSM_Taylor_expansion;
    break;
  case evaluation:
    gsm = GSM_evaluation;
    break;
  default:
    ERROR("Something wrong setting parameter: grid_setup_method ");
  }

  // ----

  int antisymmetric_lapse=0,
    averaged_lapse=0,
    pmn_lapse=0,
    brownsville_lapse=0;
  switch (params_get_int("initial_lapse")) {
  case antisymmetric:
    antisymmetric_lapse = 1;
    break;
  case averaged:
    averaged_lapse = 1;
    break;
  case psin:
    pmn_lapse = 1;
    break;
  case brownsville:
    brownsville_lapse = 1;
    break;
  }

  const double initial_lapse_psi_exponent = params_get_real("initial_lapse_psi_exponent");

  if (verbose) {
    if (pmn_lapse)
      printf("Setting initial lapse to psi^%f profile.",
	     initial_lapse_psi_exponent);
    if (brownsville_lapse)
      printf("Setting initial lapse to a Brownsville-style profile "
	     "with exp %f.",
	     initial_lapse_psi_exponent);
  }

  //SB: conformal_state  = 0, 1,2,3 -> NONSTATIC, FACTOR0,FACTOR1,FACTOR2
  //SB: check what to use here !!!
  const int conformal_state = params_get_int("conformal_state");


  const double TP_epsilon = params_get_real("TP_epsilon");
  const double TP_Tiny = params_get_real("TP_Tiny");
  const double TP_Extend_Radius = params_get_real("TP_Extend_Radius");

  const int multiply_old_lapse = params_get_int("multiply_old_lapse");
  const int swap_xz = params_get_int("swap_xz");

  if (verbose) printf ("Interpolating result");

  // Cartesian grid parameters

  const int xxx = nxyz[0];
  const int xxxyyy = nxyz[0]*nxyz[1];

#ifdef TP_OMP
#pragma omp parallel for
#endif
  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; ++i) {

        //const int ind = GFINDEX3D (i, j, k, nshift, mshift);
	const int ind = i + xxx * j + xxxyyy* k;

        double xx, yy, zz;
        xx = x[i] - center_offset[0];
        yy = y[j] - center_offset[1];
        zz = z[k] - center_offset[2];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (xx, zz);
        }

        double r_plus
          = sqrt(pow(xx - par_b, 2) + pow(yy, 2) + pow(zz, 2));
        double r_minus
          = sqrt(pow(xx + par_b, 2) + pow(yy, 2) + pow(zz, 2));

        double U;
        switch (gsm) {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition(0, nvar, n1, n2, n3, v, xx, yy, zz);
          break;
        case GSM_evaluation:
          U = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v,xx, yy, zz);
          break;
        default:
          assert (0);
        }

        r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
        r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
        if (r_plus < TP_Tiny)
	  r_plus = TP_Tiny;
        if (r_minus < TP_Tiny)
	  r_minus = TP_Tiny;
	
        double psi1 = 1
          + 0.5 * mp / r_plus
          + 0.5 * mm / r_minus + U;

#define EXTEND(M,r)                                             \
	( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) -    \
	       5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) +    \
	       15./8 / TP_Extend_Radius))
        if (r_plus < TP_Extend_Radius) {
          psi1 = 1
	    + 0.5 * EXTEND(mp,r_plus)
	    + 0.5 * mm / r_minus + U;
        }
        if (r_minus < TP_Extend_Radius) {
          psi1 = 1
	    + 0.5 * EXTEND(mm,r_minus)
	    + 0.5 * mp / r_plus + U;
        }
        double static_psi = 1;

        double Aij[3][3];
        BY_Aijofxyz (xx, yy, zz, Aij);

        double old_alp=1.0;
        if (multiply_old_lapse)
	  old_alp = alp[ind];

        if ((conformal_state > 0) || (pmn_lapse) || (brownsville_lapse)) {

          double xp, yp, zp, rp, ir;
          double s1, s3, s5;
          double p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = xx - par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
            rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* mp *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = xx + par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
            rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* mm *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (conformal_state >= 1) {
            static_psi = p;
            psi[ind] = static_psi;
          }
          if (conformal_state >= 2) {
            psix[ind] = px / static_psi;
            psiy[ind] = py / static_psi;
            psiz[ind] = pz / static_psi;
          }
          if (conformal_state >= 3) {
            psixx[ind] = pxx / static_psi;
            psixy[ind] = pxy / static_psi;
            psixz[ind] = pxz / static_psi;
            psiyy[ind] = pyy / static_psi;
            psiyz[ind] = pyz / static_psi;
            psizz[ind] = pzz / static_psi;
          }

          if (pmn_lapse)
            alp[ind] = pow(p, initial_lapse_psi_exponent);
          if (brownsville_lapse)
            alp[ind] = 2.0/(1.0+pow(p, initial_lapse_psi_exponent));

        } /* if conformal-state > 0 */

        //puncture_u[ind] = U;

        gxx[ind] = pow (psi1 / static_psi, 4);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow (psi1 / static_psi, 4);
        gyz[ind] = 0;
        gzz[ind] = pow (psi1 / static_psi, 4);

        kxx[ind] = Aij[0][0] / pow(psi1, 2);
        kxy[ind] = Aij[0][1] / pow(psi1, 2);
        kxz[ind] = Aij[0][2] / pow(psi1, 2);
        kyy[ind] = Aij[1][1] / pow(psi1, 2);
        kyz[ind] = Aij[1][2] / pow(psi1, 2);
        kzz[ind] = Aij[2][2] / pow(psi1, 2);

        // printf("(ind, gxx)=(%d, %4.8f)\n", ind, gxx[ind]);
        // printf("\n(ind, psi)=(%d, %4.8f)\n", ind, psi[ind]);

        if (antisymmetric_lapse || averaged_lapse) {
          alp[ind] =
            ((1.0 -0.5* mp /r_plus -0.5* mm/r_minus)
             /(1.0 +0.5* mp /r_plus +0.5* mm/r_minus));

          if (r_plus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(mp, r_plus) -0.5* mm/r_minus)
               /(1.0 +0.5*EXTEND(mp, r_plus) +0.5* mm/r_minus));
          }
          if (r_minus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(mm, r_minus) -0.5* mp/r_plus)
               /(1.0 +0.5*EXTEND(mp, r_minus) +0.5* mp/r_plus));
          }

          if (averaged_lapse) {
            alp[ind] = 0.5 * (1.0 + alp[ind]);
          }
        }


        if (multiply_old_lapse)
          alp[ind] *= old_alp;

        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          if (conformal_state >= 2) {
            SWAP (psix[ind], psiz[ind]);
          }
          if (conformal_state >= 3) {
            SWAP (psixx[ind], psizz[ind]);
            SWAP (psixy[ind], psiyz[ind]);
          }
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);
        } /* if swap_xz */

      } /* for i */
    }   /* for j */
  }     /* for k */

}

void TwoPunctures_finalise(ini_data *data){
  /*
    Clean up all internally allocated objects.
  */
  _dealloc_params_mem_if_req();

  if (data) {
    if (data->F) free(data->F);
    if (data->u) free_derivs(data->u);
    if (data->v) free_derivs(data->v);
    if (data->cf_v) free_derivs(data->cf_v);
    free(data);
  }
}



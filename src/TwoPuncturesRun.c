// Minimal working example for making use of this library

// TODO: valgrind still reveals a mem-leak internally

#include "TwoPunctures.h"

int main(int argc, char* argv[]) {

  char * inputfile = NULL;
  if(argc == 2){
    inputfile = argv[1];
    TwoPunctures_params_set_inputfile(inputfile);
  } else {
    // revert to defaults if no input file passed
    TwoPunctures_params_set_default();  // must be set initially
  }

  printf("Input file: %s \n", inputfile);

  /*
    Test parameter injection.
  */
  TwoPunctures_params_set_Boolean("verbose", true);
  TwoPunctures_params_set_Real("par_b", 1.0);
  TwoPunctures_params_set_Real("par_m_plus", 2.0);

  /*
    Sans parameter injection the following is equivalent to prior all NULL call
  */
  ini_data data = TwoPunctures_make_initial_data();

  // inspect struct
  int ix=0;
  for(ix=0; ix < 20; ix++){
    printf("data.F[%d]: %lf\n", ix, data.F[ix]);
  }

  // example of interpolating gxx and psi to a few nodes
  int imin[3] = {0, 0, 0};
  int imax[3] = {2, 2, 2};
  int n[3] = {2, 2, 2};

  int sz = (n[0]) * (n[1]) * (n[2]);

  double *gxx = (double *) malloc(sz * sizeof(double));
  double *psi = (double *) malloc(sz * sizeof(double));
  double *tmp = (double *) malloc(sz * sizeof(double));

  double x[2] = {-18.333333333333333, -15.0};
  double y[2] = {-15.0, -5.0};
  double z[2] = {-15.0, -5.0};


  TwoPunctures_Cartesian_interpolation
    (data,     // struct containing the previously calculated solution
     imin,         // min, max idxs of Cartesian Grid in the three directions
     imax,         // in the three dirs
     n,
     x,         // Cartesian coordinates
     y,
     z,
     tmp,       // lapse
     psi,       // conformal factor and derivatives
     NULL,      // psix
     NULL,      // psiy
     NULL,      // psiz
     NULL,      // psixx
     NULL,      // psixy
     NULL,      // psixz
     NULL,      // psiyy
     NULL,      // psiyz
     NULL,      // psizz
     gxx,       // metric components
     tmp,       // gxy
     tmp,       // gxz
     tmp,       // gyy
     tmp,       // gyz
     tmp,       // gzz
     tmp,       // kxx extrinsic curvature components
     tmp,       // kxy
     tmp,       // kxz
     tmp,       // kyy
     tmp,       // kyz
     tmp);      // kzz

  int flat_ix = 0;
  for(flat_ix=0; flat_ix<sz; flat_ix++){
    double psi4 = pow(psi[flat_ix], 4);
    printf("\ngxx[%d]=%4.8f\n", flat_ix, psi4 * gxx[flat_ix]);
    printf("\npsi[%d]=%4.8f\n", flat_ix, psi[flat_ix]);
  }

  // make sure to take care of any internal memory that must be freed!
  TwoPunctures_finalise(data);

  return 0;
}

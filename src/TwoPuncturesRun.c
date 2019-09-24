// Minimal caller

#include "TwoPunctures.h"

#define CARTESIAN_INTERP (0)

int main(int argc, char* argv[]) {

  char * inputfile = NULL;
  if (argc == 2) inputfile = argv[1];

  printf("Input file: %s \n", inputfile);

#if(CARTESIAN_INTERP)

  // example how to call it, does not work here
  TwoPunctures (inputfile, // file to set input parameters
                imin,imax, // min/max indexes of Cartesian grid in the 3 direc.
                nxyz,
                x,   y,   z, // Catersian coords
                alp, // lapse
                psi, // conf factor, and drvts:
                psix,   psiy,   psiz,
                psixx,   psixy,   psixz,
                psiyy,   psiyz,   psizz,
                // metric
                gxx,   gxy,   gxz,
                gyy,   gyz,   gzz,
                // curv
                kxx,   kxy,   kxz,
                kyy,   kyz,   kzz);

#else

  // interpolation is turned off, see TwoPunctures.h
  TwoPunctures (inputfile, // file to set input parameters
                NULL, NULL,
                NULL,
                NULL, NULL, NULL,
                // lapse
                NULL,
                // conf factor
                NULL,
                NULL, NULL, NULL,
                NULL, NULL, NULL,
                NULL, NULL, NULL,
                // metric
                NULL, NULL, NULL,
                NULL, NULL, NULL,
                // curv
                NULL, NULL, NULL,
                NULL, NULL, NULL
                );

#endif



  return 0;
}

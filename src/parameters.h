/* Parameters are temporary set ths way */

#ifndef TP_PARAMS_H
#define TP_PARAMS_H

// ---------------------------------------------
// TwoPunctures parameters
// ---------------------------------------------

// mass of the m+ puncture
#define par_m_plus (1.0)

// mass of the m- puncture
#define par_m_minus (1.0)

// x coordinate of the m+ puncture
#define par_b (1.0) 

// target ADM mass for m+
#define target_M_plus (0.5)

// target ADM mass for m-
#define target_M_minus (0.5)

// momentum of the m+ puncture
#define par_P_plus1 (0.)
#define par_P_plus2 (0.)
#define par_P_plus3 (0.)
extern double par_P_plus[3]; 

// momentum of the m- puncture
#define par_P_minus1 (0.)
#define par_P_minus2 (0.)
#define par_P_minus3 (0.)
extern double par_P_minus[3]; 

// spin of the m+ puncture
#define par_S_plus1 (0.)
#define par_S_plus2 (0.)
#define par_S_plus3 (0.)
extern double par_S_plus[3];

// spin of the m- puncture
#define par_S_minus1 (0.)
#define par_S_minus2 (0.)
#define par_S_minus3 (0.)
extern double par_S_minus[3];

// offset b=0 to position (x,y,z)
#define center_offset1 (0.)
#define center_offset2 (0.)
#define center_offset3 (0.)
extern double center_offset[3];

// User provides bare masses rather than target ADM masses
#define give_bare_mass (1)  // 1/0 = yes/no

// Tolerance of ADM masses when give_bare_mass=no
#define adm_tol (1e-10) 

// How to fill the 3D grid from the spectral grid ? 
enum{
  taylor_expansion,// use a Taylor expansion about the nearest collocation point (fast, but might be inaccurate)
  evaluation, // evaluate using all spectral coefficients (slow)
};
#define grid_setup_method (taylor_expansion) 

// How to set lapse 
enum{
  antisymmetric, // antisymmetric lapse for two puncture black holes, -1 <= alpha <= +1
  averaged, // averaged lapse for two puncture black holes, 0 <= alpha <= +1"
  psin, // Based on the initial conformal factor
  brownsville, // See Phys. Rev. D 74, 041501 (2006)
};
#define initial_lapse (psin)

// Exponent n for psi^-n initial lapse profile (<0)
#define initial_lapse_psi_exponent (-2.0)

// Number of coefficients in the compactified radial direction
#define npoints_A (30)

// Number of coefficients in the angular direction
#define npoints_B (3)

// Number of coefficients in the phi direction
#define npoints_phi (16)

// Tolerance for Newton solver
#define Newton_tol (1e-10)

// Maximum number of Newton iterations
#define Newton_maxit (5)

// A small number to smooth out singularities at the puncture locations
#define TP_epsilon (0.0)

// Tiny number to avoid nans near or at the pucture locations 
#define TP_Tiny (0.0) 

// Radius of an extended spacetime instead of the puncture
#define TP_Extend_Radius (0.0)

// Swap x and z coordinates when interpolating,
// so that the black holes are separated in the z direction
#define swap_xz (0) // 1/0 = yes/no

// Use sources?
#define use_sources (0) // 1/0 = yes/no

// "If sources are used - rescale them after solving?
#define rescale_sources (1) // 1/0 = yes/no

// Set initial guess by external function?
#define use_external_initial_guess (0) // 1/0 = yes/no

// Output debug information about the residuum
#define do_residuum_debug_output (0)

// Output debug information about initial guess
#define do_initial_debug_output (0) 

// Multiply the old lapse with the new one
#define multiply_old_lapse (0)

// Solve for momentum constraint?
#define solve_momentum_constraint (0)

//
#define verbose (0)

// ---------------------------------------------
// Additional parameters
// ---------------------------------------------

// Want to test the interpolation ?
#define TEST_CARTESIAN_INTERP (0)

// Cartesian grid parameters
#define npointsx (50)
#define npointsy (50)
#define npointsz (50)
#define GFINDEX3D (i, j, k) ( i + npointsx*j + npointsx*npointsz*k )

#endif

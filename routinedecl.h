
#include <complex.h>

#include <fftw3.h>

#include <f2c.h>
#include <clapack.h>
#include <blaswrap.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>  
#include <gsl/gsl_poly.h>
#include <gsl/gsl_rng.h>

gsl_complex root[Nx][Ny][Dim];
gsl_vector_complex * evector[Nx][Ny][Dim] ;
gsl_rng *gnalea_r3,*gnalea_r31,*gnalea_r32,*gnalea_r33,*gnalea_r4,*gnalea_r5; 

#include <time.h>



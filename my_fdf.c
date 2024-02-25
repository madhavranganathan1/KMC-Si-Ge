 
#include "variables.h"

void my_fdf (const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
     {
       *f = my_f(x, params);
       my_df(x, params, df);
     }

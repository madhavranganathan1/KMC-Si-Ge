 
#include "variables.h"

void my_fdfthomas (const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
     {
       *f = my_fthomas(x, params);
       my_dfthomas(x, params, df);
     }

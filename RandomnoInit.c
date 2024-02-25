#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "KMCheteroepitaxy.h"
#include "funcdecl.h"

/* This function initializes a Random number generator using the appropriate Seed. It returns a number one greater than the chosen Seed and
 * and this number can be used to seed the RNG the following time. The value of gnalea_r is set to the state of the RNG  */
 
unsigned long int RandomnoInit(gsl_rng *gnalea_r,unsigned long int Seed)
{
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	gsl_rng_default_seed=abs(Seed);
	T = gsl_rng_mt19937; /*gsl_rng_default;*/
	gnalea_r = gsl_rng_alloc(T);
	return Seed+1;
}

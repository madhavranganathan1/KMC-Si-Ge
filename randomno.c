#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "KMCheteroepitaxy.h"
#include "funcdecl.h"
double randomno(gsl_rng *gnalea_r)
{
	double x;
	
	x=gsl_rng_uniform_pos(gnalea_r);	
	
	return x;
}

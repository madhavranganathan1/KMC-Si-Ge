#include "variables.h"

void 1_initialconditions()
{
    int In,J,K,Zlimit,Stepseparation,Terraceheight  ; 
    unsigned long int Seed1;
    Seed1 = time(NULL);
    gsl_rng *gnalea_r ;
    const gsl_rng_type *T;
    T = gsl_rng_mt19937;
	gsl_rng_env_setup();
	gsl_rng_default_seed=abs(Seed1)+rank;
	gnalea_r = gsl_rng_alloc(T);

	  for (In= 0;In < Nx; In++) {
	     for (J =0;J < Ny; J++) {
	         Height[In][J] = DiscreteAlayers-1+Film;}}
    
    for (In=0;In < Nx; In++){
      for (J=0; J < Ny; J++) {
    	for (K=0; K < Nz; K++) {
           if ((In == Nx/2) && (J==Ny/2) && (K==0)) {
                       	    displacementType[In][J][K].xdisp = 0.0 ;
    			    displacementType[In][J][K].ydisp = 0.0 ;
    			    displacementType[In][J][K].zdisp = 0.0 ;}
     else if (K <= Height[In][J]) {
                            displacementType[In][J][K].xdisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.5) ;
    		            displacementType[In][J][K].ydisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.5) ;
    		            displacementType[In][J][K].zdisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.5) ;}
     else { displacementType[In][J][K].xdisp = 0.0 ;
    		            displacementType[In][J][K].ydisp = 0.0 ;
    	                    displacementType[In][J][K].zdisp = 0.0 ;}}}}     
   
   for (In=0;In < Nx; In++){
     for (J=0; J < Ny; J++) {
    	Zlimit = Height[In][J] ;
    	for (K=0; K < Nz ; K++) {   
    	 if (K <= Zlimit) {
    	   if (K < DiscreteAlayers) {displacementType[In][J][K].Type = 1 ;}
    	     else {displacementType[In][J][K].Type = 2 ;}}
    		else {displacementType[In][J][K].Type = 0;}}}}

}

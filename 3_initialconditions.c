#include "variables.h"

void initialconditions()
{
    int In,J,K;
    unsigned long int Seed1;
    Seed1 = time(NULL);
    gsl_rng *gnalea_r ;
    const gsl_rng_type *T;
    T = gsl_rng_mt19937;

        gsl_rng_env_setup();
        gsl_rng_default_seed=abs(Seed1)+rank;
        gnalea_r = gsl_rng_alloc(T);

    char dTfinput[256];
    sprintf(dTfinput,"500K20500x0",rank);

  FILE *fp;
  fp=fopen(dTfinput,"rb");
  fread(displacementType,sizeof(struct dT), (Nx*Ny*Nz),fp); 

for (In=0;In<Nx; In++){
     for (J=0;J<Nx;J++){
       for(K=0;K<Nz;K++){
          if(displacementType[In][J][K].Type==0){ Height[In][J]=K-1;
                                                   break;}}}}
for (In=0;In < Nx; In++){
      for (J=0; J < Ny; J++) {
        for (K=0; K < Nz; K++) {
           if ((In == Nx/2) && (J==Ny/2) && (K==0)) {
                            displacementType[In][J][K].xdisp = 0.0 ;
                            displacementType[In][J][K].ydisp = 0.0 ;
                            displacementType[In][J][K].zdisp = 0.0 ;}
     else if (K <= Height[In][J]) {
                            displacementType[In][J][K].xdisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.05) ;
                            displacementType[In][J][K].ydisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.05) ;
                            displacementType[In][J][K].zdisp = 0.01*(gsl_rng_uniform_pos(gnalea_r) - 0.05) ;}
     else { displacementType[In][J][K].xdisp = 0.0 ;
                            displacementType[In][J][K].ydisp = 0.0 ;
                            displacementType[In][J][K].zdisp = 0.0 ;}}}}
fclose(fp);
}

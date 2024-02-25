#include "variables.h"
int Add_atom_optimization(void)
{
 unsigned long int Seed2;
  Seed2 = time(NULL);
  gsl_rng *gnalea_r1, *gnalea_r2;
  const gsl_rng_type *T;
  T = gsl_rng_mt19937;
  gsl_rng_env_setup();
  gsl_rng_default_seed=abs(Seed2)+rank;
  gnalea_r1 = gsl_rng_alloc(T);
  Seed2++;

   gsl_rng_env_setup();
   gsl_rng_default_seed=abs(Seed2)+rank;
   gnalea_r2 = gsl_rng_alloc(T);

int Nomins,Minimizetest,In,J,K;
       Nomins = 0;
        do {
          Minimizetest = energyminimization();
          Nomins +=1;}
         while ((Minimizetest != 0) && (Nomins < Maxnomintry)) ;
       
//perturbe
 if (Nomins >= 1*Maxnomintry ) { printf("failed - perturbing, In the Add_atom_optimization\n");
  ss=ss+0.01,tg=tg+0.1,gt=gt+0.1;
         do {
                 displacementType[Nx/2][Ny/2][Height[Nx/2][Ny/2]-1].xdisp += gsl_rng_uniform_pos(gnalea_r1)*mismatch;                                                          Minimizetest = energyminimization();
                 Nomins += 1;}
         while( ( Minimizetest !=0) && ( Nomins < 2*Maxnomintry)) ;}

//perturbe
 if (Nomins >= 2*Maxnomintry ) { printf("failed - perturbing,  In the Add_atom_optimization\n");
  ss=ss+0.01,tg=tg+0.1,gt=gt+0.1;
         do {
                 displacementType[Nx/2][Ny/2][Height[Nx/2][Ny/2]-1].xdisp += gsl_rng_uniform_pos(gnalea_r2)*mismatch;                                                          Minimizetest = energyminimization();
                 Nomins += 1;}
         while( ( Minimizetest !=0) && ( Nomins < 3*Maxnomintry)) ;}

if (Nomins >= 3*Maxnomintry ) {
    printf("Optimization failed Even After All Try\n");
    printf("Setting All displacements to Zero\n");

    for (In=0;In < Nx; In++){
      for (J=0; J < Ny; J++) {
        for (K=0; K < Nz; K++) {
                            displacementType[In][J][K].xdisp = 0.0 ;
                            displacementType[In][J][K].ydisp = 0.0 ;
                            displacementType[In][J][K].zdisp = 0.0 ;}}}
                            }
//if(coverag%1638==0){printf("Done with Energy minimization %lf\n",energyvalue);}
//if(rank==0){printf("Add_atom_optimization %lf\n",energyvalue);}
 ss=0.01,tg=0.1,gt=0.1;

return 0;
}

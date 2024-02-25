//#include "variables.h"
//#include"energyminimization.c"

int pswap_optimization(int a,int b,int c,int d,int e,int f)
{
  unsigned long int Seed1;
  Seed1 = time(NULL);
  gsl_rng *gnalea_r1, *gnalea_r2;
  const gsl_rng_type *T;
  T = gsl_rng_mt19937;
  gsl_rng_env_setup();
  gsl_rng_default_seed=abs(Seed1);//+rank;
  gnalea_r1 = gsl_rng_alloc(T);
  Seed1++;

   gsl_rng_env_setup();
   gsl_rng_default_seed=abs(Seed1);//+rank;
   gnalea_r2 = gsl_rng_alloc(T); 
   displacementType[d][e][f].Type =2;
   displacementType[a][b][c].Type =1;


int Nomins,Minimizetest;
        Nomins = 0;
        do {
          Minimizetest = energyminimization();
          Nomins +=1;}
         while ((Minimizetest != 0) && (Nomins < Maxnomintry)) ;
       
//perturbe
 if (Nomins >= 1*Maxnomintry ) { printf("failed - perturbing,  In the optimization \n");
         do {
                 displacementType[Nx/2][Ny/2][Height[Nx/2][Ny/2]-1].xdisp += gsl_rng_uniform_pos(gnalea_r1)*mismatch;                                                          Minimizetest = energyminimization();
                 Nomins += 1;}
         while( ( Minimizetest !=0) && ( Nomins < 2*Maxnomintry)) ;}

//perturbe
 if (Nomins >= 2*Maxnomintry ) { printf("failed - perturbing,  In the optimization \n");
         do {
                 displacementType[Nx/2][Ny/2][Height[Nx/2][Ny/2]-1].xdisp += gsl_rng_uniform_pos(gnalea_r2)*mismatch;                                                          Minimizetest = energyminimization();
                 Nomins += 1;}
         while( ( Minimizetest !=0) && ( Nomins < 3*Maxnomintry)) ;}

if (Nomins >= 3*Maxnomintry ) {
    printf("Optimization failed Even After All Try\n");
     FILE *tmpout;
     char tmpfile[256];
   /*   sprintf(tmpfile,"tmpdisptype%dx%d%xd",Nx,Ny,Nz); */                                  
     sprintf(tmpfile,"tmpdisptype%dx%d%xd%xd",Nx,Ny,Nz,rank);
     tmpout = fopen(tmpfile,"wb");
    fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),tmpout);
    fclose(tmpout);
    exit(0);
      }

//if(rank==0){printf("optimization %lf\n",energyvalue);}
return 0;
}

#include <omp.h>
#include "header.h"
//#include <mpi.h>
#include<string.h>

extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz];
 
int main()//( int argc, char *argv[])
{
omp_set_num_threads(32);
/*MPI_Init(&argc,&argv);
MPI_Comm_size( MPI_COMM_WORLD, &size );
MPI_Comm_rank( MPI_COMM_WORLD, &rank );*/
 ss=0.01,tg=0.1,gt=0.1;
   
   clock_t start, end;
   double cpu_time_used;
   start = clock();
   time_t t1=time(NULL);

/*      time_t tim=time(NULL);
             char *s=ctime(&tim);
             s[strlen(s)-1]=0;   
             if(rank==0){printf("it is %s now.\n", s);}*/

     int  In,J, Printstep,AddatomNo,AddatomX,AddatomY;
    
     double fluxtester ;  /* For atom deposited, this decides the type of atom */
     double runtime;
     unsigned long int Seed;
     unsigned long int Sum=0L;
     gsl_rng *gnalea_r1,*gnalea_r2  ;

     
     /* Calculate all the spring lengths here  */

     lespringDAA = lespringLAA * pow(2.0,0.5) ;
     lespringLBB = lespringLAA * (1.0 + mismatch) ;
     lespringDBB = lespringLBB * pow(2.0,0.5) ;
     lespringLAB = (lespringLAA *  relativefluxA) + (lespringLBB * (1.0 - relativefluxA)) ;
     lespringDAB = (lespringDAA *  relativefluxA) + (lespringDBB * (1.0 - relativefluxA)) ;

         arefx = lespringLAA ;
	 arefy = lespringLAA ;   /* Mechanical equilibrium of completed layers shows relaxation in the Z direction only  */
	 arefzA = lespringLAA ;


	 /* arefzB is taken as a linear combination of the Smereka length (linearized with common spring constants )
	 * and the pure A spacing as there is also a flux of A from the vapor */

     arefzB = lespringLAA * (1.0 + ( mismatch *(kspringLBB+(4.0 * kspringDBB))/ (kspringLBB+(2.0 * kspringDBB)) ) );
     arefzAB = lespringLAA * (1.0 + ( mismatch *(kspringLBB+(4.0 * kspringDBB))/ (kspringLBB+(2.0 * kspringDBB)) )) ;

     Eigenvectors();

    initialconditions();
     optimization();
      
  // minimizes energy and generates displacements and a value of energy 

if(rank==0){printf("Done with Energy minimization %lf\n",energyvalue);
     printf("Number of Atoms going to be Added %d Number of processors in use %d \n",FLUXend,size);}

     FILE *Typefile;
     
     Seed = time(NULL);

     const gsl_rng_type *T;
     T = gsl_rng_mt19937;

	 gsl_rng_env_setup();
	 gsl_rng_default_seed=abs(Seed)+rank;
	 gnalea_r1 = gsl_rng_alloc(T);  //This RNG is used to decide where the Addatom should be DEPOSITED 
	 Seed+=1;

	 gsl_rng_env_setup();
	 gsl_rng_default_seed=abs(Seed)+rank;
	 gnalea_r2 = gsl_rng_alloc(T); // This RNG is used to determine the TYPE of addatom added by the external flux 
	 Seed+=1;

	 gsl_rng_env_setup();
	 gsl_rng_default_seed=abs(Seed)+rank;
	 gnalea_r3 = gsl_rng_alloc(T); // This RNG is used in kmcstep to determine the addatom chosen to ATTEMPT a hop 
	 Seed+=1;

          gsl_rng_env_setup();
         gsl_rng_default_seed=abs(Seed)+rank;
         gnalea_r31 = gsl_rng_alloc(T); // This RNG is used in kmcstep to determine the addatom chosen to ATTEMPT a hop 
         Seed+=1;

          gsl_rng_env_setup();
         gsl_rng_default_seed=abs(Seed)+rank;
         gnalea_r32 = gsl_rng_alloc(T); // This RNG is used in kmcstep to determine the addatom chosen to ATTEMPT a hop 
         Seed+=1;

          gsl_rng_env_setup();
         gsl_rng_default_seed=abs(Seed)+rank;
         gnalea_r33 = gsl_rng_alloc(T); // This RNG is used in kmcstep to determine the addatom chosen to ATTEMPT a hop 
         Seed+=1;


	 gsl_rng_env_setup();
	 gsl_rng_default_seed=abs(Seed)+rank;
	 gnalea_r4 = gsl_rng_alloc(T); // This RNG is used to decide the DIRECTION of the hop of an addatom that has been approved 
	 Seed+=1;

	 gsl_rng_env_setup();
	 gsl_rng_default_seed=abs(Seed)+rank;
	 gnalea_r5 = gsl_rng_alloc(T); // This RNG is used in the METROPOLIS transition probability 
	 Seed+=1;

     B1=0,B2=0,B3=0;//KMC Accepted Move counter

     int FLUXlength=FLUXend;
     int TR=Nx*Ny*ML;

//     for  (J =1; J <= TR; J++){
  //               if(J%1==0){intermixing();}
   //  for  (In =1; In <=Totalsweeps; In++)
   //  {Sum=Sum+1;
/*        if(((In-1) % Naddatom) == 0 && FLUXlength>0) {//1st if statement
        if(rank==0){printf("# atom added %d \n",(FLUXend-FLUXlength));}
             coverag=(FLUXend-FLUXlength);
             FLUXlength--;  
             AddatomNo = (int) (Nx*Ny*gsl_rng_uniform_pos(gnalea_r1)) ;

         	  AddatomX = (int) ((int) (AddatomNo) / Ny );
 	          AddatomY = ((int) (AddatomNo)) % Ny ;
         	  Height[AddatomX][AddatomY] += 1;  


         	  fluxtester = gsl_rng_uniform_pos(gnalea_r2);

         	  if (fluxtester < relativefluxA) {//2nd if statement
         	  	    displacementType[AddatomX][AddatomY][Height[AddatomX][AddatomY]].Type = 1 ;
                            Add_atom_optimization(); 

         	                                    }//end of 2nd if statement
         	  else {
             	  	displacementType[AddatomX][AddatomY][Height[AddatomX][AddatomY]].Type = 2 ;
                        Add_atom_optimization(); 
                       }
         }*/  //End of 1st if statement

     //    kmcstep(In);  
       //  }//In loop end
//  if(rank==0){
         if ((J % (prnt)==0)) {  //give which coverages to print= 500/(100*100)=.05
                                            int DA=J/(prnt);
                                            char dTfname[256];
                                            sprintf(dTfname,"700HH%d",J+52224);
                                            Typefile = fopen(dTfname,"w");
                                            fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),Typefile);
                                            fclose(Typefile);
                              }
  //           }
//     if(rank==0){if(J%prnt==0){printf("TS %lu B1 %d B2 %d B3 %d\n",Sum,B1,B2,B3);}}
     }//J loop end


   FILE *Inputfp;
   char dTInput[256];
   sprintf(dTInput,"Input1x%d",rank);
   Inputfp = fopen(dTInput,"w");
   fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),Inputfp);
   fclose(Inputfp);

  time_t t2=time(NULL);
   end = clock();
   cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
   printf(" Actual Time elapsed %f min \n",(cpu_time_used/60.));
   double et=cpu_time_used/60.;
   printf("e=%e\tf=%e\tg=%e\n",ssizeastep, tolaglobal,gtolagrad);


/*   for(In=0;In<3;In++)

 {
   free(I_co[In]);}*/
   
/*   FILE *E;
   char ET[256];
   sprintf(ET,"Totaltime");
   E=fopen(ET,"a");
   fprintf(E," elapsed time for %dx%dx%d_%drank system is %f min \n",Nx,Ny,Nz,rank,et);
   fclose(E);*/

  /*   time_t tim1=time(NULL);
             char *s1=ctime(&tim1);
             s1[strlen(s)-1]=0;
             if(rank==0){printf("it is %s now.\n", s1);}*/

 //  system("cp  ./Input1x* ./pt2ML");
 //  printf("This program is completed \n");
//   MPI_Finalize();


return 0;
}// end of main()
/* ---------------------END OF MAIN PROGRAM --------------------------------*/

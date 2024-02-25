#include "variables.h"
extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

void poorman_KMC(int Timestep) 
{
 	int NNbondsAA,NNNbondsAA,NNbondsBB,NNNbondsBB,NNbondsAB,NNNbondsAB,Atomnumber, AtomX, AtomY,AtomZ,Nomin,Mintest;
 	double initialenergy,energybond,energyvalue;
 	
 	Atomnumber = (int) (Nx*Ny*gsl_rng_uniform_pos(gnalea_r3) ) ;
 	AtomX = (int) (Atomnumber/Ny) ;
 	AtomY = Atomnumber % Ny ; 
 	if ( (AtomX < 0) || (AtomX >= Nx ) || (AtomY < 0) || (AtomY >= Ny) ) {
 		printf("Problem with Random Atom choice\n");
 		exit(0);}
 	
 	AtomZ = Height[AtomX][AtomY];
 	
 	if ((displacementType[AtomX][AtomY][AtomZ].Type) == 0) {
 	    printf("Something wrong with adatom choice; its type is zero %d %d %d\n",AtomX,AtomY,AtomZ);
 	    exit(0);}
 	
 	if (AtomZ == 0) {
 		printf ("Reached continuous substrate\n");
 		exit(0);}
 	
 	if ((displacementType[AtomX][AtomY][AtomZ+1].Type) != 0) {
 	     printf("Problem with Random Atom choice; there is an atom above it!!! %d %d %d\n",AtomX,AtomY,AtomZ);
 		exit(0);}
 	
 	/* initialize all bonds to zero */
 	NNbondsAA = 0;
 	NNNbondsAA = 0;
 	NNbondsBB = 0;
 	NNNbondsBB = 0;
 	NNbondsAB = 0;
 	NNNbondsAB = 0;

 int	LAAH = 0, LBBH = 0, LABH = 0, DAAH = 0, DBBH = 0, DABH = 0, DAAHm = 0, DBBHm = 0, DABHm = 0, DAAHp = 0, DBBHp = 0, DABHp = 0;

   int typ=displacementType[AtomX][AtomY][AtomZ].Type;

        LBH(typ,AtomZ,AtomX,AtomY,&LAAH,&LBBH,&LABH);
        DBH(typ,AtomZ,AtomX,AtomY,&DAAH,&DBBH,&DABH);
        DBHm(typ,AtomZ,AtomX,AtomY,&DAAHm,&DBBHm,&DABHm);
        DBHp(typ,AtomZ,AtomX,AtomY,&DAAHp,&DBBHp,&DABHp);

 	NNbondsAA = LAAH;
 	NNNbondsAA = DAAH+DAAHm+DAAHp;
 	NNbondsBB = LBBH;
 	NNNbondsBB = DBBH+DBBHm+DBBHp;
 	NNbondsAB = LABH;
 	NNNbondsAB = DABH+DABHm+DABHp;
 	
 	/* That ends the calculation of the number of bonds broken */
 	energybond = (NNbondsAA * gammaNNAA) + (NNbondsBB * gammaNNBB) + (NNbondsAB * gammaNNAB) + (NNNbondsAA * gammaNNNAA) + (NNNbondsBB * gammaNNNBB) + (NNNbondsAB * gammaNNNAB) ;

          
	 			
 	if (KMCtype == 1) {             /* Smereka-Russo type KMC */

 	      int No_of_Bobd=Totalbond(AtomX,AtomY);
                    //printf("bond %d \n",No_of_Bobd);
                 if(No_of_Bobd<=5){ 
                    update(); 
                        } 

 		 else if (KMCSmereka(energybond,AtomX,AtomY,AtomZ) == 1) {  /*move accepted */
                    //printf("bond %d \n",No_of_Bobd);
                         update(AtomX,AtomY);
                                                                  }
 		
 	}// end of if 
 	else if (KMCtype == 2) {         /* Accurate but slower KMC */
 	
 	/*	for (I = 0; I < Nx; I ++){
 		    for (J = 0; J < Ny; J ++){
 			    for (K=0; K < Nz; K++) {
 				   initialdisplacementType[I][J][K].xdisp = displacementType[I][J][K].xdisp ;
 				   initialdisplacementType[I][J][K].ydisp = displacementType[I][J][K].ydisp ;
 				   initialdisplacementType[I][J][K].zdisp = displacementType[I][J][K].zdisp ;
 				   initialdisplacementType[I][J][K].Type = displacementType[I][J][K].Type ;
 			    }
 		    }
 	    }
 	*/
 	
 	    initialenergy = energyvalue ;
 	    KMCcomparer(energybond, AtomX,AtomY,AtomZ,initialenergy) ; /* Move accepted  or rejected inside KMCComparer */
 	    	
 		
 	}

 	if ((Timestep % Nhopsperequilibrium ) == 0) {
                    kmc_optimization();}
 	
}//end of kmcstep
 
/* ---------------------END OF KMCSTEP ------------------------------------------------*/
 
 

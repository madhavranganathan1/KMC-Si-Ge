
#include "variables.h"

extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

    /* May be used to seed the random number generator   */

 /*        kmcstep           */
 /* This program implements the entire KMC operation. It first selects an atom at random to displace. It consists of a KMCstepper and a KMCcomparer For a rejection 
  * based algorithm, there is only one step that the KMC needs to consider that is complete dissociation  and move 
  * the atom up to KMCstepunit steps in any direction */       

void kmcstep(int Timestep) 
{
 	int NNbondsAA,NNNbondsAA,NNbondsBB,NNNbondsBB,NNbondsAB,NNNbondsAB,Atomnumber, AtomX, AtomY,AtomZ,Nomin,Mintest;
 	double initialenergy,energybond,energyvalue,toss;
 	/* Copy the displacementType structure into an initial displacementType structure */
 	
 	
 	
 	/* Now we need to select an atom to hop */ 
 	
 	
 	Atomnumber = (int) (Nx*Ny*gsl_rng_uniform_pos(gnalea_r3) ) ;
 	AtomX = (int) (Atomnumber/Ny) ;
 	AtomY = Atomnumber % Ny ; 
 	if ( (AtomX < 0) || (AtomX >= Nx ) || (AtomY < 0) || (AtomY >= Ny) ) {
 		printf("Problem with Random Atom choice\n");
 		exit(0);
 	}
 	
 		
 	/* To calculate NNbonds = number of nearest neighbour bonds 
 	 * Need Displacement function in the vicinity of the atom */
 	
 	AtomZ = Height[AtomX][AtomY];
 	/*printf("Atom Choice %d %d %d\n",AtomX,AtomY,AtomZ); */
 	
 	/* Check that you have chosen a top layer atom that is not the bottom of the simulation layers */
 	if ((displacementType[AtomX][AtomY][AtomZ].Type) == 0) {
 		
 	    printf("Something wrong with adatom choice; its type is zero %d %d %d\n",AtomX,AtomY,AtomZ);
 	    exit(0);
 	}
 	
 	if (AtomZ == 0) {
 		printf ("Reached continuous substrate\n");
 		exit(0);
 	}
 	
 	
 	if ((displacementType[AtomX][AtomY][AtomZ+1].Type) != 0) {
 	     printf("Problem with Random Atom choice; there is an atom above it!!! %d %d %d\n",AtomX,AtomY,AtomZ);
 		exit(0);
 	}
 	
 	/* initialize all bonds to zero */
 	NNbondsAA = 0;
 	NNNbondsAA = 0;
 	NNbondsBB = 0;
 	NNNbondsBB = 0;
 	NNbondsAB = 0;
 	NNNbondsAB = 0;
 	
 	/* Now we need to calculate all the bonds broken */
 	
 	if (displacementType[AtomX][AtomY][AtomZ].Type == 1) {
 		
 		/* First we look at nearest neighbours */
 		
 		/* Neighbour Z-1 */
 		if  (displacementType[AtomX][AtomY][AtomZ -1].Type == 1){
 			NNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][AtomY][AtomZ- 1].Type == 2){
 			NNbondsAB +=1;
 		}
 		
 		/* Neighbour X-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ].Type == 1){
 			NNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ].Type == 2){
 			NNbondsAB +=1;
 		} 		
 		
 		/* Neighbour X+1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ].Type == 1){
 			NNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ].Type == 2){
 			NNbondsAB +=1;
 		} 
 		
 		/* Neighbour Y-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNbondsAB +=1;
 		} 		
 		
 		/* Neighbour Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNbondsAB +=1;
 		} 
 		
 		
 		/* Now we look at Next Nearest Neighbours */

		/* Neighbour X-1 Y-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNNbondsAB +=1;
 		} 
 
        /* Neighbour X-1 Y+1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour X+1 Y-1 */
 		if  (displacementType[(AtomX+1) % Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour X+1 Y+1 */
 		if  (displacementType[(AtomX+1) % Nx][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour Y-1 Z-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ-1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ-1].Type == 2){
 			NNNbondsAB +=1;
 		} 
 
        /* Neighbour Z-1 Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ-1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ-1].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour Z+1 Y-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ+1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ+1].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour Z+1 Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ+1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ+1].Type == 2){
 			NNNbondsAB +=1;
 		} 

       /* Neighbour X-1 Z-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ-1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ-1].Type == 2){
 			NNNbondsAB +=1;
 		} 
 
        /* Neighbour X-1 Z+1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ+1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ+1].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour X+1 Z-1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ-1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ-1].Type == 2){
 			NNNbondsAB +=1;
 		} 

        /* Neighbour X+1 Z+1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ+1].Type == 1){
 			NNNbondsAA +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ+1].Type == 2){
 			NNNbondsAB +=1;
 		} 

 	} /* End of loop over Type == 1 */
 	
 	else if (displacementType[AtomX][AtomY][AtomZ].Type == 2) {
 		
 		/* First we look at nearest neighbours */
 		
 		/* Neighbour Z-1 */
 		if  (displacementType[AtomX][AtomY][AtomZ-1].Type == 1){
 			NNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][AtomY][AtomZ-1].Type == 2){
 			NNbondsBB +=1;
 		}
 		
 		/* Neighbour X-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ].Type == 1){
 			NNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ].Type == 2){
 			NNbondsBB +=1;
 		} 		
 		
 		/* Neighbour X+1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ].Type == 1){
 			NNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ].Type == 2){
 			NNbondsBB +=1;
 		} 
 		
 		/* Neighbour Y-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNbondsBB +=1;
 		} 		
 		
 		/* Neighbour Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNbondsBB +=1;
 		} 
 		
 		
 		/* Now we look at Next Nearest Neighbours */
 		
 		/* Neighbour X-1 Y-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNNbondsBB +=1;
 		} 
 
        /* Neighbour X-1 Y+1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour X+1 Y-1 */
 		if  (displacementType[(AtomX+1) % Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][(AtomY-1+Ny)%Ny][AtomZ].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour X+1 Y+1 */
 		if  (displacementType[(AtomX+1) % Nx][(AtomY+1)% Ny][AtomZ].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][(AtomY+1)% Ny][AtomZ].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour Y-1 Z-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ-1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ-1].Type == 2){
 			NNNbondsBB +=1;
 		} 
 
        /* Neighbour Z-1 Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ-1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ-1].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour Z+1 Y-1 */
 		if  (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ+1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY-1+Ny)%Ny][AtomZ+1].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour Z+1 Y+1 */
 		if  (displacementType[AtomX][(AtomY+1)% Ny][AtomZ+1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[AtomX][(AtomY+1)% Ny][AtomZ+1].Type == 2){
 			NNNbondsBB +=1;
 		} 

       /* Neighbour X-1 Z-1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ-1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ-1].Type == 2){
 			NNNbondsBB +=1;
 		} 
 
        /* Neighbour X-1 Z+1 */
 		if  (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ+1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX-1+Nx)%Nx][AtomY][AtomZ+1].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour X+1 Z-1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ-1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ-1].Type == 2){
 			NNNbondsBB +=1;
 		} 

        /* Neighbour X+1 Z+1 */
 		if  (displacementType[(AtomX+1) % Nx][AtomY][AtomZ+1].Type == 1){
 			NNNbondsAB +=1;
 		}
 		else if (displacementType[(AtomX+1) % Nx][AtomY][AtomZ+1].Type == 2){
 			NNNbondsBB +=1;
 		} 
 		 		
 	} /* End of loop over Type == 2 */
 	
    
 	
 	/* That ends the calculation of the number of bonds broken */
 	energybond = (NNbondsAA * gammaNNAA) + (NNbondsBB * gammaNNBB) + (NNbondsAB * gammaNNAB) + (NNNbondsAA * gammaNNNAA) + (NNNbondsBB * gammaNNNBB) + (NNNbondsAB * gammaNNNAB) ;
 	 	 			
 	if (KMCtype == 1) {             /* Smereka-Russo type KMC */
 
 		if (KMCSmereka(energybond,AtomX,AtomY,AtomZ) == 1) {  /*move accepted */
 			
                /* Now we randomly want to move the atom chosen to one of the neighbouring sites */
                
                     
                
                toss = gsl_rng_uniform_pos(gnalea_r4);
                
                if (toss < 0.25) {
                	
                                     	
                     displacementType[(AtomX -1 + Nx)% Nx][AtomY][(Height[(AtomX -1 + Nx)% Nx][AtomY]) + 1].Type = displacementType[AtomX][AtomY][AtomZ].Type ;                	
                	 Height[(AtomX -1 + Nx)% Nx][AtomY] += 1 ;
                }
                else if (toss < 0.5) {
                	 
                	 displacementType[(AtomX+1) % Nx][AtomY][(Height[(AtomX+1) % Nx][AtomY]) + 1].Type = displacementType[AtomX][AtomY][AtomZ].Type ;                	
                	 Height[(AtomX +1)% Nx][AtomY] += 1;
                }
                else if (toss < 0.75) {
                	 
                	 displacementType[AtomX][(AtomY-1+Ny) % Ny][(Height[AtomX][(AtomY-1+Ny) % Ny]) + 1].Type = displacementType[AtomX][AtomY][AtomZ].Type ;                	
                	 Height[AtomX][(AtomY - 1 + Ny)%Ny] += 1 ;
                }
                else if (toss < 1.0001) {
                	 displacementType[AtomX][(AtomY + 1) % Ny][(Height[AtomX][(AtomY + 1) % Ny]) + 1].Type = displacementType[AtomX][AtomY][AtomZ].Type ;                	
                	 Height[AtomX][(AtomY +1)% Ny] += 1;
                }
                
                displacementType[AtomX][AtomY][AtomZ].Type = 0 ;
 			    Height[AtomX][AtomY] -= 1;
                
 			    
 		}
 		
 	}
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
 			    	 
 			    	  Nomin = 0;
 			    	  do {
 			    	       Mintest = energyminimization();
 			    	       Nomin += 1;
 			    	       
 			    	  }
 			    	  while( ( Mintest !=0) && ( Nomin < Maxnomintry)) ;
 			    	  if (Nomin >= Maxnomintry ) {
 			    	        printf("Energyminimization failed - Trying a little perturbation\n");
 			    	        do {
 			    	        	  displacementType[Nx/2][Ny/2][Height[Nx/2][Ny/2]-1].xdisp += gsl_rng_uniform_pos(gnalea_r4)*mismatch; /* Represents a small perturbation */
 			    	              Mintest = energyminimization();
 			    	              Nomin += 1;
 			    	       
 			    	        }
 			    	        while( ( Mintest !=0) && ( Nomin < 2*Maxnomintry)) ;
 			    	        if (Nomin >= (2*Maxnomintry) ) {
 			    	                printf("Energyminimization failed - Even after perturbation\n");
 			    	                exit(0);
 			    	        }
 			    	        
 			    	  }
    }
 	
}
 
/* ---------------------END OF KMCSTEP ------------------------------------------------*/
 
 

#include "variables.h"

double Harmonic(int Atomx,int Atomy, int Atomz) 
{
/*
    int neighbours=0,In;
    double HE[51];
      for(In=0;In<=50;In++){ HE[In]=0.00;}
 */
    double harmonicenergy,extensionsq,extension,kmcfactor,arefzup, arefzdown;
           harmonicenergy = 0.0 ;
     
           if (displacementType[Atomx][Atomy][Atomz].Type == 2) harmonicenergy -= relativehoppingenergy ; /* Biases B to hop a 
           factor of exp (-relativehoppingenergy) faster than A */
		/* First we create list of neighbors taking periodic boundary conditions and the presence
				 * of the continuum substrate into account */
			
	
	
		/* We have to list all the 18 neighbour interactions carefully. I will do it using if statements,
			 * which will make it long but make sure that everything is correct
			 * 
			 * Neighbours 1-6 represent lateral neighbours and neighbours 7-18 represent diagonal neighbours
			 * 
			 *   Neighbours 1-2  - Lateral X
			 *   Neighbours 3-4  - Lateral Y
			 *   Neighbours 5-6  - Lateral Z  - Neighbour 6 absent since this is a surface atom
			 * 
			 *   Neighbours 7-10   Diagonal YZ plane
			 *   Neighbours 11-14  Diagonal XZ plane
			 *   Neighbours 15-18  Diagonal XZ plane
			 * 
			 *  */
			 
			/* For neighbour1 */ 
			
			if(Atomz < DiscreteAlayers-1) {
				arefzup = arefzA;
				arefzdown=arefzA;
			}
			else if (Atomz > (DiscreteAlayers)) {
				arefzup = arefzB;
				arefzdown=arefzB;
			}
			if(Atomz == (DiscreteAlayers-1)) {
				arefzup = arefzAB;
				arefzdown = arefzA;
			}
			else {
				arefzup = arefzB;
				arefzdown = arefzAB;
			}
			
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+ Nx) % Nx][Atomy][Atomz].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) ) 				;
                extension = pow(extensionsq,0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension) * (lespringLAA -extension) ;
                //HE[0] = kspringLAA * (lespringLAA - extension) * (lespringLAA -extension) ;
                //neighbours++;
//printf("1 harmonicenergy1 %lf\n",harmonicenergy1);
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension) * (lespringLBB -extension) ;
                //HE[1] = kspringLBB * (lespringLBB - extension) * (lespringLBB -extension) ;
                //neighbours++;
//printf("2 harmonicenergy2 %lf \n",harmonicenergy2);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB - extension) * (lespringLAB -extension) ;
                //HE[2] = kspringLAB * (lespringLAB - extension) * (lespringLAB -extension) ;
                //neighbours++;
//printf("3 harmonicenergy3 %lf \n",harmonicenergy3);
			
			}        
			
			/* For neighbour2 */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 1)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) ) ;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension ) * ( lespringLAA - extension)  ;
                //HE[3] = kspringLAA * (lespringLAA - extension ) * ( lespringLAA - extension)  ;
                //neighbours++;
			
//printf("4 harmonicenergy4 %lf \n",harmonicenergy4);
			}
			
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 2)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) ) ;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension ) * (lespringLBB - extension ) ;
                //HE[4] = kspringLBB * (lespringLBB - extension ) * (lespringLBB - extension ) ;
                //neighbours++;
			
//printf("5 harmonicenergy5 %lf \n",harmonicenergy5);
			}
			
			else if (((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 2))  || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 1)) ) {
			
			    extensionsq = ((-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp)*(-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp))
			                + ((displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp)*(displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp))
			                + ((displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp)*(displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp)) ;
			    
			    extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB - extension) * (lespringLAB - extension) ;
                //HE[5]= kspringLAB * (lespringLAB - extension) * (lespringLAB - extension) ;
                //neighbours++;
//printf("6 harmonicenergy6 %lf \n",harmonicenergy6);
			}  
		 
			 /* For neighbour 3 */
			 if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension ) * (lespringLAA - extension ) ;
                //HE[6]= kspringLAA * (lespringLAA - extension ) * (lespringLAA - extension ) ;
                //neighbours++;
//printf("7 harmonicenergy7 %lf \n",harmonicenergy7);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension ) * (lespringLBB -extension ) ;
                //HE[7] = kspringLBB * (lespringLBB - extension ) * (lespringLBB -extension ) ;
                //neighbours++;
//printf("8  harmonicenergy8 %lf \n",harmonicenergy8);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -extension)*(lespringLAB -extension) ;
                //HE[8] = kspringLAB * (lespringLAB -extension)*(lespringLAB -extension) ;
                //neighbours++;
//printf("9 harmonicenergy9 %lf\n",harmonicenergy9);
			
			}  
			
			
			/* For neighbour 4 */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                //HE[9] = kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                //neighbours++;
//printf("10 harmonicenergy10 %lf \n",harmonicenergy10);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension )* (lespringLBB - extension) ;
                //HE[10] = kspringLBB * (lespringLBB - extension )* (lespringLBB - extension) ;
                //neighbours++;
//printf("11 harmonicenergy11 %lf \n",harmonicenergy11);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                //HE[11] = kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                //neighbours++;
//printf("12  harmonicenergy12 %lf \n",harmonicenergy12);
			
			}        
			 
			
		    /* For neighbour 5 */
		    if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                //HE[12] = kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                //neighbours++;
//printf("13 harmonicenergy13 %lf \n",harmonicenergy13);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
                //HE[13] = kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
                //neighbours++;
//printf("14 harmonicenergy14 %lf \n",harmonicenergy14);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                //HE[14] = kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                //neighbours++;
//printf("15 harmonicenergy15 %lf \n",harmonicenergy15);
			
			}  
			
			/* Neighbour 6 is absent since it is a surface atom */
			/* Completes nearest neighbours 
			 * 
			 * 
			 * Now we start adding the next nearest neighbours */
			 
			 
			/* Neighbour 7 */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[15] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("16 harmonicenergy16 %lf \n",harmonicenergy16);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[16] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("17 harmonicenergy17 %lf \n",harmonicenergy17);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[17] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("18 harmonicenergy18 %lf \n",harmonicenergy18);
			
			}  
			 
			/* Neighbour 8   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[18]= kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("19 harmonicenergy19 %lf \n",harmonicenergy19);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[19] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("20 harmonicenergy20 %lf \n",harmonicenergy20);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[20] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("21 harmonicenergy21 %lf \n",harmonicenergy21);
			
			} 
			
			/* Neighbour 9   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[21] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("22 harmonicenergy22 %lf \n",harmonicenergy22);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[22] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("23 harmonicenergy23 %lf \n",harmonicenergy23);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (- arefy+ displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[23] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("24 harmonicenergy24 %lf \n",harmonicenergy24);
			
} 
			
			/* Neighbour 10   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (-arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[24] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //eighbours++;
//printf("25 harmonicenergy25 %lf \n",harmonicenergy25);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[25] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("26 harmonicenergy26 %lf \n",harmonicenergy26);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[26] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("27 harmonicenergy27 %lf \n",harmonicenergy27);
			} 
			
			/* Neighbour 11  */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[27] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("28 harmonicenergy28 %lf \n",harmonicenergy28);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[28] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("29 harmonicenergy29 %lf \n",harmonicenergy29);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[29] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("30 harmonicenergy30 %lf \n",harmonicenergy30);
			
			}   
			
			/* Neighbour 12  */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[30] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("31 harmonicenergy31 %lf \n",harmonicenergy31);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[31] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("32 harmonicenergy32 %lf \n",harmonicenergy32);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[32]= kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("33 harmonicenergy33 %lf \n",harmonicenergy33);
			
			} 
			
			/* Neighbour 13   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[33] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("34 harmonicenergy34 %lf \n",harmonicenergy34);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (-arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[34] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("35 harmonicenergy35 %lf \n",harmonicenergy35);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[35] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("36 harmonicenergy36 %lf \n",harmonicenergy36);
			
			}   

            /* Neighbour 14   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[36]= kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("37 harmonicenergy37 %lf \n",harmonicenergy37);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[37] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("38 harmonicenergy38 %lf \n",harmonicenergy38);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[38] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("39 harmonicenergy39 %lf \n",harmonicenergy39);
			
			}   				
			
			/* Neighbour 15   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[39] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("40 harmonicenergy40 %lf \n",harmonicenergy40);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[40] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("41 harmonicenergy41 %lf \n",harmonicenergy41);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[41] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("42 harmonicenergy42 %lf \n",harmonicenergy42);
			
			}  
			
			/* Neighbour 16   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[42] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("43 harmonicenergy43 %lf \n",harmonicenergy43);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy+displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[43] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("44 harmonicenergy44 %lf extension %lf lespringDBB %lf \n",harmonicenergy44,extension,lespringDBB);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy+displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[44] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("45 harmonicenergy45 %lf \n",harmonicenergy45);
			
			} 
			
			/* Neighbour 17   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[45] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("46 harmonicenergy46 %lf \n",harmonicenergy46);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[46] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("47 harmonicenergy47 %lf \n",harmonicenergy47);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[47] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("48 harmonicenergy48 %lf \n",harmonicenergy48);
			
			} 
			
			/* Neighbour 18   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //HE[48] = kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                //neighbours++;
//printf("49 harmonicenergy49 %lf \n",harmonicenergy49);
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //HE[49] = kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                //neighbours++;
//printf("50 harmonicenergy50 %lf \n",harmonicenergy50);
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //HE[50] = kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                //neighbours++;
//printf("51 harmonicenergy51 %lf \n",harmonicenergy51);
			
			} 
/*
  error_Harmonic(&HE, neighbours, BONDXY);
  printf("\n\n");
  for(In=0;In<=50;In++){
   if(HE[In]>0.00)printf("%0.10lf \n",HE[In]);
  }
 */
	
return harmonicenergy;
			
}

#include "variables.h"

double Harmonic(int Atomx,int Atomy, int Atomz) 
{
    double harmonicenergy,extensionsq,extension,kmcfactor,arefzup, arefzdown;
           harmonicenergy = 0.0 ;
           spring_Bond=0;
     
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
                spring_Bond++;
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension) * (lespringLBB -extension) ;
                spring_Bond++;
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB - extension) * (lespringLAB -extension) ;
                spring_Bond++;
			}        
			
			/* For neighbour2 */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 1)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) ) ;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension ) * ( lespringLAA - extension)  ;
                spring_Bond++;
			}
			
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 2)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp) ) ;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension ) * (lespringLBB - extension ) ;
                spring_Bond++;
			}
			
			else if (((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 2))  || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz].Type == 1)) ) {
			
			    extensionsq = ((-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp)*(-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].xdisp))
			                + ((displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp)*(displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].ydisp))
			                + ((displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp)*(displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz].zdisp)) ;
			    
			    extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB - extension) * (lespringLAB - extension) ;
                spring_Bond++;
			}  
		 
			 /* For neighbour 3 */
			 if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension ) * (lespringLAA - extension ) ;
                spring_Bond++;
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension ) * (lespringLBB -extension ) ;
                spring_Bond++;
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -extension)*(lespringLAB -extension) ;
                spring_Bond++;
			}  
			
			
			/* For neighbour 4 */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                spring_Bond++;
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension )* (lespringLBB - extension) ;
                spring_Bond++;
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                spring_Bond++;
			}        
			 
			
		    /* For neighbour 5 */
		    if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
                spring_Bond++;
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
                spring_Bond++;
			
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
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			}  
			 
			/* Neighbour 8   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy-1+Ny)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			} 
			
			/* Neighbour 9   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) * (- arefy+ displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].ydisp) )
                            + ( (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) * (arefzdown + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
} 
			
			/* Neighbour 10   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 1)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (-arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 2)) {
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) * (displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].xdisp) )
                            + ( (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) * (- arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[Atomx][(Atomy+1)%Ny][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			} 
			
			/* Neighbour 11  */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			}   
			
			/* Neighbour 12  */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 1)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 2)) {
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].Type == 1)) ){
				
				extensionsq = ( (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) * (- arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].ydisp) )
                            + ( (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) * (arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz-1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			} 
			
			/* Neighbour 13   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (-arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			}   

            /* Neighbour 14   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].xdisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) * (displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].ydisp) )
                            + ( (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) * (- arefzup + displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][Atomy][Atomz+1].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			}   				
			
			/* Neighbour 15   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			}  
			
			/* Neighbour 16   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy+displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy+displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx-1+Nx) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			} 
			
			/* Neighbour 17   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].xdisp) )
                            + ( (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) * (arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy-1+Ny)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			} 
			
			/* Neighbour 18   */
			if ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
                spring_Bond++;
			
			}
			else if  ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) {
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
                spring_Bond++;
			
			} 
			else if ( ((displacementType[Atomx][Atomy][Atomz].Type == 1) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 2)) || ((displacementType[Atomx][Atomy][Atomz].Type == 2) && (displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].Type == 1)) ){
				
				extensionsq = ( (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) * (-arefx + displacementType[Atomx][Atomy][Atomz].xdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].xdisp) )
                            + ( (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) * (-arefy + displacementType[Atomx][Atomy][Atomz].ydisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].ydisp) )
                            + ( (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) * (displacementType[Atomx][Atomy][Atomz].zdisp - displacementType[(Atomx+1) % Nx][(Atomy+1)%Ny][Atomz].zdisp) ) 				;
                extension = pow(extensionsq, 0.5) ;
                harmonicenergy += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
                spring_Bond++;
			
			} 
	
return harmonicenergy;
			
}

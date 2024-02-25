/* Energycalculations       
 * This program will be inside the energy minimization suite  - This program might not be used because of the my_f that needs to be
 * used in the energy minimization*/
/* This program calculates the energy given a displacementType structure. This should be fairly simple based on the
 * equilibrium bond lengths which can be calculated and the given spring constants.  It needs to find an approximation
 * for the substrate that will depend on how it is being treated. For now, we are assuming that the forces on the substrate atoms 
 * are calculated in the energyminimization suite and entered to the program. This file takes in the displacementType structure 
 * and the fsurface structure  */

#include "variables.h"

extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

double energycalculation()
{
	int In,J,K,Iprev,Inext,Jprev,Jnext,Kprev,Knext ;
	double energydiscrete,energycontinuum,energy,extensionsq, extension , arefzup, arefzdown ;
	/*struct fsurface{double fsx; double fsy;} finterface[Nx][Ny];*/
	double xbelow[Nx][Ny],ybelow[Nx][Ny],zbelow[Nx][Ny],datax[Nx*Ny],datay[Nx*Ny],dataz[Nx*Ny];
	   
	energydiscrete = 0.0; 
	energycontinuum = 0.0 ;
		
		
		
		/* First we create list of neighbors taking periodic boundary conditions and the presence
				 * of the continuum substrate into account */
			
	
	for (In = 0; In < Nx; In++ ) {
		if (In ==0){
			Iprev = Nx-1 ;
			Inext = 1;
		}
		else if (In == (Nx-1)){
			Iprev = Nx - 2;
			Inext = 0 ;
		}
		else {
			Iprev = In - 1;
			Inext = In + 1;
		}
		
		for (J=0; J < Ny ; J++) {
			if (J ==0){
				Jprev = Ny-1 ;
				Jnext = 1;
			}
			else if (J == (Ny-1)){
				Jprev = Ny - 2;
				Jnext = 0;
			}
			else {
				Jprev = J - 1;
				Jnext = J + 1;
			}
			
			for (K=0; K < Nz ; K ++) {
				
				if (displacementType[In][J][K].Type == 0) break;
				
				if (K == 0){
					Kprev = Nz -1 ;   /* This ensures that for the bottommost discrete layer, we only look
					at the discrete substrate above. The part below is treated by energycontinuum. It is necessary to ensure
					that the top layer is empty */ 
					Knext = 1 ;
				}
				else 
				{
					Kprev = K -1 ;
					Knext = K + 1;
				}
				
				/* We need to define the reference configuration with respect to which all displacements are 
				 * calculated. There is a choice here, but I will use the simple reference configuration in which 
				 * the lateral lengths in the X and Y directions correspond to A lattice spacing and the Z direction
				 * given by a simple Smereka type relation, without involving square roots  */
				 
				if (K < DiscreteAlayers) {
					arefzup = lespringLAA ;
					arefzdown = lespringLAA;
				}
				else if (K == DiscreteAlayers){
					arefzup = arefzAB;
					arefzdown = lespringLAA;
				}
				else if (K == (DiscreteAlayers +1)) {
					arefzup = arefzB ;
					arefzdown = arefzAB ;
				}
				else{
	     			arefzup = arefzB;
					arefzdown = arefzB;
				}
				
				/* We have to list all the 18 neighbour interactions carefully. I will do it using if statements,
				 * which will make it long but make sure that everything is correct. To avoid double counting one must 
				 * only consider 9 neighbours in the energy calculation
				 * 
				 * Neighbours 1-6 represent lateral neighbours and neighbours 7-18 represent diagonal neighbours
				 * 
				 *   Neighbours 1-2  - Lateral X
				 *   Neighbours 3-4  - Lateral Y
				 *   Neighbours 5-6  - Lateral Z
				 * 
				 *   Neighbours 7-10   Diagonal YZ plane
				 *   Neighbours 11-14  Diagonal XZ plane
				 *   Neighbours 15-18  Diagonal XY plane
				 * 
				 * 
				 * To avoid double counting, I will use Neighbours 1,3,5,7,8,11,13,15,16
				 *    
				 */
				 
				/* For neighbour1 */ 
				if ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][K].Type == 1)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq,0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA -extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][K].Type == 2)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension) * (lespringLBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][K].Type == 1)) ){
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][K].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB - extension) * (lespringLAB -extension) ;
				
				}        
				
				/* For neighbour2 */
	/*			if ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][K].Type == 1)) {
					
					extensionsq = ( (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp) * (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp) ) ;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension ) * ( lespringLAA - extension)  ;
				
				}
				
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][K].Type == 2)) {
					
					extensionsq = ( (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp) * (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp) ) ;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension ) * (lespringLBB - extension ) ;
				
				}
				
				else if (((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][K].Type == 2))  || ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][K].Type == 1)) ) {
				
				    extensionsq = ((-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp)*(-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][K].xdisp))
				                + ((displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp)*(displacementType[In][J][K].ydisp - displacementType[Inext][J][K].ydisp))
				                + ((displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp)*(displacementType[In][J][K].zdisp - displacementType[Inext][J][K].zdisp)) ;
				    
				    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB - extension) * (lespringLAB - extension) ;
				
				}  
		*/
			 
				 /* For neighbour 3 */
				 if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][K].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension ) * (lespringLAA - extension ) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][K].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension ) * (lespringLBB -extension ) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -extension)*(lespringLAB -extension) ;
				
				}  
				
				
				/* For neighbour 4 */
		/*		if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][K].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][K].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension )* (lespringLBB - extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][K].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}        
			*/	 
				
			    /* For neighbour 5 */
			    if ((displacementType[In][J][K].Type == 1) && (displacementType[In][J][Kprev].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][J][Kprev].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][J][Kprev].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}  
				printf("Energy after Neighbor 5 %lf\n",energydiscrete);
				
				/* For neighbour 6 */
		/*	    if ((displacementType[In][J][K].Type == 1) && (displacementType[In][J][Knext].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][J][Knext].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][J][Knext].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][J][Knext].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[In][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}  
				 
			*/	 
				/* Completes nearest neighbours 
				 * 
				 * 
				 * Now we start adding the next nearest neighbours */
				 
				 
				/* Neighbour 7 */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][Kprev].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][Kprev].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][Kprev].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][Kprev].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}  
				 
				/* Neighbour 8   */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][Knext].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][Knext].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jprev][Knext].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jprev][Knext].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jprev][Knext].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
				
				/* Neighbour 9   */
	/*			if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][Kprev].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][Kprev].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][Kprev].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][Kprev].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) * (- arefy+ displacementType[In][J][K].ydisp - displacementType[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				} 
				
			*/	/* Neighbour 10   */
		/*		if ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][Knext].Type == 1)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) * (-arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][Knext].Type == 2)) {
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[In][Jnext][Knext].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[In][Jnext][Knext].Type == 1)) ){
					
					extensionsq = ( (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) * (displacementType[In][J][K].xdisp - displacementType[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) * (- arefy + displacementType[In][J][K].ydisp - displacementType[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				} 
			*/	
				/* Neighbour 11  */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][Kprev].Type == 1)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) * (arefx + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][Kprev].Type == 2)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][Kprev].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   
				
				/* Neighbour 12  */
	/*			if ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][Kprev].Type == 1)) {
					
					extensionsq = ( (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) * (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][Kprev].Type == 2)) {
					
					extensionsq = ( (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) * (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][Kprev].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) * (- arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Kprev].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) * (arefzdown + displacementType[In][J][K].zdisp - displacementType[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
		*/		
				/* Neighbour 13   */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][Knext].Type == 1)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][Knext].Type == 2)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) )
                                + ( (-arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][J][Knext].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][J][Knext].Type == 1)) ){
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Iprev][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   

                /* Neighbour 14   */
		/*		if ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][Knext].Type == 1)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][Knext].Type == 2)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][J][Knext].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][J][Knext].Type == 1)) ){
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][J][Knext].xdisp) )
                                + ( (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) * (displacementType[In][J][K].ydisp - displacementType[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) * (- arefzup + displacementType[In][J][K].zdisp - displacementType[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   				
		*/		
				/* Neighbour 15   */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][Jprev][K].Type == 1)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][Jprev][K].Type == 2)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][Jprev][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}  
				
				/* Neighbour 16   */
				if ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][Jnext][K].Type == 1)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][Jnext][K].Type == 2)) {
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy +displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Iprev][Jnext][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Iprev][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) * (arefx + displacementType[In][J][K].xdisp - displacementType[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy +displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Iprev][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
				printf("Energy after Neighbor 16 %lf   K %d\n",energydiscrete,K);
				
				/* Neighbour 17   */
			/*	if ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][Jprev][K].Type == 1)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][Jprev][K].Type == 2)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][Jprev][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jprev][K].xdisp) )
                                + ( (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) * (arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jprev][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
				
				 Neighbour 18   */
		/*		if ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][Jnext][K].Type == 1)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][Jnext][K].Type == 2)) {
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((displacementType[In][J][K].Type == 1) && (displacementType[Inext][Jnext][K].Type == 2)) || ((displacementType[In][J][K].Type == 2) && (displacementType[Inext][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) * (-arefx + displacementType[In][J][K].xdisp - displacementType[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) * (-arefy + displacementType[In][J][K].ydisp - displacementType[Inext][Jnext][K].ydisp) )
                                + ( (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) * (displacementType[In][J][K].zdisp - displacementType[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}
				
             End of calculations of neighbours */ 
			}
		}
	}

                /* -----------------------------*/
				/* End of discrete energy calculation  */
				
/* Here we need to calculate the interface forces. We will adapt the procedure of Smereka and Russo to our problem */

/* We will use the variables datax,datay,dataz for the fft routines */


    printf("energydiscrete %lf\n",energydiscrete);
    for (In=0; In< Nx; In ++) {
    	for (J=0; J< Ny; J ++) {
    	        datax[(In*Ny)+J] = displacementType[In][J][0].xdisp ;
		    	datay[(In*Ny)+J] = displacementType[In][J][0].ydisp ;
		    	dataz[(In*Ny)+J] = displacementType[In][J][0].zdisp ;
		   	
		 }
    }    				
    
    
	/* Now we start the fft part of the routine */
	
	fftw_complex *outx,*outy,*outz ;
	fftw_plan px,py,pz,pxb,pyb,pzb;
	double complex c1,c2,c3,rvec[Dim][Dim],det;
	int Dim1,Dim2 ;
	 
	outx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*((Ny/2)+1)) ) ;
	px = fftw_plan_dft_r2c_2d(Nx,Ny,datax,outx,FFTW_ESTIMATE);
	outy = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*((Ny/2)+1)) ) ;
	py = fftw_plan_dft_r2c_2d(Nx,Ny,datay,outy,FFTW_ESTIMATE);
	outz = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*((Ny/2)+1)) ) ;
	pz = fftw_plan_dft_r2c_2d(Nx,Ny,dataz,outz,FFTW_ESTIMATE);
	
	fftw_execute(px);
	fftw_execute(py);
	fftw_execute(pz);
	
	/*for (In =0 ; In< 3*Nx ; In++) printf("Outx %lf  %g  %gi\n",datax[In],creal(outx[In]),cimag(outx[In])); */
	fftw_destroy_plan(px);
	fftw_destroy_plan(py);
	fftw_destroy_plan(pz);

	for (In =0;In < Nx; In++)  {
	     	for (J=0; J<= Ny/2; J++) {
	     		 
	     		 
	     		 for (Dim1=0; Dim1<Dim ;Dim1++) {
	     		 	for (Dim2=0 ; Dim2 < Dim ;Dim2++) {
	     		 	/*	rvec[Dim1][Dim2][0] = GSL_REAL(gsl_vector_complex_get(evector[In][J][Dim1],Dim2));
	                	rvec[Dim1][Dim2][1] = GSL_IMAG(gsl_vector_complex_get(evector[In][J][Dim1],Dim2));   */
	                  /*  GSL_SET_COMPLEX(&rvec[Dim1][Dim2],GSL_REAL(gsl_vector_complex_get(evector[In][J][Dim1],Dim2)),GSL_IMAG(gsl_vector_complex_get(evector[In][J][Dim1],Dim2)));
                        printf("GSL Complex Test evec  %lf %lf  rvec %lf %lf\n",GSL_REAL(rvec[Dim1][Dim2]),GSL_IMAG(rvec[Dim1][Dim2]),GSL_REAL(gsl_vector_complex_get(evector[In][J][Dim1],Dim2)),GSL_IMAG(gsl_vector_complex_get(evector[In][J][Dim1],Dim2)));*/
	     		 	    
	     		 	    rvec[Dim1][Dim2] = GSL_REAL(gsl_vector_complex_get(evector[In][J][Dim1],Dim2)) + (I*GSL_IMAG(gsl_vector_complex_get(evector[In][J][Dim1],Dim2))) ;
	     		 	 /*   printf("rvalue %d %d %g %gi\n",Dim1,Dim2,creal(rvec[Dim1][Dim2]),cimag(rvec[Dim1][Dim2])); */
	     		 	}
	     		 }
	             
	             /* The determinant of the Matrix of eigenvectors is here */
	             
	             det = (rvec[0][0] * ((rvec[1][1]*rvec[2][2]) -(rvec[1][2]*rvec[2][1])))+ (rvec[0][1] * ((rvec[1][2]*rvec[2][0]) -(rvec[1][0]*rvec[2][2])))+ (rvec[0][2] * ((rvec[1][0]*rvec[2][1]) -(rvec[1][1]*rvec[2][0]))) ;
	             
	             /* Now we calculate the coefficients c1,c2,c3. Special situation is when kx=Nx/4,ky=Ny/4*/
	            
	            
	             c1= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))/det
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))/det
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2])))/det ;
	                    
	             c2= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))/det
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))/det
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0])))/det ;
	                       
	             c3= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))/det
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))/det
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1])))/det ;
	            
	        /*    
	            if(J <= Ny/2) {
	              
	                     c1= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2]))) ;
	                    
	                     c2= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0]))) ;
	                       
	                     c3= (outx[(In*((Ny/2)+1))+J] * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))
	                        +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))
	                        +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1]))) ;
	                      
	             }
	             else {
	             	     
	             	     c1= (conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))
	                        +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))
	                        +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2]))) ;
	                    
	                     c2= (conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))
	                        +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))
	                        +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0]))) ;
	                       
	                     c3= (conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))
	                        +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))
	                        +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1]))) ;
	             	     
	             }
	          */   
	             /* printf("%d   %d    %g   %gi    %g   %gi    %g   %gi\n ", In, J, creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
	                printf("%d   %d    %lf   %lf i    %lf   %lf i    %lf   %lf i\n ",In,J,GSL_REAL(root[In][J][0]),GSL_IMAG(root[In][J][0]),GSL_REAL(root[In][J][1]),GSL_IMAG(root[In][J][1]),GSL_REAL(root[In][J][2]),GSL_IMAG(root[In][J][2]));
	           */  /* Now we calculate the Fourier transform of the positions of the layer below */
	             
	             if ((In==0)&&(J==0))   {
                          outx[(In*((Ny/2)+1))+J]=0.0 + I * 0.0;
                          outy[(In*((Ny/2)+1))+J]=0.0 + I * 0.0;
                          outz[(In*((Ny/2)+1))+J]=0.0 + I * 0.0 ;
	             }
	             else if ((In == Nx/4)&&(J == Ny/4)) {
	                      
	                      outx[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[0][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[0][1]) ; 
	                      outy[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[1][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[1][1]) ; 
	                      outz[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[2][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[2][1]) ; 
	                      	     	             	   
	             }
	             else {
	                      
	                      outx[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[0][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[0][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[0][2]) ; 
	                      outy[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[1][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[1][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[1][2]) ; 
	                      outz[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[2][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[2][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[2][2]) ; 
	       	             	
	             }
	             /*printf("%g  %gi \n",creal(outx[(In*((Ny/2)+1))+J]),cimag(outx[(In*((Ny/2)+1))+J])); */
	             
	        } /* end of loop over J*/
	
	} /* End of loop over In */
	
	pxb = fftw_plan_dft_c2r_2d(Nx,(Ny/2)+1,outx,datax,FFTW_MEASURE);
	pyb = fftw_plan_dft_c2r_2d(Nx,(Ny/2)+1,outy,datay,FFTW_MEASURE);
	pzb = fftw_plan_dft_c2r_2d(Nx,(Ny/2)+1,outz,dataz,FFTW_MEASURE);
	
	fftw_execute(pxb);
	fftw_execute(pyb);
	fftw_execute(pzb);
	
	fftw_destroy_plan(pxb);
	fftw_destroy_plan(pyb);
	fftw_destroy_plan(pzb);
	
	fftw_free(outx);
	fftw_free(outy);
	fftw_free(outz);
	/*printf("made it here at end of backward fourier loop \n");*/
	
	
	energycontinuum = 0.0 ;

	for (In = 0; In < Nx; In ++) {
		for (J =0; J < Ny; J++) {
			
		    xbelow[In][J] = datax[(In*Ny)+J]/(Nx*Ny) ;
		    ybelow[In][J] = datay[(In*Ny)+J]/(Nx*Ny) ;
		    zbelow[In][J] = dataz[(In*Ny)+J]/(Nx*Ny) ;
		/*    printf("J  %d below %lf  %lf  %lf  energycontinuum %lf\n",J,xbelow[In][J],ybelow[In][J],zbelow[In][J],energycontinuum);
		*/	
		   	
		}
	}
	for (In = 0; In < Nx; In ++) {
		for (J =0; J < Ny; J++) {
			
		    
	        energycontinuum += (0.5 * kspringLAA * (displacementType[In][J][0].zdisp -zbelow[In][J])*(displacementType[In][J][0].zdisp -zbelow[In][J]) )
                        + (0.25 * kspringDAA * (-displacementType[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+ displacementType[In][J][0].zdisp - zbelow[(In+1)%Nx][J])*(-displacementType[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+ displacementType[In][J][0].zdisp - zbelow[(In+1)%Nx][J]) )
                        + (0.25 * kspringDAA * (-displacementType[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ displacementType[In][J][0].zdisp - zbelow[In][(J+1)%Ny])*(-displacementType[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ displacementType[In][J][0].zdisp - zbelow[In][(J+1)%Ny]) )
                        + (0.25 * kspringDAA * (displacementType[In][J][0].xdisp - xbelow[(In-1+Nx)% Nx][J]+ displacementType[In][J][0].zdisp - zbelow[(In-1+Nx)%Nx][J])*(-displacementType[In][J][0].xdisp - xbelow[(In+Nx-1)% Nx][J]+ displacementType[In][J][0].zdisp - zbelow[(In+Nx-1)%Nx][J]) )
                        + (0.25 * kspringDAA * (displacementType[In][J][0].ydisp - ybelow[In][(J+Ny-1)% Ny]+ displacementType[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny])*(-displacementType[In][J][0].ydisp + ybelow[In][(J+Ny-1)% Ny]+ displacementType[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny]) );   
			
			             
		/*	printf("J  %d below %lf  %lf  %lf  energycontinuum %lf\n",J,xbelow[In][J],ybelow[In][J],zbelow[In][J],energycontinuum);   */
			
		}
	}
	printf("Econtinuum  %lf  Ediscrete %lf\n",energycontinuum, energydiscrete);  
	energy = (0.5*energydiscrete) + energycontinuum ;				
	printf("Energy %lf\n",energy);
	
	return energy;
}/* End of energycalculation */

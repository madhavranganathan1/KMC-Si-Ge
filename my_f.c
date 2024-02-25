/* My_f.c is a copy of energycalculation. It basically sets takes a vector v of length 3 * Nx * Ny * Nz and converts
it into a set of displacementTypes and calculates the energy using exactly the procedure in  the energy calculation function */

#include "variables.h"

struct dT {double xdisp; 
           double ydisp; 
           double zdisp; 
           int Type:3;} displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 
	
double my_f(const gsl_vector *v, void *params)
{       
        double START1 = omp_get_wtime();
        int SIZE=omp_get_num_threads();	
	int In,J,K,Iprev,Inext,Jprev,Jnext,Kprev,Knext,i;
	double energydiscrete,energycontinuum,energy,extensionsq,extension,arefzup,arefzdown ;
	/*struct fsurface{double fsx; double fsy;} finterface[Nx][Ny];*/
	double xbelow[Nx][Ny],ybelow[Nx][Ny],zbelow[Nx][Ny],datax[Nx*Ny],datay[Nx*Ny],dataz[Nx*Ny];
	double *dp = (double *)params;
        double t=0;;
	   
//	energydiscrete = 0.0; 
	energycontinuum = 0.0 ;
	
        double x = gsl_vector_get(v, 0);
        double y = gsl_vector_get(v, 1);
//       printf("v%.10f\t%.10f\n",x,y);
	
	/* I need to assign disTy to the appropriate elements of v */
		
		
	for (In=0;In< Nx; In++) {
		for (J=0; J < Ny; J++) {
			for (K=0; K< Nz; K++) {
	              disTy[In][J][K].xdisp = gsl_vector_get(v,3*(K + (J*Nz) + (In*Ny*Nz))) ;
	              disTy[In][J][K].ydisp = gsl_vector_get(v,1+3*(K + J*Nz + In*Ny*Nz)) ;
	              disTy[In][J][K].zdisp = gsl_vector_get(v,2+3*(K + J*Nz + In*Ny*Nz)) ;
	              disTy[In][J][K].Type = dp[K + J*Nz + In*Ny*Nz] ;
	                                      }}}
		
		
		
		/* First we create list of neighbors taking periodic boundary conditions and the presence
				 * of the continuum substrate into account */
int ID,y_start,x_start,y_end,x_end;

#pragma omp parallel private(ID,i,In,J,K,y_start,x_start,y_end,x_end,Iprev,Inext,Jprev,Jnext,Kprev,Knext,arefzup,arefzdown,extensionsq,extension,energydiscrete) shared(SIZE,disTy) reduction(+:t)
{ 
//int A= omp_get_thread_num();
//    printf("%d\n",A);

energydiscrete=0;
#pragma parallel for nowait 
for (i =0; i <SIZE;i++){
ID= omp_get_thread_num();
/*    x_start = (ID)/2 * 32;
    y_start = (ID)%2 * 32;
    x_end = x_start+32;
    y_end = y_start+32;*/
//printf("Ran %d\n", ID);
    x_start =(ID)%8 * 16;
    y_start =(ID)/8 * 32;
    x_end = x_start+16;
    y_end = y_start+32;


	for (In = x_start; In <x_end; In++ ) {
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
		
		for (J=y_start; J <y_end ; J++) {
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
				
				if (disTy[In][J][K].Type == 0) break;
				
				if (K == 0) {
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
				 
				if (K < DiscreteAlayers-1) {
					arefzup = lespringLAA ;
					arefzdown = lespringLAA;
				}
				else if (K == DiscreteAlayers-1){
					arefzup = arefzAB;
					arefzdown = lespringLAA;
				}
				else if (K == (DiscreteAlayers )) {
					arefzup = arefzB ;
					arefzdown = arefzAB ;
				}
				else {
	     			arefzup = arefzB;
					arefzdown = arefzB;
				}
				
			/*	printf("arefx arefy arefzup arefzdown %lf %lf %lf %lf\n",arefx,arefy,arefzup,arefzdown);	*/
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
				/*printf("arefx %lf arefy %lf arefzup %lf arefzdown %lf\n",arefx,arefy,arefzup,arefzdown);*/
			/*	printf("kspringLAB %lf kspringLBB %lf kspringDAB %lf kspringDBB %lf\n",kspringLAB, kspringLBB,kspringDAB,kspringDBB); 
				printf("lespringLAB %lf lespringLBB %lf lespringDAB %lf lespringDBB %lf\n",lespringLAB, lespringLBB,lespringDAB,lespringDBB); */
				
				if ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][K].Type == 1)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq,0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA -extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][K].Type == 2)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension) * (lespringLBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][K].Type == 1)) ){
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][K].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][J][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB - extension) * (lespringLAB -extension) ;
				
				}        
			/*	printf("Energy after Neighbor 1 %lf\n",energydiscrete); */
				
				/* For neighbour2 */
	/*			if ((disTy[In][J][K].Type == 1) && (disTy[Inext][J][K].Type == 1)) {
					
					extensionsq = ( (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp) * (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp) ) ;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension ) * ( lespringLAA - extension)  ;
				
				}
				
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][K].Type == 2)) {
					
					extensionsq = ( (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp) * (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp) ) ;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension ) * (lespringLBB - extension ) ;
				
				}
				
				else if (((disTy[In][J][K].Type == 1) && (disTy[Inext][J][K].Type == 2))  || ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][K].Type == 1)) ) {
				
				    extensionsq = ((-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp)*(-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][K].xdisp))
				                + ((disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp)*(disTy[In][J][K].ydisp - disTy[Inext][J][K].ydisp))
				                + ((disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp)*(disTy[In][J][K].zdisp - disTy[Inext][J][K].zdisp)) ;
				    
				    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB - extension) * (lespringLAB - extension) ;
				
				}  
		*/
			 
				 /* For neighbour 3 */
				 if ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][K].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension ) * (lespringLAA - extension ) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][K].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension ) * (lespringLBB -extension ) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -extension)*(lespringLAB -extension) ;
				
				}  
			/*	printf("Energy after Neighbor 3 %lf\n",energydiscrete); */
				
				
				/* For neighbour 4 */
		/*		if ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][K].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][K].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension )* (lespringLBB - extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][K].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[In][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}        
			*/	 
				
			    /* For neighbour 5 */
			    if ((disTy[In][J][K].Type == 1) && (disTy[In][J][Kprev].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][J][Kprev].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][J][Kprev].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}  
		 /*		printf("Energy after Neighbor 5 %lf\n",energydiscrete); */
				
				/* For neighbour 6 */
		/*	    if ((disTy[In][J][K].Type == 1) && (disTy[In][J][Knext].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAA * (lespringLAA - extension) * (lespringLAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][J][Knext].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLBB * (lespringLBB - extension)*(lespringLBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][J][Knext].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][J][Knext].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[In][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringLAB * (lespringLAB -  extension ) * (lespringLAB - extension) ;
				
				}  
				 
			*/	 
				/* Completes nearest neighbours 
				 * 
				 * 
				 * Now we start adding the next nearest neighbours */
				 
				 
				/* Neighbour 7 */
				if ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][Kprev].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][Kprev].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][Kprev].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][Kprev].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Kprev].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jprev][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}  
			/*	printf("Energy after Neighbor 7 %lf\n",energydiscrete); */
				
				/* Neighbour 8   */
				if ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][Knext].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][Knext].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jprev][Knext].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jprev][Knext].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jprev][Knext].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[In][Jprev][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jprev][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
				/* printf("Energy after Neighbor 8 %lf\n",energydiscrete);*/
				
				/* Neighbour 9   */
	/*			if ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][Kprev].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][Kprev].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][Kprev].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][Kprev].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Kprev].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) * (- arefy+ disTy[In][J][K].ydisp - disTy[In][Jnext][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[In][Jnext][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				} 
				
				*//* Neighbour 10   */
		/*		if ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][Knext].Type == 1)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) * (-arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][Knext].Type == 2)) {
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[In][Jnext][Knext].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[In][Jnext][Knext].Type == 1)) ){
					
					extensionsq = ( (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) * (disTy[In][J][K].xdisp - disTy[In][Jnext][Knext].xdisp) )
                                + ( (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) * (- arefy + disTy[In][J][K].ydisp - disTy[In][Jnext][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[In][Jnext][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				} 
			*/	
				/* Neighbour 11  */
				if ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][Kprev].Type == 1)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][Kprev].Type == 2)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][Kprev].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Iprev][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   
			/*	printf("Energy after Neighbor 11 %lf\n",energydiscrete);*/
				
				/* Neighbour 12  */
	/*			if ((disTy[In][J][K].Type == 1) && (disTy[Inext][J][Kprev].Type == 1)) {
					
					extensionsq = ( (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) * (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][Kprev].Type == 2)) {
					
					extensionsq = ( (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) * (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Inext][J][Kprev].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][Kprev].Type == 1)) ){
					
					extensionsq = ( (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) * (- arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Kprev].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Kprev].ydisp) )
                                + ( (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) * (arefzdown + disTy[In][J][K].zdisp - disTy[Inext][J][Kprev].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
		*/		
				/* Neighbour 13   */
				if ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][Knext].Type == 1)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][Knext].Type == 2)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) )
                                + ( (-arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Iprev][J][Knext].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Iprev][J][Knext].Type == 1)) ){
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Iprev][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Iprev][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   
                /*printf("Energy after Neighbor 13 %lf\n",energydiscrete); */
				
                /* Neighbour 14   */
		/*		if ((disTy[In][J][K].Type == 1) && (disTy[Inext][J][Knext].Type == 1)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][Knext].Type == 2)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Inext][J][Knext].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Inext][J][Knext].Type == 1)) ){
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][J][Knext].xdisp) )
                                + ( (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) * (disTy[In][J][K].ydisp - disTy[Inext][J][Knext].ydisp) )
                                + ( (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) * (- arefzup + disTy[In][J][K].zdisp - disTy[Inext][J][Knext].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}   				
		*/		
				/* Neighbour 15   */
				if ((disTy[In][J][K].Type == 1) && (disTy[Iprev][Jprev][K].Type == 1)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Iprev][Jprev][K].Type == 2)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Iprev][Jprev][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Iprev][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}  
			/*	printf("Energy after Neighbor 15 %lf\n",energydiscrete); */
				
				/* Neighbour 16   */
				if ((disTy[In][J][K].Type == 1) && (disTy[Iprev][Jnext][K].Type == 1)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Iprev][Jnext][K].Type == 2)) {
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy +disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Iprev][Jnext][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Iprev][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) * (arefx + disTy[In][J][K].xdisp - disTy[Iprev][Jnext][K].xdisp) )
                                + ( (-arefy +disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Iprev][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Iprev][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB - extension) ;
				
				} 
				/*printf("Energy after Neighbor 16 %lf  arefx %lf  arefy %lf\n",energydiscrete,arefx,arefy); */
				
				/* Neighbour 17   */
			/*	if ((disTy[In][J][K].Type == 1) && (disTy[Inext][Jprev][K].Type == 1)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Inext][Jprev][K].Type == 2)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Inext][Jprev][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Inext][Jprev][K].Type == 1)) ){
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jprev][K].xdisp) )
                                + ( (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) * (arefy + disTy[In][J][K].ydisp - disTy[Inext][Jprev][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jprev][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				} 
				
				 Neighbour 18   */
		/*		if ((disTy[In][J][K].Type == 1) && (disTy[Inext][Jnext][K].Type == 1)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAA * (lespringDAA - extension) * (lespringDAA - extension) ;
				
				}
				else if  ((disTy[In][J][K].Type == 2) && (disTy[Inext][Jnext][K].Type == 2)) {
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDBB * (lespringDBB - extension)*(lespringDBB -extension) ;
				
				} 
				else if ( ((disTy[In][J][K].Type == 1) && (disTy[Inext][Jnext][K].Type == 2)) || ((disTy[In][J][K].Type == 2) && (disTy[Inext][Jnext][K].Type == 1)) ){
					
					extensionsq = ( (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) * (-arefx + disTy[In][J][K].xdisp - disTy[Inext][Jnext][K].xdisp) )
                                + ( (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) * (-arefy + disTy[In][J][K].ydisp - disTy[Inext][Jnext][K].ydisp) )
                                + ( (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) * (disTy[In][J][K].zdisp - disTy[Inext][Jnext][K].zdisp) ) 				;
                    extension = pow(extensionsq, 0.5) ;
                    energydiscrete += kspringDAB * (lespringDAB - extension)*(lespringDAB-extension) ;
				
				}
				
             End of calculations of neighbours */ 
		      // } //if loop
		}
	}
}
//printf("E=%lf_%d\n",energydiscrete,ID);
}//closed pragma for
t=energydiscrete;
//printf("t=%lf_%d\n",energydiscrete,ID);
}
                /* -----------------------------*/
				/* End of discrete energy calculation  */
				
/* Here we need to calculate the interface forces. We will adapt the procedure of Smereka and Russo to our problem */

/* We will use the variables datax,datay,dataz for the fft routines */

    for (In=0; In< Nx; In ++) {
    	for (J=0; J< Ny; J ++) {
    	                datax[(In*Ny)+J] = disTy[In][J][0].xdisp ;
		    	datay[(In*Ny)+J] = disTy[In][J][0].ydisp ;
		    	dataz[(In*Ny)+J] = disTy[In][J][0].zdisp ; 
		   	        }}    				
    
    
	/* Now we start the fft part of the routine */
	
	fftw_complex *outx,*outy,*outz ;
	fftw_plan px,py,pz,pxb,pyb,pzb;
	double complex c1,c2,c3,rvec[Dim][Dim],det;
	int Dim1,Dim2 ;
	 
	outx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*Ny) ) ;
	px = fftw_plan_dft_r2c_2d(Nx,Ny,datax,outx,FFTW_ESTIMATE);
	outy = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*Ny) ) ;
	py = fftw_plan_dft_r2c_2d(Nx,Ny,datay,outy,FFTW_ESTIMATE);
	outz = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *(Nx*Ny) ) ;
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
	            if ((fabs(creal(det)) < 0.0001) &&(fabs(cimag(det)) < 0.0001) )  {
	            	printf("Determinant zero for In  J  %d  %d  det   %lf + I %lf\n",In,J,creal(det),cimag(det));
	            	exit(0);
	            }
	         /*   if ((det*det) <= 0.000001) printf("In  J  K  %d  %d  %d  det   %lf\n",In,J,K,det); */
	             /* Now we calculate the coefficients c1,c2,c3. Special situation is when kx=Nx/4,ky=Ny/4*/
	            
	            
	            c1= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2]))))/det ;
	                    
	            c2= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0]))))/det ;
	                       
	            c3= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1]))))/det ;
	           
	            
	       /*     if(J <= Ny/2) {
	              
	                     c1= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2]))))/det ;
	                    
	                     c2= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0]))))/det ;
	                       
	                     c3= ((outx[(In*((Ny/2)+1))+J] * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))
	                         +(outy[(In*((Ny/2)+1))+J] * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))
	                         +(outz[(In*((Ny/2)+1))+J] * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1]))))/det ;
	                      
	             }
	             else {
	             	     
	             	     c1= ((conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][1]*rvec[2][2])-(rvec[2][1]*rvec[1][2])))
	                         +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][1]*rvec[0][2])-(rvec[0][1]*rvec[2][2])))
	                         +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][1]*rvec[1][2])-(rvec[1][1]*rvec[0][2]))))/det ;
	                    
	                     c2= ((conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][2]*rvec[2][0])-(rvec[2][2]*rvec[1][0])))
	                         +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][2]*rvec[0][0])-(rvec[0][2]*rvec[2][0])))
	                         +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][2]*rvec[1][0])-(rvec[1][2]*rvec[0][0]))))/det ;
	                       
	                     c3= ((conj(outx[(In*((Ny/2)+1))+Ny-J]) * ((rvec[1][0]*rvec[2][1])-(rvec[2][0]*rvec[1][1])))
	                         +(conj(outy[(In*((Ny/2)+1))+Ny-J]) * ((rvec[2][0]*rvec[0][1])-(rvec[0][0]*rvec[2][1])))
	                         +(conj(outz[(In*((Ny/2)+1))+Ny-J]) * ((rvec[0][0]*rvec[1][1])-(rvec[1][0]*rvec[0][1]))))/det ;
	             	     
	             } */
	             
	             /* printf("%d   %d    %g   %gi    %g   %gi    %g   %gi\n ", In, J, creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
	                printf("%d   %d    %lf   %lf i    %lf   %lf i    %lf   %lf i\n ",In,J,GSL_REAL(root[In][J][0]),GSL_IMAG(root[In][J][0]),GSL_REAL(root[In][J][1]),GSL_IMAG(root[In][J][1]),GSL_REAL(root[In][J][2]),GSL_IMAG(root[In][J][2]));
	         */    /* Now we calculate the Fourier transform of the positions of the layer below */
	             
	             if ((In==0)&&(J==0))   {
                          outx[(In*((Ny/2)+1))+J]=0.0 + I * 0.0;
                          outy[(In*((Ny/2)+1))+J]=0.0 + I * 0.0;
                          outz[(In*((Ny/2)+1))+J]=0.0 + I * 0.0 ;
	             }
	             /*else if ((In == Nx/4)&&(J == Ny/4)) */ 
	             else if (gsl_complex_abs(root[In][J][2]) < 0.0001) {
	                      
	                      outx[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[0][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[0][1]) ; 
	                      outy[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[1][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[1][1]) ; 
	                      outz[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[2][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[2][1]) ; 
	                      	     	             	   
	             }
	          /*   else if (J <= Ny/2) { */
	             
	             else   {
	                      
	                      outx[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[0][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[0][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[0][2]) ; 
	                      outy[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[1][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[1][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[1][2]) ; 
	                      outz[(In*((Ny/2)+1))+J] = (c1/(GSL_REAL(root[In][J][0]) + I*GSL_IMAG(root[In][J][0])) * rvec[2][0]) + (c2/(GSL_REAL(root[In][J][1]) + I*GSL_IMAG(root[In][J][1])) * rvec[2][1]) + (c3/(GSL_REAL(root[In][J][2]) + I*GSL_IMAG(root[In][J][2])) * rvec[2][2]) ; 
	       	             	
	             }
	             /*printf("%g  %gi \n",creal(outx[(In*((Ny/2)+1))+J]),cimag(outx[(In*((Ny/2)+1))+J])); */
	             
	        } /* end of loop over J*/
	
	} /* End of loop over In */
	
	pxb = fftw_plan_dft_c2r_2d(Nx,Ny,outx,datax,FFTW_MEASURE);
	pyb = fftw_plan_dft_c2r_2d(Nx,Ny,outy,datay,FFTW_MEASURE);
	pzb = fftw_plan_dft_c2r_2d(Nx,Ny,outz,dataz,FFTW_MEASURE);
	
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
/*
	for (In = 0; In < Nx; In ++) {
		for (J =0; J < Ny; J++) {
			
		    
	        energycontinuum += (0.5 * kspringLAA * (disTy[In][J][0].zdisp -zbelow[In][J])*(disTy[In][J][0].zdisp -zbelow[In][J]) )
                        + (0.25 * kspringDAA * (-disTy[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+1)%Nx][J])*(-disTy[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+1)%Nx][J]) )
                        + (0.25 * kspringDAA * (-disTy[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+1)%Ny])*(-disTy[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+1)%Ny]) )
                        + (0.25 * kspringDAA * (disTy[In][J][0].xdisp - xbelow[(In-1+Nx)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In-1+Nx)%Nx][J])*(-disTy[In][J][0].xdisp - xbelow[(In+Nx-1)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+Nx-1)%Nx][J]) )
                        + (0.25 * kspringDAA * (disTy[In][J][0].ydisp - ybelow[In][(J+Ny-1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny])*(-disTy[In][J][0].ydisp + ybelow[In][(J+Ny-1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny]) );   
			
		}
	}
*/
	/*printf("Econtinuum  %lf  Ediscrete %lf\n",energycontinuum, energydiscrete); */

for (In = 0; In < Nx; In ++) {
                for (J =0; J < Ny; J++) {
  energycontinuum += (0.5 * kspringLAA * (disTy[In][J][0].zdisp - zbelow[In][J])*(disTy[In][J][0].zdisp -zbelow[In][J]) )
         + (0.25 * kspringDAA * (-disTy[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+1)%Nx][J]) *(-disTy[In][J][0].xdisp + xbelow[(In+1)% Nx][J]+  disTy[In][J][0].zdisp - zbelow[(In+1)%Nx][J]) )
         + (0.25 * kspringDAA * (-disTy[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+1)%Ny]) * (-disTy[In][J][0].ydisp + ybelow[In][(J+1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+1)%Ny]) )
         + (0.25 * kspringDAA * (disTy[In][J][0].xdisp - xbelow[(In+Nx-1)% Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+Nx-1)%Nx][J]) * (disTy[In][J][0].xdisp - xbelow[(In+Nx-1)%Nx][J]+ disTy[In][J][0].zdisp - zbelow[(In+Nx-1)%Nx][J]) )
         + (0.25 * kspringDAA * (disTy[In][J][0].ydisp - ybelow[In][(J+Ny-1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny]) * (disTy[In][J][0].ydisp - ybelow[In][(J+Ny-1)% Ny]+ disTy[In][J][0].zdisp - zbelow[In][(J+Ny-1)%Ny]) );
                }
        }

//printf("Edis=%lf\n",t);
//printf("Econ=%lf\n",energycontinuum);
	energy = (0.5*t) + energycontinuum ;	/* Testing the discrete part for minimization */

//printf("E=%lf\n",energy);
	return energy;
}

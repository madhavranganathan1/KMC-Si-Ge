
#include "variables.h"

/* This program constructs Eigenvectors which will be used by the energycalculation routine to get the substrate energy
It is going to use a singular value decomposition technique to construct the eigenvectors in case eigenvalues are the same.
Actually due to the symmetry of the problems the fourier transforms of one half are related to the other half. This means that
the first half of eigenvalues and eigenvectors are related to the other half. However, we will calculate all of the eigenvalues
and eigenvectors here. This can help to check the code if necessary. */

void Eigenvectors()
{

	    gsl_complex evecelement[Dim] ;
     	gsl_complex elemomega11,elemomega12,elemomega13,elemomega22,elemomega23,elemomega33 ;
        gsl_complex det,det0,det1,det2;  /*determinant of the eigenvectors */
     	/* Here we define a bunch of quantities for the CLAPACK svd routine */

     	doublecomplex omegalapack[Dim*Dim],ulapack[Dim*Dim],vt[Dim*Dim],work[201]; /* Lapack definitions */
     	integer Mlapack,Nlapack,LDU,LDA,LDVT;

        Mlapack=Nlapack=LDU=LDA=LDVT=Dim;
        char JOBU, JOBVT;
        integer LWORK=201,INFO;
        doublereal slapack[Dim],rwork[201];


     	double rroot,iroot; /* real and imaginary parts of roots */
     	int In,J,K,L,Kcount  ;
        double ksi,eta,ck,sk,ce,se,a[(2*Dim) + 1],z[4*Dim],ba[(2*Dim) - 1],kd,kl ; /* ksi is the wavevector angle, s is sin(ksi), c is cos(ksi) and a[7] corresponds to the
        7 coefficients for the 6th order equation, z represents the roots of the polynomial equation with real and imaginary parts */

   /*     FILE *evecprint;  */



        JOBU = 'N';
     	JOBVT = 'A';

       /* evecprint = fopen("Eigenvectors.res","w");
	    fclose(evecprint);  */

        for (In= 0;In < Nx;In++) {

        	for (J=0;J < Ny;J++)      {

        		for (K=0; K < Dim ; K++)   {

		        		if ((evector[In][J][K] = gsl_vector_complex_alloc(Dim))==NULL)  {
		        			printf("vector allocation failed\n") ;
		        			exit(0) ;
		        		}

        		}

       	    }

        }

        kd = kspringDAA ;
        kl = kspringLAA ;


     	for (In=0;In < Nx;In++) {

        	for (J=0; J < Ny; J++) {

               if ((In==0)&&(J==0))  {

               	    GSL_SET_COMPLEX(&root[0][0][0],1.0,0.0);
	     	        GSL_SET_COMPLEX(&root[0][0][1],1.0,0.0);
	     	        GSL_SET_COMPLEX(&root[0][0][2],1.0,0.0);
	     	        gsl_vector_complex_set(evector[In][J][0],0,GSL_COMPLEX_ONE);
	     	        gsl_vector_complex_set(evector[In][J][0],1,GSL_COMPLEX_ZERO);
	     	        gsl_vector_complex_set(evector[In][J][0],2,GSL_COMPLEX_ZERO);
	     	        gsl_vector_complex_set(evector[In][J][1],0,GSL_COMPLEX_ZERO);
	     	        gsl_vector_complex_set(evector[In][J][1],1,GSL_COMPLEX_ONE);
               	    gsl_vector_complex_set(evector[In][J][1],2,GSL_COMPLEX_ZERO);
               	    gsl_vector_complex_set(evector[In][J][2],0,GSL_COMPLEX_ZERO);
	     	        gsl_vector_complex_set(evector[In][J][2],1,GSL_COMPLEX_ZERO);
	     	        gsl_vector_complex_set(evector[In][J][2],2,GSL_COMPLEX_ONE);


	    	/*        evecprint = fopen("Eigenvectors.res","a");

	     	        for (K=0;K< Dim; K++)  {
	     	              fprintf(evecprint,"In %d  J% d  K %d  Root %lf+%lf I\n",In,J,K,GSL_REAL(root[In][J][K]),GSL_IMAG(root[In][J][K]));
	                      gsl_vector_complex_fprintf(evecprint,evector[In][J][K],"%lf");

	                      fprintf(evecprint,"\n");

	                }
              */
               	    det0 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],0) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],1),gsl_vector_complex_get(evector[In][J][2],2)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],1),gsl_vector_complex_get(evector[In][J][1],2)) )  );
	                det1 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],1) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],2),gsl_vector_complex_get(evector[In][J][2],0)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],2),gsl_vector_complex_get(evector[In][J][1],0)) )  );
	                det2 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],2) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],0),gsl_vector_complex_get(evector[In][J][2],1)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],0),gsl_vector_complex_get(evector[In][J][1],1)) )  );
	                det = gsl_complex_add(gsl_complex_add(det0,det1),det2);
               	    if (gsl_complex_abs(det) < 0.0001) {
               	    	printf("det=0 for In = %d  J=%d\n",In,J);
               	        exit(0);
               	    }
               	  /*  fprintf(evecprint,"Determinant  %lf  +I %lf\n  \n",GSL_REAL(det), GSL_IMAG(det));
               	    fclose(evecprint);
	    */

               }
               else {

		              if(In<Nx/2){
                               ksi = 2.0 * (float)(In) * pi /(float) (Nx) ;}
                               else {ksi = 2.0 * (float)(In-Nx) * pi /(float) (Nx) ;}
                               if(J<Ny/2){
                               eta = 2.0 * (float)(J) * pi /(float) (Ny) ;}
                                else{ eta = 2.0 * (float)(J-Ny) * pi /(float) (Ny) ;}
		               ck = cos(ksi);
		               sk = sin(ksi);
		               ce = cos(eta);
		               se = sin(eta);

		               a[0] = kd * kd *( (kd *ck) + (ce *((kl * ck )+kd))) ;
		               a[6] = a[0];
		               a[1] = 2.0*kd * ( (-3.0*kd*kd) - (2.0*kd*kl) - (kl*kl*ck) + (kd *cos(2.0*eta)*((kl*ck) + (kd*cos(2.0*ksi))))
		                               + (ce*( (ck*(-(4.0*kd*kd)-(3.0*kd*kl)+(2.0*kl*kl)))-(kl*kl)+(kl*kd*cos(2.0*ksi)))) );
		               a[5] = a[1];
		               a[2] = (2.0*kl*((5.0*kd*kd)+(6.0*kd*kl)+(2.0*kl*kl))) + (kd*kd*kd*cos(3.0*eta))
		                    + (((22.0*kd*kd*kd) + (26.0*kd*kd*kl)+(4.0*kd*kl*kl)-(4.0*kl*kl*kl)+(kd*kd*kl*cos(3.0*eta)))*ck)
		                    - (2.0*kd*kl*kl*cos(2.0*ksi))
		                    - (2.0*kd*cos(2.0*eta)*((kl*kl)+(ck*((4.0*kd*kd)+(3.0*kd*kl)-(2.0*kl*kl)))-(kd*kl*cos(2.0*ksi))))
		                    + (kd*kd*kd*cos(3.0*ksi))
		                    + (ce*((22.0*kd*kd*kd) + (26.0*kd*kd*kl) + (4.0*kd*kl*kl) - (4.0*kl*kl*kl)+(kl*ck*(-(39.0*kd*kd)-(24.0*kl*kd)+(4.0*kl*kl))) + (kd*cos(2.0*ksi)*((-8.0*kd*kd)-(6.0*kd*kl)+(4.0*kl*kl))) + (kd*kd*kl*cos(3.0*ksi)) ) ) ;
		               a[4] = a[2];
		               a[3] = 2.0*( (-4.0*((2.0*kd)+kl)*((2.0*kd)+kl)*((2.0*kd)+kl)) +(((2.0*kl*((5.0*kd*kd) +(6.0*kd*kl)+ (2.0*kl*kl))) + (kd*kd*kd*cos(3.0*eta)))*ck)
		                          - (2.0*kd*cos(2.0*eta)*((3.0*kd*kd) + (2.0*kd*kl) + (kl*kl*ck)))
		                          - (2.0*kd*kd*cos(2.0*ksi)*((3.0*kd)+(2.0*kl)))
		                          + (ce* ((ck*((22.0*kd*kd*kd) + (26.0*kd*kd*kl)+(4.0*kd*kl*kl)-(4.0*kl*kl*kl)))+(2.0*kl*((5.0*kd*kd)+(6.0*kd*kl)+(2.0*kl*kl)-(kd*kl*cos(2.0*ksi))))+(kd*kd*kd*cos(3.0*ksi))) )
		                          ) ;
		         /*      printf("Coeffs: %lf   %lf   %lf    %lf\n",a[0],a[1],a[2],a[3]);*/

		               if (fabs(a[6]) < 0.00001) {
		            	    ba[0]=a[1];
		            	    ba[1]=a[2];
		            	    ba[2]=a[3];
		            	    ba[3]=a[2];
		            	    ba[4]=a[1];

		            	    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (5);
		                    gsl_poly_complex_solve (ba,5,w,z);
		                    gsl_poly_complex_workspace_free (w);
		                    z[8]=0.0;
		                    z[9]=0.0;
		                    z[10]=0.0;
		                    z[11]=0.0;
		               }
		               else {
		                    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (7);
		                    gsl_poly_complex_solve (a,7,w,z);
		                    gsl_poly_complex_workspace_free (w);
		               }

		               Kcount=0;
		       /*        printf("In %d\n",In);            */

		               for (L=0;L<(2*Dim);L++) {
		                /*  printf("z %lf z %lf  ",z[2*J],z[2*J+1]);  */

		     	    	  if(fabs(a[6]) < 0.00001) {
		     	    	  	    GSL_SET_COMPLEX(&root[In][J][2],0.0,0.0);
		     	    	  } /* Sets the third root to zero */

		     	    	  if ( (pow(z[2*L],2.0)+pow(z[(2*L)+1],2.0)) > 1.0 )  {
		     	    	    	Kcount++ ;
		     	    	 	    if (Kcount ==1) {
		     	    	 		    GSL_SET_COMPLEX(&root[In][J][0],z[2*L],z[(2*L)+1]) ;
		     	    /*	 		    printf("Root%d   %lf + I %lf\n",Kcount,GSL_REAL(root[In][J][0]),GSL_IMAG(root[In][J][0]));    */
		     	    	        }
		     	    	        else if (Kcount ==2) {
		     	    	            GSL_SET_COMPLEX(&root[In][J][1],z[2*L],z[(2*L)+1]) ;
		     	  /*  	 		    printf("Root%d   %lf + I %lf\n",Kcount,GSL_REAL(root[In][J][1]),GSL_IMAG(root[In][J][1]));      */

		     	    	        }
		     	    	        else {
		     	    	 		    if ( (Kcount != 3) && (fabs(a[6]) > 0.00001)) {
		     	    	 			    printf ("incorrect number of roots Kcount %d\n",Kcount);
		     	    	 			    exit(0);
		     	    	 		    }
		     	    	 		   /*if ((In==Nx/4)&&(J==Ny/4))  {
		     	    	 		    	GSL_SET_COMPLEX(&root[In][J][2],1000.0,0.0) ;  *//* This value has arbitrarily been chosen to a very large one
		     	    	 		    	                                                - This is because you have a lower order polynomial which in this
		     	    	 		    	                                                 language corresponds to two roots; 0 and infinity. The root greater
		     	    	 		    	                                                 than 1 is infinity, so we choose a large number. The exact value does not matter
		     	    	 		    	                                                 because we will explicitly put in the eigenvector (-1,1,0) and we will set the coefficient
		     	    	 		    	                                                 for the layer below to zero corresponding to a root of infinity  */
		     	    	 		   /* }
		     	    	 		    */
		     	    	 		    else {
			     	    	 		    GSL_SET_COMPLEX(&root[In][J][2],z[2*L],z[(2*L)+1]) ;
			     	/*    	 		    printf("Root%d   %lf + I %lf\n",Kcount,GSL_REAL(root[In][J][2]),GSL_IMAG(root[In][J][2]));  */
			     	    	 		}
		     	    	        }

		     	    	   }
		     	        }  /* end of loop over L - calculation of eigenvalues*/


		     /*	    	 		    printf("Root%d   %lf + I %lf\n",Kcount,GSL_REAL(root[In][J][2]),GSL_IMAG(root[In][J][2]));

		     	        } */

		           /* This completes Eigenvalue calculation. We have to only treat the cases of fabs(a[0]) < 0.00001 exceptionally */

		       /*         evecprint = fopen("Eigenvectors.res","a");   */

		     	        for (K=0; K< Dim;K++)  {

			     	    	  rroot = GSL_REAL(root[In][J][K]);
			                  iroot = GSL_IMAG(root[In][J][K]);

			                  /*   printf("root %d: %lf +I %lf\n",K,rroot,iroot); */

			                  GSL_SET_COMPLEX(&elemomega11,(2.0*kl*(ck-1)*rroot)+(kd*((ck*(pow(rroot,2.0)- pow(iroot,2.0)+1.0))+(2*rroot*((ck*ce)-2.0)))),(2*kl*(ck-1)*iroot)+ (kd*((ck*2.0*iroot*rroot)+(2*iroot*((ck*ce)-2.0)))) );
			     	    	  GSL_SET_COMPLEX(&elemomega22,(2.0*kl*(ce-1)*rroot)+(kd*((ce*(pow(rroot,2.0)- pow(iroot,2.0)+1.0))+(2*rroot*((ck*ce)-2.0)))),(2*kl*(ce-1)*iroot)+ (kd*((ce*2.0*iroot*rroot)+(2*iroot*((ck*ce)-2.0)))) );
			     	    	  GSL_SET_COMPLEX(&elemomega12,-2.0*rroot*kd*se*sk,-2.0*iroot*kd*se*sk);
			     	    	  GSL_SET_COMPLEX(&elemomega13,-2.0*kd*sk*rroot*iroot,kd*sk*( pow(rroot,2.0)- pow(iroot,2.0)-1.0 ));
			     	    	  GSL_SET_COMPLEX(&elemomega23,-2.0*kd*se*rroot*iroot,kd*se*( pow(rroot,2.0)- pow(iroot,2.0)-1.0 ));
			     	    	  GSL_SET_COMPLEX(&elemomega33,(kl*(pow(rroot,2.0)-pow(iroot,2.0)-(2.0*rroot)+1.0)) + (kd*( ((ce+ck)*(pow(rroot,2.0)-pow(iroot,2.0)+1)) -(4.0*rroot))) , (kl*(2.0*iroot*(rroot-1.0))) + (kd*((2.0*iroot*rroot*(ck+ce))-(4.0*iroot) ) ) );

		     	    	      /* printf("ce   %lf  Omega22real %lf\n",ce,GSL_REAL(elemomega22));*/
			     	    	  /* Now we do the SVD in CLAPACK */
			     	    	  /* First we need to define the matrix omegalapack */


			     	          if ((K>0) &&( fabs(rroot -GSL_REAL(root[In][J][K-1])) < 0.001) && (fabs(iroot -GSL_IMAG(root[In][J][K-1])) < 0.001))  {

		              	          if ( (K==2) && (fabs(rroot -GSL_REAL(root[In][J][0])) < 0.001) && (fabs(iroot -GSL_IMAG(root[In][J][0])) < 0.001) ) {   /* All three roots equal */

		              	          	  GSL_SET_COMPLEX(&evecelement[0],vt[0].r,vt[0].i);
					     	    	  GSL_SET_COMPLEX(&evecelement[1],vt[3].r,vt[3].i);
					     	    	  GSL_SET_COMPLEX(&evecelement[2],vt[6].r,vt[6].i);

		              	          }
		              	          else {

			              	          GSL_SET_COMPLEX(&evecelement[0],vt[1].r,vt[1].i);
					     	    	  GSL_SET_COMPLEX(&evecelement[1],vt[4].r,vt[4].i);
					     	    	  GSL_SET_COMPLEX(&evecelement[2],vt[7].r,vt[7].i);


					           /*     printf ("VT   %d   %lf +I %lf   %lf +I %lf   %lf +I %lf\n         %lf +I %lf  %lf +I %lf   %lf +I %lf\n         %lf +I %lf  %lf +I %lf   %lf +I %lf\n",K, vt[0].r,vt[0].i,vt[1].r,vt[1].i,vt[2].r,vt[2].i,vt[3].r,vt[3].i,vt[4].r,vt[4].i,vt[5].r,vt[5].i,vt[6].r,vt[6].i,vt[7].r,vt[7].i,vt[8].r,vt[8].i);
					     	    	  printf ("U   %d   %lf +I %lf   %lf +I %lf\n          %lf +I %lf   %lf +I %lf\n",K, ulapack[0].r,ulapack[0].i,ulapack[1].r,ulapack[1].i,ulapack[2].r,ulapack[2].i,ulapack[3].r,ulapack[3].i);
					     	      	  printf ("S   In=%d  J=%d  K=%d   %lf    %lf     %lf\n",In,J,K,slapack[0],slapack[1],slapack[2]);
			                    */
		              	          }/* End of Clause for all three roots equal */
			                  	  gsl_vector_complex_set(evector[In][J][K],0,evecelement[0]);
					     	      gsl_vector_complex_set(evector[In][J][K],1,evecelement[1]);
					     	      gsl_vector_complex_set(evector[In][J][K],2,evecelement[2]);

			     	          }	 /* clause for equal roots */

			     	          else {

						     	   /*  if ((In==Nx/4) && (J==Ny/4) &&(K==2)) {

				     	    	  	      gsl_vector_complex_set(evector[In][J][K],0,gsl_complex_sub(GSL_COMPLEX_ZERO,GSL_COMPLEX_ONE));
				     	    	          gsl_vector_complex_set(evector[In][J][K],1,GSL_COMPLEX_ONE);
				     	    	          gsl_vector_complex_set(evector[In][J][K],2,GSL_COMPLEX_ZERO);

				     	    	     }
						     	     else {
						     	     */
						     	     	  omegalapack[0].r = GSL_REAL(elemomega11);
						     	          omegalapack[0].i = GSL_IMAG(elemomega11);
						     	          omegalapack[1].r = GSL_REAL(elemomega12);
						     	          omegalapack[1].i = GSL_IMAG(elemomega12);
						     	    	  omegalapack[2].r = GSL_REAL(elemomega13);
						     	          omegalapack[2].i = GSL_IMAG(elemomega13);
						     	    	  omegalapack[3].r = GSL_REAL(elemomega12);
						     	          omegalapack[3].i = GSL_IMAG(elemomega12);
						     	          omegalapack[4].r = GSL_REAL(elemomega22);
						     	          omegalapack[4].i = GSL_IMAG(elemomega22);
						     	          omegalapack[5].r = GSL_REAL(elemomega23);
						     	          omegalapack[5].i = GSL_IMAG(elemomega23);
						     	    	  omegalapack[6].r = GSL_REAL(elemomega13);
						     	          omegalapack[6].i = GSL_IMAG(elemomega13);
						     	    	  omegalapack[7].r = GSL_REAL(elemomega23);
						     	          omegalapack[7].i = GSL_IMAG(elemomega23);
						     	          omegalapack[8].r = GSL_REAL(elemomega33);
						     	          omegalapack[8].i = GSL_IMAG(elemomega33);


						     	    	/* This completes the definition of omegalapack */
						     	    	/* Now we can start  the SVD */
						     	   /*      printf ("Omega   %d   %lf +I %lf   %lf +I %lf  %lf +I %lf \n  %lf +I %lf  %lf +I %lf   %lf +I %lf \n %lf +I %lf   %lf +I %lf  %lf +I %lf\n",K, omegalapack[0].r,omegalapack[0].i,omegalapack[1].r,omegalapack[1].i,omegalapack[2].r,omegalapack[2].i,omegalapack[3].r,omegalapack[3].i,omegalapack[4].r,omegalapack[4].i,omegalapack[5].r,omegalapack[5].i,omegalapack[6].r,omegalapack[6].i,omegalapack[7].r,omegalapack[7].i,omegalapack[8].r,omegalapack[8].i);
						     	   */


						     	    	  zgesvd_(&JOBU,&JOBVT,&Mlapack, &Nlapack, omegalapack, &LDA, slapack, ulapack, &LDU, vt, &LDVT, work, &LWORK, rwork,&INFO);
						     	          GSL_SET_COMPLEX(&evecelement[0],vt[2].r,vt[2].i);
					     	    	      GSL_SET_COMPLEX(&evecelement[1],vt[5].r,vt[5].i);
						     	    	  GSL_SET_COMPLEX(&evecelement[2],vt[8].r,vt[8].i);

						     	    	  gsl_vector_complex_set(evector[In][J][K],0,evecelement[0]);
						     	    	  gsl_vector_complex_set(evector[In][J][K],1,evecelement[1]);
						     	    	  gsl_vector_complex_set(evector[In][J][K],2,evecelement[2]);


						     	    /* } */ /* End of else clause for In==Nx/4, J==Ny/4 */




			     	          } /* End of else clause for identical eigenvalues */

			     	      /*     printf ("VT   %d   %lf +I %lf   %lf +I %lf   %lf +I %lf\n         %lf +I %lf  %lf +I %lf   %lf +I %lf\n         %lf +I %lf  %lf +I %lf   %lf +I %lf\n",K, vt[0].r,vt[0].i,vt[1].r,vt[1].i,vt[2].r,vt[2].i,vt[3].r,vt[3].i,vt[4].r,vt[4].i,vt[5].r,vt[5].i,vt[6].r,vt[6].i,vt[7].r,vt[7].i,vt[8].r,vt[8].i);
					     	   printf ("U   %d   %lf +I %lf   %lf +I %lf\n          %lf +I %lf   %lf +I %lf\n",K, ulapack[0].r,ulapack[0].i,ulapack[1].r,ulapack[1].i,ulapack[2].r,ulapack[2].i,ulapack[3].r,ulapack[3].i);
					     	   printf ("S   I=%d  J=%d  K=%d   %lf    %lf     %lf\n",In,J,K,slapack[0],slapack[1],slapack[2]);
					     	*/

		         /*             fprintf(evecprint,"In %d  J% d  K %d  Root %lf+%lf I\n",In,J,K,GSL_REAL(root[In][J][K]),GSL_IMAG(root[In][J][K]));
			                  gsl_vector_complex_fprintf(evecprint,evector[In][J][K],"%lf");
		                      fprintf(evecprint,"\n");
			     	*/

			     	    }/* End of loop over K - for each root.  */

		     	        det0 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],0) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],1),gsl_vector_complex_get(evector[In][J][2],2)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],1),gsl_vector_complex_get(evector[In][J][1],2)) )  );
		                det1 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],1) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],2),gsl_vector_complex_get(evector[In][J][2],0)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],2),gsl_vector_complex_get(evector[In][J][1],0)) )  );
		                det2 = gsl_complex_mul( gsl_vector_complex_get(evector[In][J][0],2) , gsl_complex_sub (gsl_complex_mul(gsl_vector_complex_get(evector[In][J][1],0),gsl_vector_complex_get(evector[In][J][2],1)) , gsl_complex_mul(gsl_vector_complex_get(evector[In][J][2],0),gsl_vector_complex_get(evector[In][J][1],1)) )  );
		                det = gsl_complex_add(gsl_complex_add(det0,det1),det2);
	               	    if (gsl_complex_abs(det) < 0.0001) {
               	        	printf("det=0 for In = %d  J=%d\n",In,J);
               	        	printf("Corresponding Eigenvectors are  given below\n");
               	        	gsl_vector_complex_fprintf(stdout,evector[In][J][0],"%lf");
               	            gsl_vector_complex_fprintf(stdout,evector[In][J][1],"%lf");
               	            gsl_vector_complex_fprintf(stdout,evector[In][J][2],"%lf");
               	            printf("Corresponding Eigenvalues are  %lf %lf %lf %lf %lf %lf\n",GSL_REAL(root[In][J][0]),GSL_IMAG(root[In][J][0]),GSL_REAL(root[In][J][1]),GSL_IMAG(root[In][J][1]),GSL_REAL(root[In][J][2]),GSL_IMAG(root[In][J][2]));
               	            printf("Coeffs: %lf   %lf   %lf    %lf\n",a[0],a[1],a[2],a[3]);
               	        	exit(0);
               	        }


	               	    /*fprintf(evecprint,"Determinant  %lf  +I %lf\n  \n",GSL_REAL(det), GSL_IMAG(det));
	               	    fclose(evecprint);
		     	       */
	     	   }  /* End of else clause of If In=J=0.  */

     		}/* end of loop over J */

     	} /* This represents the loop over In */



}/* End of eigenvectors */

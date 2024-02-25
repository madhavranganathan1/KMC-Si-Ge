int energyminimization()
{	 
   const size_t Natoms = Nx*Ny*Nz ; /* This has to be the product of the top two numbers */
   const size_t N = Dim*Nx*Ny*Nz ; /* this is twice the above number with one atom fixed*/
   size_t iter = 0;
   int status,In,J,K; 
       
   const gsl_multimin_fdfminimizer_type *T;
         gsl_multimin_fdfminimizer *s;   
     
       
        double par[Natoms] ;
        gsl_vector *x ;
        x = gsl_vector_alloc (N);
        gsl_multimin_function_fdf my_func;   
        
        if(SpringType==1) /* 1 = Harmonic Springs with Square roots*/
        {        
             my_func.f = &my_f;
             my_func.df = &my_df;
             my_func.fdf = &my_fdf;   
        }
       // if(SpringType==2) /* 2 = Linear Springs with No Square roots*/
      /*  {        
             my_func.f = &my_fthomas;
             my_func.df = &my_dfthomas;
             my_func.fdf = &my_fdfthomas;   
        }*/
        my_func.n = N;
        my_func.params = &par;

        
        for (K = 0; K < Natoms; K++) {
        	gsl_vector_set(x,3*K,displacementType[K/(Ny*Nz)][(K%(Ny*Nz))/Nz][(K%(Ny*Nz))%Nz].xdisp);
        	gsl_vector_set(x,(3*K)+1,displacementType[K/(Ny*Nz)][(K%(Ny*Nz))/Nz][(K%(Ny*Nz))%Nz].ydisp);
        	gsl_vector_set(x,(3*K)+2,displacementType[K/(Ny*Nz)][(K%(Ny*Nz))/Nz][(K%(Ny*Nz))%Nz].zdisp);
        	par[K] = displacementType[K/(Ny*Nz)][(K%(Ny*Nz))/Nz][(K%(Ny*Nz))%Nz].Type;
                                     }
                //printf (" %0.10f %0.10lf\n\n", gsl_vector_get (x, 0),displacementType[0][0][0].xdisp);

      //  if(optimization_type==0){T = gsl_multimin_fdfminimizer_conjugate_fr;}
      //  if(optimization_type==1){T = gsl_multimin_fdfminimizer_conjugate_pr;}
      //  if(optimization_type==2){T = gsl_multimin_fdfminimizer_vector_bfgs2;}
      //  if(optimization_type==3){T = gsl_multimin_fdfminimizer_vector_bfgs;}
        if(optimization_type==4){T = gsl_multimin_fdfminimizer_steepest_descent;}

        s = gsl_multimin_fdfminimizer_alloc (T, N);
        gsl_multimin_fdfminimizer_set (s, &my_func, x,  ssizeastep, tolaglobal);
       
 do
     {
       iter++;
       status = gsl_multimin_fdfminimizer_iterate (s);
       if (status)
         break;
           status = gsl_multimin_test_gradient (s->gradient, gtolagrad);
                    //printf ("%5d %.10f %10.5f\n", iter, gsl_vector_get (s->x, 0),s->f);
       if (status == GSL_SUCCESS) {energyvalue = s->f;
                                    if(energyvalue<0.000001){ printf("\n There is some problem in energyminimization \n");
                                      exit(0);}

                                  }

       //printf ("Minimum found at:\n"); 
       //printf ("%5d %.10f %10.5f\n", iter, gsl_vector_get (s->x, 0),s->f);
     }
  while (status == GSL_CONTINUE && iter < 1000);
//printf ("Minimum not found at:\n");

        for (In = 0; In < Nx; In++) {
                for (J=0; J < Ny; J++) {
                        for (K = 0; K < Nz; K++) {
                                displacementType[In][J][K].xdisp = gsl_vector_get(s->x,3*(K + (J*Nz) + (In*Ny*Nz)));
                                displacementType[In][J][K].ydisp = gsl_vector_get(s->x,1+3*(K + (J*Nz) + (In*Ny*Nz)));
                                displacementType[In][J][K].zdisp = gsl_vector_get(s->x,2+3*(K + (J*Nz) + (In*Ny*Nz)));
                        }
                }
        }
  //printf ("%5d %.10f %0.10f %10.5f\n", iter, gsl_vector_get (s->x, 0),displacementType[0][0][0].xdisp, s->f);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);


return status;
	
}

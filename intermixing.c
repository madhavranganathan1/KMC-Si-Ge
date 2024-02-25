#include"chemical_bond.c"
#include"swap.c"
#include"valuemin.c"
#include"pswap_optimization.c"
#include"swap_update.c"
#include"swap_optimization.c"
void intermixing()

{

int i,j,k,xp,yp,xm,ym,atomx,atomy,atomz,Atomx,Atomy,Atomz,N_int=0,TYPe,x,y,z,n=0,c=0;
double a,b=0,xl=0,d=0,e=0,f=0,g=0,h,l,m,o,p,q,ru,energydiff,tos,exchange_factor;
for(k=0;k<24;k++)
{

  for(i=0;i<Nx;i++)
  {
   for(j=0;j<Ny;j++)
    {
      if(displacementType[i][j][k].Type==1)
        {
          	 xp=(i+1)%Nx, xm=(i-1+Nx)%Nx, yp=(j+1)%Nx, ym=(j-1+Nx)%Nx;
          if((displacementType[xp][j][k].Type==2)||(displacementType[i][yp][k].Type==2)||(displacementType[xm][j][k].Type==2)||(displacementType[i][ym][k].Type==2)||(displacementType[i][j][k-1].Type==2)||(displacementType[i][j][k+1].Type==2))
	  {
		  N_int=N_int+1;

	  }
	  
             }}}}	 printf("interface_atoms=%d\n",N_int);

           for(i=0;i<3;i++){
                   I_co[i]=(int *)calloc(N_int,sizeof(int));
                   if((I_co[i]=(int *)calloc(N_int,sizeof(int)))==NULL){printf("no space avaliable\n");
                                                                             exit;}
                  }
           
           for(k=0;k<24;k++)
   {       
        
  for(i=0;i<Nx;i++)
  {       
   for(j=0;j<Ny;j++)
    {      
      if(displacementType[i][j][k].Type==1)
        {     
                 xp=(i+1)%Nx, xm=(i-1+Nx)%Nx, yp=(j+1)%Nx, ym=(j-1+Nx)%Nx;
          if((displacementType[xp][j][k].Type==2)||(displacementType[i][yp][k].Type==2)||(displacementType[xm][j][k].Type==2)||(displacementType[i][ym][k].Type==2)||(displacementType[i][j][k-1].Type==2)||(displacementType[i][j][k+1].Type==2))
            {
		    *(I_co[0])=i;
		      I_co[0]++;
		      *(I_co[1])=j;
		        I_co[1]++;
			*(I_co[2])=k;
			  I_co[2]++;
			  
	    }		  
       }}}}
          
           
        for(i=0;i<3;i++){I_co[i]=I_co[i]-N_int;}
/*        for(i=0;i<N_int;i++){
            printf("i=%d\tj=%d\tk=%d\n",*(I_co[0]+i),*(I_co[1]+i),*(I_co[2]+i));}*/

 
        int r=  (int) (N_int*gsl_rng_uniform_pos(gnalea_r31) ) ;
	printf("r=%d\n",r);
	atomx=*(I_co[0]+r);
	atomy=*(I_co[1]+r);
	atomz=*(I_co[2]+r);
	int type=displacementType[atomx][atomy][atomz].Type;
	if(type==1){a=chemical_bond(type,atomx,atomy,atomz);// start of first if statement
		    printf("atomx=%d\tatomy=%d\tatomz=%d\n",atomx,atomy,atomz);
		    int typ=displacementType[(atomx+1)%Nx][atomy][atomz].Type;
	            if(typ==2){b=chemical_bond(typ,((atomx+1)%Nx),atomy,atomz);}
                    int Typ=displacementType[atomx][(atomy+1)%Nx][atomz].Type;
                    if(Typ==2){xl=chemical_bond(Typ,atomx,(atomy+1)%Nx,atomz);}
	            int TYp=displacementType[(atomx-1+Nx)%Nx][atomy][atomz].Type;		    
                    if(TYp==2){d=chemical_bond(TYp,((atomx-1+Nx)%Nx),atomy,atomz);}
		    int tYp=displacementType[atomx][(atomy-1+Ny)%Ny][atomz].Type;
		    if(tYp==2){e=chemical_bond(tYp,atomx,((atomy-1+Ny)%Ny),atomz);}
		    int tyP=displacementType[atomx][atomy][atomz-1].Type;	    
                    if(tyP==2){f=chemical_bond(tyP,atomx,atomy,(atomz-1));}
		    int tYpe=displacementType[atomx][atomy][atomz+1].Type;
		    if(tYpe==2){g=chemical_bond(tYpe,atomx,atomy,(atomz+1));}
	             h= valuemin(b,xl,d,e,f,g);
                     if(h==b){Atomx=(atomx+1)%Nx,Atomy=atomy,Atomz=atomz,TYPe=typ;}
		     if(h==xl){Atomx=atomx,Atomy=(atomy+1)%Ny,Atomz=atomz,TYPe=Typ;}
		     if(h==d){Atomx=(atomx-1+Nx)%Nx,Atomy=atomy,Atomz=atomz,TYPe=TYp;}
		     if(h==e){Atomx=atomx,Atomy=(atomy-1+Ny)%Ny,Atomz=atomz,TYPe=tYp;}
		     if(h==f){Atomx=atomx,Atomy=atomy,Atomz=atomz-1,TYPe=tyP;}
		     if(h==g){Atomx=atomx,Atomy=atomy,Atomz=atomz+1,TYPe=tYpe;}
		     if(displacementType[Atomx][Atomy][Atomz].Type==2){//START OF SECOND IF STATEMENT	
		     q=energyvalue;
	             m=(a+h);
                   //  l=chemical_bond(type,Atomx,Atomy,Atomz);
		     p=swap(atomx,atomy,atomz,Atomx,Atomy,Atomz);
		     pswap_optimization(Atomx,Atomy,Atomz,atomx,atomy,atomz);
		     ru=energyvalue;
		     o= p;
		    
		     energydiff=m-o+ru-q+0.4;
                       printf("q=%lf\tru=%lf\nm=%lf\to=%lf\tenergydiff=%lf\n",q,ru,m,o,energydiff);

		     if(h==b||h==xl||h==d||h==e){printf("swap is on the plane\n\n");}
		     if(energydiff<0){displacementType[atomx][atomy][atomz].Type=2;
		                      displacementType[Atomx][Atomy][Atomz].Type=1;
				      swap_update(N_int,Atomx,Atomy,Atomz,r);
				      swap_optimization(); 
		                        n++;  }//move accepted
		      else if(energydiff>0){
                                              tos=gsl_rng_uniform_pos(gnalea_r32);
					      printf("tos=%lf\n",tos);
					       exchange_factor=exp(-energydiff/(kB * temperature));

			                      if(tos<=exchange_factor){displacementType[atomx][atomy][atomz].Type=2;
                                                                  displacementType[Atomx][Atomy][Atomz].Type=1;
								  swap_update(N_int,Atomx,Atomy,Atomz,r);
								  swap_optimization();
				                                           n++;	      }//move accepted
		                              else{displacementType[atomx][atomy][atomz].Type=1;
                                                   displacementType[Atomx][Atomy][Atomz].Type=2;
					         c++;       }// move rejected  
		      }
		      printf("Atomx=%d\tAtomy=%d\tAtomz=%d\n\n",Atomx,Atomy,Atomz);}//END OF SECOND IF STATEMENT
		     
			

	
	


      }	//end of first if 		    


     printf("accepted_moves=%d\trejected_moves=%d\n",n,c);
		     

/*for(i=0;i<3;i++)
 {
   free(I_co[i]);}*/




                    		    
                    
 }
         	

void swap_update(int n,int Atomx,int Atomy,int Atomz,int l)
{     
	int i,j,k,N_intnew=0,xm,ym,xp,yp;	
	printf("l=%d\n",l);
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
                  N_intnew=N_intnew+1;
          }
     }}}}
     printf("interfacial_Atom=%d\n",N_intnew);


	if(N_intnew==n){  // start of first if statement
                      *(I_co[0]+l)=Atomx;
                      *(I_co[1]+l)=Atomy;
                      *(I_co[2]+l)=Atomz;} //end of first if statement
         if(N_intnew<n|| N_intnew>n)//start of second if statement
                        {/*int iy=1; for(i=0;i<3;i++){  
                            I_co[i]=(int*)realloc(I_co[i],N_intnew*sizeof(int));
                            if(I_co[i]>0 && I_co[i]==NULL)
                              {
                                printf("memory reallocation failed in initial_class\n");
                               }}*/ /**(I_co[0]+l)=Atomx;
                      *(I_co[1]+l)=Atomy;
                      *(I_co[2]+l)=Atomz;*/
                          /*  for(j=0;j<N_int;j++){x=*(I_co[0]+j),y=*(I_co[1]+j),z=*(I_co[2]+j);
                                xp=(x+1)%Nx,yp=(y+1)%Ny,xm=(x-1+Nx)%Nx,ym=(y-1+Ny)%Ny;    
                 if((displacementType[xp][y][z].Type!=2)&&(displacementType[x][yp][z].Type!=2)&&(displacementType[xm][y][z].Type!=2)&&(displacementType[x][ym][z].Type!=2)&&(displacementType[x][y][z-1].Type!=2)&&(displacementType[x][y][z+1].Type!=2))                 {    *(I_co[0]+j)=*(I_co[0]+N_int-iy);
                                                                    *(I_co[0]+j)=*(I_co[0]+N_int-iy);
                                                                    *(I_co[0]+j)=*(I_co[0]+N_int-iy);
                                                                               iy=iy+1;  }}
                                                             for(i=0;i<3;i++){  
                                                             I_co[i]=(int*)realloc(I_co[i],N_intnew*sizeof(int));
                                                             if(I_co[i]>0 && I_co[i]==NULL)
                                                             {
                                                                printf("memory reallocation failed in initial_class\n");}}
                                                                                                                           }//end of if
         
                        if(N_intnew>N_int){
                                 *(I_co[0]+r)=Atomx;
                                 *(I_co[1]+r)=Atomy;
                                 *(I_co[2]+r)=Atomz;
                                 printf("x=%d\ty=%d\tz=%d\n\n",*(I_co[0]+r),*(I_co[1]+r),*(I_co[2]+r));*/

                                 for(i=0;i<3;i++){
                            I_co[i]=(int*)realloc(I_co[i],N_intnew*sizeof(int));
                           if( I_co[i]==NULL)
                              {
                                printf("memory reallocation failed in initial_class\n");}}

                        for(k=0;k<24;k++){
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
			
   	 for(i=0;i<3;i++){I_co[i]=I_co[i]-N_intnew;}
/*	 for(i=0;i<N_intnew;i++){
            printf("i=%d\tj=%d\tk=%d\n",*(I_co[0]+i),*(I_co[1]+i),*(I_co[2]+i));}*/


			}// end of second if statement


}                                                           

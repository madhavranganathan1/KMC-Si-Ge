
extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

int EnBH(int t,int h,int x,int y,double *AAH,double *BBH,double *ABH)//Leteral Bond at height H
{unsigned  xp,xm,yp,ym;
       xp=(x+1)%Nx, xm=(x-1+Nx)%Nx, yp=(y+1)%Nx, ym=(y-1+Nx)%Nx;
  
  unsigned int TA[4];//Atom type
  unsigned int MD[4];
  unsigned int BA[3][4]={0,0,0,0,0,0,0,0,0,0,0,0};//Atom bond
  
  
    //double delta=0.25;  
double rE;
//if(h>2){rE=.02;}//(RE/(double)(h-/*+*/1));}
    //else{rE=.23}
    //printf("%lf_%d\n",rE,h);
if(h>2){rE=(RE/(double)(h-1));}

else{rE=RE;}
//rE=RE;
  TA[0] =displacementType[xp][y][h].Type;
  TA[1] =displacementType[xm][y][h].Type;
  TA[2] =displacementType[x][yp][h].Type;
  TA[3] =displacementType[x][ym][h].Type;
  
  MD[0]=t&TA[0];
  MD[1]=t&TA[1];
  MD[2]=t&TA[2];
  MD[3]=t&TA[3];

  if(MD[0]==1){ BA[0][0]=1;}
  else if(MD[0]==2){ BA[1][0]=1;}
  else{ MD[0]=t|TA[0];
           if(MD[0]==3){BA[2][0]=1;}
           else{BA[2][0]=0;}}

  if(MD[1]==1){ BA[0][1]=1;}
  else if(MD[1]==2){ BA[1][1]=1;}
  else{ MD[1]=t|TA[1];
           if(MD[1]==3){BA[2][1]=1;}
           else{BA[2][1]=0;}}

  if(MD[2]==1){ BA[0][2]=1;}
  else if(MD[2]==2){ BA[1][2]=1;}
  else{ MD[2]=t|TA[2];
           if(MD[2]==3){BA[2][2]=1;}
           else{BA[2][2]=0;}}

  if(MD[3]==1){ BA[0][3]=1;}
  else if(MD[3]==2){ BA[1][3]=1;}
  else{ MD[3]=t|TA[3];
           if(MD[3]==3){BA[2][3]=1;}
           else{BA[2][3]=0;}}
//============================================
//printf("BA[0][0]=%u, BA[0][1]=%u, BA[0][2]=%u, BA[0][3]=%u\n",BA[0][0], BA[0][1], BA[0][2], BA[0][3]);
//printf("BA[1][1]=%u, BA[1][2]=%u, BA[1][3]=%u, BA[0][4]=%u\n",BA[1][0], BA[1][1], BA[1][2], BA[1][3]);
//printf("BA[2][1]=%u, BA[2][2]=%u, BA[2][3]=%u, BA[0][4]=%u\n",BA[2][0], BA[2][1], BA[2][2], BA[2][3]);
//===========================================


     if(h%2==0){if((h-Height[x][yp])>1 || (h-Height[x][ym])>1){rE=0;}
     if(x%2==0){
           *AAH= ((double)BA[0][0]*(gammaNNAA+rE))+((double)BA[0][1]*(gammaNNAA-rE))+((double)BA[0][2]+(double)BA[0][3])*gammaNNAA;
           *BBH= ((double)BA[1][0]*(gammaNNBB+rE))+((double)BA[1][1]*(gammaNNBB-rE))+((double)BA[1][2]+(double)BA[1][3])*gammaNNBB;
           *ABH= ((double)BA[2][0]*(gammaNNAB+rE))+((double)BA[2][1]*(gammaNNAB-rE))+((double)BA[2][2]+(double)BA[2][3])*gammaNNAB;
 	       }

else if(x%2==1){
           *AAH= ((double)BA[0][0]*(gammaNNAA-rE))+((double)BA[0][1]*(gammaNNAA+rE))+((double)BA[0][2]+(double)BA[0][3])*gammaNNAA;
           *BBH= ((double)BA[1][0]*(gammaNNBB-rE))+((double)BA[1][1]*(gammaNNBB+rE))+((double)BA[1][2]+(double)BA[1][3])*gammaNNBB;
           *ABH= ((double)BA[2][0]*(gammaNNAB-rE))+((double)BA[2][1]*(gammaNNAB+rE))+((double)BA[2][2]+(double)BA[2][3])*gammaNNAB;
               }
               }
else if(h%2==1){if((h-Height[xp][y])>1 || (h-Height[xm][y])>1){rE=0;}
     if(y%2==0){
           *AAH= ((double)BA[0][0]+(double)BA[0][1])*gammaNNAA+((double)BA[0][2]*(gammaNNAA+rE))+((double)BA[0][3]*(gammaNNAA-rE));
           *BBH= ((double)BA[1][0]+(double)BA[1][1])*gammaNNBB+((double)BA[1][2]*(gammaNNBB+rE))+((double)BA[1][3]*(gammaNNBB-rE));
           *ABH= ((double)BA[2][0]+(double)BA[2][1])*gammaNNAB+((double)BA[2][2]*(gammaNNAB+rE))+((double)BA[2][3]*(gammaNNAB-rE));
               }
else if(y%2==1){
           *AAH= ((double)BA[0][0]+(double)BA[0][1])*gammaNNAA+((double)BA[0][2]*(gammaNNAA-rE))+((double)BA[0][3]*(gammaNNAA+rE));
           *BBH= ((double)BA[1][0]+(double)BA[1][1])*gammaNNBB+((double)BA[1][2]*(gammaNNBB-rE))+((double)BA[1][3]*(gammaNNBB+rE));
           *ABH= ((double)BA[2][0]+(double)BA[2][1])*gammaNNAB+((double)BA[2][2]*(gammaNNAB-rE))+((double)BA[2][3]*(gammaNNAB+rE));
               }
               }


return 0;
}


extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	
int DBH(int t,int h,int x,int y,int *DAAH,int *DBBH,int *DABH)
{unsigned  xp,xm,yp,ym;
       xp=(x+1)%Nx, xm=(x-1+Nx)%Nx, yp=(y+1)%Nx, ym=(y-1+Nx)%Nx;
   
   unsigned int TA[4],
                MD[4],
             BA[3][4]={0,0,0,0,0,0,0,0,0,0,0,0};

  TA[0] =displacementType[xp][yp][h].Type;
  TA[1] =displacementType[xp][ym][h].Type;
  TA[2] =displacementType[xm][yp][h].Type;
  TA[3] =displacementType[xm][ym][h].Type;
  
  MD[0] =t&TA[0];
  MD[1] =t&TA[1];
  MD[2] =t&TA[2];
  MD[3] =t&TA[3];

  
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
  
  *DAAH=BA[0][0]+BA[0][1]+BA[0][2]+BA[0][3];
  *DBBH=BA[1][0]+BA[1][1]+BA[1][2]+BA[1][3];
  *DABH=BA[2][0]+BA[2][1]+BA[2][2]+BA[2][3];
return 0;
}

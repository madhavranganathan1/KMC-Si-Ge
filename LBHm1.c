
extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

int LBHm1(int t,int h,int x,int y,int *LAAHm,int *LBBHm,int *LABHm)//Leteral Bond at height H
{
  unsigned int TA[1];//Atom type
  unsigned int MD[1];
  unsigned int BA[3][1]={0,0,0};//Atom bond

  TA[0] =displacementType[x][y][h-1].Type;
  
  MD[0]=t&TA[0];

  if(MD[0]==1){ BA[0][0]=1;}
  else if(MD[0]==2){ BA[1][0]=1;}
  else{ MD[0]=t|TA[0];
           if(MD[0]==3){BA[2][0]=1;}
           else{BA[2][0]=0;}}

           *LAAHm= BA[0][0];
           *LBBHm= BA[1][0];
           *LABHm= BA[2][0];

 	
return 0;
}


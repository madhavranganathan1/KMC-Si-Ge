
extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

int LBHm(int t,int h,int x,int y,int *LAAHm,int *LBBHm,int *LABHm,int *laahm ,int *lbbhm ,int *labhm)//Leteral Bond at height H
{
  unsigned int TA[2];//Atom type
  unsigned int MD[2];
  unsigned int BA[3][2]={0,0,0,0,0,0};//Atom bond

  TA[0] =displacementType[x][y][h-1].Type;
  TA[1] =displacementType[x][y][h+1].Type;

  MD[0]=t&TA[0];
  MD[1]=t&TA[1];

  

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



           *LAAHm= BA[0][0];
           *LBBHm= BA[1][0];
           *LABHm= BA[2][0];
           *laahm=BA[0][1];

           *lbbhm=BA[1][1];

           *labhm=BA[2][1];

return 0;
}


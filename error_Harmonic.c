int error_Harmonic(double *HE,int neighbours, int BONDXY)
{
  int error;
  int In;
  if(BONDXY==1 && neighbours==5){error=0;}
  if(BONDXY==2 && neighbours==5){error=0;}
  if(BONDXY==3 && neighbours==5){error=0;}
  if(BONDXY==4 && neighbours==5){error=0;}
  if(BONDXY==5 && neighbours==5){error=0;}





  for(In=0;In<=50;In++){
   if(HE[In]>0.0)printf("  %0.10lf \n",HE[In]);
  }

return error;
}

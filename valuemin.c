//#include<stdio.h>
double valuemin(double q,double r,double s,double t,double u,double v)
{
   double minimum,z[6];
   int i;
  
   z[0]=q,z[1]=r,z[2]=s,z[3]=t,z[4]=u,z[5]=v;
   for(i=0;i<6;i++){if(z[i]!=0 ){minimum=z[i];}}
  
   for(i=0;i<6;i++){
	   if(z[i]<minimum && z[i]!=0){minimum=z[i];}}
   printf("minimum value =%lf\n",minimum);
   return minimum;
}


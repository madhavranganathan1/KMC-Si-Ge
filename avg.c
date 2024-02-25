#include <stdio.h>

int main(void)
{
   FILE *fp;
   float dummy,avg;
   avg=0.;
  
  fp = fopen("avg.dat","r");
   while(fscanf(fp,"%f",&dummy)!=EOF)
       avg=avg+dummy;
   fclose(fp);
  avg=avg/20.0;
    printf("%f\n",avg);
          
fclose(fp);
}

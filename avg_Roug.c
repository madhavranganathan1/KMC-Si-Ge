#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main ()
{
int i,k;
double a=0;
int M=0;
double j,b,d1,d2;
//==============================
  FILE *ptr1;
  ptr1=fopen("count","r");
  while (fscanf(ptr1,"%lf\n",&d1)!=EOF){M++;}
      fclose(ptr1);
     // printf("Trajectry%d\n",M);
//==============================
        double data1[M], data2[M];
        FILE *ptr;
        ptr = fopen("count","r");
        for(i=0;i<M;i++){fscanf(ptr,"%lf\n", &data1[i]);}

       
        for(i=0;i<M;i++){
                         a=a+data1[i];
                         }
                        b =a/(double)M;
printf("%lf\n",b);

                        
return 0;
}


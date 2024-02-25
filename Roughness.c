#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

#define Nx 128
#define Ny 128
#define Nz 17
double H[Nx][Ny];
struct DT {
 double xdisp;
 double ydisp;
 double zdisp;
 int Type;
   } dT[Nx][Ny][Nz];

int main(int argc, char *argv[])
{
  int I,J,K;
  char dTfname[256];
  sprintf(dTfname,"CONFIG");
  FILE *fp;
  fp=fopen(dTfname,"rb");
  fread(dT,sizeof(struct DT), (Nx*Ny*Nz),fp);

  for (I=0;I<Nx;I++){
    for (J=0;J<Nx;J++){
       for(K=0;K<Nz;K++){
          if(dT[I][J][K].Type==0){ H[I][J]=(double)K-1;
                                       break;
                                     }}}}
double avgH=0.0;
    for (I=0;I<Nx;I++){
     for (J=0;J<Ny;J++){avgH=avgH+H[I][J];}}

     avgH=avgH/(double)(Nx*Ny);

double Rough=0.0;
    for (I=0;I<Nx;I++){
     for (J=0;J<Nx;J++){
         Rough=Rough+((H[I][J]-avgH)*(H[I][J]-avgH));
           }}
Rough=Rough/(double)(Nx*Ny);
Rough=pow(Rough,0.5);


printf("%f\n",Rough);
//printf("%f  %f \n",(double)atom/Nx*NY, Rough);
fclose(fp);
return 0;
}


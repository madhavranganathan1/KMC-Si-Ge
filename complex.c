#include<stdio.h>
#include<gsl/gsl_complex_math.h>

typedef struct
{
  double dat[2];
  gsl_complex;
}
int main()
{
  double x ,y, z;
  printf("enter the value of x and y\n");
  scanf("%lf %lf ",&x, &y);
  GSL_SET_COMPLEX(&z,x, y);
  return(z);
}

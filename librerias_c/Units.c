#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G  6.6740831e-11
double units(double uM, double uL, double uT);

int main(void){

  printf("Hola mundo \n");
  printf("uT = %lf \n",units(1.98e30 , 150.e9 , -1.0)/86400);
  return 0;
}

double units(double uM, double uL, double uT){
  if(uM == -1.0){
    return uL*uL*uL/(uT*uT*G);
  }

  else if(uL == -1.0){
    return pow(G*uM*uT*uT,1.0/3.0);
  }

  else if(uT == -1.0){
    return pow(uL*uL*uL/(G*uM),0.5);
  }

  return 0;
}



#include "variables.h"

int KMCSmereka2(double bondenergy,int Atomx,int Atomy,int Atomz) 
{
           
 double extensionsq,extension,energychange,toss,kmcfactor,harmonicenergy,E1,E2;;
 int Returnvalue;

double ref_Energy=.68;

//if(NoOfBond==1){energychange=-1.0;}
//else{
//  harmonicenergy=2.5*Harmonic(Atomx,Atomy,Atomz);
  energychange =/* - 0.500* harmonicenergy +*/ bondenergy - ref_Energy;
  //  }
  
//if(spring_Bond==5){harmonicenergy=0.00;}
//		      energychange = - 0.500* harmonicenergy + bondenergy - ref_Energy;
//  if(spring_Bond<5){energychange=-1.0;}


                        if (energychange<=0.0) {Returnvalue = 1;
                       
                                                }
                        else if (energychange > 0.0) {
                                toss = gsl_rng_uniform_pos(gnalea_r5);
                                kmcfactor = exp(-energychange/(kB * temperature));
                                if (toss <= kmcfactor) {
                                        Returnvalue = 1;

                                }
                                else
                                {
                                        Returnvalue = 0; // Move rejected 
                          
                                }
                        }

//if(Returnvalue==1 && (bondenergy-ref_Energy)>0 && spring_Bond==6){printf("%lf_%lf_bond_%d\tenergy=%lf\n",bondenergy-ref_Energy,.500* harmonicenergy,spring_Bond,energychange);}
return Returnvalue;
}

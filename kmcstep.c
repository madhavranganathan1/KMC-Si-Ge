#include "variables.h"
extern struct dT displacementType[Nx][Ny][Nz],disTy[Nx][Ny][Nz]; 	

void kmcstep(int Timestep) 
{

 	int NNbondsAA,NNNbondsAA,NNbondsBB,NNNbondsBB,NNbondsAB,NNNbondsAB,Atomnumber, AtomX, AtomY,AtomZ,Nomin,Mintest;
 	double initialenergy,energybond,energyvalue;
 	
 	Atomnumber = (int) (Nx*Ny*gsl_rng_uniform_pos(gnalea_r3) ) ;
 	AtomX = (int) (Atomnumber/Ny) ;
 	AtomY = (int) (Atomnumber % Ny) ;
 	AtomZ = Height[AtomX][AtomY];
          while(AtomZ==0){Atomnumber = (int) (Nx*Ny*gsl_rng_uniform_pos(gnalea_r3) ) ;
                AtomX = (int) (Atomnumber/Ny) ;
                AtomY = (int) (Atomnumber % Ny) ;
                AtomZ = Height[AtomX][AtomY];
               }


 	if ( (AtomX < 0) || (AtomX >= Nx ) || (AtomY < 0) || (AtomY >= Ny) ) {
 		printf("Problem with Random Atom choice\n");
                maxminHeight();
                FILE *tmpout;
                char tmpfile[256];
                sprintf(tmpfile,"tmpdisptype%dx%dx%dx%d",Nx,Ny,Nz,rank);
                tmpout = fopen(tmpfile,"wb");
                fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),tmpout);
                fclose(tmpout);
//                MPI_Finalize();
 		exit(0);}
 	
 	if ((displacementType[AtomX][AtomY][AtomZ].Type) == 0) {
 	    printf("Something wrong with adatom choice; its type is zero %d %d %d _rank%d\n",AtomX,AtomY,AtomZ,rank);
            FILE *tmpout;
                char tmpfile[256];
                sprintf(tmpfile,"tmpdisptype%dx%dx%dx%d",Nx,Ny,Nz,rank);
                tmpout = fopen(tmpfile,"wb");
                fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),tmpout);
                 fclose(tmpout);
//                MPI_Finalize();
 	         exit(0);}
 	
 	if (AtomZ == -1) {
 		printf ("Reached continuous substrate\n");
                maxminHeight();
                FILE *tmpout;
                char tmpfile[256];
                sprintf(tmpfile,"tmpdisptype%dx%dx%dx%d",Nx,Ny,Nz,rank);
                tmpout = fopen(tmpfile,"wb");
                fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),tmpout);
                fclose(tmpout);
//                MPI_Finalize();
 		exit(0);}
 	
 	if ((displacementType[AtomX][AtomY][AtomZ+1].Type) != 0) {
 	     printf("Problem with Random Atom choice; there is an atom above it!!! %d %d %d\n",AtomX,AtomY,AtomZ);
             maxminHeight();
                FILE *tmpout;
                char tmpfile[256];
                sprintf(tmpfile,"tmpdisptype%dx%dx%dx%d",Nx,Ny,Nz,rank);
                tmpout = fopen(tmpfile,"wb");
                fwrite(displacementType,sizeof(struct dT), (Nx*Ny*Nz),tmpout);
                fclose(tmpout);
//                MPI_Finalize();
 		exit(0);}
	
 	NNbondsAA = 0, NNbondsBB = 0, NNbondsAB = 0;
        int LAAHm=0, LBBHm=0, LABHm=0,
            LAAH=0,  LBBH=0, LABH=0,
            DAAH=0,  DBBH=0,  DABH=0,
            DAAHm=0, DBBHm=0, DABHm=0,
            DAAHp=0, DBBHp=0, DABHp=0;

        double AAH=0, BBH=0, ABH=0;
        
   int typ=displacementType[AtomX][AtomY][AtomZ].Type;

        LBHm1(typ,AtomZ,AtomX,AtomY,&LAAHm,&LBBHm,&LABHm);
//        LBH(typ,AtomZ,AtomX,AtomY,&LAAH,&LBBH,&LABH);
        EnBH(typ,AtomZ,AtomX,AtomY,&AAH,&BBH,&ABH);
        DBH (typ,AtomZ,AtomX,AtomY,&DAAH,&DBBH,&DABH);
        DBHm(typ,AtomZ,AtomX,AtomY,&DAAHm,&DBBHm,&DABHm);
        DBHp(typ,AtomZ,AtomX,AtomY,&DAAHp,&DBBHp,&DABHp);


// 	NNbondsAA = LAAH;
// 	NNbondsBB = LBBH;
// 	NNbondsAB = LABH;

 	
 	/* That ends the calculation of the number of bonds broken */
 //	double NNenergybond   = (((double)NNbondsAA) * gammaNNAA)  +(((double)NNbondsBB)  * gammaNNBB)  + (((double)NNbondsAB)  * gammaNNAB);
        double NNenergybond    = AAH+BBH+ABH;

        double NNenergybondHm  = (((double)LAAHm)     * gammaNNAAHm)+(((double)LBBHm)      * gammaNNBBHm)+ (((double)LABHm)    * gammaNNABHm);
        double NNNenergybond   = (((double)DAAH)      * gammaNNNAA)+(((double)DBBH)        * gammaNNNBB)+ (((double)DABH)      * gammaNNNAB);
        double NNNenergybondHm = ((((double)DAAHm)    * gammaNNNAAHm)+(((double)DBBHm)     * gammaNNNBBHm)+ (((double)DABHm)   * gammaNNNABHm));
        double NNNenergybondHp = (((double)DAAHp)     * gammaNNNAA)+(((double)DBBHp)       * gammaNNNBB)+ (((double)DABHp)     * gammaNNNAB);


 energybond=NNenergybond+NNenergybondHm + NNNenergybond+ NNNenergybondHm + NNNenergybondHp; 
int BONDXY=1+Totalbond(AtomX,AtomY);
//========================================================================================


       if ((BONDXY > 1) && (KMCSmereka1(energybond,AtomX,AtomY,AtomZ) == 1)) {
                                           update(AtomX,AtomY);
                                   if(CHK==10000){kmc_optimization();
                                              CHK=0;}
                                   else {CHK=CHK+1;}
                                              B3++;}

          else if((BONDXY==1)/* && Height[AtomX][AtomY]>DiscreteAlayers*/ &&(KMCSmereka2(energybond,AtomX,AtomY,AtomZ,BONDXY)==1)){
                          int BondAU=update(AtomX,AtomY);//Bond of atom XY Afer Update
                            if(BONDXY+BondAU>2){if(CHK==10000){kmc_optimization();
                                                         CHK=0;}
                                               else {CHK=CHK+1;}
                                               }
                                                 B2++;}


             //else if(BONDXY==1 && Height[AtomX][AtomY]==DiscreteAlayers){
               //            int BondAU=update(AtomX,AtomY);//Bond of atom XY Afer Update
               //               if(BONDXY+BondAU>2) {kmc_optimization();}
               //                                    B1++;}


//==========================================================================
/*            if(BONDXY==1 && Height[AtomX][AtomY]==(DiscreteAlayers)){update(AtomX,AtomY);}

            else if( (BONDXY==1) && Height[AtomX][AtomY]>(DiscreteAlayers) && (KMCSmereka(energybond,AtomX,AtomY,AtomZ) == 1) ){
                   update(AtomX,AtomY);//Bond of atom XY Afer Update
                         }

            else if ((BONDXY > 1) && (KMCSmereka(energybond,AtomX,AtomY,AtomZ) == 1)) {update(AtomX,AtomY);}
            if((Timestep%100000)==0){kmc_optimization();}
*/
//if(KMCSmereka(energybond,AtomX,AtomY,AtomZ)== 1){update(AtomX,AtomY);}
}
/* ---------------------END OF KMCSTEP ------------------------------------------------*/

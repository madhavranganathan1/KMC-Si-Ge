#include "variables.h"

int Totalbond(int AddatomX, int AddatomY)
{
        
      int NNbondsAA,NNNbondsAA,NNbondsBB,NNNbondsBB,NNbondsAB,NNNbondsAB;
      int AddatomZ=Height[AddatomX][AddatomY];

        NNbondsAA = 0;
        NNNbondsAA = 0;
        NNbondsBB = 0;
        NNNbondsBB = 0;
        NNbondsAB = 0;
        NNNbondsAB = 0;

 //int    LAAH = 0, LBBH = 0, LABH = 0, DAAH = 0, DBBH = 0, DABH = 0, DAAHm = 0, DBBHm = 0, DABHm = 0, DAAHp = 0, DBBHp = 0, DABHp = 0;
 int    LAAH = 0, LBBH = 0, LABH = 0;

   int typ=displacementType[AddatomX][AddatomY][AddatomZ].Type;

        LBH(typ,AddatomZ,AddatomX,AddatomY,&LAAH,&LBBH,&LABH);
        //DBH(typ,AddatomZ,AddatomX,AddatomY,&DAAH,&DBBH,&DABH);
        //DBHm(typ,AddatomZ,AddatomX,AddatomY,&DAAHm,&DBBHm,&DABHm);
        //DBHp(typ,AddatomZ,AddatomX,AddatomY,&DAAHp,&DBBHp,&DABHp);

        NNbondsAA = LAAH;
        //NNNbondsAA = DAAH+DAAHm+DAAHp;
        NNbondsBB = LBBH;
        //NNNbondsBB = DBBH+DBBHm+DBBHp;
        NNbondsAB = LABH;
        //NNNbondsAB = DABH+DABHm+DABHp;


//return (NNbondsAA + NNNbondsAA + NNbondsBB + NNNbondsBB + NNbondsAB + NNNbondsAB);
return ( NNbondsAA + NNbondsBB +  NNbondsAB );
}

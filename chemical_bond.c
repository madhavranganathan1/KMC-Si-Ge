
double chemical_bond(int s,int x,int y,int z)
{
       double energybond;
       int LAAHm=0, LBBHm=0, LABHm=0,
            laahm=0,lbbhm=0,labhm=0,
            LAAH=0,  LBBH=0, LABH=0,
            DAAH=0,  DBBH=0,  DABH=0,
            DAAHm=0, DBBHm=0, DABHm=0,
            DAAHp=0, DBBHp=0, DABHp=0;

        double AAH=0, BBH=0, ABH=0;

   int typ=s;      //displacementType[x][y][z].Type;

        LBHm(typ,z,x,y,&LAAHm,&LBBHm,&LABHm,&laahm,&lbbhm,&labhm);
//        LBH(typ,AtomZ,AtomX,AtomY,&LAAH,&LBBH,&LABH);
        EnBH(typ,z,x,y,&AAH,&BBH,&ABH);
        DBH (typ,z,x,y,&DAAH,&DBBH,&DABH);
        DBHm(typ,z,x,y,&DAAHm,&DBBHm,&DABHm);
        DBHp(typ,z,x,y,&DAAHp,&DBBHp,&DABHp);


//      NNbondsAA = LAAH;
//      NNbondsBB = LBBH;
//      NNbondsAB = LABH;


        /* That ends the calculation of the number of bonds broken */
 //     double NNenergybond   = (((double)NNbondsAA) * gammaNNAA)  +(((double)NNbondsBB)  * gammaNNBB)  + (((double)NNbondsAB)  * gammaNNAB);
        double NNenergybond    = AAH+BBH+ABH;

        double NNenergybondHm  = (((double)LAAHm)     * gammaNNAAHm)+(((double)LBBHm)      * gammaNNBBHm)+ (((double)LABHm)    * gammaNNABHm)+(((double)laahm)* gammaNNAAHm)+(((double)lbbhm)* gammaNNBBHm)+ (((double)labhm)    * gammaNNABHm);

        double NNNenergybond   = (((double)DAAH)      * gammaNNNAA)+(((double)DBBH)        * gammaNNNBB)+ (((double)DABH)      * gammaNNNAB);
        double NNNenergybondHm = ((((double)DAAHm)    * gammaNNNAAHm)+(((double)DBBHm)     * gammaNNNBBHm)+ (((double)DABHm)   * gammaNNNABHm));
        double NNNenergybondHp = (((double)DAAHp)     * gammaNNNAA)+(((double)DBBHp)       * gammaNNNBB)+ (((double)DABHp)     * gammaNNNAB);

        
        energybond=NNenergybond+NNenergybondHm + NNNenergybond+ NNNenergybondHm + NNNenergybondHp;
        int BONDXY=LAAHm+LBBHm+LABHm+laahm+labhm+lbbhm+DAAH+DBBH+DABH+DAAHm+DBBHm+DABHm+DAAHp+DBBHp+DABHp;
     //   printf("bonds=%d\tnnh=%lf\tnnv=%lf\tnnnh=%lf\tnnnhm=%lf\tnnnhp=%lf\n",BONDXY,NNenergybond,NNenergybondHm, NNNenergybond,NNNenergybondHm, NNNenergybondHp );

        return energybond;
}

#define Dim 3
#define Nx  128 
#define Ny  128
#define Nz  24
#define Nz_old  16

#define Totalsweeps   10676000//1ML/aec
#define Naddatom      10676000//5871655  
#define ML 6
#define FLUXend (int) (ML*((double)Nx)*((double)Ny))

#define Printperlayer  256       /* The results are printed out Printperlayer times 
each completed sweep */
#define prnt  256    //Cov=prnt*(Totalsweeps/Naddatom)

#define mismatch  0.0430//0.0835 /* Lattice mismatch */
#define relativefluxA 0.0 /* Relative flux of A with respect to B in the film. Thus 
the film has form A_xB_1-x */
#define relativehoppingenergy 0.0 


#define Initialsurface 1 /* This is 1 for flat initial surface or uniform vicinal 
surface, 2 for given initial surface input from a file KMCInitialCondition. If the 
initial surface is vicinal, then one has to define the number of steps on the initial
 surface. If this number is zero, we recover the flat initial surface. */

#define Nsteps 0    /* This is zero for a flat surface and nonzero for vicinal surfaces. 
                                                */

#define KMCtype 1 

#define Nhopsperequilibrium   15165
#define Nlighteratomratio 1
#define Maxnomintry 12/* Maximum number of minimizations attempted before the program gives up and says "Energy Minimization failed" */

#define Energytype 1 /* This defines the type of Energy function used. We will use 1
 to denote harmonic springs */

#define kspringLAA  .73  
#define kspringDAA  .96 
#define kspringLBB  .73 
#define kspringDBB  .96 
#define kspringLAB  .73 
#define kspringDAB  .96 


#define lespringLAA 2.7155  /*  2.71 angstrom unit spring equilibrium length for A-A lateral */


/* Bond Energies are described by gamma variables */
#define gammaNNAA 0.24     
#define gammaNNBB 0.19
#define gammaNNAB 0.21

#define gammaNNNAA 0.03
#define gammaNNNBB 0.03
#define gammaNNNAB 0.03

#define gammaNNNAAHm 0.10
#define gammaNNNBBHm 0.10
#define gammaNNNABHm 0.10

#define gammaNNAAHm 0.47//18.55 //These are strength of Bonds At the height (H-1)
#define gammaNNBBHm 0.28 
#define gammaNNABHm 0.40
//=============================
#define RE 0.14
int coverag;
//=============================
#define DiscreteAlayers 2 
#define Film 0
#define SpringType 1 /* SpringType 1 corresponds to Harmnic springs with a square root and SpringType 2 
                        corresponds to linearized Harmonic springs */

#define temperature .0604////K
#define kB  1//8.617e-5 
#define pi 3.141592653 

double ss,tg,gt;
 /* GSL Conjugate gradient parameters */
#define ssizeastep ss     /*step size */
#define tolaglobal tg    /* global tolerance */
#define gtolagrad  gt    /* Gradient tolerance */
#define optimization_type  4 //[0 conjugate_fr] [1 conjugate_pr]  [2 vector_bfgs2]  [3 vector_bfgs] [4 steepest_descent] 

int Height[Nx][Ny],spring_Bond;
double energyvalue;
double lespringDAA,lespringLBB,lespringDBB,lespringLAB,lespringDAB,arefx,arefy,arefzA,arefzB,arefzAB ;   
int rank,size,rc,B1,B2,B3,*I_co[3];

#define Diff(X,Y) X>Y?(X-Y):(Y-X)

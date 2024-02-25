 rm *.o
 rm Printout
 rm  tmpdisptype***
 rm Totaltime
 clear

   echo -e '"........This iS KMC with ELASTICITY!"\n'
   echo '"......COMPILING YOUR PROGRAM please wait ...."'

export OMPI_MPICC=gcc
export CC=gcc
   gcc -I/opt/apps/openmpi/3.1.5/include -fopenmp -I/home/nidhig/gsl/include -I/home/nidhig/fftw/include -I/home/nidhig/CLAPACK-3.1.1.1/INCLUDE  -c KMCheteroepitaxy.c -o object.o

   #gcc -I/opt/intel/openmpi-1.4.4/include/  -fopenmp -c KMCheteroepitaxy.c -o  object.o
if [ -a ./object.o ]
  then
    echo -e "......congratulation COMPILATION is SUCCESSFUL \n....NOW LINKING object files to the LIBRARY FILES, \nPlease wait for a while ..."
 else
    echo "....COMPILATION is not SUCCESSFUL "
    exit 1;
 fi
   gcc object.o -fopenmp -I/opt/apps/openmpi/3.1.5/lib -L/home/nidhig/gsl/lib -lgsl -lgslcblas -L/home/nidhig/fftw/lib -lfftw3 -L/home/nidhig/lib  -llapack -lfblaswr -lblas -lf2c -lm   -o KMCtest
#mpicc object.o                                                              -L/opt/intel/Compiler/11.1/046/lib/intel64 -limf  -fopenmp -o KMC
#mpicc  object.o -fopenmp -I/opt/intel/openmpi-1.4.4/lib -L/usr/local/lib -lgsl -lgslcblas -llapack -lblas -lF77 -lI77 -lfftw3 -limf -lm  -o KMCtest

if [ -x ./KMCtest ]; then
    echo -e '"....LINKING SUCCESSFUL"\n\n'
else
    echo ".....SORRY LINKING is not successful"
    exit 1;
 fi

#if [ -x ./KMCtest ]; then
 #   echo -e '"Want to run the program "\n\n'
#fi
 #          OPTIONS="Y N"
  #         select opt in $OPTIONS; do
   #            if [ "$opt" = "Y" ]; then

#qsub submit.sh
#mpirun -np 4 
#./KMCtest
 #exit
  #             elif [ "$opt" = "N" ]; then
   #             exit
    #           else
     #           clear
      #          echo bad option
       #        fi
        #   done

#exit 0;
-- INSERT --                                  
    




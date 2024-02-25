 rm KMCtest
 rm *.o
 rm Printout
 rm  tmpdisptype***
 rm Totaltime
 clear

   echo -e '"........This iS KMC with ELASTICITY!"\n'
   echo '"......COMPILING YOUR PROGRAM please wait ...."'

export OMPI_MPICC=gcc
export CC=gcc
# mpicc -I/home/paramita/include/ -I/home/paramita/fftw/include -I/home/paramita/LAPACK/CLAPACK-3.1.1.1/INCLUDE -fopenmp -c KMCheteroepitaxy.c -o object.o

  mpicc -I/opt/intel/openmpi-1.4.4/include/  -fopenmp -c KMCheteroepitaxy.c -o  object.o
if [ -a ./object.o ]
  then
    echo -e "......congratulation COMPILATION is SUCCESSFUL \n....NOW LINKING object files to the LIBRARY FILES, \nPlease wait for a while ..."
 else
    echo "....COMPILATION is not SUCCESSFUL "
    exit 1;
 fi
#mpicc object.o -L/usr/mpi/gcc/openmpi-1.6.5/lib64 -L/home/paramita/fftw/lib -lfftw3 -L/home/paramita/lib -lgsl -lgslcblas -llapack -lblas -lf2c -lm -fopenmp  -o KMCtest
#mpicc object.o                                                              -L/opt/intel/Compiler/11.1/046/lib/intel64 -limf  -fopenmp -o KMC
mpicc  object.o -fopenmp -I/opt/intel/openmpi-1.4.4/lib -L/usr/local/lib -lgsl -lgslcblas -llapack -lblas -lF77 -lI77 -lfftw3 -limf -lm  -o KMCtest

if [ -x ./KMCtest ]; then
    echo -e '"....LINKING SUCCESSFUL"\n\n'
else
    echo ".....SORRY LINKING is not successful"
    exit 1;
 fi

if [ -x ./KMCtest ]; then
    echo -e '"Want to run the program "\n\n'
fi
           OPTIONS="Y N"
           select opt in $OPTIONS; do
               if [ "$opt" = "Y" ]; then
               
#qsub submit.sh
#mpirun -np 4 
./KMCtest
 exit
               elif [ "$opt" = "N" ]; then
                exit
               else
                clear
                echo bad option
               fi
           done

exit 0;

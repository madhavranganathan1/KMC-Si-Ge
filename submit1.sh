#!/bin/bash
#PBS -N p22_k1.5 
#PBS -q batch
#PBS -l nodes=compute-0-2-ib:ppn=1,walltime=43200000
#PBS -j oe
cd $PBS_O_WORKDIR
#export I_MPI_FABRICS=shm:dapl
#export I_MPI_MPD_TMPDIR=/scratch/paramita
#export OMP_NUM_THREADS=8

mpirun -machinefile $PBS_NODEFILE  ./KMCtest cat>>Printout


#!/bin/bash
#PBS -N K5test_1
#PBS -q batch
# -l nodes=4:ppn=8,walltime=43200000
#PBS -l select=1:ncpus=8:mpiprocess
#PBS -o out.log
#PBS -e err.log

cd $PBS_O_WORKDIR

mpirun ./KMC > Printout
# -machinefile $PBS_NODEFILE -n ./KMC >> Printout
#-npernode 1 ./KMC >> Printout


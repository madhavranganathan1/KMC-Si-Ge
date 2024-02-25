#!/bin/bash
export OMPI_MPICC=gcc
mpicc -I/opt/intel/openmpi-1.4.4/include/  -fopenmp -c test.c -o  object.o
mpicc object.o -L/opt/intel/Compiler/11.1/046/lib/intel64 -limf  -fopenmp -o KMC
qsub  submit_test.sh


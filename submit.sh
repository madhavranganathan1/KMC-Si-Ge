#!/bin/bash
#PBS -N 700K_hanneal
#PBS -q sim
#PBS -l nodes=1:ppn=10
cd $PBS_O_WORKDIR
echo Start Job
date
LD_LIBRARY_PATH=/home/nidhig/gsl/lib
export LD_LIBRARY_PATH



(cd /home/nidhig/kmc_128/700K_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K1_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K2_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K4_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K5_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K6_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K7_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/700K8_1.7_0.3/nidhi; ./KMCtest>>printout1) &

(cd /home/nidhig/kmc_128/nidhi; ./KMCtest>>printout1) &


sleep 168h




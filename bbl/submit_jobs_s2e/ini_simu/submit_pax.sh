#!/bin/bash

#set job name:
#SBATCH -J s2e

#128 tasks on 4 nodes:
#SBATCH -N 4

#one process for each physical core:
#SBATCH -c 2

#run on the broadwell partition (pax11):
#SBATCH -p broadwell

#runtime of 20 minutes:
#SBATCH -t 47:59:59

#copy environment variables from submit environment:
#SBATCH --get-user-env
#send mail on all occasions:
#SBATCH --mail-type=ALL

cores=128

module purge
module add gnu13 impi

# impt
#----------
cd ./00_impt
genimpacttin lte.impt line
mpirun -np $cores ImpactT.exe 

# impz
#----------
cd ../01_impz
mv ../00_impt/fort.50 ./particle.in
addlineNum particle.in

genimpactzin lte.impz line
mpirun -np $cores ImpactZ.exe

# genesis-V4
#--------------------------
impz2sliceh5 gen4_dist.1003
mv gen4_sliced.h5 ../02_gen4

cd ../02_gen4
module purge && module add python/3.9 mpi/openmpi-x86_64 phdf5/1.14.4 szip/2.1.1
mpirun -np $cores genesis4 gen4.in




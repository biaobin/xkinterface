#!/bin/bash

# one process for each physical core:
#SBATCH -c 2

# reserved for debugging
##SBATCH --reservation=pitz_28  

# run on the haswell partition (pax11):
##SBATCH -p haswell
##SBATCH -p backfill
#SBATCH -p broadwell

# run on the rome partition (pax12):
##SBATCH -p rome

# runtime 
#SBATCH --time=47:59:59

# number of tasks/cores, processes
#SBATCH --ntasks=64

# Job name:
#SBATCH --job-name=test
#SBATCH --error=test-%N-%j.err
#SBATCH --output=test-%N-%j.out

# copy environment variables from submit environment
#SBATCH --get-user-env

# send mail on all occasions:
#SBATCH --mail-type=ALL

#/project/singularity/images/pax.img
#source /afs/ifh.de/group/pitz/data/lixiangk/work/apps/python3/3.9.18/bin/activate

#module purge
module add python/3.9
module add gnu14 openmpi5

# Get input from command line, has to come after setting up the partition
use_mpi=${1:-1}
if [ "$use_mpi" == "1" ]; then
  ntasks=128
elif [ "$use_mpi" == "0" ]; then
  ntasks=64
fi
echo $ntasks

jobname=${2:-test4}

script=test_pymoo_MPI2.py
direc='test'
direc=$jobname
echo $direc

mkdir -p $direc
cd $direc

cp ../$script ./$script
cp ../my_eval.py ./my_eval.py
cp ../sub.sh ./sub.sh

if [ "$use_mpi" == "1" ]; then
  mpirun python $script $use_mpi
  #srun python $script
elif [ "$use_mpi" == "0" ]; then
  python $script $use_mpi
fi

cd ..
exit

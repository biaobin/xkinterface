#!/bin/bash

# one process for each physical core:
#SBATCH -c 2

# run on the haswell partition (pax11):
##SBATCH -p haswell
#SBATCH -p broadwell

# runtime 
#SBATCH --time=47:59:59

# number of tasks/cores, processes
#SBATCH --ntasks=1

# Job name:
#SBATCH --job-name=test
##SBATCH --error=test-%N-%j.err
##SBATCH --output=test-%N-%j.out
#SBATCH --error=test.err
#SBATCH --output=test.out

# copy environment variables from submit environment
#SBATCH --get-user-env

# send mail on all occasions:
##SBATCH --mail-type=ALL

cd .
module add genesis/4.6.6
python run_job4.py 1  2  2>&1 | tee run_job4.log

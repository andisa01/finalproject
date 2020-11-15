#!/bin/bash

#SBATCH -J GBS_SNP_filter
#SBATCH -p day
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=5:00:00
#SBATCH --mail-user a.andis@yale.edu
#SBATCH --mail-type ALL

./runall.sh

mkdir outs/

mv *.out outs/

#!/bin/bash

#SBATCH --job-name=iqtree_test
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load IQ-TREE/1.6.12

iqtree -s ../../2.subsampling/YMFmain/YMFmain.subset.0.50_0.50_1000.min1.phy -nt AUTO


#!/bin/bash

#SBATCH --job-name=iqtree_test
#SBATCH --time=30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load IQ-TREE/1.6.12

iqtree -s ../../2.subsampling/YMFpilot/YMFpilot.0.5_0.5_1000.min1.phy -nt AUTO


#!/bin/bash

# Define variables
VCFIN=~/project/2019_RASY_YMF_PopGen/filtered_vcfs/filter_SNPs/YMFmain/YMFmain.subset.0.50_0.50.vcf
N_LOCI=1000

# Set up the directory
mkdir code

cp -R ../code/ .

cd code

# Run the script
# Note, eventually I will need to loop through this
./subVCF2phylip.sh $VCFIN $N_LOCI

# Move the file out to the main directory
mv *.phy ../

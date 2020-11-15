#!/bin/bash

# Define variables
VCFIN=../../../1.filtering/YMFpilot/GBS_SNP_filter/YMFpilot.0.5_0.5.vcf
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

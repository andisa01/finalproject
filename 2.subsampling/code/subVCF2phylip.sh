#!/bin/bash
# This script subsamples the SNP from a VCF and outputs a phylip file
# This script calls vcf2phylip.py from Ortiz:
# Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861

# Load modules

module load Python

# Define variables
VCFIN=$1  
N_LOCI=$2
VCFOUT=$(basename $VCFIN .vcf)
SUB=${VCFOUT}_${N_LOCI}

# Creat a subsample of the VCF
head -n 11 $VCFIN > head.vcf 

shuf -n $N_LOCI $VCFIN > body.vcf 

cat head.vcf body.vcf > $SUB.vcf 

# Convert the file to phylip
python vcf2phylip.py --input $SUB.vcf -m 1

# Remove intermediate files
rm head.vcf body.vcf $SUB.vcf

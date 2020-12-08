#!/bin/bash

module load VCFtools/

VCFINPATH=../0.filterSNPs/YMFcomb.oneSNP.vcf
#VCFIN=$(basename $VCFINPATH)

#LOCI_THRESHOLD=$(echo 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99)
LOCI_THRESHOLD=$(echo 0.34 0.36 0.38 0.42 0.44 0.46 0.48 0.52 0.54 0.56 0.58 0.62 0.64 0.66 0.68 0.72 0.74 0.76)
#LOCI_THRESHOLD=$(echo 0.1 0.9)

for i in $LOCI_THRESHOLD
do 
  OUT=$(echo YMFcomb_${i})
  vcftools --vcf $VCFINPATH --max-missing $i --missing-indv --out $OUT
done

#!/bin/bash

# Load required modules
module load R/3.5.3-foss-2018a-X11-20180131
module load PLINK/2.00-alpha2.3
module load VCFtools/0.1.15-foss-2018a-Perl-5.26.1

# Define variables
VCFINPATH=~/project/2019_RASY_YMF_PopGen/consensus_reads/YMF_main/myersPopGen_main_minsamp013_outfiles/myersPopGen_main_minsamp013.vcf
VCFIN=$(basename $VCFINPATH)

# Set up the directories
cp -r ../code/GBS_SNP_filter-master/ .

mv GBS_SNP_filter-master/ GBS_SNP_filter/

cd GBS_SNP_filter/

# Copy in the data
cp ../../data/popmap_main.txt .

mv popmap_main.txt popmap.txt

ln -s $VCFINPATH

# Make a new parameters file
printf "$VCFIN\n0.5\n0.5\n0.05\n0.5\n17\n#CHROM\n\n" > GBS_SNP_filter.txt

# Run the script
chmod +x GBS_SNP_filter.sh

./GBS_SNP_filter.sh

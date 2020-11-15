We used extremely permissive thresholds for missingness when assembling the reads with ipyrad.

In this step, we filter based thresholds for missingness of loci across individuals and missingness of individuals across loci.
- We additionally filter out multiallelic sites, retaining only biallelic sites.
- And finally, we retain only a single SNP with the greatest coverage across individuals for each locus.
- This is performed with GBS_SNP_filter script by Alexander: Alexander, A. 2018. GBS_SNP_filter v1.17. Available from https://github.com/laninsky/GBS_SNP_filter

The final outputs used in this analysis include a moderate threshold of missingess of 50% and 50%.
That is, loci must be present in at least 50% of the samples to be retained and individuals must have loci for at least 50% of sites to be retained.



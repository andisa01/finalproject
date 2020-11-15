# Missing the Point: A sensitivity analysis of missing data on tree topology inferred from genomic data

## 1. Introduction:
Understanding how the patterns of ecology drive the process of evolution is one of the most fundamental questions in biology and harkens back to debates at the very outset of the field of genetics (Wright, 1932). More recently, biologists have come to understand that evolution can happen very quickly, sometimes reciprocally driving ecological patterns (Hendry and Kinnison, 1999; Hairston et al., 2005). This may be good news for contemporary wildlife populations if they can rapidly adapt to novel selection pressures in the Anthropocene (Gonzalez et al., 2013). However, humans have not only altered the global environmental conditions--we also disrupt the landscape matrix and movement of wildlife through expanding infrastructure (Forman and Alexander, 1998). Thus, understanding how population structure helps or hinders rapid adaptation is paramount for conservation (Kinnison and Hairston, 2007).

Wood frogs are an excellent model system with which to study population structure and rapid adaptation because their natural history is a story of climate adaptation (Lee-Yaw et al., 2008) and they live in natural metapopulations. Over the past two decades, research on wood frogs at Yale Myer’s Forests (YMF) in northeastern Connecticut has intimated that populations are undergoing fine-scale, microgeographic adaptation (e.g. (Freidenburg and Skelly, 2004; Skelly, 2004; Ligon and Skelly, 2009)) in response to selection pressures of warming climate ((Arietta et al., 2020), Gahm et al. 2020). However, just how gene flow affects microgeographic divergence is still under debate (Richardson et al., 2014). Understanding the fine-scale population structure among wood frog breeding ponds at Yale Myers Forest is an important step in understanding this phenomenon.

## 1.1 Goals:
My goal is to perform a sensitivity analysis of missingness in my wood frog dataset creating trees from concatenated SNPs with varying thresholds of missing data. Missing data can be a problem for large RADseq based datasets (Leaché et al., 2015). As additional individuals are included in the analysis, the chance that data will be missing across a given locus increases. In most phylogeny building software, missing or unknown values for a given site are dropped from the site-likelihood estimation, essentially estimating only the site-likelihood for the subtree containing extant values (Minh et al., 2020). Thus, the amount of missingness can bias phylogenetic estimates (Lewis, 2001). These biases can be exacerbated by the manner in which RAD-seq data are passed to phylogenetic programs (Leaché et al., 2015). Because individual RAD loci generally contain little information--at most a few and often zero polymorphic sites--they are usually concatenated into supermatrices. This can be accomplished by including all sequence data, including SNPs and invariants sites, only polymorphic sites, or a single randomly selected SNP from each locus (e.g. (Emerson et al., 2010; Wagner et al., 2013; Yoder et al., 2013). Dropping invariant sites increases computation speeds, but results in acquisition bias and leads specious estimates of excessive divergence (and overestimates of branch lengths) (Bertels et al., 2014). Methods have been developed to account for acquisition bias of concatenated SNPs (Kuhner et al., 2000; Lewis, 2001). However, the ability to correct for acquisition bias is negatively correlated with the amount of missing data (Leaché et al., 2015).

## 2. Methods:
In 2018, I collected 277 tissue sample from 41 wood frog populations at YMF and three outgroup populations from New Haven Co. and Guilford Co. Connecticut. The samples from YMF represent approximately 10% of each pond population. To avoid collecting full-siblings, embryos were collected from individual egg masses shortly after oviposition and reared to hatching. We sequenced the samples using double-digest RAD-seq (Peterson et al., 2012).

An initial set of 16 individuals representing the geographic and environmental range of populations at YMF was sequenced at high depth up to 4M reads per sample. This read depth yielded 12,378 loci. These pilot data were then used to determine the necessary sequencing effort for the full dataset. The remaining samples were genotyped using the same ddRAD protocols yielding ~120,000 reads per sample. Base calls were made with pyRAD (Eaton, 2014).

Initially, I will test for the best phylogenetic model using IQ-TREE (Minh et al., 2020) on a concatenation of full sites (SNPs and invariant sites) from my pilot dataset of 16 representative individuals. These data contain low levels of missingness in general and I will exclude any sites missing in more than 1 individual. Then, using the best fitting model, I will assess the effect sample-wise and site-wise missingness on tree topology using concatenated SNPs. First, I will assess the effect of removing individuals based on percentage of missing sites from 10% to 90% by increments of 10%. Second, I will assess the effect of removing sites based on percentage of missing individuals with thresholds of less than 1%, 5%, 10%, and 50% missingness. Finally, I will test the difference in topologies of trees generated from each dataset (I am still figuring out how to do this!)
  
## 3. Results

## 3.1 MVA
Much of my effort has been in processing the data upstream of this project. This involved performing a de novo assembly of the RADseq data in iPyrad (Eaton ). The assembly parameters used in iPyrad, especially the clustering threshold can greatly impact results and must be empirically tuned (McCartney-Melstad et al. 2019). An inappropriate threshold will result in over- or under-splitting loci, which can be especially problematic in amphibians which can have high levels of paralogous regions. Following the pipeline suggested in McCartney-Melstad et al. (2019), I chose an optimal clustering threshold and assembled the sequence data.

Although the current analysis only involved 16 samples sequenced to high depth as part of a pilot sequencing project, I also assembled the full metapopulation-wide dataset that includes 249 samples with multiple individuals per pond population and included outgroups from other sites in the state of Connecticut. I made two assemblies of this full dataset—a de novo assembly and a reference-based assembly—both with the same parameters as established for the pilot dataset. I used the draft genome for the North American bullfrog (Rana catesbeiana) for the reference (Hammond et al. 2017). Although both wood frogs and bull frogs are Ranids, they are not particularly closely related and there is contention over the placement of wood frogs in the Ranid phylogeny (Hillis & Wilcox 2005), this is the closest genome available.

Following assembly, I filtered each dataset to retain only biallelic, variable markers and only the single SNP with greatest depth across all samples for each locus. In order to reduce the sparseness of the resulting data matrices and reduce the dataset size, I additionally filtered the datasets to retain only loci represented in at least 50% of individuals and individuals with at least 50% of loci. This filtering step was performed with GBS_SNP_filter (Alexander 2018). I used vcf2phylip (Ortiz 2019) to covert the filtered vcf files to phylip format for further phylogenetic analysis.

For my initial minimum viable analysis, I estimated phylogenetic relationships and maximum-likelihood trees with IQ-TREE (Nguyen et al. 2015).

## 3.2 Issues:
In addition to these phylogenetic analysis, I also performed a principle component analysis on the full metapopulation dataset (Figure). There seem to be two distinct clusters among samples. Interestingly, all pond populations have representatives in both clusters. I performed many diagnostics to confirm that this is a biological signal and not methodological. These clusters persist even when using the reference mapped data. 

![Figure. PCA of full dataset with geographic locations](/figs/PCA_full.png)

I have two theories on these clusters that I am currently looking into. First, the clusters may represent clusters. Second, the clusters may represent a larger inversion or other structural polymorphism in the population.

In any case, I am worried that this signal will drive the phylogenetic estimates. Looking at the tree inferred for the full dataset, it certainly seems that there is a very week relationship between clades and pond populations as I would expect.

While looking into these issues, I’ve also realized that these type of fine-scale, intra-specific data violate all of the assumption of the IQ-TREE and most (all?) of the phylogenetic inference programs. Given the likely extreme level of incomplete lineage sorting from interbreeding populations, I’m not even sure if it makes sense to use these tools at all.


## 4. Discussion

These results indicate...

The biggest difficulty in implementing these analyses was...

If I did these analyses again, I would...


## References:
- Alexander, A. 2018. GBS_SNP_filter v1.17. Available from https://github.com/laninsky/GBS_SNP_filter
- Arietta, A. Z. A., Freidenburg, L. K., Urban, M. C., Rodrigues, S. B., Rubinstein, A., and Skelly, D. K. (2020). Phenological delay despite warming in wood frog Rana sylvatica reproductive timing: a 20‐year study. Ecography 52, 27. doi:10.1111/ecog.05297.
- Bertels, F., Silander, O. K., Pachkov, M., Rainey, P. B., and van Nimwegen, E. (2014). Automated reconstruction of whole-genome phylogenies from short-sequence reads. Mol. Biol. Evol. 31, 1077–1088. doi:10.1093/molbev/msu088.
- Eaton, D. A. R. (2014). PyRAD: assembly of de novo RADseq loci for phylogenetic analyses. Bioinformatics 30, 1844–1849. Available at: https://academic.oup.com/bioinformatics/article-abstract/30/13/1844/2422183.
- Emerson, K. J., Merz, C. R., Catchen, J. M., Hohenlohe, P. A., Cresko, W. A., Bradshaw, W. E., et al. (2010). Resolving postglacial phylogeography using high-throughput sequencing. Proc. Natl. Acad. Sci. U. S. A. 107, 16196–16200. doi:10.1073/pnas.1006538107.
- Forman, R. T. T., and Alexander, L. E. (1998). Roads and Their Major Ecological Effects. Annu. Rev. Ecol. Syst. 29, 231, C2. Available at: http://www.jstor.org/stable/221707.
- Freidenburg, L. K., and Skelly, D. K. (2004). Microgeographical variation in thermal preference by an amphibian. Ecol. Lett. 7, 369–373. doi:10.1111/j.1461-0248.2004.00587.x.
- Gahm, K., Arietta, A. Z. A., & Skelly, D. K. (2020). Temperature‐mediated tradeoff between development and performance in larval wood frogs (Rana sylvatica). Journal of Experimental Zoology A.
- Gonzalez, A., Ronce, O., Ferriere, R., and Hochberg, M. E. (2013). Evolutionary rescue: an emerging focus at the intersection between ecology and evolution. Philos. Trans. R. Soc. Lond. B Biol. Sci. 368, 20120404. doi:10.1098/rstb.2012.0404.
- Hairston, N. G., Ellner, S. P., Geber, M. A., Yoshida, T., and Fox, J. A. (2005). Rapid evolution and the convergence of ecological and evolutionary time. Ecol. Lett. 8, 1114–1127. doi:10.1111/j.1461-0248.2005.00812.x.
- Hammond, S. A., Warren, R. L., Vandervalk, B. P., Kucuk, E., Khan, H., Gibb, E. A., … Birol, I. (2017). The North American bullfrog draft genome provides insight into hormonal regulation of long noncoding RNA. Nature Communications, 8(1), 1433.
- Hendry, A. P., and Kinnison, M. T. (1999). Perspective: The pace of modern fife: measuring rates of contemporary microevolution. Evolution 53, 1637–1653. doi:10.2307/2640428.
- Hillis, D. M., & Wilcox, T. P. (2005). Phylogeny of the New World true frogs (Rana). Molecular Phylogenetics and Evolution, 34(2), 299–314.
- Kinnison, M. T., and Hairston, N. G. (2007). Eco-evolutionary conservation biology: contemporary evolution and the dynamics of persistence. Functional Ecology 21, 444–454. doi:10.1111/j.1365-2435.2007.01278.x.
- Kuhner, M. K., Beerli, P., Yamato, J., and Felsenstein, J. (2000). Usefulness of single nucleotide polymorphism data for estimating population parameters. Genetics 156, 439–447. Available at: https://www.ncbi.nlm.nih.gov/pubmed/10978306.
- Leaché, A. D., Banbury, B. L., Felsenstein, J., de Oca, A. N.-M., and Stamatakis, A. (2015). Short Tree, Long Tree, Right Tree, Wrong Tree: New Acquisition Bias Corrections for Inferring SNP Phylogenies. Syst. Biol. 64, 1032–1047. doi:10.1093/sysbio/syv053.
- Lee-Yaw, J. A., Irwin, J. T., and Green, D. M. (2008). Postglacial range expansion from northern refugia by the wood frog, Rana sylvatica. Mol. Ecol. 17, 867–884. doi:10.1111/j.1365-294X.2007.03611.x.
- Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete morphological character data. Syst. Biol. 50, 913–925. doi:10.1080/106351501753462876.
- Ligon, N. F., and Skelly, D. K. (2009). Cryptic divergence: countergradient variation in the wood frog. Evol. Ecol. Res. 11, 1099–1109.
- McCartney-Melstad, E., Gidiş, M., & Shaffer, H. B. (2019). An empirical pipeline for choosing the optimal clustering threshold in RADseq studies. Molecular Ecology Resources, 19(5), 1195–1204.
- Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., et al. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37, 1530–1534. doi:10.1093/molbev/msaa015.
- Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268–274.
- Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861
- Peterson, B. K., Weber, J. N., Kay, E. H., Fisher, H. S., and Hoekstra, H. E. (2012). Double digest RADseq: an inexpensive method for de novo SNP discovery and genotyping in model and non-model species. PLoS One 7, e37135. doi:10.1371/journal.pone.0037135.
- Richardson, J. L., Urban, M. C., Bolnick, D. I., and Skelly, D. K. (2014). Microgeographic adaptation and the spatial scale of evolution. Trends Ecol. Evol. 29, 165–176. doi:10.1016/j.tree.2014.01.002.
- Skelly, D. K. (2004). Microgeographic countergradient variation in the wood frog, Rana sylvatica. Evolution 58, 160–165. doi:10.1111/j.0014-3820.2004.tb01582.x.
- Wagner, C. E., Keller, I., Wittwer, S., Selz, O. M., Mwaiko, S., Greuter, L., et al. (2013). Genome-wide RAD sequence data provide unprecedented resolution of species boundaries and relationships in the Lake Victoria cichlid adaptive radiation. Mol. Ecol. 22, 787–798. doi:10.1111/mec.12023.
- Wright, S. (1932). The roles of mutation, inbreeding, crossbreeding, and selection in evolution. in Proceedings of the Sixth International Congress on Genetics, 355–366. Available at: http://www.blackwellpublishing.com/ridley/classictexts/wright.pdf.
- Yoder, J. B., Briskine, R., Mudge, J., Farmer, A., Paape, T., Steele, K., et al. (2013). Phylogenetic signal variation in the genomes of Medicago (Fabaceae). Syst. Biol. 62, 424–438. doi:10.1093/sysbio/syt009.



- For information on formatting text files with markdown, see https://guides.github.com/features/mastering-markdown/ . You can use markdown to include images in this document by linking to files in the repository, eg `![GitHub Logo](/images/logo.png)`.



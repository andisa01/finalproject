# Missing the Point: A sensitivity analysis of missing data on tree topology inferred from genomic data

## 1. Introduction:
Understanding how the patterns of ecology drive the process of evolution is one of the most fundamental questions in biology and harkens back to debates at the very outset of the field of genetics (Wright, 1932). More recently, biologists have come to understand that evolution can happen very quickly, sometimes reciprocally driving ecological patterns (Hendry and Kinnison, 1999; Hairston et al., 2005). This may be good news for contemporary wildlife populations if they can rapidly adapt to novel selection pressures in the Anthropocene (Gonzalez et al., 2013). However, humans have not only altered the global environmental conditions--we also disrupt the landscape matrix and movement of wildlife through expanding infrastructure (Forman and Alexander, 1998). Thus, understanding how population structure helps or hinders rapid adaptation is paramount for conservation (Kinnison and Hairston, 2007).

Wood frogs are an excellent model system with which to study population structure and rapid adaptation because their natural history is a story of climate adaptation (Lee-Yaw et al., 2008) and they live in natural metapopulations. Over the past two decades, research on wood frogs at Yale Myer’s Forests (YMF) in northeastern Connecticut has intimated that populations are undergoing fine-scale, microgeographic adaptation (e.g. (Freidenburg and Skelly, 2004; Skelly, 2004; Ligon and Skelly, 2009)) in response to selection pressures of warming climate ((Arietta et al., 2020), Gahm et al. 2020). However, just how gene flow affects microgeographic divergence is still under debate (Richardson et al., 2014). Understanding the fine-scale population structure among wood frog breeding ponds at Yale Myers Forest is an important step in understanding this phenomenon.

Traditionally, isolation-by-distance (IBD) is considered the standard null model for tests of spatial genetic relatedness (Jenkins et al., 2010). Genetic relatedness that correlates strongly with environment or ecology, termed isolation-by-environment (or -ecology) (IBE) (Wang & Bradburd, 2014) in opposition or excess of IBD is often considered evidence for adaptive divergence (Shafer & Wolf, 2013). However, it is important to recognize that signatures of IBE are difficult to determine if adaptation is driving divergence or if ecological partitioning is driving divergence (Räsänen & Hendry, 2008). To further complicate inferences, colonization/extinction dynamics may also drive differentiation between populations, especially in species with metapopulation structure as exhibited by wood frogs wherein founder effects of colonizers drive differentiation (isolation-by-colonization: IBC) (Nadeau et al. 2016; Orsini et al. 2013).

Missing data can be a problem for large RADseq based datasets (Leaché et al., 2015). As additional individuals are included in the analysis, the chance that data will be missing across a given locus increases. In most phylogeny building software, missing or unknown values for a given site are dropped from the site-likelihood estimation, essentially estimating only the site-likelihood for the subtree containing extant values (Minh et al., 2020). Thus, the amount of missingness can bias phylogenetic estimates (Lewis, 2001). These biases can be exacerbated by the manner in which RAD-seq data are passed to phylogenetic programs (Leaché et al., 2015). Because individual RAD loci generally contain little information--at most a few and often zero polymorphic sites--they are usually concatenated into supermatrices. This can be accomplished by including all sequence data, including SNPs and invariants sites, only polymorphic sites, or a single randomly selected SNP from each locus (e.g. (Emerson et al., 2010; Wagner et al., 2013; Yoder et al., 2013). Dropping invariant sites increases computation speeds, but results in acquisition bias and leads specious estimates of excessive divergence (and overestimates of branch lengths) (Bertels et al., 2014). Methods have been developed to account for acquisition bias of concatenated SNPs (Kuhner et al., 2000; Lewis, 2001). However, the ability to correct for acquisition bias is negatively correlated with the amount of missing data (Leaché et al., 2015).

RADseq, in particular, is susceptible to biased missingness due to mutations within the restriction enzyme cut sites (Cariou et al. 2013). Mutations of this sort lead to systematic loss of phylogenetic information since the greater the evolutionary divergence between samples, the higher the likelihood of mutations within cut sites.
A second problem is that, unlike sampling strategies that sequence a single gene deeply for each sample, most RADseq strategies seek to maximize the number of samples and genomic coverage, a tradeoff against depth at each locus. Thus, RADseq datasets generally have lots of data, but also quite a bit of missingness for individual samples or loci. Insufficient, uneven sequencing efforts. This is compounded by the fact that the only way to increase average sequencing depth for a given locus without increasing total sequencing is to include more cutters to reduce the number of loci recovered—but this simply increases the chances of allelic-dropout due to mutations at the site (Huang & Knowles, 2016; Buerkle & Gompert, 2013).
Third, RADseq is commonly used in applications where DNA is degraded. This is because RADseq relies on shearing DNA and sequencing only short fragments. This lends it to applications in ancient DNA or museum specimens where DNA is already degraded (Smith et al., 2020).
Fourth, the multi-part processing of RADseq, in laboratory and bioinformatics procedures, allows multiple points where technical errors can introduce non-random variation in locus recovery (Escudero et al., 2014). This is especially an issue in computationally assessing homologous versus paralogous loci in de novo sequence alignment (Eaton, 2014). 
The purposes of this study are two-fold. First, I determine the magnitude and impact of missing data in a RADseq dataset generated from wood frog tissue samples from Yale Myers Forest. Second, I consider the utility of phylogenetic methods (tree based) in contrast to population genetics methods (cluster based) in assessing the relative roles of IBD, IBE, and IBC driving population structure in wood frogs.
## 2. Methods:
In 2018, I collected 277 tissue sample from 41 wood frog populations at YMF and three outgroup populations from New Haven Co. and Guilford Co. Connecticut. The samples from YMF represent approximately 10% of each pond population. To avoid collecting full-siblings, embryos were collected from individual egg masses shortly after oviposition and reared to hatching. We sequenced the samples using double-digest RAD-seq (Peterson et al., 2012).

An initial set of 16 individuals representing the geographic and environmental range of populations at YMF was sequenced at high depth up to 4M reads per sample. This read depth yielded 12,378 loci. These pilot data were then used to determine the necessary sequencing effort for the full dataset. The remaining samples were genotyped using the same ddRAD protocols yielding ~120,000 reads per sample.

I performed a de novo assembly of the full set of samples in iPyRad (Eaton, 2014). The assembly parameters used in iPyrad, especially the clustering threshold can greatly impact results and must be empirically tuned (McCartney-Melstad et al. 2019). An inappropriate threshold will result in over- or under-splitting loci, which can be especially problematic in amphibians which can have high levels of paralogous regions. Following the pipeline suggested in McCartney-Melstad et al. (2019), I chose an optimal clustering threshold and assembled the sequence data.

Missingness
Following assembly, I filtered each dataset to retain only biallelic, variable markers and only the single SNP with greatest depth across all samples for each locus with GBS_SNP_filter (Alexander, 2018). In order to quantify missingness of sites within samples and samples across sites, I used VCFtools to filter out individual samples looping through thresholds from 1-99% missingness and output a report for each detailing the number of nonmissing loci for retained samples. I then looped through loci missingness for each set in R and visualized combinations of thresholds as heat maps. To determine the most appropriate joint thresholding values, I computed a weighted score that maximized samples and total SNPs retained while weighting against missingness ((Standardized number of samples + standardized number of SNPs) - (Locus-wise missingness + sample-wise missingness)). I used these values to output a moderately filtered dataset and strict (< 10% missingness) and lenient (< 90% missingness) datasets for further comparison.

Phylogenetic relationships
I tested for phylogenetic signal in my dataset by fitting a phylogenetic model with IQ-TREE (Minh et al., 2020) on a concatenation of full sites (SNPs and invariant sites).

PCA
As an alternative to model-based phylogenetic analysis, which is likely inappropriate for these data given high levels of ongoing gene flow, I tested for population structure using PCA as non-model-based clustering method.

Population-wise test for IBD versus IBC
While the prior analyses test for relationships between individual samples, structure can also be tested at the population-wise level using differentiation metrics, such as Fst. I tested for IBD versus IBC by conducting a network analysis using Fst values as edge weights. In IBD, Fst values should be greater between pairs of populations separated by the greatest geographic distance whereas with IBC Fst values will be greatest between population pairs with more dissimilar harmonic mean population sizes.

## 3. Results
As expected, more extreme filtering thresholds drastically reduce the number of samples and SNPs in the dataset.
![Fig 1. Samples retained filtering on maximum % missing loci or minimum % missingness among samples](/figs/MissingSamps.png)
![Fig 2. Total SNPs retained (x 10^6) filtering on maximum % missing loci or minimum % missingness among samples](/figs/MissingSNPs.png)

Weighted scoring of missingess against total data indicates that retaining samples with up to 50% missingness and loci present in at least 58% of samples results in the largest, most complete matrix and includes 265 of 276 samples.
![Fig 3. Weighted scores of missingess against total data](/figs/MissingScores.png)

Model-based cladogram does not recover pond populations as unique clades.
![Fig 4. Cladogram of full metapopulation data](/figs/YMFmain_phylo.png)

Rather, clustering-based PCA indicates that genetic relatedness is associated primarily with two distinct clusters that are not associated with geography.
![Fig 5. PCA of full dataset with geographic locations](/figs/PCA_full.png)

The SNPs driving these associations are not related to levels of missingness in the dataset.
![Fig 6. PC loadings of individual SNPs are not associated with locus frequency or missingness in the dataset](/figs/VariantPCloadings.png)

Instead of individual-based measures of differentiation, I analyzed population-wise levels of average differentiation with Wright’s Fst. Using Fst as edge weights, I considered if the disparity in harmonic population size or geographic distance associated with differentiation. While Fst was not strongly associated with geography, ponds with more similar harmonic mean population sizes tended to have lower pairwise genetic differentiation.
![Fig 7. The relationship between genetic similarity and geography among pond populations. Denser edge weights indicate lower Fst (i.e. lower differentiation). Point size indicates harmonic mean population size and color indicates extinction probability (red = high probability) based on census studies. To prevent overplotting, only edges in the lowest 2.5th, 5th, 10th and 20th cumulative percentiles are shown.](/figs/Network.png)
![Fig 8. Network diagram of populations. The spatial distribution of populations is based on network stress--population with lowest differentiation to between more pairs cluster closer to the center. Only edges in the lowest 20th percentiles are included. Populations without pairwise edges in the lowest 20th percentile are shown on the left.](/figs/VariantPCloadings.png)


## 4. Discussion

These results indicate that colonization/extinction dynamics likely play a larger role than geographic distance in structuring populations of wood frogs as Yale Myers Forest. More work is needed to compare ecological differences between population nodes to determine if IBC is more prevalent than IBE.

An important conclusion from these analyses is that model-based phylogenetic estimates of population structure are inappropriate for populations with ongoing gene flow. At the spatial scale of the data in this study, even individual based clustering that are unreliant on evolutionary models are prone to genetic structure that may not be related to population structure. For example, the two clusters uncovered in the PCA analysis may be defined by sex determining regions of the genome. In closely related populations with high gene flow, variance between populations may be obscure by variance between sexes. Ongoing efforts to sequence the wood frog genome should help in understanding the clusters seen in these data.


## References:
- Alexander, A. 2018. GBS_SNP_filter v1.17. Available from https://github.com/laninsky/GBS_SNP_filter
- Arietta, A. Z. A., Freidenburg, L. K., Urban, M. C., Rodrigues, S. B., Rubinstein, A., and Skelly, D. K. (2020). Phenological delay despite warming in wood frog Rana sylvatica reproductive timing: a 20‐year study. Ecography 52, 27. doi:10.1111/ecog.05297.
- Bertels, F., Silander, O. K., Pachkov, M., Rainey, P. B., and van Nimwegen, E. (2014). Automated reconstruction of whole-genome phylogenies from short-sequence reads. Mol. Biol. Evol. 31, 1077–1088. doi:10.1093/molbev/msu088.
- Buerkle A.C. Gompert Z. 2013. Population genomics based on low coverage sequencing: how low should we go? Mol. Ecol. 22:3028–3035.
- Cariou M. Duret L. Charlat S. 2013. Is RAD-seq suitable for phylogenetic inference? An in silico assessment and optimization. Ecol. Evol. 3:846–852.
- Eaton, D. A. R. (2014). PyRAD: assembly of de novo RADseq loci for phylogenetic analyses. Bioinformatics 30, 1844–1849. Available at: https://academic.oup.com/bioinformatics/article-abstract/30/13/1844/2422183.
- Emerson, K. J., Merz, C. R., Catchen, J. M., Hohenlohe, P. A., Cresko, W. A., Bradshaw, W. E., et al. (2010). Resolving postglacial phylogeography using high-throughput sequencing. Proc. Natl. Acad. Sci. U. S. A. 107, 16196–16200. doi:10.1073/pnas.1006538107.
- Escudero M. Eaton D.A. Hahn M. Hipp A.L. 2014. Genotyping-by-sequencing as a tool to infer phylogeny and ancestral hybridization: A case study in carex (cyperaceae). Mol. Phylogenet. Evol. 79:359–367.
- Forman, R. T. T., and Alexander, L. E. (1998). Roads and Their Major Ecological Effects. Annu. Rev. Ecol. Syst. 29, 231, C2. Available at: http://www.jstor.org/stable/221707.
- Freidenburg, L. K., and Skelly, D. K. (2004). Microgeographical variation in thermal preference by an amphibian. Ecol. Lett. 7, 369–373. doi:10.1111/j.1461-0248.2004.00587.x.
- Gahm, K., Arietta, A. Z. A., & Skelly, D. K. (2020). Temperature‐mediated tradeoff between development and performance in larval wood frogs (Rana sylvatica). Journal of Experimental Zoology A.
- Gonzalez, A., Ronce, O., Ferriere, R., and Hochberg, M. E. (2013). Evolutionary rescue: an emerging focus at the intersection between ecology and evolution. Philos. Trans. R. Soc. Lond. B Biol. Sci. 368, 20120404. doi:10.1098/rstb.2012.0404.
- Hairston, N. G., Ellner, S. P., Geber, M. A., Yoshida, T., and Fox, J. A. (2005). Rapid evolution and the convergence of ecological and evolutionary time. Ecol. Lett. 8, 1114–1127. doi:10.1111/j.1461-0248.2005.00812.x.
- Hammond, S. A., Warren, R. L., Vandervalk, B. P., Kucuk, E., Khan, H., Gibb, E. A., … Birol, I. (2017). The North American bullfrog draft genome provides insight into hormonal regulation of long noncoding RNA. Nature Communications, 8(1), 1433.
- Hendry, A. P., and Kinnison, M. T. (1999). Perspective: The pace of modern fife: measuring rates of contemporary microevolution. Evolution 53, 1637–1653. doi:10.2307/2640428.
- Hillis, D. M., & Wilcox, T. P. (2005). Phylogeny of the New World true frogs (Rana). Molecular Phylogenetics and Evolution, 34(2), 299–314.
- Huang H. Knowles L.L. 2016. Unforeseen Consequences of Excluding Missing Data from Next-Generation Sequences: Simulation study of RAD Sequences. Syst. Biol. 65:357–365.
- Jenkins, D. G., Carey, M., Czerniewska, J., Fletcher, J., Hether, T., Jones, A., … Tursi, R. (2010). A meta-analysis of isolation by distance: relic or reference standard for landscape genetics? Ecography, 33(2), 315–320.
- Kinnison, M. T., and Hairston, N. G. (2007). Eco-evolutionary conservation biology: contemporary evolution and the dynamics of persistence. Functional Ecology 21, 444–454. doi:10.1111/j.1365-2435.2007.01278.x.
- Kuhner, M. K., Beerli, P., Yamato, J., and Felsenstein, J. (2000). Usefulness of single nucleotide polymorphism data for estimating population parameters. Genetics 156, 439–447. Available at: https://www.ncbi.nlm.nih.gov/pubmed/10978306.
- Leaché, A. D., Banbury, B. L., Felsenstein, J., de Oca, A. N.-M., and Stamatakis, A. (2015). Short Tree, Long Tree, Right Tree, Wrong Tree: New Acquisition Bias Corrections for Inferring SNP Phylogenies. Syst. Biol. 64, 1032–1047. doi:10.1093/sysbio/syv053.
- Lee-Yaw, J. A., Irwin, J. T., and Green, D. M. (2008). Postglacial range expansion from northern refugia by the wood frog, Rana sylvatica. Mol. Ecol. 17, 867–884. doi:10.1111/j.1365-294X.2007.03611.x.
- Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete morphological character data. Syst. Biol. 50, 913–925. doi:10.1080/106351501753462876.
- Ligon, N. F., and Skelly, D. K. (2009). Cryptic divergence: countergradient variation in the wood frog. Evol. Ecol. Res. 11, 1099–1109.
- McCartney-Melstad, E., Gidiş, M., & Shaffer, H. B. (2019). An empirical pipeline for choosing the optimal clustering threshold in RADseq studies. Molecular Ecology Resources, 19(5), 1195–1204.
- Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., et al. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. Evol. 37, 1530–1534. doi:10.1093/molbev/msaa015.
- Nadeau, S., Meirmans, P. G., Aitken, S. N., Ritland, K., & Isabel, N. (2016). The challenge of separating signatures of local adaptation from those of isolation by distance and colonization history: The case of two white pines. Ecology and Evolution, 6(24), 8649–8664.
- Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268–274.
- Orsini, L., Vanoverbeke, J., Swillen, I., Mergeay, J., & De Meester, L. (2013). Drivers of population genetic differentiation in the wild: isolation by dispersal limitation, isolation by adaptation and isolation by colonization. Molecular Ecology, 22(24), 5983–5999.
- Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861
- Peterson, B. K., Weber, J. N., Kay, E. H., Fisher, H. S., and Hoekstra, H. E. (2012). Double digest RADseq: an inexpensive method for de novo SNP discovery and genotyping in model and non-model species. PLoS One 7, e37135. doi:10.1371/journal.pone.0037135.
- Räsänen, K., & Hendry, A. P. (2008). Disentangling interactions between adaptive divergence and gene flow when ecology drives diversification. Ecology Letters, 11(6), 624–636.
- Richardson, J. L., Urban, M. C., Bolnick, D. I., and Skelly, D. K. (2014). Microgeographic adaptation and the spatial scale of evolution. Trends Ecol. Evol. 29, 165–176. doi:10.1016/j.tree.2014.01.002.
- Shafer, A., & Wolf, J. B. W. (2013). Widespread evidence for incipient ecological speciation: a meta-analysis of isolation-by-ecology. Ecology Letters, 16(7), 940–950.
- Skelly, D. K. (2004). Microgeographic countergradient variation in the wood frog, Rana sylvatica. Evolution 58, 160–165. doi:10.1111/j.0014-3820.2004.tb01582.x.
- Smith, B. T., Mauck, W. M., Benz, B. W., & Andersen, M. J. (2020). Uneven Missing Data Skew Phylogenomic Relationships within the Lories and Lorikeets. Genome Biology and Evolution, 12(7), 1131–1147.
- Wagner, C. E., Keller, I., Wittwer, S., Selz, O. M., Mwaiko, S., Greuter, L., et al. (2013). Genome-wide RAD sequence data provide unprecedented resolution of species boundaries and relationships in the Lake Victoria cichlid adaptive radiation. Mol. Ecol. 22, 787–798. doi:10.1111/mec.12023.
- Wang, I. J., & Bradburd, G. S. (2014). Isolation by environment. Molecular Ecology, 23(23), 5649–5662.
- Wright, S. (1932). The roles of mutation, inbreeding, crossbreeding, and selection in evolution. in Proceedings of the Sixth International Congress on Genetics, 355–366. Available at: http://www.blackwellpublishing.com/ridley/classictexts/wright.pdf.
- Yoder, J. B., Briskine, R., Mudge, J., Farmer, A., Paape, T., Steele, K., et al. (2013). Phylogenetic signal variation in the genomes of Medicago (Fabaceae). Syst. Biol. 62, 424–438. doi:10.1093/sysbio/syt009.





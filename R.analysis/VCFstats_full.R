library(tidyverse)

var_qual <- read_delim("./data/min013.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
  # This doesn't work because ipyrad changes all of the quality scores to 13.

var_depth <- read_delim("./data/min013.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
summary(var_depth$mean_depth) # This is obviously also not correct because I specified a minimum depth of 6 for a locus to be retained. So, this must be pulling the depth output from ipyrad which is including all variants.

var_miss <- read_delim("./data/min013.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
  # Many sites have lots of missing data, but many also are represented across all individuals. I specified in ipyrad that we would only retain loci present in 10% of the population, or 13 individuals, so there should be no variants with less than 10% missingness.
summary(var_miss$fmiss)
  # We have about 65% missingness on average.

var_freq <- read_delim("./data/min013.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
summary(var_freq$maf)
  # Most variants have low frequencies


# Mean depth per individual
ind_depth <- read_delim("./data/min013.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

# Proportion of missing data per individual
ind_miss  <- read_delim("./data/min013.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()


# Heterozygosity and inbreeding coefficient per individual
ind_het <- read_delim("./data/min013.het", delim = "\t",
           col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
  # It seems like we have an excess of heterozygosity overall



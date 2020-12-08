# Redoing the analysis with PCA computed after removing sample 105

library(tidyverse)
library(cowplot)
library(ggalt)

### Read in the data from PLINK

PCs <- read_delim("./data/processed/YMFmain.subset.l50i50.eigenvec", delim = "\t", col_names = TRUE, trim_ws = TRUE) %>%
  select(-'#FID') %>%
  rename(ID = IID)
eigenval <- scan("./data/processed/YMFmain.subset.l50i50.eigenval")
sample_ID <- read.csv("./data/processed/YMF2018_RASYpopgen_sampleinfo.csv") %>% select(-X)

# Convert to percentage variance explained by the PCs
pve <- data.frame(PC = 1:10, eig = eigenval)

pve %>%
  ggplot(aes(x = PC, y = eig)) + 
  geom_bar(stat = "identity") +
  labs(y = "Eigenvalues") +
  scale_x_continuous(breaks=seq(1, 10, 1))

# Read in the VCF stats
readVCFstats <- function(LOCI, INDV){
# Mean read depth among individuals
  # Note that this might not be correct since ipyrad already did some filtering
assign(paste0("var_idepth_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".idepth"), delim = "\t", 
                                                           col_names = c("chr", "pos", "qual"), skip = 1),
       envir = parent.frame())

# Mean depth of coverage per locus (essentially, number of reads mapped to each position)
assign(paste0("var_ldepth_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".ldepth.mean"), delim = "\t", 
                                                           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1),
       envir = parent.frame())

# Mean depth of coverage per individual
assign(paste0("var_idepth_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".idepth"), delim = "\t", 
                                                           col_names = c("ind", "nsites", "depth"), skip = 1),
       envir = parent.frame())

# Missingness by locus (how many individuals lack a genotype call at this site)
assign(paste0("var_lmiss_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".lmiss"), delim = "\t", 
                                                          col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1),
       envir = parent.frame())

# Missingness by individual
assign(paste0("var_imiss_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".imiss"), delim = "\t", 
                                                          col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1),
       envir = parent.frame())

# Minor allele frequency
assign(paste0("var_maf_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".frq"), delim = "\t", 
                                                         col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1),
       envir = parent.frame())

# Minor allele frequency
assign(paste0("var_het_l", LOCI, "i", INDV), read_delim(paste0("./data/processed/YMF_l", LOCI, "i", INDV, ".het"), delim = "\t", 
                                                        col_names = c("ind","ho", "he", "nsites", "f"), skip = 1),
       envir = parent.frame())
}
readVCFstats("50", "50")

# Combine all of the VCF statistics
VCFstats <- var_idepth_l50i50 %>% rename(nsites.depth = nsites) %>% 
  left_join(var_imiss_l50i50, by = "ind") %>%
  left_join(var_het_l50i50 %>% rename(nsites.het = nsites), by = "ind") %>%
  rename(ID = ind) %>%
  separate(ID, c("tmpPop", "tmpID", "Run1", "Run2")) %>% 
  mutate(Pop = ifelse(tmpPop == "YNL", tmpID, tmpPop), ID = ifelse(tmpPop == "YNL", Run1, tmpID)) %>%
  mutate(ID = as.integer(ID)) %>%
  select(-tmpPop, -tmpID, -Run1, -Run2)

### Sort out the population IDs
## Combine the dataset

# PCdat <- PCs %>% mutate(FileName = ID) %>%
#   separate(ID, c("tmpPop", "tmpID", "Run1", "Run2")) %>% 
#   mutate(Pop = ifelse(tmpPop == "YNL", tmpID, tmpPop), ID = ifelse(tmpPop == "YNL", Run1, tmpID)) %>%
#   mutate(ID = as.integer(ID)) %>%
#   left_join(sample_ID, by = c("ID" = "Extract")) %>%
#   select(-Pop.x, -tmpPop, -tmpID)


PCdat <- PCs %>% mutate(FileName = ID) %>%
  separate(ID, c("tmpPop", "tmpID", "Run1", "Run2")) %>% 
  mutate(Pop = ifelse(tmpPop == "YNL", tmpID, tmpPop), ID = ifelse(tmpPop == "YNL", Run1, tmpID)) %>%
  mutate(ID = as.integer(ID)) %>%
  select(-tmpPop, -tmpID, -Run2) %>%
  left_join(VCFstats, by = c("Pop", "ID")) %>%
  left_join(sample_ID, by = c("ID" = "Extract"))



PC.1.2 <- PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = Pond, label = Pond)) +
  geom_text() +
  theme_minimal() +
  labs(title = "SNPs")

PC.1.3 <- PCdat %>% 
  ggplot(aes(x = PC1, y = PC3, col = Pond, label = Pond)) +
  geom_text() +
  theme_minimal() +
  labs(title = "SNPs")

PC.2.3 <- PCdat %>% 
  ggplot(aes(x = PC2, y = PC3, col = Pond, label = Pond)) +
  geom_text() +
  theme_minimal() +
  labs(title = "SNPs")

Geo <- PCdat %>% 
  ggplot(aes(x = Long, y = Lat, col = Pond, label = Pond)) +
  geom_text() +
  theme_minimal() +
  labs(title = "Geography")

plot_grid(PC.1.2 + theme(legend.position = "none"),
          PC.1.3 + theme(legend.position = "none"),
          PC.2.3 + theme(legend.position = "none"),
             Geo + theme(legend.position = "none"))

# Restricting to only sites with 10 or more samples. 
PCdat %>% 
  inner_join(PCdat %>% group_by(Pond) %>% tally() %>% arrange(desc(n)) %>% filter(n > 9), by = "Pond") %>%
  ggplot(aes(x = PC1, y = PC2, col = Pond, label = Pond)) +
  geom_text() +
  stat_ellipse(aes(x = PC1, y = PC2, col = Pond)) +
  theme_minimal() +
  labs(title = "SNPs")


# Checking correlation with other variables
PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = sample_coverage)) +
  geom_text(aes(, label = Pond), cex = 4) +
  theme_minimal() +
  scale_color_gradient(high = 'red', low = 'blue') +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = Age)) +
  geom_text(aes(, label = Pond), cex = 4) +
  theme_minimal() +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = fmiss)) +
  geom_text(aes(, label = Pond), cex = 4) +
  theme_minimal() +
  scale_color_gradient(high = 'red', low = 'blue') +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = nsites.het)) +
  geom_text(aes(, label = Pond), cex = 4) +
  theme_minimal() +
  scale_color_gradient(high = 'red', low = 'blue') +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = he)) +
  geom_text(aes(label = Pond), cex = 4) +
  theme_minimal() +
  scale_color_gradient(high = 'red', low = 'blue') +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = f)) +
  geom_text(aes(label = Pond), cex = 4) +
  theme_minimal() +
  scale_color_gradient(high = 'red', low = 'blue') +
  labs(title = "SNPs")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = Conc.Qubit.ng_ul)) +
  geom_text(aes(, label = Pond), cex = 4) +
  scale_color_gradient(high = 'red', low = 'blue') +
  theme_minimal() +
  labs(title = "SNPs", col = "Het")

PCdat %>% 
  ggplot(aes(x = PC1, y = PC2, col = Extract.Date)) +
  geom_text(aes(, label = Pond), cex = 4) +
  theme_minimal() +
  labs(title = "SNPs")



PCdat %>%
  ggplot(aes(x = fmiss, y = f)) +
  geom_jitter() +
  labs(y = "Heterozygosity", x = "% sites missing")

PCdat %>%
  ggplot(aes(x = Conc.Qubit.ng_ul, y = fmiss)) +
  geom_jitter() +
  labs(y = "% sites missing", x = "DNA extract quant.")

PCdat %>%
  ggplot(aes(x = sample_coverage, y = fmiss)) +
  geom_jitter() +
  labs(y = "% sites missing", x = "Sample coverage")

PCdat %>%
  ggplot(aes(x = Conc.Qubit.ng_ul, y = sample_coverage)) +
  geom_jitter() +
  labs(y = "Sample coverage", x = "DNA extract quant.")

PCdat %>%
  ggplot(aes(x = Age, y = fmiss)) +
  geom_jitter() +
  labs(y = "% sites missing", x = "Tissue type")

PCdat %>%
  ggplot(aes(x = Extract.Date, y = fmiss, col = Age)) +
  geom_jitter() +
  labs(y = "% sites missing", x = "Extraction Date")

### Figuring out cluster members

PCdat %>% 
  mutate(Group = ifelse(PC1 <= 0 & PC2 >= -0.02, 1, 0)) %>%
  mutate(Group = ifelse(PC1 >= 0 & PC2 <= -0.02, 2, Group)) %>%
  ggplot(aes(x = PC1, y = PC2, label = Pond, col = Group)) +
  geom_point() +
  theme_minimal() +
  labs(title = "SNPs") +
  theme(legend.position = "none")

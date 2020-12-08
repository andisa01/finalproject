
library(tidyverse)

PCs <- read_delim("./data/raw/YMFmain-RCref.subset.l01i99.eigenvec.var", delim = "\t", col_names = TRUE, trim_ws = TRUE) %>%
  rename(CHROM = '#CHROM')

PCs %>% view()


PCs %>% arrange(desc(abs(PC1)))

length(PCs$CHROM)
length(unique(PCs$CHROM))

PCs %>% group_by(CHROM) %>% tally() %>% arrange(desc(n))
# At maximum we have 25 SNPs per chromosome

PCs %>% group_by(CHROM) %>% tally() %>% 
  ggplot(aes(x = n)) +
    geom_histogram(binwidth = 1)

PCs %>% left_join(PCs %>% group_by(CHROM) %>% tally(), by = c("CHROM")) %>%
  #filter(n >= 4) %>%
  mutate(CHROM = fct_reorder(CHROM, desc(n))) %>%
  ggplot(aes(x = CHROM, y = abs(PC1), col = n)) +
    geom_point() +
    scale_color_gradient(high = "firebrick", low = "grey50") +
    labs(y = "PC1", color = "per chromosome") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
#ggsave("figs/VariantPC1loadings.pdf", width = 9, height = 5, dpi = 300, device = cairo_pdf)

PCs %>% left_join(PCs %>% group_by(CHROM) %>% tally(), by = c("CHROM")) %>%
  #filter(n >= 4) %>%
  mutate(CHROM = fct_reorder(CHROM, desc(n))) %>%
  ggplot(aes(x = CHROM, y = abs(PC2), col = n)) +
    geom_point() +
    scale_color_gradient(high = "firebrick", low = "grey50") +
    labs(y = "PC2", color = "per chromosome") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
#ggsave("figs/VariantPC2loadings.pdf", width = 9, height = 5, dpi = 300, device = cairo_pdf)

PCs_byTally <- PCs %>% left_join(PCs %>% group_by(CHROM) %>% tally(), by = c("CHROM")) %>%
  group_by(n) %>%
  summarise(medianPC1 = median(abs(PC1)), meanPC1 = mean(abs(PC1)), medianPC2 = median(abs(PC2)), meanPC2 = mean(abs(PC2)))

PCs_byTally %>%
  ggplot(aes(x = n, y = meanPC1)) +
  geom_point() +
  geom_smooth(col = 1)

PCs_byTally %>%
  ggplot(aes(x = n, y = meanPC2)) +
  geom_point() +
  geom_smooth(col = 1)

PCs_byTally %>%
  ggplot(aes(x = n, y = medianPC1)) +
  geom_point() +
  geom_smooth(col = 1)

PCs_byTally %>%
  ggplot(aes(x = n, y = medianPC2)) +
  geom_point() +
  geom_smooth(col = 1)



PCs %>% mutate(jointLoad = sqrt(abs(PC1)^2 + abs(PC2)^2)) %>%
  ggplot(aes(x = jointLoad)) +
    geom_density()

PCs %>% mutate(jointLoad = sqrt(abs(PC1)^2 + abs(PC2)^2)) %>%
  arrange(desc(jointLoad)) %>%
  print(n = 40)

PCs %>% left_join(PCs %>% group_by(CHROM) %>% tally(), by = c("CHROM")) %>%
  mutate(jointLoad = sqrt(abs(PC1)^2 + abs(PC2)^2)) %>%
  filter(n >= 4) %>%
  mutate(CHROM = fct_reorder(CHROM, desc(n))) %>%
  ggplot(aes(x = CHROM, y = jointLoad, col = n)) +
    geom_point() +
    scale_color_gradient(high = "firebrick", low = "grey50") +
    labs(y = "jointLoad", color = "per chromosome") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

PCs %>% mutate(jointLoad = sqrt(abs(PC1)^2 + abs(PC2)^2)) %>%
  filter(jointLoad >= 2) %>% View()


### PC by loc, not CHROM

PCs <- read_delim("./data/processed/YMFmain.l50i50.eigenvec.var", delim = "\t", col_names = TRUE, trim_ws = TRUE) %>%
  rename(chr = '#CHROM')

PCs %>% arrange(desc(abs(PC1)))

# Read in VCF stats ====
# Mean depth of coverage per locus (essentially, number of reads mapped to each position)
var_ldepth_l50i50

var_lmiss_l50i50

VarDat <- PCs %>%
  mutate(jointLoad = sqrt(abs(PC1)^2 + abs(PC2)^2)) %>%
  select(chr:PC2) %>%
  inner_join(var_ldepth_l50i50, by = c("chr")) %>%
  inner_join(var_lmiss_l50i50, by = c("chr"))

VarDat %>%
  ggplot(aes(x = mean_depth, y = PC1)) +
  geom_point()

VarDat %>%
  ggplot(aes(x = mean_depth, y = PC2)) +
  geom_point()

VarDat %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

PCdep <- VarDat %>% arrange(mean_depth) %>%
  ggplot(aes(x = PC1, y = PC2, col = mean_depth)) +
  geom_point(pch = 16, size = 2) +
  scale_color_gradient(high = "firebrick", low = "grey80") +
  labs(title = "Mean Locus Depth") +
  theme(legend.position = "bottom")

PCdpevar <- VarDat %>% arrange(var_depth) %>%
  ggplot(aes(x = PC1, y = PC2, col = var_depth)) +
  geom_point(pch = 16, size = 2) +
  scale_color_gradient(high = "firebrick", low = "grey80") +
  labs(title = "Variance in Locus Depth") +
  theme(legend.position = "bottom")

PCmis <- VarDat %>% arrange(fmiss) %>%
  ggplot(aes(x = PC1, y = PC2, col = fmiss)) +
  geom_point(pch = 16, size = 2) +
  scale_color_gradient(high = "firebrick", low = "grey80") +
  labs(title = "Locus missingness (freq.)") +
  theme(legend.position = "bottom")

plot_grid(PCdep, PCdpevar, PCmis, ncol = 3)
#ggsave("figs/VariantPCloadings.pdf", width = 10, height = 4.25, dpi = 300, device = cairo_pdf)

VarDat %>% arrange(abs(PC1)) %>%
  ggplot(aes(x = log(var_depth), y = fmiss, col = abs(PC1))) +
  geom_point(pch = 16, size = 2) +
  scale_color_gradient(high = "firebrick", low = "grey80") +
  labs(title = "Variance in Locus Depth") +
  theme(legend.position = "bottom")

VarDat %>% arrange(abs(PC2)) %>%
  ggplot(aes(x = log(var_depth), y = fmiss, col = abs(PC2))) +
  geom_point(pch = 16, size = 2) +
  scale_color_gradient(high = "firebrick", low = "grey80") +
  labs(title = "Variance in Locus Depth") +
  theme(legend.position = "bottom")

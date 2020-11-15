
library(ape)
library(tidyverse)
library(ggtree)

tree_YMFpilot1 <- read.tree("./data/YMFpilot.0.01_0.90_100.min1.phy.treefile")

tree_YMFpilot1$tip.label <- c("BS", "MI", "SH", "W1", "KE", "E8", "WF", "BPS", "GB", "LO", "BO", "FHE", "LA", "WM", "PB", "WP")

plot(tree_YMFpilot1)



tree_YMFmain <- read.tree("./2.subsampling/YMFmain/YMFmain.subset.0.50_0.50_1000.min1.phy.treefile")

popmap_YMFmain <- read_table2("1.filtering/YMFmain/GBS_SNP_filter/popmap.txt", col_names = c("Sample", "Pop")) %>%
  separate(Pop, into = c("Locality", "Population")) %>%
  separate(Sample, into = c("LocalityYear", "Num", "Lib")) %>%
  filter(LocalityYear == "YMF2018") %>%
  select(Num, Population) %>%
  group_by(Num) %>%
  summarise(Population = first(Population))

YMFmain_labels <- as.data.frame(tree_YMFmain$tip.label) %>% rename(Sample = "tree_YMFmain$tip.label") %>% separate(Sample, into = c("Sample", "Num", "Lib", "ext")) %>%
  select(Num) %>% 
  left_join(popmap_YMFmain, by = "Num") %>%
  rename(tip.label = Num) %>%
  select(tip.label, Population)

tree_YMFmain$tip.label <- YMFmain_labels$tip.label

plot(tree_YMFmain)

class(tree_YMFmain)

p <- ggtree(tree_YMFmain, layout = 'circular', branch.length = 'none') %<+% YMFmain_labels

p + geom_tiplab(aes(label = Population, col = Population), size = 3) +
  theme(legend.position = 'none')
ggsave("./figs/YMFmain_phylo.png", width = 8, height = 8)

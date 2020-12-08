library(tidyverse)

X <- c()
for (i in list.files("./data/processed/", pattern = "YMFcomb_*")) {
  DF <- read_table2(paste0("data/processed/",i)) %>% 
    mutate(fLOCI = as.numeric(str_extract(i, "\\d+\\.*\\d*")))
  X <- bind_rows(X, DF)
  }

X %>% group_by(fLOCI) %>%
  summarise(miss = mean(F_MISS)) %>%
  ggplot(aes(x = fLOCI, y = miss)) +
  geom_line()

X2 <- c()
for (i in unique(X$fLOCI)) {
  DF <- X %>% group_by(fLOCI) %>%
    filter(F_MISS <= i) %>%
    tally() %>%
    mutate(fINDV = i)
  X2 <- bind_rows(X2, DF)
}

# X2 %>% pivot_wider(names_from = fINDV, values_from = n) %>%
#   replace(is.na(.), 0) %>%
#   arrange(fLOCI)

X3 <- c()
for (i in unique(X$fLOCI)) {
  DF <- X %>% group_by(fLOCI) %>%
    mutate(N_LOCI = N_DATA - N_MISS) %>%
    filter(F_MISS <= i) %>%
    summarise(L_DEPTH = mean(N_LOCI)) %>%
    mutate(fINDV = i)
  X3 <- bind_rows(X3, DF)
}

# X3 %>% pivot_wider(names_from = fINDV, values_from = L_DEPTH) %>%
#   replace(is.na(.), 0) %>%
#   arrange(desc(fLOCI))

Xall <- as_tibble(expand.grid(fLOCI = unique(X$fLOCI), fINDV = unique(X$fLOCI))) %>%
  left_join(X2, by = c("fLOCI", "fINDV")) %>%
  left_join(X3, by = c("fLOCI", "fINDV")) %>%
  mutate(SNPs = n * L_DEPTH) %>%
  replace(is.na(.), 0)

Xall <- Xall %>% mutate(zn = (n - mean(Xall$n))/sd(n), 
                zL_DEPTH = (L_DEPTH - mean(Xall$L_DEPTH))/sd(L_DEPTH),
                zSNPs = (SNPs - mean(Xall$SNPs))/sd(Xall$SNPs))

# Rough-grain figures

Xall %>% 
  filter(fLOCI %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1)) & fINDV %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1))) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = heat.colors(2)[1], high = heat.colors(2)[2]) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()

Xall %>% 
  filter(fLOCI %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1)) & fINDV %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1))) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = L_DEPTH)) +
    geom_tile() +
    geom_text(aes(label = round(L_DEPTH/1000, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = heat.colors(2)[1], high = heat.colors(2)[2]) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Loci Depth") +
    theme_minimal()
    
Xall %>%
  filter(fLOCI %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1)) & fINDV %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1))) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = SNPs)) +
    geom_tile() +
    geom_text(aes(label = round(SNPs/1000000, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = heat.colors(2)[1], high = heat.colors(2)[2]) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "SNPs") +
    theme_minimal()

Xall %>%
  filter(fLOCI %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1)) & fINDV %in% c(0.01, 0.99, seq(0.1, 0.9, by = 0.1))) %>%
  mutate(Z = (zSNPs + zL_DEPTH + zn*2)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = Z)) +
    geom_tile() +
    geom_text(aes(label = round(Z, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Z", subtitle = "Z score = standardized: number of samples*2 + mean loci freq across samples + total SNPs") +
    theme_minimal()

# Fine-grain figures

Xall %>% 
  filter(fLOCI %in% seq(0.34, 0.76, by = 0.02) & fINDV %in% seq(0.34, 0.76, by = 0.02)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 230) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()

Xall %>% 
  filter(fLOCI %in% seq(0.34, 0.76, by = 0.02) & fINDV %in% seq(0.34, 0.76, by = 0.02)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = L_DEPTH)) +
    geom_tile() +
    geom_text(aes(label = round(L_DEPTH/1000, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = mean(Xall$L_DEPTH)) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()
    
Xall %>%
  filter(fLOCI %in% seq(0.34, 0.76, by = 0.02) & fINDV %in% seq(0.34, 0.76, by = 0.02)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = SNPs)) +
    geom_tile() +
    geom_text(aes(label = round(SNPs/1000000, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = mean(Xall$SNPs)) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()

Xall %>%
  filter(fLOCI %in% seq(0.34, 0.76, by = 0.02) & fINDV %in% seq(0.34, 0.76, by = 0.02)) %>%
  mutate(Z = (zSNPs + zL_DEPTH + zn*2)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = Z)) +
    geom_tile() +
    geom_text(aes(label = round(Z, 2))) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Z", subtitle = "Z score = standardized: number of samples*2 + mean loci freq across samples + total SNPs") +
    theme_minimal()

### All

# Fine-grain figures

Xall %>% 
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), size = 2) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 230) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()

Xall %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = L_DEPTH)) +
    geom_tile() +
    geom_text(aes(label = round(L_DEPTH/1000, 2)), size = 2) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = mean(Xall$L_DEPTH)) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()
    
Xall %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = SNPs)) +
    geom_tile() +
    geom_text(aes(label = round(SNPs/1000000, 2)), size = 2) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = mean(Xall$SNPs)) +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Samples") +
    theme_minimal()

Xall %>%
  mutate(Z = (zSNPs + zL_DEPTH + zn)) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = Z)) +
    geom_tile() +
    #geom_text(aes(label = round(Z, 2)), size = 2) +
    geom_text(aes(label = n, col = ifelse(n >= 250, "black", "grey70")), size = 3) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0) +
    scale_color_identity() +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Z", subtitle = "Z score = standardized: number of samples*2 + mean loci freq across samples + total SNPs") +
    theme_minimal()

Xall %>%
  mutate(Z = (SNPs/(max(SNPs)-min(SNPs))) + (L_DEPTH/(max(L_DEPTH)-min(L_DEPTH))) + (n/(max(n)-min(n))) / 3) %>%
  mutate(penalty = ((1-fLOCI)/0.98 + fINDV/0.98)/2) %>%
  mutate(score = Z - penalty) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = score)) +
    geom_tile() +
    #geom_text(aes(label = round(Z, 2)), size = 2) +
    geom_text(aes(label = n, col = ifelse(score >= 0.5 & n >= 250, "black", "grey70")), size = 3) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0) +
    scale_color_identity() +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Z", subtitle = "Z score = standardized: number of samples*2 + mean loci freq across samples + total SNPs") +
    theme_minimal()

Xall %>%
  mutate(Z = (SNPs/(max(SNPs)-min(SNPs))) + (n/(max(n)-min(n))) / 2) %>%
  mutate(penalty = ((1-fLOCI)/0.98 + fINDV/0.98)/2) %>%
  mutate(score = Z - penalty) %>%
  ggplot(aes(x = as.factor(fLOCI), y = as.factor(fINDV), fill = score)) +
    geom_tile() +
    #geom_text(aes(label = round(Z, 2)), size = 2) +
    geom_text(aes(label = n, col = ifelse(score >= 0.5 & n >= 260, "black", "grey70")), size = 3) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0) +
    scale_color_identity() +
    labs(x = "Loci (Prop. samples with locus)", y = "Sample missingness threshold", fill = "Z", subtitle = "Z score = standardized (#SNPs + #Indv) - ((1-Loci threshold) + (Sample threshold))") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


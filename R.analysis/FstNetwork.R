
library(ggraph)
library(tidyverse)
library(graphlayouts)
library(igraph)
library(cowplot)

### Useful tutorials:
# Intro to networks: https://www.jessesadler.com/post/network-analysis-with-r/
# GGraph intro: https://www.data-imaginist.com/2017/ggraph-introduction-layouts/
# GGraph edges: https://www.data-imaginist.com/2017/ggraph-introduction-edges/
# Tidygraph intro: https://www.data-imaginist.com/2017/introducing-tidygraph/
# Tidygraph update: https://www.data-imaginist.com/2018/tidygraph-1-1-a-tidy-hope/
# Tidygraph example: https://rviews.rstudio.com/2019/03/06/intro-to-graph-analysis/
# Network visualizations: http://mr.schochastics.net/netVizR.html


### Read in the data ====
FstEdges <- read_table2("data/processed/FstEdges.txt") # These are the Fst values estimated from a VCF of the YMF ponds filtered to loci represented in at least 50% of samples and samples with at least 50% coverage. These estimates are pond-wise average Fst.

PondID <- read.csv("C:\\Users\\Andis\\Google Drive\\YMFdata\\YMFpondID.csv", header = TRUE)

OutgroupCoords <- data.frame(Pond = c("Concord", "Colonial", "BethBog", "ScranProp", "YP05", "YP16", "NAT", "CPN", "CPS", "RichProp"), 
           Lat = c(41.344701, 41.326933, 41.4261, 41.271603, 41.319019, 41.320662, 41.9677148, 41.949582, 41.949156,  41.265173), 
           Long = c(-72.614925, -72.633661, -72.9928, -72.730821, -72.990415, -72.992365, -72.1320442, -72.120016, -72.119785, -72.735702)) # This is actually not necessary for this dataset.

# The oviposition data need a bit of cleaning because CPN and CPS have been considered the same pond in some years during surveys and separate in others. For the egg mass count data, I just collapsed all the counts into a single pond "CP", but for the network graph, I need to split them bach into two ponds. I've done this by just assuming that half the population came from CPN and half from CPS.
OviData_CPS <- read.csv("data/raw/YMF_EggSurveyData_00to19_CLEANED.csv") %>% select(-X) %>%
  filter(Pond == "CP") %>%
  mutate(Pond = ifelse(Pond == "CP", "CPS", as.character(Pond))) # This creates a dataset for CPS data which we will bind to the main dataset.

OviData <- read.csv("data/raw/YMF_EggSurveyData_00to19_CLEANED.csv") %>% select(-X) %>%
  mutate(Pond = ifelse(Pond == "LP", "LO", as.character(Pond))) %>%
  mutate(Pond = ifelse(Pond == "CP", "CPN", as.character(Pond))) %>%
  bind_rows(OviData_CPS) %>% # Bind the CP duplicates
  mutate(Count = ifelse(Pond == "CPS" | Pond == "CPN", Count/2, Count)) %>% # Halve their population size
  select(-DOY) %>%
  filter(!is.na(DOYobs))
#write.csv(OviData, "./OvipositonData.csv")

Ovi_allyrs <- OviData %>% group_by(Pond) %>%
  summarise(Pop_allyrs = mean(Count)) # Estimate population size across all years

Ovi_2018 <- OviData %>% filter(Year == 2018) %>%
  group_by(Pond) %>%
  summarise(Pop_2018 = mean(Count)) # Estimate populations size in 2018 when we collected the data

Ovi_extprob <- OviData %>% group_by(Pond) %>%
  summarise(Obs = n(), Extinct = sum(Count == 0), ExtProb = Extinct/Obs) %>%
  select(-Obs, -Extinct) # Estimate the extinction probability for each pond.

Ovi_allyrs_harmonic <- OviData %>% mutate(Count2 = 1/Count) %>%
  group_by(Pond) %>%
  summarise(Pop_harmonic = 1/mean(Count2))

PondInfo <- read.csv("C:/Users/Andis/Google Drive/YMFdata/YM_Pond_Coords.csv", header = TRUE, na.strings = "") %>% mutate(AKA = as.character(Name)) %>%
  left_join(PondID, by = "AKA") %>% 
  select(Pond, Lat, Long) %>%
  mutate(Pond = ifelse(Pond == "LP", "LO", as.character(Pond))) %>% # Fix old ID code for Long Pond to new code "LO"
  bind_rows(OutgroupCoords) %>%
  mutate(name = ifelse(!grepl("YP|Concord|Colonial|BethBog|Prop", Pond), paste0("YMF_", Pond), Pond)) %>% # Reformat pond names as in the genetic dataset
  group_by(Pond) %>%
  summarise_all(list(first)) %>% # For some reason NAT is duplicated, so this removes all duplicates
  ungroup() %>%
  left_join(Ovi_allyrs, by = "Pond") %>%
  left_join(Ovi_2018, by = "Pond") %>%
  left_join(Ovi_extprob, by = "Pond") %>%
  left_join(Ovi_allyrs_harmonic, by = "Pond")
  
anti_join(FstEdges, PondInfo, by = c("from" = "name")) # Check to make sure all node names are accounted for
anti_join(FstEdges, PondInfo, by = c("to" = "name")) # Check to make sure all node names are accounted for

### Make network figures ====

# First, need to make a graph object and scale the variables
FstGraph <- as_tbl_graph(FstEdges) %>% # Make the Fst data into a graph object
  activate(nodes) %>%
    left_join(PondInfo) %>% # Add the pond-wise variables to the nodes
    group_by() %>%
    mutate(sPop_allyrs = Pop_allyrs/max(Pop_allyrs),
           sPop_2018 = Pop_2018/max(Pop_2018),
           sPop_harmonic = Pop_harmonic/max(Pop_harmonic)) %>% # Scale the population sizes
  activate(edges) %>%
    mutate(Fst = weight) %>%
    mutate(weight = log(Fst)/max(log(Fst))) # Estimate a log scaled Fst (basically, as Fst decreases i.e. more connectivity, this value increases exponentially. We will use this as the alpha scale for the edges

FstEdges %>% mutate(Fst = weight) %>%
  mutate(weight = log(Fst)/max(log(Fst))) %>%
  ggplot(aes(x = Fst, y = weight)) +
  geom_point() # This graph just shows how the alpha for the edge weighing compares to Fst. This makes the edges with lower Fst exponentially darker.

# Estimate the percentiles for Fst values
FstGraph %>%
  activate(edges) %>%
  as_tibble() %>%
  group_by() %>%
  summarise(Fst = quantile(Fst, c(0.1, 0.25, 0.5, 0.75, 0.9, 1)), q = c(0.1, 0.25, 0.5, 0.75, 0.9, 1))

### Plots ====
## . Geographic networks ====
## . . Geographic network plots by Fst percentile ====
Fst_plot20 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.2)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_link0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 20th percentile") + 
  coord_fixed(ratio = 1)
  
Fst_plot10 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_link0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 10th percentile") + 
  coord_fixed(ratio = 1)

Fst_plot05 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.05)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_link0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 5th percentile") + 
  coord_fixed(ratio = 1)

Fst_plot02 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.025)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_link0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 2.5nd percentile") + 
  coord_fixed(ratio = 1)

Geo_net <- plot_grid(Fst_plot02 + theme(legend.position = "none"),
                     Fst_plot05 + theme(legend.position = "none"), 
                     Fst_plot10 + theme(legend.position = "none"),
                     Fst_plot20 + theme(legend.position = "none"),
                     nrow = 1,
                     rel_widths = c(1,1,1,1))
Geo_net
#ggsave("figs/FstNetwork_line.pdf", width = 14, height = 6, dpi = 300, device = cairo_pdf)

## . . (diagonal) Geographic network plots by Fst percentile ====
Fst_plot20 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.2)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 20th percentile") + 
  coord_fixed(ratio = 1)
  
Fst_plot10 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 10st percentile") + 
  coord_fixed(ratio = 1)

Fst_plot05 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.05)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 5th percentile") + 
  coord_fixed(ratio = 1)

Fst_plot02 <- FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.025)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Geographic", subtitle = "Fst 2.5nd percentile") + 
  coord_fixed(ratio = 1)

Geo_net <- plot_grid(Fst_plot02 + theme(legend.position = "none"),
                     Fst_plot05 + theme(legend.position = "none"), 
                     Fst_plot10 + theme(legend.position = "none"),
                     Fst_plot20 + theme(legend.position = "none"),
                     nrow = 1,
                     rel_widths = c(1,1,1,1))
Geo_net
#ggsave("figs/FstNetwork_diag.pdf", width = 14, height = 6, dpi = 300, device = cairo_pdf)


# . . Geographic based on centrality ====
FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.05)) %>%
  ggraph(x = Long, y = Lat) +
  geom_edge_diagonal(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = centrality, fill = centrality), pch = 21, col = "black") +
  scale_size(range = c(3, 10)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Fst 5th percentile") +
  coord_fixed(ratio = 1)

# . Chord diagram ====
FstChord02 <- FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.025)) %>%
  as.igraph() %>%
  create_layout(layout = 'linear', circular = TRUE) %>%
  ggraph() +
  geom_edge_arc(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Chord diagram", subtitle = "Fst 2.5nd percentile") +
  coord_fixed()

FstChord05 <- FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.05)) %>%
  as.igraph() %>%
  create_layout(layout = 'linear', circular = TRUE) %>%
  ggraph() +
  geom_edge_arc(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Chord diagram", subtitle = "Fst 5th percentile") +
  coord_fixed()

FstChord10 <- FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  as.igraph() %>%
  create_layout(layout = 'linear', circular = TRUE) %>%
  ggraph() +
  geom_edge_arc(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Chord diagram", subtitle = "Fst 10th percentile") +
  coord_fixed()

FstChord20 <- FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.20)) %>%
  as.igraph() %>%
  create_layout(layout = 'linear', circular = TRUE) %>%
  ggraph() +
  geom_edge_arc(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Chord diagram", subtitle = "Fst 20th percentile") +
  coord_fixed()
Chord_net <- plot_grid(FstChord02 + theme(legend.position = "none"),
                    FstChord05 + theme(legend.position = "none"), 
                     FstChord10 + theme(legend.position = "none"),
                     FstChord20 + theme(legend.position = "none"),
                     nrow = 1,
                     rel_widths = c(1,1,1,1))
Chord_net
#ggsave("figs/FstNetwork_chord.pdf", width = 14, height = 6, dpi = 300, device = cairo_pdf)

plot_grid(Geo_net, Chord_net, cols = 1)
#ggsave("figs/FstNetwork_diag_chord.pdf", width = 14, height = 10, dpi = 300, device = cairo_pdf)

# . MDS diagram ====
FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  as.igraph() %>%
  create_layout(layout = 'mds', circular = FALSE) %>%
  ggraph() +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Multi-Dim Scaling", subtitle = "Fst 10th percentile")

# . Fruchterman-Reingold layout ====
FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  as.igraph() %>%
  create_layout(layout = 'fr', circular = FALSE) %>%
  ggraph() +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Fruchterman-Reingold layout", subtitle = "Fst 10th percentile")

# . Stress layout ====
FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.50)) %>%
  as.igraph() %>%
  create_layout(layout = 'stress', circular = FALSE) %>%
  ggraph() +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Stress assorted layout", subtitle = "Fst median")
#ggsave("figs/FstNetwork_stress50.pdf", width = 18, height = 6, dpi = 300, device = cairo_pdf)

# . Stress arangement layout
FstGraph %>%
  activate(nodes) %>%
  mutate(centrality = centrality_authority()) %>%
  activate(edges) %>% 
  filter(Fst <= quantile(FstEdges$weight, 0.20)) %>%
  as.igraph() %>%
  create_layout(layout = 'stress', circular = FALSE) %>%
  ggraph() +
  geom_edge_diagonal0(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_harmonic, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Stress assorted layout", subtitle = "Fst 20th percentile")
#ggsave("figs/FstNetwork_stress20_harmonic.pdf", width = 18, height = 6, dpi = 300, device = cairo_pdf)

# . Kamada-Kawai layout ====
FstGraph %>%
  activate(edges) %>%
  filter(Fst <= quantile(FstEdges$weight, 0.10)) %>%
  ggraph() +
  geom_edge_diagonal(aes(alpha = weight), width = 1, colour = 'black') +
  geom_node_point(aes(size = sPop_allyrs, fill = ExtProb), pch = 21, col = "black") +
  scale_size(range = c(2, 15)) +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = median(PondInfo$ExtProb, na.rm = TRUE)) +
  geom_node_text(aes(label = Pond), size = 3) +
  theme_graph() +
  labs(title = "Kamada-Kawai layout", subtitle = "Fst 10th percentile")

OviData %>% 
  left_join(OviData %>% group_by(Year) %>% summarise(meanDOY = mean(DOYobs)) %>% select(Year, meanDOY), by = c("Year")) %>%
  mutate(lag = DOYobs - meanDOY) %>%
  filter(Pond %in% c("DT", "LO", "SD")) %>%
  ggplot(aes(x = Pond, y = lag)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_boxplot(width = 0.1) +
    theme_minimal()
  

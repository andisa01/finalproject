
library(ape)

tree_YMFpilot1 <- read.tree("./data/YMFpilot.0.01_0.90_100.min1.phy.treefile")

tree_YMFpilot1$tip.label <- c("BS", "MI", "SH", "W1", "KE", "E8", "WF", "BPS", "GB", "LO", "BO", "FHE", "LA", "WM", "PB", "WP")

plot(tree_YMFpilot1)


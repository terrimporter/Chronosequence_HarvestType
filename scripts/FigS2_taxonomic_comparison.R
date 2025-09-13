# Teresita M. Porter, June 30, 2023
# This file generates the tax_meta file that combines taxonomy with metadata in later scripts
# Script to summarize taxonomy NO rarefaction
# For reads and ESVs
# For Disturbance and DevStage

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(scales)
library(viridis)
library(cowplot)
library(treemap)
library(grid) #veiwport

# Read infile with meta ----

# read in taxonomy & read data from combined file
# denoised to remove ESVs with < 0.01% reads
AM <- read.csv("outfiles/tax_meta2.csv", header = TRUE, stringsAsFactors = FALSE)

# Relative abundance ----

# Calculate heatmaps using read counts
# pivot table
A.1 <- data.frame(AM %>% group_by(Disturbance, DevStage, Phylum, Class, Amplicon) %>% 
                    dplyr::summarize(sum(ESVsize)))
names(A.1)[6] <- "Reads"

# convert to relative abundance within each Disturbance-DevStage
A.2 <- A.1 %>% group_by(Disturbance,DevStage) %>% dplyr::summarize(sum=sum(Reads))

# add sum column to A.1
A.3 <- merge(A.1, A.2, by=c("Disturbance","DevStage"), all.x = TRUE)

# calc relative read abundance
A.3$relabund <- A.3$Reads/A.3$sum*100

# create new field to help organize taxonomy
A.3$tax <- paste(A.3$Disturbance, A.3$DevStage, A.3$Phylum, A.3$Class, A.3$Disturbance, sep="; ")

# set log scale breaks in semi-automated manner
# https://stackoverflow.com/questions/55113333/ggplot2-introduce-breaks-on-a-x-log-scale
# fill_breaks = 10^pretty(log10(A.3$relabund))

# create factors
A.3$Amplicon <- factor(A.3$Amplicon, levels=c("BE", "F230", "16S", "ITS"))
A.3$Disturbance <- factor(A.3$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
                          labels=c("Wildfire", "Clearcut", "Salvage"))
A.3$DevStage <- factor(A.3$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                       labels=c("E", "CC", "EST", "LST", "M"))

# plot taxonomy for each marker separately
A_16S <- A.3[A.3$Amplicon=="16S",]
A_ITS <- A.3[A.3$Amplicon=="ITS",]
# fix formatting at class rank
A_ITS$Class <- gsub("c__", "", A_ITS$Class)
A_ITS$Class <- gsub("p__", "", A_ITS$Class)
A_ITS$Class <- gsub("k__Fungi_unidentified_unidentified", "Fungi unidentified", A_ITS$Class)

A_BE <- A.3[A.3$Amplicon=="BE",]
A_F230 <- A.3[A.3$Amplicon=="F230",]
A_COI <- rbind(A_BE, A_F230)
  
  
p1 <- ggplot(A_16S, aes(x=Amplicon, y=Class, fill=relabund)) + 
  geom_tile() +
  ggtitle("Bacteria") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#F8766D') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "Relative read abundance (%)", 
                                       title.position = "top")))
p1

p2 <- ggplot(A_ITS, aes(x=Amplicon, y=Class, fill=relabund)) + 
  geom_tile() +
  ggtitle("Fungi") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#00BA38') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "Relative read abundance (%)", 
                                       title.position = "top")))
p2

p3 <- ggplot(A_COI, aes(x=Amplicon, y=Class, fill=relabund)) + 
  geom_tile() +
  ggtitle("Arthropoda") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#619CFF') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "Relative read abundance (%)", title.position = "top")))
p3

l <- plot_grid(p1, nrow=1, labels=c("A"))
l
r <- plot_grid(p2, p3, ncol=1, rel_heights = c(0.7,0.3), labels=c("B","C"))
r
g <- plot_grid(l, r, nrow=1, rel_widths = c(1,1))
g
ggsave("outfiles/FigS2_relabund.pdf", g, height = 10, width = 8)

# ESV Richness ----

# Calculate heatmaps using ESV counts
# pivot table
A.4 <- data.frame(AM %>% 
                    group_by(Disturbance, DevStage, Phylum, Class, Amplicon) %>% 
                    dplyr::summarize(length(unique(GlobalESV))))
names(A.4)[6] <- "ESVs"

# convert ESVsize from integer to numeric
A.4$ESVs <- as.numeric(A.4$ESVs)

# create new field to help organize taxonomy
A.4$tax <- paste(A.4$Disturbance, A.4$DevStage, A.4$Phylum, A.4$Class, sep="; ")

# create factors
A.4$Amplicon <- factor(A.4$Amplicon, levels=c("BE", "F230", "16S", "ITS"))
A.4$Disturbance <- factor(A.4$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
                          labels=c("Wildfire", "Clearcut", "Salvage"))
A.4$DevStage <- factor(A.4$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                       labels=c("E", "CC", "EST", "LST", "M"))

# plot taxonomy for each marker separately
A_16S2 <- A.4[A.4$Amplicon=="16S",]
A_ITS2 <- A.4[A.4$Amplicon=="ITS",]
# fix formatting at class rank
A_ITS2$Class <- gsub("c__", "", A_ITS2$Class)
A_ITS2$Class <- gsub("p__", "", A_ITS2$Class)
A_ITS2$Class <- gsub("k__Fungi_unidentified_unidentified", "Fungi unidentified", A_ITS2$Class)

A_BE2 <- A.4[A.4$Amplicon=="BE",]
A_F2302 <- A.4[A.4$Amplicon=="F230",]
A_COI2 <- rbind(A_BE2, A_F2302)


p4 <- ggplot(A_16S2, aes(x=Amplicon, y=Class, fill=ESVs)) + 
  geom_tile() +
  ggtitle("Bacteria") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#F8766D') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "ESV Richness", 
                                       title.position = "top")))
p4

p5 <- ggplot(A_ITS2, aes(x=Amplicon, y=Class, fill=ESVs)) + 
  geom_tile() +
  ggtitle("Fungi") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#00BA38') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "ESV richness", 
                                       title.position = "top")))
p5

p6 <- ggplot(A_COI2, aes(x=Amplicon, y=Class, fill=ESVs)) + 
  geom_tile() +
  ggtitle("Arthropoda") +
  labs(x="", y="") +
  facet_grid(Disturbance ~ DevStage, switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "white", high = '#619CFF') +
  theme_bw() +
  theme(
    title = element_text(size=7),
    axis.title = element_text(size=7),
    axis.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7)) +
  guides(guides(fill = guide_colourbar(title = "ESV Richness", title.position = "top")))
p6

l2 <- plot_grid(p4, nrow=1, labels=c("A"))
l2
r2 <- plot_grid(p5, p6, ncol=1, rel_heights = c(0.7,0.3), labels=c("B","C"))
r2
g <- plot_grid(l2, r2, nrow=1, rel_widths = c(1,1))
g
ggsave("outfiles/FigS2_ESVrichness.pdf", g, height = 10, width = 8)






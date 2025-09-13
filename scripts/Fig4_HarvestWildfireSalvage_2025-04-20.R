# Teresita M. Porter, Apr 20, 2025

# Create plate for Harvest Type
# Richness, PERMANOVA, NMDS, RDA?

# Reference condition is wildfire 
## compare harvest and salvage to wildfire at each dev stage, E and CC

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(gridExtra) # grid.arrange
library(cowplot) # get_legend
library(ggpubr) # normality
library(dplyr) # summarize_all
library(RColorBrewer)

# Richness ----
## Bacteria ----


# consider working with non-rarefied data (just drop poorly sequenced samples) since they are all saturated anyways?

# read in merged taxonomy and metadata from tax_summary.R
a <- read.csv("outfiles/tax_meta2.csv", header=TRUE, stringsAsFactors = FALSE)
# 106634     53

# pool results from F230 and BE
a$Amplicon <- gsub("F230", "COI", a$Amplicon)
a$Amplicon <- gsub("BE", "COI", a$Amplicon)
a$GlobalESV <- gsub("F230", "COI", a$GlobalESV)
a$GlobalESV <- gsub("BE", "COI", a$GlobalESV)

# Create new column for dcast, pool across layers
a$Sample <- paste(a$SiteCode, a$Disturbance, a$DevStage, a$Replicate, sep="_")

# just bacteria
b <- a[grepl("16S_", a$GlobalESV),]

# only keep E and CC
b <- b[(b$DevStage=="Establishment" |
          b$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    54 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 54 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 480

# vegan
richness <- specnumber(esv.notnull2)

# create df
b.df <- data.frame(sample=names(richness), richness=richness)

# visual normality check
ggqqplot(b.df$richness)

#Shapiro-Wilk test of normality
shapiro.test(b.df$richness)
# FTR null hypothesis of normality (normal)

# Get separate method and siterep cols
setDT(b.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(b.df)[colnames(b.df)=="S1"] <- "SiteCode"
colnames(b.df)[colnames(b.df)=="S2"] <- "Disturbance"
colnames(b.df)[colnames(b.df)=="S3"] <- "DevStage"
colnames(b.df)[colnames(b.df)=="S4"] <- "Replicate"

# create factors
b.df$Disturbance <- factor(b.df$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
                           labels=c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
b.df$DevStage <- factor(b.df$DevStage, levels=c("Establishment", "CrownClosure"),
                        labels=c("Establishment", "Crown Closure"))

# blues <- brewer.pal(9,"Blues")
# # blues <- blues[c(8,4)]
# reds <- brewer.pal(9,"Reds")
# # reds <- reds[c(8,4)]
# cols <- c(reds[8], blues[8], "burlywood3")
# # cols for disturbance type
# Star Trek color palette lol (red, blue, amber)
cols <- c("#df0000", "#0099f6", "#f2c300")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# [6] "#4292C6" "#2171B5" "#08519C" "#08306B"

# blues <- brewer.pal(9,"Blues")
# cols <- cols <- c("red", "darkblue", "chocolate4")

p1.tmp <- ggplot(b.df, aes(x=Disturbance, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE, outlier.size = 0.5) +
  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  ggtitle("Bacteria") +
  labs(x="", y="Sequence Variant Richness") +
  facet_wrap(~DevStage) +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme( text = element_text(size=8),
         plot.title = element_text(size=8, hjust=0.01),
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(),
         axis.title.x = element_text(vjust = -0.75),
         legend.title = element_blank(),
         legend.text = element_text(size=8),
         legend.position = "top",
         strip.background = element_blank(),
         strip.text.x = element_text(angle = 0, hjust = 0.5))
p1.tmp

l <- get_legend(p1.tmp)

my_comparisons <- list(c("Wildfire", "Clearcut\nHarvest"), c("Clearcut\nHarvest", "Salvage\nHarvest"), c("Wildfire", "Salvage\nHarvest"),
                       c("Wildfire", "Clearcut\nHarvest"), c("Clearcut\nHarvest", "Salvage\nHarvest"), c("Wildfire", "Salvage\nHarvest"))

my_y <- c(420, 450, 480)

p1 <- ggplot(b.df, aes(x=Disturbance, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE, outlier.size = 0.5) +
  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  ggtitle("Bacteria") +
  labs(x="", y="Sequence Variant Richness") +
  facet_wrap(~DevStage) +
  ylim(0,500) +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme( text = element_text(size=8),
         plot.title = element_text(size=8, hjust=0.01),
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(),
         axis.title.x = element_text(vjust = -0.75),
         legend.title = element_blank(),
         legend.text = element_text(size=8),
         legend.position = "none",
         strip.background = element_blank(),
         strip.text.x = element_text(angle = 0, hjust = 0.5)) +
  stat_compare_means( 
    comparisons = my_comparisons, 
    label.y=my_y,
    method = "t.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"), 
    show.legend=FALSE)
p1

# inspired by https://github.com/kassambara/ggpubr/issues/65
## build plot
p1.1 <- ggplot_build(p1)

## look at data and replace the label with the properly calculated p.adj
#p1.1$data[[3]]$label <- p1.1$data[[3]]$p.adj

# convert col from chr to numeric
p1.1$data[[3]]$annotation <- as.numeric(unlist(p1.1$data[[3]]$annotation))


# reformat to use star notation
p1.1$data[[3]]$annotation <- ifelse(p1.1$data[[3]]$annotation < 0.001, "***",
                                    ifelse(p1.1$data[[3]]$annotation < 0.01, "**",
                                           ifelse(p1.1$data[[3]]$annotation < 0.05, "*", "")))
## plot anew
plot(ggplot_gtable(p1.1))


## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# only keep E and CC
f <- f[(f$DevStage=="Establishment" |
          f$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    54 680

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 679

# vegan
richness <- specnumber(esv.notnull2)

# create df
f.df <- data.frame(sample=names(richness), richness=richness)

# visual normality check
ggqqplot(f.df$richness)

#Shapiro-Wilk test of normality
shapiro.test(f.df$richness)
# FTR null hypothesis of normality (normal)

# Get separate method and siterep cols
setDT(f.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(f.df)[colnames(f.df)=="S1"] <- "SiteCode"
colnames(f.df)[colnames(f.df)=="S2"] <- "Disturbance"
colnames(f.df)[colnames(f.df)=="S3"] <- "DevStage"
colnames(f.df)[colnames(f.df)=="S4"] <- "Replicate"

# create factors
f.df$Disturbance <- factor(f.df$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
                           labels=c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
f.df$DevStage <- factor(f.df$DevStage, levels=c("Establishment", "CrownClosure"),
                        labels=c("Establishment", "Crown Closure"))

# my_comparisons <- list(c("Wildfire", "Harvest"), c("Wildfire", "Salvage"),
                       # c("Wildfire", "Harvest"), c("Wildfire", "Salvage"))

my_y <- c(330, 360, 390)

p2 <- ggplot(f.df, aes(x=Disturbance, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE, outlier.size = 0.5) +
  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  ggtitle("Fungi") +
  labs(x="", y="Sequence Variant Richness") +
  facet_wrap(~DevStage) +
  ylim(0,500) +
  theme(legend.title=element_blank()) +
  # scale_color_viridis_d(begin=0, end=1)+
  scale_color_manual(values = cols) +
  # scale_fill_viridis_d(begin=0, end=1)+
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme( text = element_text(size=8),
         plot.title = element_text(size=8, hjust=0.01),
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(),
         axis.title.x = element_text(vjust = -0.75),
         legend.title = element_blank(),
         legend.text = element_text(size=8),
         legend.position = "none",
         strip.background = element_blank(),
         strip.text.x = element_text(angle = 0, hjust = 0.5)) +
  stat_compare_means( 
    comparisons = my_comparisons, 
    label.y=my_y,
    method = "t.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"), 
    show.legend=FALSE)
p2

# inspired by https://github.com/kassambara/ggpubr/issues/65
## build plot
p2.1 <- ggplot_build(p2)

## look at data and replace the label with the properly calculated p.adj
#p2.1$data[[3]]$label <- p2.1$data[[3]]$p.adj

# convert col from chr to numeric
p2.1$data[[3]]$annotation <- as.numeric(unlist(p2.1$data[[3]]$annotation))

# reformat to use star notation
p2.1$data[[3]]$annotation <- ifelse(p2.1$data[[3]]$annotation < 0.001, "***",
                                    ifelse(p2.1$data[[3]]$annotation < 0.01, "**",
                                           ifelse(p2.1$data[[3]]$annotation < 0.05, "*", "")))
## plot anew
plot(ggplot_gtable(p2.1))


## Arthropoda ----


# just arthropoda
c <- a[grepl("COI_", a$GlobalESV),]

# keep E and CC
c <- c[(c$DevStage=="Establishment" |
          c$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]     54 458

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 457

# vegan
richness <- specnumber(esv.notnull2)

# create df
c.df <- data.frame(sample=names(richness), richness=richness)

# visual normality check
ggqqplot(c.df$richness)

#Shapiro-Wilk test of normality
shapiro.test(c.df$richness)
# FTR null hypothesis of normality (normal)

# Get separate method and siterep cols
setDT(c.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(c.df)[colnames(c.df)=="S1"] <- "SiteCode"
colnames(c.df)[colnames(c.df)=="S2"] <- "Disturbance"
colnames(c.df)[colnames(c.df)=="S3"] <- "DevStage"
colnames(c.df)[colnames(c.df)=="S4"] <- "Replicate"

# create factors
c.df$Disturbance <- factor(c.df$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
                           labels=c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
c.df$DevStage <- factor(c.df$DevStage, levels=c("Establishment", "CrownClosure"),
                        labels=c("Establishment", "Crown Closure"))


# my_comparisons <- list(c("Wildfire", "Harvest"), c("Wildfire", "Salvage"),
                       # c("Wildfire", "Harvest"), c("Wildfire", "Salvage"))

my_y <- c(180, 210, 240)

p3 <- ggplot(c.df, aes(x=Disturbance, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE, outlier.size = 0.5) +
  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  ggtitle("Arthropoda") +
  labs(x="", y="Sequence Variant Richness") +
  facet_wrap(~DevStage) +
  ylim(0,500) +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme( text = element_text(size=8),
         plot.title = element_text(size=8, hjust=0.01),
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(),
         axis.title.x = element_text(vjust = -0.75),
         legend.title = element_blank(),
         legend.text = element_text(size=8),
         legend.position = "none",
         strip.background = element_blank(),
         strip.text.x = element_text(angle = 0, hjust = 0.5)) +
  stat_compare_means( 
    comparisons = my_comparisons, 
    label.y=my_y,
    method = "t.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"), 
    show.legend=FALSE)
p3

# inspired by https://github.com/kassambara/ggpubr/issues/65
## build plot
p3.1 <- ggplot_build(p3)

## look at data and replace the label with the properly calculated p.adj
#p3.1$data[[3]]$label <- p3.1$data[[3]]$p.adj

# convert col from chr to numeric
p3.1$data[[3]]$annotation <- as.numeric(unlist(p3.1$data[[3]]$annotation))

# reformat to use star notation
p3.1$data[[3]]$annotation <- ifelse(p3.1$data[[3]]$annotation < 0.001, "***",
                                    ifelse(p3.1$data[[3]]$annotation < 0.01, "**",
                                           ifelse(p3.1$data[[3]]$annotation < 0.05, "*", "")))
## plot anew
plot(ggplot_gtable(p3.1))








library(data.table) # setDT
library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplotf
library(goeveg) # scree
library(plyr) # ddply
library(gridExtra) #grid.arrange
library(otuSummary) #matrixConvert
library(usedist) #dist_subset
library(dplyr)

#library(devtools)
#Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# NMDS ----

## Bacteria ----


# read in merged taxonomy and metadata from tax_summary.R
a <- read.csv("outfiles/tax_meta2.csv", header=TRUE, stringsAsFactors = FALSE)
# 106634     53

# pool results from F230 and BE
a$Amplicon <- gsub("F230", "COI", a$Amplicon)
a$Amplicon <- gsub("BE", "COI", a$Amplicon)
a$GlobalESV <- gsub("F230", "COI", a$GlobalESV)
a$GlobalESV <- gsub("BE", "COI", a$GlobalESV)

# Create new column for dcast, pool across layers
a$Sample <- paste(a$SiteCode, a$Disturbance, a$DevStage, a$Replicate, sep="_")

# just bacteria
b <- a[grepl("16S_", a$GlobalESV),]

# keep E and CC
b <- b[(b$DevStage=="Establishment" |
          b$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    54 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 54 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 480

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# Do 3 dimensional NMDS
# Default uses Bray Curtis dissimilarity
# When using presence-absence data ~ Sorensen dissimilarity
nmds3_b <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.07396399 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Fig4_stressplot_bacteria.pdf")
# stressplot(nmds3_b)
# gof <-goodness(nmds3_b)
# gof
# plot(nmds3_b, display = "sites", type="n", main="SSU")
# points(nmds3_b, display="sites",cex=2*gof/mean(gof))
# dev.off()
#linear R2 = 0.972


# Create grouping matrix for samples by grabbing row names from above matrix
names_b <- data.frame(row.names(esv.notnull2), stringsAsFactors = FALSE)
# Rename the column
names(names_b) <- "sample"
# Copy column to row names
row.names(names_b) <- names_b$sample
# Split first column into their own fields
names_b.1 <- data.frame(names_b, do.call(rbind, strsplit(names_b$sample,'_')), stringsAsFactors = FALSE)
names(names_b.1)[2:5]<-c("SiteCode", "Disturbance", "DevStage", "Replicate")
# Remove first column
names_b.1 <- names_b.1[,-1]
# Grab sites/species scores from NMDS output
df.b <- data.frame(scores(nmds3_b, display = "sites"))
# Put it all in one df for ggplot
gg_b <- merge(df.b, names_b.1, by="row.names")

# create factors
gg_b$Disturbance <- factor(gg_b$Disturbance,
                           levels = c("Fire", "Harvest", "Salvage"),
                           labels = c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
gg_b$DevStage <- factor(gg_b$DevStage,
                        levels = c("Establishment","CrownClosure"),
                        labels = c("E", "CC"))

# Create metadata from rownames 'sample'
env_b <- gg_b[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_b <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
#bd.layer <- betadisper(sor_b, as.factor(env_b$Layer))
bd.disturbance <- betadisper(sor_b, as.factor(env_b$Disturbance))
bd.devstage <- betadisper(sor_b, as.factor(env_b$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
#anova(bd.layer) # 0.041 *
anova(bd.disturbance) # n/s
anova(bd.devstage) # n/s
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# # plot distance to centroid for fire and harvest separately for devstage
# fire_filter <- grepl("Fire", names(bd.devstage$distances))
# distances <- bd.devstage$distances[fire_filter]
# group <- bd.devstage$group[fire_filter]
# b.bd.devstage.f <- data.frame(distances=distances, group=group)
# b.bd.devstage.f$disturbance <- "Fire"
# 
# harvest_filter <- grepl("Harvest", names(bd.devstage$distances))
# distances <- bd.devstage$distances[harvest_filter]
# group <- bd.devstage$group[harvest_filter]
# b.bd.devstage.h <- data.frame(distances=distances, group=group)
# b.bd.devstage.h$disturbance <- "Harvest"
# 
# salvage_filter <- grepl("Salvage", names(bd.devstage$distances))
# distances <- bd.devstage$distances[salvage_filter]
# group <- bd.devstage$group[salvage_filter]
# b.bd.devstage.s <- data.frame(distances=distances, group=group)
# b.bd.devstage.s$disturbance <- "Salvage"
# 
# # put it together for ggplot
# b.bd.devstage <- rbind(b.bd.devstage.f, b.bd.devstage.h, b.bd.devstage.s)
# b.bd.devstage$amplicon <- "16S"

### Pairwise PERMANOVA ----

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
b.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_b.e <- vegdist(b.e, "bray", binary=TRUE)
env_b.e <- env_b[grepl("E", env_b$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_b.e, env_b.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.e <- betadisper(sor_b.e, as.factor(env_b.e$Disturbance))
bd.disturbance.b.e <- data.frame(distances=bd.disturbance.e$distances, group=bd.disturbance.e$group)
bd.disturbance.b.e$amplicon <- "16S"
bd.disturbance.b.e$DevStage <- "E"

b.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_b.cc <- vegdist(b.cc, "bray", binary=TRUE)
env_b.cc <- env_b[grepl("CC", env_b$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_b.cc, env_b.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.cc <- betadisper(sor_b.cc, as.factor(env_b.cc$Disturbance))
bd.disturbance.b.cc <- data.frame(distances=bd.disturbance.cc$distances, group=bd.disturbance.cc$group)
bd.disturbance.b.cc$amplicon <- "16S"
bd.disturbance.b.cc$DevStage <- "CC"

# put it together
comp.all <- rbind(comp.e, comp.cc)
comp.all$taxon <- "Bacteria"


# Fit environmental variables
# Read in metadata

m <- read.csv(file='infiles/metadata_2024-11-10_DM.csv', head=TRUE)

# create sample col to match sample from matrix
# 12H_13_B_Harvest_Establishment_1
m$SampleName <- paste(m$SiteCode, m$Disturbance, m$DevStage, m$Replicate, sep="_")

# remove samples that are missing from the nmds plot
missing <- setdiff(m$SampleName, env_b$Row.names)

# reorder columns
m.b <- m[!m$SampleName %in% missing,]
m.b <- m.b[,c(21,1:20)]
#rownames(m.b) <- m.b$SampleName

# average values across soil layers where needed
m.b.2 <- data.frame(m.b %>% 
                      group_by(SampleName) %>% 
                      dplyr::summarize(
                                       pH_mean=mean(pH, na.rm = TRUE), 
                                       TOC_mean=mean(TOC_gkg, na.rm = TRUE), 
                                       TN_mean=mean(TN_gkg, na.rm = TRUE)))
rownames(m.b.2) <- m.b.2$SampleName
m.b.2$SampleName <- NULL

# reorder rows based on rownames(b)
m.b.2 <- m.b.2[match(rownames(esv.notnull2), rownames(m.b.2)), ]
# keep continuous data, focus on "F" part
#m.b <- m.b[,c(22:24,26,28)]
names(m.b.2) <- c( "pH", "TOC", "TN")

# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_b, m.b.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.b <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.b$x <- 0
fit.df.sig.b$y <- 0


# plot disturbance * devstage separately

# color by development stage
gg_b.e <- gg_b[gg_b$DevStage=="E",]
chulls.disturbance.b.e <- ddply(gg_b.e, .(Disturbance), function(gg_b.e) gg_b.e[chull(gg_b.e$NMDS1, gg_b.e$NMDS2), ])

# NMDS plot, add sig env vars
b_fitted_disturbance12_e.tmp <- ggplot(data=gg_b.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.b.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.e, size = 2.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.b, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Bacteria - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
b_fitted_disturbance12_e.tmp

l2 <- get_legend(b_fitted_disturbance12_e.tmp)

# NMDS plot, add sig env vars
b_fitted_disturbance12_e <- ggplot(data=gg_b.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.b.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.b, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-0.9,0.8) +
  ylim(-1,1) +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Bacteria - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
b_fitted_disturbance12_e

# color by development stage
gg_b.cc <- gg_b[gg_b$DevStage=="CC",]
chulls.disturbance.b.cc <- ddply(gg_b.cc, .(Disturbance), function(gg_b.cc) gg_b.cc[chull(gg_b.cc$NMDS1, gg_b.cc$NMDS2), ])

# NMDS plot, add sig env vars
b_fitted_disturbance12_cc <- ggplot(data=gg_b.cc, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.b.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.b, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-0.9,0.8) +
  ylim(-1,1) +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Bacteria - Crown Closure") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
b_fitted_disturbance12_cc

g1 <- grid.arrange(b_fitted_disturbance12_e, b_fitted_disturbance12_cc, ncol=2)


## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# keep E and CC
f <- f[(f$DevStage=="Establishment" |
          f$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    54 680

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 679

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# Do 3 dimensional NMDS
nmds3_f <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.1111676

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Fig4_stressplot_fungi.pdf")
# stressplot(nmds3_f)
# gof <-goodness(nmds3_f)
# gof
# plot(nmds3_f, display = "sites", type="n", main="ITS")
# points(nmds3_f, display="sites",cex=2*gof/mean(gof))
# dev.off()
#linear R2 = 0.926




# Create grouping matrix for samples by grabbing row names from above matrix
names_f <- data.frame(row.names(esv.notnull2), stringsAsFactors = FALSE)
# Rename the column
names(names_f) <- "sample"
# Copy column to row names
row.names(names_f) <- names_f$sample
# Split first column into their own fields
names_f.1<-data.frame(names_f, do.call(rbind, strsplit(names_f$sample,'_')), stringsAsFactors = FALSE)
names(names_f.1)[2:5]<-c("SiteCode", "Disturbance", "DevStage", "Replicate")
# Remove first column
names_f.1 <- names_f.1[,-1]
# Grab sites/species scores from NMDS output
df.f <- data.frame(scores(nmds3_f, display = "sites"))
# Put it all in one df for ggplot
gg_f <- merge(df.f, names_f.1, by="row.names")

# create factors
gg_f$Disturbance <- factor(gg_f$Disturbance,
                           levels = c("Fire", "Harvest", "Salvage"),
                           labels = c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
gg_f$DevStage <- factor(gg_f$DevStage,
                        levels = c("Establishment","CrownClosure"),
                        labels = c("E", "CC"))

# Create metadata from rownames 'sample'
env_f <- gg_f[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_f <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
# bd.layer <- betadisper(sor_f, as.factor(env_f$Layer))
bd.disturbance <- betadisper(sor_f, as.factor(env_f$Disturbance))
bd.devstage <- betadisper(sor_f, as.factor(env_f$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
# anova(bd.layer) # n/s
anova(bd.disturbance) # 0.0005438 ***
anova(bd.devstage) # n/s
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# # plot distance to centroid for fire and harvest separately for devstage
# fire_filter <- grepl("Fire", names(bd.devstage$distances))
# distances <- bd.devstage$distances[fire_filter]
# group <- bd.devstage$group[fire_filter]
# f.bd.devstage.f <- data.frame(distances=distances, group=group)
# f.bd.devstage.f$disturbance <- "Fire"
# 
# harvest_filter <- grepl("Harvest", names(bd.devstage$distances))
# distances <- bd.devstage$distances[harvest_filter]
# group <- bd.devstage$group[harvest_filter]
# f.bd.devstage.h <- data.frame(distances=distances, group=group)
# f.bd.devstage.h$disturbance <- "Harvest"
# 
# salvage_filter <- grepl("Salvage", names(bd.devstage$distances))
# distances <- bd.devstage$distances[salvage_filter]
# group <- bd.devstage$group[salvage_filter]
# f.bd.devstage.s <- data.frame(distances=distances, group=group)
# f.bd.devstage.s$disturbance <- "Salvage"
# 
# # put it together for ggplot
# f.bd.devstage <- rbind(f.bd.devstage.f, b.bd.devstage.h, f.bd.devstage.s)
# f.bd.devstage$amplicon <- "ITS"

### Pairwise PERMANOVA ----

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
f.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_f.e <- vegdist(f.e, "bray", binary=TRUE)
env_f.e <- env_f[grepl("E", env_f$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_f.e, env_f.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.e <- betadisper(sor_f.e, as.factor(env_f.e$Disturbance))
bd.disturbance.f.e <- data.frame(distances=bd.disturbance.e$distances, group=bd.disturbance.e$group)
bd.disturbance.f.e$amplicon <- "ITS"
bd.disturbance.f.e$DevStage <- "E"

f.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_f.cc <- vegdist(f.cc, "bray", binary=TRUE)
env_f.cc <- env_f[grepl("CC", env_f$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_f.cc, env_f.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.cc <- betadisper(sor_f.cc, as.factor(env_f.cc$Disturbance))
bd.disturbance.f.cc <- data.frame(distances=bd.disturbance.cc$distances, group=bd.disturbance.cc$group)
bd.disturbance.f.cc$amplicon <- "ITS"
bd.disturbance.f.cc$DevStage <- "CC"

# put it together
comp.all2 <- rbind(comp.e, comp.cc)
comp.all2$taxon <- "Fungi"


# Fit environmental variables 
# Read in metadata

#m <- read.csv(file='LV_metadata_110321.csv', head=TRUE)

# create sample col to match sample from matrix
# 12H_13_B_Harvest_Establishment_1
#m$SampleName <- paste(m$SiteCode, m$Disturbance, m$DevStage, m$Replicate, sep="_")

# remove samples that are missing from the nmds plot
missing <- setdiff(m$SampleName, env_f$Row.names)

# reorder columns
m.f <- m[!m$SampleName %in% missing,]
m.f <- m.f[,c(21,1:20)]
#rownames(m.f) <- m.f$SampleName

# average values across soil layers where needed
m.f.2 <- data.frame(m.f %>% 
                      group_by(SampleName) %>% 
                      dplyr::summarize(
                                       pH_mean=mean(pH, na.rm = TRUE), 
                                       TOC_mean=mean(TOC_gkg, na.rm = TRUE), 
                                       TN_mean=mean(TN_gkg, na.rm = TRUE)))
rownames(m.f.2) <- m.f.2$SampleName
m.f.2$SampleName <- NULL

# reorder rows based on rownames(b)
m.f.2 <- m.f.2[match(rownames(esv.notnull2), rownames(m.f.2)), ]
# keep continuous data, focus on "F" part
#m.f <- m.f[,c(22:24,26,28)]
names(m.f.2) <- c( "pH", "TOC", "TN")

# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_f, m.f.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.f <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.f$x <- 0
fit.df.sig.f$y <- 0


# plot disturbance * devstage separately

# color by development stage
gg_f.e <- gg_f[gg_f$DevStage=="E",]
chulls.disturbance.f.e <- ddply(gg_f.e, .(Disturbance), function(gg_f.e) gg_f.e[chull(gg_f.e$NMDS1, gg_f.e$NMDS2), ])

# NMDS plot, add sig env vars
f_fitted_disturbance12_e.tmp <- ggplot(data=gg_f.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.f.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.e, size = 2.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.f, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Fungi - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
f_fitted_disturbance12_e.tmp

# l <- get_legend(f_fitted_disturbance12_e.tmp)

# NMDS plot, add sig env vars
f_fitted_disturbance12_e <- ggplot(data=gg_f.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.f.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.f, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-1.1,1.2) +
  #ylim() +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Fungi - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
f_fitted_disturbance12_e

# color by development stage
gg_f.cc <- gg_f[gg_f$DevStage=="CC",]
chulls.disturbance.f.cc <- ddply(gg_f.cc, .(Disturbance), function(gg_f.cc) gg_f.cc[chull(gg_f.cc$NMDS1, gg_f.cc$NMDS2), ])

# NMDS plot, add sig env vars
f_fitted_disturbance12_cc <- ggplot(data=gg_f.cc, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.f.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.f, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-1.1,1.2) +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Fungi - Crown Closure") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
f_fitted_disturbance12_cc

g2 <- grid.arrange(f_fitted_disturbance12_e, f_fitted_disturbance12_cc, ncol=2)


## Arthropoda ----


# just arthropoda
c <- a[grepl("COI_", a$GlobalESV),]

# keep E and CC
c <- c[(c$DevStage=="Establishment" |
          c$DevStage=="CrownClosure"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    54 458

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 54 457

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# Do 3 dimensional NMDS
nmds3_a <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.125591


# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Fig4_stressplot_arthropoda.pdf")
# stressplot(nmds3_a)
# gof <-goodness(nmds3_a)
# gof
# plot(nmds3_a, display = "sites", type="n", main="COI")
# points(nmds3_a, display="sites",cex=2*gof/mean(gof))
# dev.off()
# #linear R2 = 0.901


# Create grouping matrix for samples by grabbing row names from above matrix
names_a <- data.frame(row.names(esv.notnull2), stringsAsFactors = FALSE)
# Rename the column
names(names_a) <- "sample"
# Copy column to row names
row.names(names_a) <- names_a$sample
# Split first column into their own fields
names_a.1<-data.frame(names_a, do.call(rbind, strsplit(names_a$sample,'_')), stringsAsFactors = FALSE)
names(names_a.1)[2:5]<-c("SiteCode", "Disturbance", "DevStage", "Replicate")
# Remove first column
names_a.1 <- names_a.1[,-1]
# Grab sites/species scores from NMDS output
df.a <- data.frame(scores(nmds3_a, display = "sites"))
# Put it all in one df for ggplot
gg_a <- merge(df.a, names_a.1, by="row.names")

# create factors
gg_a$Disturbance <- factor(gg_a$Disturbance,
                           levels = c("Fire", "Harvest", "Salvage"),
                           labels = c("Wildfire", "Clearcut\nHarvest", "Salvage\nHarvest"))
gg_a$DevStage <- factor(gg_a$DevStage,
                        levels = c("Establishment","CrownClosure"),
                        labels = c("E", "CC"))

# Create metadata from rownames 'sample'
env_a <- gg_a[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_a <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
# bd.layer <- betadisper(sor_a, as.factor(env_a$Layer))
bd.disturbance <- betadisper(sor_a, as.factor(env_a$Disturbance))
bd.devstage <- betadisper(sor_a, as.factor(env_a$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
# anova(bd.layer) # n/s
anova(bd.disturbance) # 1.243e-07 ***
anova(bd.devstage) # 0.0454 *
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# # plot distance to centroid for fire and harvest separately for devstage
# fire_filter <- grepl("Fire", names(bd.devstage$distances))
# distances <- bd.devstage$distances[fire_filter]
# group <- bd.devstage$group[fire_filter]
# a.bd.devstage.f <- data.frame(distances=distances, group=group)
# a.bd.devstage.f$disturbance <- "Fire"
# 
# harvest_filter <- grepl("Harvest", names(bd.devstage$distances))
# distances <- bd.devstage$distances[harvest_filter]
# group <- bd.devstage$group[harvest_filter]
# a.bd.devstage.h <- data.frame(distances=distances, group=group)
# a.bd.devstage.h$disturbance <- "Harvest"
# 
# salvage_filter <- grepl("Salvage", names(bd.devstage$distances))
# distances <- bd.devstage$distances[salvage_filter]
# group <- bd.devstage$group[salvage_filter]
# a.bd.devstage.s <- data.frame(distances=distances, group=group)
# a.bd.devstage.s$disturbance <- "Salvage"
# 
# # put it together for ggplot
# a.bd.devstage <- rbind(a.bd.devstage.f, a.bd.devstage.h, a.bd.devstage.s)
# a.bd.devstage$amplicon <- "COI"

### Pairwise PERMANOVA ----

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
a.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_a.e <- vegdist(a.e, "bray", binary=TRUE)
env_a.e <- env_a[grepl("E", env_a$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_a.e, env_a.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.e <- betadisper(sor_a.e, as.factor(env_a.e$Disturbance))
bd.disturbance.a.e <- data.frame(distances=bd.disturbance.e$distances, group=bd.disturbance.e$group)
bd.disturbance.a.e$amplicon <- "COI"
bd.disturbance.a.e$DevStage <- "E"

a.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_a.cc <- vegdist(a.cc, "bray", binary=TRUE)
env_a.cc <- env_a[grepl("CC", env_a$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_a.cc, env_a.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"
# compare distances to centroid for each disturbance at each dev stage
bd.disturbance.cc <- betadisper(sor_a.cc, as.factor(env_a.cc$Disturbance))
bd.disturbance.a.cc <- data.frame(distances=bd.disturbance.cc$distances, group=bd.disturbance.cc$group)
bd.disturbance.a.cc$amplicon <- "COI"
bd.disturbance.a.cc$DevStage <- "CC"

# put it together
comp.all3 <- rbind(comp.e, comp.cc)
comp.all3$taxon <- "Arthropoda"


# Fit environmental variables 
# Read in metadata

#m <- read.csv(file='LV_metadata_110321.csv', head=TRUE)

# create sample col to match sample from matrix
# 12H_13_B_Harvest_Establishment_1
#m$SampleName <- paste(m$SiteCode, m$Disturbance, m$DevStage, m$Replicate, sep="_")

# remove samples that are missing from the nmds plot
missing <- setdiff(m$SampleName, env_a$Row.names)

# reorder columns
m.a <- m[!m$SampleName %in% missing,]
m.a <- m.a[,c(21,1:20)]

# average values across soil layers where needed
m.a.2 <- data.frame(m.a %>% 
                      group_by(SampleName) %>% 
                      dplyr::summarize(
                                       pH_mean=mean(pH, na.rm = TRUE), 
                                       TOC_mean=mean(TOC_gkg, na.rm = TRUE), 
                                       TN_mean=mean(TN_gkg, na.rm = TRUE)))
rownames(m.a.2) <- m.a.2$SampleName
m.a.2$SampleName <- NULL

# reorder rows based on rownames(b)
m.a.2 <- m.a.2[match(rownames(esv.notnull2), rownames(m.a.2)), ]
# keep continuous data, focus on "F" part
names(m.a.2) <- c( "pH", "TOC", "TN")

# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_a, m.a.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.a <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.a$x <- 0
fit.df.sig.a$y <- 0


# plot disturbance * devstage separately

# color by development stage
gg_a.e <- gg_a[gg_a$DevStage=="E",]
chulls.disturbance.a.e <- ddply(gg_a.e, .(Disturbance), function(gg_a.e) gg_a.e[chull(gg_a.e$NMDS1, gg_a.e$NMDS2), ])

# NMDS plot, add sig env vars
a_fitted_disturbance12_e.tmp <- ggplot(data=gg_a.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.a.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_a.e, size = 2.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.a, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.a,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.a)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Arthropod - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
a_fitted_disturbance12_e.tmp

# l <- get_legend(a_fitted_disturbance12_e.tmp)

# NMDS plot, add sig env vars
a_fitted_disturbance12_e <- ggplot(data=gg_a.e, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.a.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_a.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.a, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-1.1, 1.2) +
  ylim(-0.5, 1.0) +
  geom_text(fit.df.sig.a,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.a)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Arthropod - Establishment") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
a_fitted_disturbance12_e

# color by development stage
gg_a.cc <- gg_a[gg_a$DevStage=="CC",]
chulls.disturbance.a.cc <- ddply(gg_a.cc, .(Disturbance), function(gg_a.cc) gg_a.cc[chull(gg_a.cc$NMDS1, gg_a.cc$NMDS2), ])

# NMDS plot, add sig env vars
a_fitted_disturbance12_cc <- ggplot(data=gg_a.cc, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.disturbance.a.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_a.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  geom_segment(fit.df.sig.a, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  xlim(-1.1, 1.2) +
  ylim(-0.5, 1.0) +
  geom_text(fit.df.sig.a,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.a)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  ggtitle("Arthropod - Crown Closure") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(
    #color="none",
    #fill="none",
    #    shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)
  ) 
a_fitted_disturbance12_cc

g3 <- grid.arrange(a_fitted_disturbance12_e, a_fitted_disturbance12_cc, ncol=2)

# Pairwise PERMANOVA ----

# summarize results for pairwise PERMANOVA
comp.all.taxa <- rbind(comp.all, comp.all2, comp.all3)

# rename pairs for easier comparisons
comp.all.taxa$pairs <- gsub("Clearcut\nHarvest vs Wildfire", "Wildfire vs Clearcut\nHarvest", comp.all.taxa$pairs)
#comp.all.taxa$pairs <- gsub("Clearcut Harvest vs Wildfire", "Wildfire vs Clearcut Harvest", comp.all.taxa$pairs)
#comp.all.taxa$pairs <- gsub("Wildfire vs Clearcut Harvest", "Wildfire vs Clearcut Harvest", comp.all.taxa$pairs)
#comp.all.taxa$pairs <- gsub("Wildfire vs Salvage Harvest", "Wildfire vs Salvage Harvest", comp.all.taxa$pairs)
comp.all.taxa$pairs <- gsub("Salvage\nHarvest vs Clearcut\nHarvest", "Clearcut\nHarvest vs Salvage\nHarvest", comp.all.taxa$pairs)

# only keep comparisons of interest
# comp.all.taxa <- comp.all.taxa[(comp.all.taxa$pairs=="Wildfire vs Harvest" |
                                  # comp.all.taxa$pairs=="Wildfire vs Salvage"), ]



# Create factors
comp.all.taxa$pairs <- factor(comp.all.taxa$pairs, levels=rev(c("Wildfire vs Clearcut\nHarvest", "Wildfire vs Salvage\nHarvest", "Clearcut\nHarvest vs Salvage\nHarvest")),
                              labels=rev(c("W vs CCH", "W vs SH", "CCH vs SH")))

comp.all.taxa$taxon <- factor(comp.all.taxa$taxon, levels=rev(c("Bacteria", "Fungi", "Arthropoda")),
                              labels=rev(c("B", "F", "A")))

comp.all.taxa$DevStage <- factor(comp.all.taxa$DevStage, levels=c("Establishment", "CrownClosure"),
                                 labels=c("Establishment", "Crown Closure"))

p1 <- ggplot(comp.all.taxa, aes(x=p.adjusted, y=pairs, shape=taxon)) + 
  geom_text(aes(label=taxon), size=2.5, position = position_jitterdodge(seed=3)) +
  labs(x = "Holm adjusted p-values", y = "Stand Development Stages") +
  geom_hline(yintercept=0.05, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  scale_shape_manual(values = c("B", "F", "A")) +
  facet_wrap(~DevStage) +
  theme_bw() + 
  theme(
    text = element_text(size=8),
    plot.title = element_text(size=8, hjust=0.01),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0.5)) 
p1


left <- plot_grid(ggplot_gtable(p1.1), ggplot_gtable(p2.1), ggplot_gtable(p3.1), ncol=1)

#right_top <- p1

right_top <- plot_grid(g1, g2, g3, l2, ncol=1, rel_heights=c(1,1,1,0.2))

right_bottom <- p1

right <- plot_grid(right_top, right_bottom, ncol=1, labels=c("(b)", "(c)"), rel_heights=c(1, 0.4))

plate <- plot_grid(left, right, nrow=1, labels=c("(a)", ""))


ggsave("outfiles/Fig4_DisturbanceTypePlate.pdf", plate, width = 8, height = 10)

















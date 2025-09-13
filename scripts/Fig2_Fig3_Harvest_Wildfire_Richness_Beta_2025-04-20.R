# Teresita M. Porter, Apr. 20, 2024
# Q2 - fire vs harvest at each dev stage
# - harvest vs refcond (fire, mature)

# Richness convergence ----

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(gridExtra) # grid.arrange
library(cowplot) # get_legend
library(ggpubr) # normality, stat_compare_means
library(dplyr) # summarize_all
library(RColorBrewer)

# read in merged taxonomy and metadata from tax_summary.R
a <- read.csv("outfiles/tax_meta2.csv", header=TRUE, stringsAsFactors = FALSE)
# 106,634     53

# pool results from F230 and BE
a$Amplicon <- gsub("F230", "COI", a$Amplicon)
a$Amplicon <- gsub("BE", "COI", a$Amplicon)
a$GlobalESV <- gsub("F230", "COI", a$GlobalESV)
a$GlobalESV <- gsub("BE", "COI", a$GlobalESV)

# Create new column for dcast, pool across layers
a$Sample <- paste(a$SiteCode, a$Disturbance, a$DevStage, a$Replicate, sep="_")



## Bacteria ----


# Now do pairwise comparisons for each devstage
# just bacteria
b <- a[grepl("16S_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
b <- b[(b$Disturbance=="Harvest") |
         (b$Disturbance=="Fire"),]
b <- b[!grepl("Fire_Mature", b$Sample),]

# pivot to make esv matrix
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 66 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 66 480

# vegan SV richness
richness <- specnumber(esv.notnull2)

# create df
b.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(b.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(b.df)[colnames(b.df)=="S1"] <- "SiteCode"
colnames(b.df)[colnames(b.df)=="S2"] <- "Disturbance"
colnames(b.df)[colnames(b.df)=="S3"] <- "DevStage"
colnames(b.df)[colnames(b.df)=="S4"] <- "Replicate"

# create factors
b.df$Disturbance <- factor(b.df$Disturbance, levels=c("Harvest", "Fire"), 
                           labels=c("Clear Cut Harvest", "Wildfire"))
b.df$DevStage <- factor(b.df$DevStage, 
                        levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning"),
                        labels=c("E", "CC", "EST", "LST"))

# create col Disturbance + DevStage
b.df$label <- paste(b.df$Disturbance, b.df$DevStage, sep="_")

b.df$taxon <- "Bacteria"


## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
f <- f[(f$Disturbance=="Harvest") |
         (f$Disturbance=="Fire"),]
f <- f[!grepl("Fire_Mature", f$Sample),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 676

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
#  66 675

# vegan
richness <- specnumber(esv.notnull2)

# create df
f.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(f.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(f.df)[colnames(f.df)=="S1"] <- "SiteCode"
colnames(f.df)[colnames(f.df)=="S2"] <- "Disturbance"
colnames(f.df)[colnames(f.df)=="S3"] <- "DevStage"
colnames(f.df)[colnames(f.df)=="S4"] <- "Replicate"

# create factors
f.df$Disturbance <- factor(f.df$Disturbance, levels=c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
f.df$DevStage <- factor(f.df$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning"),
                        labels=c("E", "CC", "EST", "LST"))

# create Disturbance+DevStage field
f.df$label <- paste(f.df$Disturbance, f.df$DevStage, sep="_")

f.df$taxon <- "Fungi"




## Arthropoda ----


# just arthropoda
c <- a[grepl("COI_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
c <- c[(c$Disturbance=="Harvest") |
         (c$Disturbance=="Fire"),]
c <- c[!grepl("Fire_Mature", c$Sample),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 475

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
#  66 474

# vegan
richness <- specnumber(esv.notnull2)

# create df
c.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(c.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(c.df)[colnames(c.df)=="S1"] <- "SiteCode"
colnames(c.df)[colnames(c.df)=="S2"] <- "Disturbance"
colnames(c.df)[colnames(c.df)=="S3"] <- "DevStage"
colnames(c.df)[colnames(c.df)=="S4"] <- "Replicate"

# create factors
c.df$Disturbance <- factor(c.df$Disturbance, levels=c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
c.df$DevStage <- factor(c.df$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning"),
                        labels=c("E", "CC", "EST","LST"))

# create Disturbance + DevStage field
c.df$label <- paste(c.df$Disturbance, c.df$DevStage, sep="_")

c.df$taxon <- "Arthropoda"



# visual normality check
ggqqplot(b.df$richness)
#  normal
ggqqplot(f.df$richness)
# normal
ggqqplot(c.df$richness)
# normal

#Shapiro-Wilk test of normality
shapiro.test(b.df$richness)
# Reject null hypothesis of normality (not normal)
shapiro.test(f.df$richness)
# Reject null hypothesis of normality (not normal)
shapiro.test(c.df$richness)
# FTR null hypothesis of normality (normal)



# put these together into a single df

# Run tests below twice: do Wilcox for bacteria and fungi, t-test for arthropods 
# edit plots in Inkscape to get it to show up right



# put it together
all.df <- rbind(b.df, f.df, c.df)

my_comparisons <- list(c("Wildfire_E", "Clear Cut Harvest_E"), c("Wildfire_CC", "Clear Cut Harvest_CC"), c("Wildfire_EST", "Clear Cut Harvest_EST"), c("Wildfire_LST", "Clear Cut Harvest_LST")) 

# blues <- brewer.pal(9,"Blues")
# blues <- blues[c(8,4)]
# reds <- brewer.pal(9,"Reds")
# reds <- reds[c(8,4)]
cols <- c("#0099f6", "#df0000")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# [6] "#4292C6" "#2171B5" "#08519C" "#08306B"

# create factors
all.df$taxon <- factor(all.df$taxon, levels=c("Bacteria", "Fungi", "Arthropoda"))
all.df$Disturbance <- factor(all.df$Disturbance, levels=c("Clear Cut Harvest", "Wildfire"))

p1.tmp <- ggplot(all.df, aes(x=DevStage, y=richness, color=Disturbance, group=label)) +
  geom_boxplot(position=position_dodge(width = 0.8), outlier.size = 0.5) +
  # geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  labs(x="", y="ESV Richness") +
  facet_wrap(~taxon, ncol=1, scales = "free") +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  # scale_fill_manual(values = blues) +
  theme_bw() + 
  theme(
    text = element_text(size=8),
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
    strip.text.x = element_text(angle = 0, hjust = 0)) +
  # bacteria & fungi (wilcox test)
  stat_compare_means(
    method = "wilcox.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"),
    label.y=200,
    show.legend=FALSE, size=3)
  # arthropods (t-test)
  # stat_compare_means( 
  #   method = "t.test", paired=F, tip.length = 0,
  #   method.args=list(p.adjust.method = "holm"), 
  #   label.y=200, 
  #   show.legend=FALSE, size=3)
p1.tmp

l1 <- get_legend(p1.tmp)

p1 <- ggplot(all.df, aes(x=DevStage, y=richness, color=Disturbance, group=label)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE, outlier.size = 0.5) +
  # geom_dotplot(aes(fill=DevStage), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  labs(x="", y="Richness") +
  facet_wrap(~taxon, ncol=1, scales = "free") +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme(
    text = element_text(size=12),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0)) +
  # # bacteria & fungi (wilcox test)
  stat_compare_means(
    method = "wilcox.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"),
    label.y=200,
    show.legend=FALSE, size=3)
  # arthropods (t-test)
  # stat_compare_means( 
  #   method = "t.test", paired=F, tip.length = 0,
  #   method.args=list(p.adjust.method = "holm"), 
  #   label.y=200, 
  #   show.legend=FALSE, size=3)
p1

# inspired by https://github.com/kassambara/ggpubr/issues/65
## build plot
partA <- ggplot_build(p1)

## look at data and replace the label with the properly calcualted p.adj
partA$data[[2]]$label <- partA$data[[2]]$p.adj

# reformat to use star notation
partA$data[[2]]$label <- ifelse(partA$data[[2]]$label < 0.001, "***",
                                ifelse(partA$data[[2]]$label < 0.01, "**",
                                       ifelse(partA$data[[2]]$label < 0.05, "*", "")))


g1 <- plot_grid(ggplot_gtable(partA))
g1




# Richness recovery ----
# Compare each harvest devstage with refcond fire-mature
## Bacteria ----
b <- a[grepl("16S_", a$GlobalESV),]

# keep all Harvest
# keep Fire if Mature
b <- b[(b$Disturbance=="Harvest") |
         (b$Disturbance=="Fire" & b$DevStage=="Mature"),]

# pivot to make esv matrix
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 42 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 480

# vegan
richness <- specnumber(esv.notnull2)

# create df
b.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(b.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(b.df)[colnames(b.df)=="S1"] <- "SiteCode"
colnames(b.df)[colnames(b.df)=="S2"] <- "Disturbance"
colnames(b.df)[colnames(b.df)=="S3"] <- "DevStage"
colnames(b.df)[colnames(b.df)=="S4"] <- "Replicate"

# create factors
b.df$Disturbance <- factor(b.df$Disturbance, levels=c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
b.df$DevStage <- factor(b.df$DevStage, 
                        levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                        labels=c("E", "CC","EST","LST","M"))

b.df$taxon <- "Bacteria"


## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
f <- f[(f$Disturbance=="Harvest") |
         (f$Disturbance=="Fire" & f$DevStage=="Mature"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 634

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 633

# vegan
richness <- specnumber(esv.notnull2)

# create df
f.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(f.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(f.df)[colnames(f.df)=="S1"] <- "SiteCode"
colnames(f.df)[colnames(f.df)=="S2"] <- "Disturbance"
colnames(f.df)[colnames(f.df)=="S3"] <- "DevStage"
colnames(f.df)[colnames(f.df)=="S4"] <- "Replicate"

# create factors
f.df$Disturbance <- factor(f.df$Disturbance, levels=c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
f.df$DevStage <- factor(f.df$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                        labels=c("E","CC","EST","LST","M"))

f.df$taxon <- "Fungi"




## Arthropoda ----


# just arthropoda
c <- a[grepl("COI_", a$GlobalESV),]

# keep all Harvest 
# keep Fire if Mature
c <- c[(c$Disturbance=="Harvest") |
         (c$Disturbance=="Fire" & c$DevStage=="Mature"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 445

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 444

# vegan
richness <- specnumber(esv.notnull2)

# create df
c.df <- data.frame(sample=names(richness), richness=richness)

# Get separate method and siterep cols
setDT(c.df)[, paste0("S", 1:4) := tstrsplit(sample, "_")]
colnames(c.df)[colnames(c.df)=="S1"] <- "SiteCode"
colnames(c.df)[colnames(c.df)=="S2"] <- "Disturbance"
colnames(c.df)[colnames(c.df)=="S3"] <- "DevStage"
colnames(c.df)[colnames(c.df)=="S4"] <- "Replicate"

# create factors
c.df$Disturbance <- factor(c.df$Disturbance, levels=c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
c.df$DevStage <- factor(c.df$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                        labels=c("E","CC","EST","LST","M"))

c.df$taxon <- "Arthropoda"




# check for sig diff richness between disturbance types, for each development stage (check for convergence)

# visual normality check
ggqqplot(b.df$richness)
# not normal
ggqqplot(f.df$richness)
# normal
ggqqplot(c.df$richness)
# normal

#Shapiro-Wilk test of normality
shapiro.test(b.df$richness)
# Reject null hypothesis of normality (not normal)
shapiro.test(f.df$richness)
# FTR null hypothesis of normality (normal)
shapiro.test(c.df$richness)
# FTR null hypothesis of normality (normal)



# put these together into a single df

# Run tests below twice: do Wilcox for bacteria, t-tests for fungi and arthropods, 
# edit plots in Inkscape to get it to show up right


# combine the plots

all.df <- rbind(b.df, f.df, c.df)

# create factors
all.df$taxon <- factor(all.df$taxon, levels=c("Bacteria", "Fungi", "Arthropoda"))

#my_comparisons <- list(c("E", "M"), c("CC", "M"), c("EST", "M"), c("LST", "M")) 

cols <- c("#0099f6", "#df0000")
blues <- brewer.pal(9,"Blues")
reds <- brewer.pal(9,"Reds")
cols2 <- c(blues[c(4,5,6,7)], "#df0000")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# [6] "#4292C6" "#2171B5" "#08519C" "#08306B"


# create factors
all.df$DevStage <- factor(all.df$DevStage, levels=c("E", "CC", "EST", "LST", "M"),
                          labels=c("E", "CC", "EST", "LST", "M"))
all.df$DevStage2 <- all.df$DevStage
all.df$DevStage2 <- factor(all.df$DevStage, levels=c("E", "CC", "EST", "LST", "M"))


p2.tmp <- ggplot(all.df, aes(x=DevStage2, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), outlier.size = 0.5) +
#  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  labs(x="", y="Fungal Richness") +
  facet_wrap(~taxon, ncol=1, scales = "free") +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
#  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    text = element_text(size=8),
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
    strip.text.x = element_text(angle = 0, hjust = 0)) +
  # bacteria - Wilcox (more conservative)
  stat_compare_means(
    ref.group = "M", method = "wilcox.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"),
    label.y=290,
    show.legend=FALSE, size=3)
  # fungi & arthropods - T.test
  # stat_compare_means(
  #   ref.group = "M", method = "t.test", paired=F, tip.length = 0,
  #   method.args=list(p.adjust.method = "holm"),
  #   label.y=290,
  #   show.legend=FALSE, size=3)
p2.tmp

l2 <- get_legend(p2.tmp)

p2 <- ggplot(all.df, aes(x=DevStage2, y=richness, color=Disturbance)) +
  geom_boxplot(position=position_dodge(width = 0.8), outlier.size = 0.5) +
#  geom_dotplot(aes(fill=Disturbance), binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize=0.5) +
  labs(x="", y="Richness") +
  facet_wrap(~taxon, ncol=1, scales = "free") +
  theme(legend.title=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() + 
  theme(
    text = element_text(size=12),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0)) +
  # bacteria - Wilcox (more conservative)
  stat_compare_means(
    ref.group = "M", method = "wilcox.test", paired=F, tip.length = 0,
    method.args=list(p.adjust.method = "holm"),
    label.y=290,
    show.legend=FALSE, size=3)
  # fungi & arthropods - T.test
  # stat_compare_means(
  #   ref.group = "M", method = "t.test", paired=F, tip.length = 0,
  #   method.args=list(p.adjust.method = "holm"),
  #   label.y=290,
  #   show.legend=FALSE, size=3)
p2

# inspired by https://github.com/kassambara/ggpubr/issues/65
## build plot
partB <- ggplot_build(p2)

## look at data and replace the label with the properly calculated p.adj
partB$data[[2]]$label <- partB$data[[2]]$p.adj

# reformat to use star notation
partB$data[[2]]$label <- ifelse(partB$data[[2]]$label < 0.001, "***",
                                ifelse(partB$data[[2]]$label < 0.01, "**",
                                       ifelse(partB$data[[2]]$label < 0.05, "*", "")))

g2 <- plot_grid(ggplot_gtable(partB))
g2








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






# NMDS convergence ----

## Bacteria ----


# just bacteria
b <- a[grepl("16S_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
b <- b[(b$Disturbance=="Harvest") |
         (b$Disturbance=="Fire" ),]
b <- b[!grepl("Fire_Mature", b$Sample),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 66 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 66 480

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
pdf("outfiles/supporting/Q2_Scree_bacteria_pairwise.pdf")
# check dims
dimcheckMDS(esv.notnull2)
dev.off()
# use k=3

# Do 3 dimensional NMDS
nmds3_b <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.07814812 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_bacteria_pairwise.pdf")
# stressplot(nmds3_b)
# gof <-goodness(nmds3_b)
# gof
# plot(nmds3_b, display = "sites", type="n", main="SSU")
# points(nmds3_b, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.979

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
                           levels = c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
gg_b$DevStage <- factor(gg_b$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning"),
                        labels = c("E", "CC", "T", "LST"))

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
anova(bd.disturbance) # 0.01734 *
anova(bd.devstage) # 7.633e-05 ***
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_bacteria_pairwise.pdf")
# par(mfrow=c(2,2))
# #boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_b ~ Disturbance*DevStage, data=env_b, permutations=999)
#                       Df SumOfSqs      R2       F Pr(>F)    
# Disturbance           1  0.04733 0.03748  4.1790  0.004 ** 
# DevStage              3  0.35791 0.28347 10.5345  0.001 ***
# Disturbance:DevStage  3  0.20052 0.15881  5.9018  0.001 ***
# Residual             58  0.65685 0.52024                   
# Total                65  1.26261 1.00000 
# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# subset the distance matrix (don't recreate it)
b.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_b.e <- vegdist(b.e, "bray", binary=TRUE)
env_b.e <- env_b[grepl("E", env_b$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_b.e, env_b.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"

b.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_b.cc <- vegdist(b.cc, "bray", binary=TRUE)
env_b.cc <- env_b[grepl("CC", env_b$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_b.cc, env_b.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"

b.t <- esv.notnull2[grepl("Thinning", rownames(esv.notnull2)),]
sor_b.t <- vegdist(b.t, "bray", binary=TRUE)
env_b.t <- env_b[grepl("T", env_b$DevStage),]
# add comparisons to table 
comp.t <- pairwise.adonis(sor_b.t, env_b.t$Disturbance, p.adjust.m="holm")
comp.t$DevStage <- "Thinning"

b.lst <- esv.notnull2[grepl("LateStageThinning", rownames(esv.notnull2)),]
sor_b.lst <- vegdist(b.lst, "bray", binary=TRUE)
env_b.lst <- env_b[grepl("LST", env_b$DevStage),]
# add comparisons to table 
comp.lst <- pairwise.adonis(sor_b.lst, env_b.lst$Disturbance, p.adjust.m="holm")
comp.lst$DevStage <- "LateStageThinning"

# put it together
comp.b <- rbind(comp.e, comp.cc, comp.t, comp.lst)

comp.b$taxon <- "Bacteria"


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
names(m.b.2) <- c("pH", "TOC", "TN")



# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_b, m.b.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.b <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.b$x <- 0
fit.df.sig.b$y <- 0

# # fit environmental variables to NMDS2 & NMDS3
# fit23 <- envfit(nmds3_b, m.b.2, perm = 999, na.rm = TRUE, choices=c(2,3))
# fit.vectors23 <- fit23[[1]]
# fit.pvals23 <- fit23[[1]]$pvals
# fit.df23 <- as.data.frame(fit.vectors23[[1]])
# fit.df23$pvals <- fit.pvals23
# 
# fit.df.sig.b23 <- fit.df23[fit.df23$pvals <0.05,]
# 
# fit.df.sig.b23$x <- 0
# fit.df.sig.b23$y <- 0
# 
# 
# 
# # fit environmental variables to NMDS3 & NMDS1
# fit31 <- envfit(nmds3_b, m.b.2, perm = 999, na.rm = TRUE, choices=c(3,1))
# fit.vectors31 <- fit31[[1]]
# fit.pvals31 <- fit31[[1]]$pvals
# fit.df31 <- as.data.frame(fit.vectors31[[1]])
# fit.df31$pvals <- fit.pvals31
# 
# fit.df.sig.b31 <- fit.df31[fit.df31$pvals <0.05,]
# 
# fit.df.sig.b31$x <- 0
# fit.df.sig.b31$y <- 0

# subset by devstage
gg_b.e <- gg_b[gg_b$DevStage=="E",]
# color by disturbance
chulls.disturbance.b.e <- ddply(gg_b.e, .(Disturbance), function(gg_b.e) gg_b.e[chull(gg_b.e$NMDS1, gg_b.e$NMDS2), ])


# # color by development stage
# chulls.devstage.b <- ddply(gg_b, .(DevStage), function(gg_b) gg_b[chull(gg_b$NMDS1, gg_b$NMDS2), ])

# blues <- brewer.pal(9,"Blues")
# blues <- blues[c(5,8)]
# 4,6,7,8 -> E, CC, T, LST

# reds <- brewer.pal(9,"Reds")
# reds <- reds[c(5:9)]
# 5,6,7,8,9 -> E, CC, T, LST,M
# [1] "#FB6A4A" "#EF3B2C" "#CB181D" "#A50F15" "#67000D"

cols <- c("#0099f6", "#df0000")

# NMDS plot, add sig env vars
b_fitted_devstage12_e <- ggplot(data=gg_b.e, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.b.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-0.7, 0.4) +
  ylim(-0.20, 0.4) +
  geom_segment(fit.df.sig.b,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               linewidth=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.3, color="black") +
  annotate(geom = 'text', label = '', x = 0, y = 0.4, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.5, hjust = 1, vjust = 1) +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

b_fitted_devstage12_e

# NMDS plot, add sig env vars
b_fitted_devstage12_e.tmp <- ggplot(data=gg_b.e, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.b.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-0.7, 0.4) +
  ylim(-0.20, 0.4) +
  annotate(geom = 'text', label = '', x = 0, y = 0.4, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = ' ', x = Inf, y = 0.5, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.b,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               linewidth=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "bottom")
b_fitted_devstage12_e.tmp

l3 <- get_legend(b_fitted_devstage12_e.tmp)

# color by devstage
gg_b.cc <- gg_b[gg_b$DevStage=="CC",]
# facet by disturbance
chulls.disturbance.b.cc <- ddply(gg_b.cc, .(Disturbance), function(gg_b.cc) gg_b.cc[chull(gg_b.cc$NMDS1, gg_b.cc$NMDS2), ])

# NMDS plot, add sig env vars
b_fitted_devstage12_cc <- ggplot(data=gg_b.cc, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.b.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-0.7, 0.4) +
  ylim(-0.20, 0.4) +
  annotate(geom = 'text', label = '', x = 0, y = 0.4, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = ' ', x = Inf, y = 0.5, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.b,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
b_fitted_devstage12_cc

# subset by devstage
gg_b.t <- gg_b[gg_b$DevStage=="T",]
# color by disturbance
chulls.disturbance.b.t <- ddply(gg_b.t, .(Disturbance), function(gg_b.t) gg_b.t[chull(gg_b.t$NMDS1, gg_b.t$NMDS2), ])

# NMDS plot, add sig env vars
b_fitted_devstage12_t <- ggplot(data=gg_b.t, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.b.t, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.t, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-0.7, 0.4) +
  ylim(-0.20, 0.4) +
  annotate(geom = 'text', label = '', x = 0, y = 0.4, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = ' ', x = Inf, y = 0.5, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.b,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
b_fitted_devstage12_t

# subset by devstage
gg_b.lst <- gg_b[gg_b$DevStage=="LST",]
# color by disturbance
chulls.disturbance.b.lst <- ddply(gg_b.lst, .(Disturbance), function(gg_b.lst) gg_b.lst[chull(gg_b.lst$NMDS1, gg_b.lst$NMDS2), ])

# NMDS plot, add sig env vars
b_fitted_devstage12_lst <- ggplot(data=gg_b.lst, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.b.lst, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b.lst, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-0.7, 0.4) +
  ylim(-0.20, 0.4) +
  annotate(geom = 'text', label = '', x = 0, y = 0.4, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = ' ', x = Inf, y = 0.5, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.b,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
b_fitted_devstage12_lst

p3 <- plot_grid(b_fitted_devstage12_e, b_fitted_devstage12_cc,
                b_fitted_devstage12_t, b_fitted_devstage12_lst, nrow=1)
p3













## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
f <- f[(f$Disturbance=="Harvest") |
         (f$Disturbance=="Fire" ),]
f <- f[!grepl("Fire_Mature", f$Sample),]
# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 676

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL
# 66 675

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 66 675

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("outfiles/supporting/Q2_Scree_fungi_pairwise.pdf")
# # check dims
# dimcheckMDS(esv.notnull2)
# dev.off()
# # use k=3

# Do 3 dimensional NMDS
nmds3_f <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.1130992 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_fungi_pairwise.pdf")
# stressplot(nmds3_f)
# gof <-goodness(nmds3_f)
# gof
# plot(nmds3_f, display = "sites", type="n", main="SSU")
# points(nmds3_f, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.941

# Create grouping matrix for samples by grabbing row names from above matrix
names_f <- data.frame(row.names(esv.notnull2), stringsAsFactors = FALSE)
# Rename the column
names(names_f) <- "sample"
# Copy column to row names
row.names(names_f) <- names_f$sample
# Split first column into their own fields
names_f.1 <- data.frame(names_f, do.call(rbind, strsplit(names_f$sample,'_')), stringsAsFactors = FALSE)
names(names_f.1)[2:5]<-c("SiteCode", "Disturbance", "DevStage", "Replicate")
# Remove first column
names_f.1 <- names_f.1[,-1]
# Grab sites/species scores from NMDS output
df.f <- data.frame(scores(nmds3_f, display = "sites"))
# Put it all in one df for ggplot
gg_f <- merge(df.f, names_f.1, by="row.names")

# create factors
gg_f$Disturbance <- factor(gg_f$Disturbance,
                           levels = c("Harvest", "Fire"),
                           labels=c("Clear Cut Harvest", "Wildfire"))
gg_f$DevStage <- factor(gg_f$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning"),
                        labels = c("E", "CC", "T", "LST"))

# Create metadata from rownames 'sample'
env_f <- gg_f[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_f <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
#bd.layer <- betadisper(sor_b, as.factor(env_b$Layer))
bd.disturbance <- betadisper(sor_f, as.factor(env_f$Disturbance))
bd.devstage <- betadisper(sor_f, as.factor(env_f$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
#anova(bd.layer) # 0.041 *
anova(bd.disturbance) # 7.118e-05 ***
anova(bd.devstage) # 1.839e-05 ***
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_fungi_pairwise.pdf")
# par(mfrow=c(2,2))
# #boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_f ~ Disturbance*DevStage, data=env_f, permutations=999)
#                       Df SumOfSqs      R2      F Pr(>F)    
# Disturbance           1   0.4981 0.05325 5.0701  0.001 ***
# DevStage              3   2.0746 0.22181 7.0395  0.001 ***
# Disturbance:DevStage  3   1.0828 0.11577 3.6741  0.001 ***
# Residual             58   5.6976 0.60917                  
# Total                65   9.3530 1.00000  

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# subset the distance matrix (don't recreate it)
f.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_f.e <- vegdist(f.e, "bray", binary=TRUE)
env_f.e <- env_f[grepl("E", env_f$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_f.e, env_f.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"

f.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_f.cc <- vegdist(f.cc, "bray", binary=TRUE)
env_f.cc <- env_f[grepl("CC", env_f$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_f.cc, env_f.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"

f.t <- esv.notnull2[grepl("Thinning", rownames(esv.notnull2)),]
sor_f.t <- vegdist(f.t, "bray", binary=TRUE)
env_f.t <- env_f[grepl("T", env_f$DevStage),]
# add comparisons to table 
comp.t <- pairwise.adonis(sor_f.t, env_f.t$Disturbance, p.adjust.m="holm")
comp.t$DevStage <- "Thinning"

f.lst <- esv.notnull2[grepl("LateStageThinning", rownames(esv.notnull2)),]
sor_f.lst <- vegdist(f.lst, "bray", binary=TRUE)
env_f.lst <- env_f[grepl("LST", env_f$DevStage),]
# add comparisons to table 
comp.lst <- pairwise.adonis(sor_f.lst, env_f.lst$Disturbance, p.adjust.m="holm")
comp.lst$DevStage <- "LateStageThinning"

# put it together
comp.f <- rbind(comp.e, comp.cc, comp.t, comp.lst)

comp.f$taxon <- "Fungi"


# Fit environmental variables 
# Read in metadata

m <- read.csv(file='infiles/metadata_2024-11-10_DM.csv', head=TRUE)

# create sample col to match sample from matrix
# 12H_13_B_Harvest_Establishment_1
m$SampleName <- paste(m$SiteCode, m$Disturbance, m$DevStage, m$Replicate, sep="_")

# remove samples that are missing from the nmds plot
missing <- setdiff(m$SampleName, env_f$Row.names)

# reorder columns
m.f <- m[!m$SampleName %in% missing,]
m.f <- m.f[,c(21,1:20)]

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
#m.b <- m.b[,c(22:24,26,28)]
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

# # fit environmental variables to NMDS2 & NMDS3
# fit23 <- envfit(nmds3_f, m.f.2, perm = 999, na.rm = TRUE, choices=c(2,3))
# fit.vectors23 <- fit23[[1]]
# fit.pvals23 <- fit23[[1]]$pvals
# fit.df23 <- as.data.frame(fit.vectors23[[1]])
# fit.df23$pvals <- fit.pvals23
# 
# fit.df.sig.f23 <- fit.df23[fit.df23$pvals <0.05,]
# 
# fit.df.sig.f23$x <- 0
# fit.df.sig.f23$y <- 0
# 
# 
# 
# # fit environmental variables to NMDS3 & NMDS1
# fit31 <- envfit(nmds3_f, m.f.2, perm = 999, na.rm = TRUE, choices=c(3,1))
# fit.vectors31 <- fit31[[1]]
# fit.pvals31 <- fit31[[1]]$pvals
# fit.df31 <- as.data.frame(fit.vectors31[[1]])
# fit.df31$pvals <- fit.pvals31
# 
# fit.df.sig.f31 <- fit.df31[fit.df31$pvals <0.05,]
# 
# fit.df.sig.f31$x <- 0
# fit.df.sig.f31$y <- 0

# subset by devstage
gg_f.e <- gg_f[gg_f$DevStage=="E",]
# color by disturbance
chulls.disturbance.f.e <- ddply(gg_f.e, .(Disturbance), function(gg_f.e) gg_f.e[chull(gg_f.e$NMDS1, gg_f.e$NMDS2), ])


# # color by development stage
# chulls.devstage.b <- ddply(gg_b, .(DevStage), function(gg_b) gg_b[chull(gg_b$NMDS1, gg_b$NMDS2), ])

# blues <- brewer.pal(9,"Blues")
# # blues <- blues[c(5,8)]
# blues <- blues[c(4,8)]
# # "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# # [6] "#4292C6" "#2171B5" "#08519C" "#08306B"

# chulls.devstage.b$refcond <- chulls.devstage.b$DevStage
# chulls.devstage.b$refcond <- gsub("E", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("CC", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("LST", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("T", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("M", "1", chulls.devstage.b$refcond)
# 
# 
# gg_b$refcond <- gg_b$DevStage
# gg_b$refcond <- gsub("E", "0", gg_b$refcond)
# gg_b$refcond <- gsub("CC", "0", gg_b$refcond)
# gg_b$refcond <- gsub("LST", "0", gg_b$refcond)
# gg_b$refcond <- gsub("T", "0", gg_b$refcond)
# gg_b$refcond <- gsub("M", "1", gg_b$refcond)

# NMDS plot, add sig env vars
f_fitted_devstage12_e <- ggplot(data=gg_f.e, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.f.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1, 0.5) +
  ylim(-0.75, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.6, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.6, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.f,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

f_fitted_devstage12_e

# subset by devstage
gg_f.cc <- gg_f[gg_f$DevStage=="CC",]
# color by disturbance
chulls.disturbance.f.cc <- ddply(gg_f.cc, .(Disturbance), function(gg_f.cc) gg_f.cc[chull(gg_f.cc$NMDS1, gg_f.cc$NMDS2), ])

# NMDS plot, add sig env vars
f_fitted_devstage12_cc <- ggplot(data=gg_f.cc, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.f.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1, 0.5) +
  ylim(-0.75, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.6, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.6, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.f,
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
f_fitted_devstage12_cc

# subset by devstage
gg_f.t <- gg_f[gg_f$DevStage=="T",]
# color by disturbance
chulls.disturbance.f.t <- ddply(gg_f.t, .(Disturbance), function(gg_f.t) gg_f.t[chull(gg_f.t$NMDS1, gg_f.t$NMDS2), ])

# NMDS plot, add sig env vars
f_fitted_devstage12_t <- ggplot(data=gg_f.t, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.f.t, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.t, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1, 0.5) +
  ylim(-0.75, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.6, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.6, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.f,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
f_fitted_devstage12_t

# subset by devstage
gg_f.lst <- gg_f[gg_f$DevStage=="LST",]
# color by disturbance
chulls.disturbance.f.lst <- ddply(gg_f.lst, .(Disturbance), function(gg_f.lst) gg_f.lst[chull(gg_f.lst$NMDS1, gg_f.lst$NMDS2), ])

# NMDS plot, add sig env vars
f_fitted_devstage12_lst <- ggplot(data=gg_f.lst, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.f.lst, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f.lst, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1, 0.5) +
  ylim(-0.75, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.6, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.6, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.f,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.f)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
f_fitted_devstage12_lst

p4 <- plot_grid(f_fitted_devstage12_e, f_fitted_devstage12_cc,
                f_fitted_devstage12_t, f_fitted_devstage12_lst, nrow=1)

p4








## Arthropod ----


# just arthropods
c <- a[grepl("COI_", a$GlobalESV),]

# keep all Harvest & Fire
# remove Fire if Mature
c <- c[(c$Disturbance=="Harvest") |
         (c$Disturbance=="Fire" ),]
c <- c[!grepl("Fire_Mature", c$Sample),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    66 475

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL
# 66 474

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 66 474

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("outfiles/supporting/Q2_Scree_arth_pairwise.pdf")
# # check dims
# dimcheckMDS(esv.notnull2)
# dev.off()
# # use k=3

# Do 3 dimensional NMDS
nmds3_c <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.1312807 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_arth_pairwise.pdf")
# stressplot(nmds3_c)
# gof <-goodness(nmds3_c)
# gof
# plot(nmds3_c, display = "sites", type="n", main="SSU")
# points(nmds3_c, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.916

# Create grouping matrix for samples by grabbing row names from above matrix
names_c <- data.frame(row.names(esv.notnull2), stringsAsFactors = FALSE)
# Rename the column
names(names_c) <- "sample"
# Copy column to row names
row.names(names_c) <- names_c$sample
# Split first column into their own fields
names_c.1 <- data.frame(names_c, do.call(rbind, strsplit(names_c$sample,'_')), stringsAsFactors = FALSE)
names(names_c.1)[2:5]<-c("SiteCode", "Disturbance", "DevStage", "Replicate")
# Remove first column
names_c.1 <- names_c.1[,-1]
# Grab sites/species scores from NMDS output
df.c <- data.frame(scores(nmds3_c, display = "sites"))
# Put it all in one df for ggplot
gg_c <- merge(df.c, names_c.1, by="row.names")

# create factors
gg_c$Disturbance <- factor(gg_c$Disturbance,
                           levels = c("Harvest", "Fire"),
                           labels = c("Clear Cut Harvest", "Wildfire"))
gg_c$DevStage <- factor(gg_c$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning"),
                        labels = c("E", "CC", "T", "LST"))

# Create metadata from rownames 'sample'
env_c <- gg_c[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_c <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
#bd.layer <- betadisper(sor_b, as.factor(env_b$Layer))
bd.disturbance <- betadisper(sor_c, as.factor(env_c$Disturbance))
bd.devstage <- betadisper(sor_c, as.factor(env_c$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
#anova(bd.layer) # 0.041 *
anova(bd.disturbance) # n/s
anova(bd.devstage) # 1.46e-07 ***
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_arth_pairwise.pdf")
# par(mfrow=c(2,2))
# #boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_c ~ Disturbance*DevStage, data=env_c, permutations=999)
#                       Df SumOfSqs      R2      F Pr(>F)    
# Disturbance           1   0.4587 0.04235 3.7460  0.001 ***
# DevStage              3   2.1167 0.19544 5.7628  0.001 ***
# Disturbance:DevStage  3   1.1538 0.10653 3.1411  0.001 ***
# Residual             58   7.1013 0.65568                  
# Total                65  10.8305 1.00000 

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# subset the distance matrix (don't recreate it)
c.e <- esv.notnull2[grepl("Establishment", rownames(esv.notnull2)),]
sor_c.e <- vegdist(c.e, "bray", binary=TRUE)
env_c.e <- env_c[grepl("E", env_c$DevStage),]
# add comparisons to table 
comp.e <- pairwise.adonis(sor_c.e, env_c.e$Disturbance, p.adjust.m="holm")
comp.e$DevStage <- "Establishment"

c.cc <- esv.notnull2[grepl("CrownClosure", rownames(esv.notnull2)),]
sor_c.cc <- vegdist(c.cc, "bray", binary=TRUE)
env_c.cc <- env_c[grepl("CC", env_c$DevStage),]
# add comparisons to table 
comp.cc <- pairwise.adonis(sor_c.cc, env_c.cc$Disturbance, p.adjust.m="holm")
comp.cc$DevStage <- "CrownClosure"

c.t <- esv.notnull2[grepl("Thinning", rownames(esv.notnull2)),]
sor_c.t <- vegdist(c.t, "bray", binary=TRUE)
env_c.t <- env_c[grepl("T", env_c$DevStage),]
# add comparisons to table 
comp.t <- pairwise.adonis(sor_c.t, env_c.t$Disturbance, p.adjust.m="holm")
comp.t$DevStage <- "Thinning"

c.lst <- esv.notnull2[grepl("LateStageThinning", rownames(esv.notnull2)),]
sor_c.lst <- vegdist(c.lst, "bray", binary=TRUE)
env_c.lst <- env_c[grepl("LST", env_c$DevStage),]
# add comparisons to table 
comp.lst <- pairwise.adonis(sor_c.lst, env_c.lst$Disturbance, p.adjust.m="holm")
comp.lst$DevStage <- "LateStageThinning"

# put it together
comp.c <- rbind(comp.e, comp.cc, comp.t, comp.lst)

comp.c$taxon <- "Arthropoda"


# Fit environmental variables 
# Read in metadata

m <- read.csv(file='infiles/metadata_2024-11-10_DM.csv', head=TRUE)

# create sample col to match sample from matrix
# 12H_13_B_Harvest_Establishment_1
m$SampleName <- paste(m$SiteCode, m$Disturbance, m$DevStage, m$Replicate, sep="_")

# remove samples that are missing from the nmds plot
missing <- setdiff(m$SampleName, env_c$Row.names)

# reorder columns
m.c <- m[!m$SampleName %in% missing,]
m.c <- m.c[,c(21,1:20)]

# average values across soil layers where needed
m.c.2 <- data.frame(m.c %>% 
                      group_by(SampleName) %>% 
                      dplyr::summarize(
                                       pH_mean=mean(pH, na.rm = TRUE), 
                                       TOC_mean=mean(TOC_gkg, na.rm = TRUE), 
                                       TN_mean=mean(TN_gkg, na.rm = TRUE)))
rownames(m.c.2) <- m.c.2$SampleName
m.c.2$SampleName <- NULL

# reorder rows based on rownames(b)
m.c.2 <- m.c.2[match(rownames(esv.notnull2), rownames(m.c.2)), ]
# keep continuous data, focus on "F" part
#m.b <- m.b[,c(22:24,26,28)]
names(m.c.2) <- c( "pH", "TOC", "TN")



# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_c, m.c.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.c <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.c$x <- 0
fit.df.sig.c$y <- 0

# # fit environmental variables to NMDS2 & NMDS3
# fit23 <- envfit(nmds3_c, m.c.2, perm = 999, na.rm = TRUE, choices=c(2,3))
# fit.vectors23 <- fit23[[1]]
# fit.pvals23 <- fit23[[1]]$pvals
# fit.df23 <- as.data.frame(fit.vectors23[[1]])
# fit.df23$pvals <- fit.pvals23
# 
# fit.df.sig.c23 <- fit.df23[fit.df23$pvals <0.05,]
# 
# fit.df.sig.c23$x <- 0
# fit.df.sig.c23$y <- 0
# 
# 
# 
# # fit environmental variables to NMDS3 & NMDS1
# fit31 <- envfit(nmds3_c, m.c.2, perm = 999, na.rm = TRUE, choices=c(3,1))
# fit.vectors31 <- fit31[[1]]
# fit.pvals31 <- fit31[[1]]$pvals
# fit.df31 <- as.data.frame(fit.vectors31[[1]])
# fit.df31$pvals <- fit.pvals31
# 
# fit.df.sig.c31 <- fit.df31[fit.df31$pvals <0.05,]
# 
# fit.df.sig.c31$x <- 0
# fit.df.sig.c31$y <- 0

# subset by devstage
gg_c.e <- gg_c[gg_c$DevStage=="E",]
# color by disturbance
chulls.disturbance.c.e <- ddply(gg_c.e, .(Disturbance), function(gg_c.e) gg_c.e[chull(gg_c.e$NMDS1, gg_c.e$NMDS2), ])


# # color by development stage
# chulls.devstage.b <- ddply(gg_b, .(DevStage), function(gg_b) gg_b[chull(gg_b$NMDS1, gg_b$NMDS2), ])

# blues <- brewer.pal(9,"Blues")
# # blues <- blues[c(5,8)]
# blues <- blues[c(4,8)]
# # "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# # [6] "#4292C6" "#2171B5" "#08519C" "#08306B"

# chulls.devstage.b$refcond <- chulls.devstage.b$DevStage
# chulls.devstage.b$refcond <- gsub("E", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("CC", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("LST", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("T", "0", chulls.devstage.b$refcond)
# chulls.devstage.b$refcond <- gsub("M", "1", chulls.devstage.b$refcond)
# 
# 
# gg_b$refcond <- gg_b$DevStage
# gg_b$refcond <- gsub("E", "0", gg_b$refcond)
# gg_b$refcond <- gsub("CC", "0", gg_b$refcond)
# gg_b$refcond <- gsub("LST", "0", gg_b$refcond)
# gg_b$refcond <- gsub("T", "0", gg_b$refcond)
# gg_b$refcond <- gsub("M", "1", gg_b$refcond)

# NMDS plot, add sig env vars
c_fitted_devstage12_e <- ggplot(data=gg_c.e, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.c.e, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_c.e, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1.25, 0.75) +
  ylim(-0.5, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.6, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.6, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.c,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.c,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.c)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

c_fitted_devstage12_e

# subset by devstage
gg_c.cc <- gg_c[gg_c$DevStage=="CC",]
# color by disturbance
chulls.disturbance.c.cc <- ddply(gg_c.cc, .(Disturbance), function(gg_c.cc) gg_c.cc[chull(gg_c.cc$NMDS1, gg_c.cc$NMDS2), ])

# NMDS plot, add sig env vars
c_fitted_devstage12_cc <- ggplot(data=gg_c.cc, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.c.cc, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_c.cc, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1.25, 0.75) +
  ylim(-0.5, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.7, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '* ', x = Inf, y = 0.7, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.c,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.c,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.c)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
c_fitted_devstage12_cc

# subset by devstage
gg_c.t <- gg_c[gg_c$DevStage=="T",]
# color by disturbance
chulls.disturbance.c.t <- ddply(gg_c.t, .(Disturbance), function(gg_c.t) gg_c.t[chull(gg_c.t$NMDS1, gg_c.t$NMDS2), ])

# NMDS plot, add sig env vars
c_fitted_devstage12_t <- ggplot(data=gg_c.t, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.c.t, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_c.t, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1.25, 0.75) +
  ylim(-0.5, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.7, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '', x = Inf, y = 0.7, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.c,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.c,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.c)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
c_fitted_devstage12_t

# subset by devstage
gg_c.lst <- gg_c[gg_c$DevStage=="LST",]
# color by disturbance
chulls.disturbance.c.lst <- ddply(gg_c.lst, .(Disturbance), function(gg_c.lst) gg_c.lst[chull(gg_c.lst$NMDS1, gg_c.lst$NMDS2), ])

# NMDS plot, add sig env vars
c_fitted_devstage12_lst <- ggplot(data=gg_c.lst, aes(x=NMDS1, y=NMDS2)) +
  geom_polygon(data=chulls.disturbance.c.lst, aes(x=NMDS1, y=NMDS2, fill=Disturbance), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_c.lst, size = 1.5, aes(color = Disturbance, shape=Disturbance, fill = Disturbance)) +
  xlim(-1.25, 0.75) +
  ylim(-0.5, 0.75) +
  annotate(geom = 'text', label = '', x = 0, y = 0.7, hjust = 0.5, vjust = 1) +
  annotate(geom = 'text', label = '', x = Inf, y = 0.7, hjust = 1, vjust = 1) +
  geom_segment(fit.df.sig.c,
               mapping=aes(x=x,y=y,xend=NMDS1*0.50, yend=NMDS2*0.50),
               arrow=arrow(angle=25, length=unit(0.10, "inches")),
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.c,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.c)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")
c_fitted_devstage12_lst

p5 <- plot_grid(c_fitted_devstage12_e, c_fitted_devstage12_cc,
                c_fitted_devstage12_t, c_fitted_devstage12_lst, nrow=1)
p5


# summarize results for pairwise comparisons of beta diversity
comp.all.taxa2 <- rbind(comp.b, comp.f, comp.c)

comp.all.taxa2$pairs <- gsub("Clear Cut Harvest vs Wildfire", "Wildfire vs Clear Cut Harvest", comp.all.taxa2$pairs) 

# Create factors
comp.all.taxa2$DevStage <- factor(comp.all.taxa2$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning"),
                                  labels=c("CCH.E vs W.E", "CCH.CC vs W.CC", "CCH.EST vs W.EST", "CCH.LST vs W.LST"))

comp.all.taxa2$taxon <- factor(comp.all.taxa2$taxon, levels=rev(c("Bacteria", "Fungi", "Arthropoda")),
                               labels=rev(c("B", "F", "A")))

g3 <- ggplot(comp.all.taxa2, aes(x=p.adjusted, y=DevStage, shape=taxon)) + 
  geom_text(aes(label=taxon), size = 3, position = position_jitterdodge(seed=3)) +
  labs(x = "Holm adjusted p-values", y = "") +
  geom_hline(yintercept=0.05, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  scale_shape_manual(values = c("B", "F", "A")) +
  theme_bw() + 
  theme(
    text = element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(hjust = 0.5, vjust = 0),
    axis.text = element_text(size=8),
    axis.title.x = element_text(vjust = -0.75),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom",
    strip.text = element_text(size = 8)) 
g3

ggsave("outfiles/Fig3_BetaDiv_Convergence_pvals.pdf", g3, height = 4, width = 4)




# NMDS recovery ----


## Bacteria ----


# just bacteria
b <- a[grepl("16S_", a$GlobalESV),]

# keep all Harvest
# keep Fire if Mature
b <- b[(b$Disturbance=="Harvest") |
         (b$Disturbance=="Fire" & b$DevStage=="Mature"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
b.esv <- reshape2::dcast(b, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 481

# move sample to rownames then delete
rownames(b.esv) <- b.esv$Sample
b.esv$Sample <- NULL
# 42 480

#remove any columns with only zeros
esv.notnull <- b.esv[,colSums(b.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 480

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("outfiles/supporting/Q2_Scree_bacteria.pdf")
# # check dims
# dimcheckMDS(esv.notnull2)
# dev.off()
# # use k=3

# Do 3 dimensional NMDS
nmds3_b <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.08029379 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_bacteria.pdf")
# stressplot(nmds3_b)
# gof <-goodness(nmds3_b)
# gof
# plot(nmds3_b, display = "sites", type="n", main="SSU")
# points(nmds3_b, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.978

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
                           levels = c("Harvest", "Fire"),
                           labels = c("Clear Cut Harvest", "Wildfire"))
gg_b$DevStage <- factor(gg_b$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning","Mature"),
                        labels = c("CCH.E", "CCH.CC", "CCH.EST", "CCH.LST","W.M"))

# Create metadata from rownames 'sample'
env_b <- gg_b[,c(1,5:8)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_b <- vegdist(esv.notnull2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.disturbance <- betadisper(sor_b, as.factor(env_b$Disturbance))
bd.devstage <- betadisper(sor_b, as.factor(env_b$DevStage))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
anova(bd.disturbance) # n/s
anova(bd.devstage) # n/s
# PERMANOVA differences due to differences between factors (not within),
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_bacteria.pdf")
# par(mfrow=c(2,2))
# #boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_b ~ Disturbance*DevStage, data=env_b, permutations=999)
#             Df SumOfSqs      R2      F Pr(>F)    
# Disturbance  1  0.03831 0.07215 4.0902  0.007 ** 
# DevStage     3  0.14609 0.27514 5.1989  0.001 ***
# Residual    37  0.34657 0.65271                  
# Total       41  0.53097 1.00000  

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# subset the distance matrix (don't recreate it)

comp <- pairwise.adonis(sor_b, env_b$DevStage, p.adjust.m="holm")
comp$taxon <- "Bacteria"


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
names(m.b.2) <- c("pH", "TOC", "TN")

# fit environmental variables to NMDS1 & NMDS2
fit <- envfit(nmds3_b, m.b.2, perm = 999, na.rm = TRUE)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.b <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.b$x <- 0
fit.df.sig.b$y <- 0

# color by development stage
chulls.devstage.b <- ddply(gg_b, .(DevStage), function(gg_b) gg_b[chull(gg_b$NMDS1, gg_b$NMDS2), ])

blues <- brewer.pal(9,"Blues")
# blues <- blues[c(5,8)]
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6"
# [6] "#4292C6" "#2171B5" "#08519C" "#08306B"

# reds <- brewer.pal(9,"Reds")
# [1] "#FFF5F0" "#FEE0D2" "#FCBBA1" "#FC9272" "#FB6A4A"
# [6] "#EF3B2C" "#CB181D" "#A50F15" "#67000D"
# reds <- reds[c(5,8)]

cols <- c(blues[c(4,5,6,7)], "#df0000")

# NMDS plot, add sig env vars
b_fitted_devstage12 <- ggplot(data=gg_b, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.devstage.b, aes(x=NMDS1, y=NMDS2, fill=DevStage), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b, size = 1.5, aes(color = DevStage, shape=DevStage, fill = DevStage)) +
  geom_segment(fit.df.sig.b, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("     Bacteria") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

b_fitted_devstage12

# NMDS plot, add sig env vars
b_fitted_devstage12.tmp <- ggplot(data=gg_b, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.devstage.b, aes(x=NMDS1, y=NMDS2, fill=DevStage), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_b, size = 1.5, aes(color = DevStage, shape=DevStage, fill = DevStage)) +
  geom_segment(fit.df.sig.b, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.5, yend=NMDS2*0.5), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.b,
            mapping=aes(x=NMDS1*0.5, y=NMDS2*0.5, label=rownames(fit.df.sig.b)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25))+
  ggtitle("     Bacteria") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "bottom")+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) 
b_fitted_devstage12.tmp

l4 <- get_legend(b_fitted_devstage12.tmp)

p6 <- b_fitted_devstage12







## Fungi ----


# just fungi
f <- a[grepl("ITS_", a$GlobalESV),]

# keep all Harvest
# keep Fire if Mature
f <- f[(f$Disturbance=="Harvest") |
         (f$Disturbance=="Fire" & f$DevStage=="Mature"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
f.esv <- reshape2::dcast(f, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 634

# move sample to rownames then delete
rownames(f.esv) <- f.esv$Sample
f.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- f.esv[,colSums(f.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 633

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("outfiles/supporting/Q2_Scree_fungi.pdf")
# # check dims
# dimcheckMDS(esv.notnull2)
# dev.off()
# # use k=3

# Do 3 dimensional NMDS
nmds3_f <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.1165377 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_fungi.pdf")
# stressplot(nmds3_f)
# gof <-goodness(nmds3_f)
# gof
# plot(nmds3_f, display = "sites", type="n", main="SSU")
# points(nmds3_f, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.926

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
                           levels = c("Harvest", "Fire"),
                           labels = c("Clear Cut Harvest", "Wildfire"))
gg_f$DevStage <- factor(gg_f$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning","Mature"),
                        labels = c("CCH.E", "CCH.CC", "CCH.EST", "CCH.LST", "W.M"))

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
anova(bd.disturbance) # n/s
anova(bd.devstage) # 0.004793 **
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_fungi.pdf")
# par(mfrow=c(2,2))
# # boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# plot distance to centroid for fire and harvest separately for devstage
fire_filter <- grepl("Fire", names(bd.devstage$distances))
distances <- bd.devstage$distances[fire_filter]
group <- bd.devstage$group[fire_filter]
f.bd.devstage.f <- data.frame(distances=distances, group=group)
f.bd.devstage.f$disturbance <- "Fire"

harvest_filter <- grepl("Harvest", names(bd.devstage$distances))
distances <- bd.devstage$distances[harvest_filter]
group <- bd.devstage$group[harvest_filter]
f.bd.devstage.h <- data.frame(distances=distances, group=group)
f.bd.devstage.h$disturbance <- "Harvest"

# put it together for ggplot
f.bd.devstage <- rbind(f.bd.devstage.f, f.bd.devstage.h)
f.bd.devstage$amplicon <- "ITS"



# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_f ~ Disturbance*DevStage, data=env_f, permutations=999)
#             Df SumOfSqs      R2      F Pr(>F)    
# Disturbance  1   0.3627 0.07512 3.7690  0.001 ***
# DevStage     3   0.9050 0.18742 3.1345  0.001 ***
# Residual    37   3.5609 0.73746                  
# Total       41   4.8286 1.00000  

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage and layer
# subset the distance matrix (don't recreate it)

comp2 <- pairwise.adonis(sor_f, env_f$DevStage, p.adjust.m="holm")
comp2$taxon <- "Fungi"


# Fit environmental variables for diatom
# Read in metadata

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
names(m.f.2) <- c("pH", "TOC", "TN")

# fit environmental variables to NMDS3 & NMDS1
fit31 <- envfit(nmds3_f, m.f.2, perm = 999, na.rm = TRUE, choices=c(3,1))
fit.vectors31 <- fit31[[1]]
fit.pvals31 <- fit31[[1]]$pvals
fit.df31 <- as.data.frame(fit.vectors31[[1]])
fit.df31$pvals <- fit.pvals31

fit.df.sig.f31 <- fit.df31[fit.df31$pvals <0.05,]

fit.df.sig.f31$x <- 0
fit.df.sig.f31$y <- 0

# color by development stage
chulls.devstage.f31 <- ddply(gg_f, .(DevStage), function(gg_f) gg_f[chull(gg_f$NMDS3, gg_f$NMDS1), ])

# NMDS plot, add sig env vars
f_fitted_devstage31 <- ggplot(data=gg_f, aes(x=NMDS3, y=NMDS1)) + 
  geom_polygon(data=chulls.devstage.f31, aes(x=NMDS3, y=NMDS1, fill=DevStage), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_f, size = 1.5, aes(color = DevStage, shape=DevStage, fill = DevStage)) +
  # xlim(-0.8,0.8) +
  # ylim(-0.8,0.8) +
  geom_segment(fit.df.sig.f31, 
               mapping=aes(x=x,y=y,xend=NMDS3*0.5, yend=NMDS1*0.5), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.f31,
            mapping=aes(x=NMDS3*0.5, y=NMDS1*0.5, label=rownames(fit.df.sig.f31)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("     Fungi") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

f_fitted_devstage31

p7 <- f_fitted_devstage31






## Arthropod ----


# just arthropoda
c <- a[grepl("COI_", a$GlobalESV),]

# keep all Harvest
# keep Fire if Mature
c <- c[(c$Disturbance=="Harvest") |
         (c$Disturbance=="Fire" & c$DevStage=="Mature"),]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
c.esv <- reshape2::dcast(c, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    42 445

# move sample to rownames then delete
rownames(c.esv) <- c.esv$Sample
c.esv$Sample <- NULL

#remove any columns with only zeros
esv.notnull <- c.esv[,colSums(c.esv) !=0]

#remove any rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]
# 42 444

# convert to presence-absence
esv.notnull2[esv.notnull2>0] <- 1

# # Scree plots to determine number of dimensions to use for NMDS
# pdf("outfiles/supporting/Q2_Scree_arthropoda.pdf")
# # check dims
# dimcheckMDS(esv.notnull2)
# dev.off()
# # use k=3

# Do 3 dimensional NMDS
nmds3_a <- metaMDS(esv.notnull2, k=3, trymax=100)
# stress = 0.1395952 

# # Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
# pdf("outfiles/supporting/Q2_stressplot_arthropoda.pdf")
# stressplot(nmds3_a)
# gof <-goodness(nmds3_a)
# gof
# plot(nmds3_a, display = "sites", type="n", main="SSU")
# points(nmds3_a, display="sites",cex=2*gof/mean(gof))
# dev.off()
# linear R2 = 0.891


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
                           levels = c("Harvest", "Fire"),
                           labels = c("Clear Cut Harvest", "Wildfire"))
gg_a$DevStage <- factor(gg_a$DevStage,
                        levels = c("Establishment","CrownClosure","Thinning","LateStageThinning", "Mature"),
                        labels = c("CCH.E", "CCH.CC", "CCH.EST", "CCH.LST", "W.M"))

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
#anova(bd.layer) # n/s
anova(bd.disturbance) # 0.003445 **
anova(bd.devstage) # 2.967e-05 ***
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

# pdf("outfiles/supporting/Q2_BetaDispersion_arthropoda.pdf")
# par(mfrow=c(2,2))
# # boxplot(bd.layer, main="Layer")
# boxplot(bd.disturbance, main="Disturbance")
# boxplot(bd.devstage, main="DevStage")
# dev.off()

# plot distance to centroid for fire and harvest separately for devstage
fire_filter <- grepl("Fire", names(bd.devstage$distances))
distances <- bd.devstage$distances[fire_filter]
group <- bd.devstage$group[fire_filter]
a.bd.devstage.f <- data.frame(distances=distances, group=group)
a.bd.devstage.f$disturbance <- "Fire"

harvest_filter <- grepl("Harvest", names(bd.devstage$distances))
distances <- bd.devstage$distances[harvest_filter]
group <- bd.devstage$group[harvest_filter]
a.bd.devstage.h <- data.frame(distances=distances, group=group)
a.bd.devstage.h$disturbance <- "Harvest"

# put it together for ggplot
a.bd.devstage <- rbind(a.bd.devstage.f, a.bd.devstage.h)
a.bd.devstage$amplicon <- "COI"




# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis2(sor_a ~ Disturbance*DevStage, data=env_a, permutations=999)
#             Df SumOfSqs      R2      F Pr(>F)    
# Disturbance  1   0.3846 0.06501 3.2909  0.001 ***
# DevStage     3   1.2075 0.20411 3.4443  0.001 ***
# Residual    37   4.3237 0.73088                  
# Total       41   5.9158 1.00000

# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage
# separate out datasets by devstage then do pairwise comparisons of disturbance
# do pairwise comparisons of disturbance and devstage and layer
# subset the distance matrix (don't recreate it)

comp3 <- pairwise.adonis(sor_a, env_a$DevStage, p.adjust.m="holm")
comp3$taxon <- "Arthropoda"


# Read in metadata

# m <- read.csv(file='infiles/metadata_2024-11-10_DM.csv', head=TRUE)

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
                        #Temp_mean=mean(MeanAnnualSeasonalTemp_C),
                         #              Precip_mean=mean(MeanAnnualSeasonalPrecip_mm),
                                       pH_mean=mean(pH, na.rm = TRUE), 
                                       TOC_mean=mean(TOC_gkg, na.rm = TRUE), 
                                       TN_mean=mean(TN_gkg, na.rm = TRUE)))
rownames(m.a.2) <- m.a.2$SampleName
m.a.2$SampleName <- NULL

# reorder rows based on rownames(b)
m.a.2 <- m.a.2[match(rownames(esv.notnull2), rownames(m.a.2)), ]
# keep continuous data, focus on "F" part
names(m.a.2) <- c("pH", "TOC", "TN")

# fit environmental variables to NMDS3 & NMDS1
fit31 <- envfit(nmds3_a, m.a.2, perm = 999, na.rm = TRUE, choices=c(3,1))
fit.vectors31 <- fit31[[1]]
fit.pvals31 <- fit31[[1]]$pvals
fit.df31 <- as.data.frame(fit.vectors31[[1]])
fit.df31$pvals <- fit.pvals31

fit.df.sig.a31 <- fit.df31[fit.df31$pvals <0.05,]

fit.df.sig.a31$x <- 0
fit.df.sig.a31$y <- 0

# color by development stage
chulls.devstage.a31 <- ddply(gg_a, .(DevStage), function(gg_a) gg_a[chull(gg_a$NMDS3, gg_a$NMDS1), ])

# NMDS plot, add sig env vars
a_fitted_devstage31 <- ggplot(data=gg_a, aes(x=NMDS3, y=NMDS1)) + 
  geom_polygon(data=chulls.devstage.a31, aes(x=NMDS3, y=NMDS1, fill=DevStage), alpha=0.25, show.legend = FALSE) +
  geom_point(data=gg_a, size = 1.5, aes(color = DevStage, shape=DevStage, fill = DevStage)) +
  geom_segment(fit.df.sig.a31, 
               mapping=aes(x=x,y=y,xend=NMDS3*0.5, yend=NMDS1*0.5), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.a31,
            mapping=aes(x=NMDS3*0.5, y=NMDS1*0.5, label=rownames(fit.df.sig.a31)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle("     Arthropoda") +
  theme_bw() +
  theme(
    text = element_text(size=11),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "none")

a_fitted_devstage31

p8 <- a_fitted_devstage31




# # keep bacteria, fungi, and arthropod pairwise permanova comparisons separate
# 
# # summarize results for pairwise comparisons of beta diversity
# #comp.all.taxa <- rbind(comp, comp2, comp3)
# 
# # Bacteria pairwise permanova comparisons (checking convergence at each stage)
# # only keep comparisons of interest
# comp <- comp[(comp$pairs=="CCH.E vs W.M" |
#               comp$pairs=="CCH.CC vs W.M" |
#               comp$pairs=="CCH.EST vs W.M" |
#               comp$pairs=="W.M vs CCH.LST"), ]
# 
# # Create factors
# comp$pairs <- factor(comp$pairs, 
#                      levels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "W.M vs CCH.LST")),
#                      labels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "CCH.LST vs W.M")))
# 
# comp$taxon <- factor(comp$taxon, 
#                      levels=rev(c("Bacteria", "Fungi", "Arthropoda")),
#                      labels=rev(c("B", "F", "A")))
# 
# g4 <- ggplot(comp, aes(x=p.adjusted, y=pairs)) + 
#   geom_point(shape = "*", size=6) +
# #  geom_text(aes(label=taxon), size = 2.5, position = position_jitterdodge(seed=3)) +
#   labs(x = "p-values", y = "") +
#   geom_hline(yintercept=0.05, linetype="dashed") +
#   geom_vline(xintercept = 0.05, linetype="dashed") +
# #  scale_shape_manual(values = c("B", "F", "A")) +
#   ggtitle("Bacteria") +
#   theme_bw() + 
#   theme(
#     text = element_text(size=11),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.text.x = element_text(hjust = 0.5, vjust = 0),
#     axis.text = element_text(size=8),
#     axis.title.x = element_text(vjust = -0.75),
#     legend.title = element_blank(),
#     legend.text = element_text(size=11),
#     legend.position = "bottom",
#     strip.text = element_text(size = 11)) 
# g4
# 
# 
# 
# # Fungi pairwise permanova comparisons (checking convergence at each stage)
# # only keep comparisons of interest
# comp2 <- comp2[(comp2$pairs=="CCH.E vs W.M" |
#                 comp2$pairs=="CCH.CC vs W.M" |
#                 comp2$pairs=="CCH.EST vs W.M" |
#                 comp2$pairs=="W.M vs CCH.LST"), ]
# 
# # Create factors
# comp2$pairs <- factor(comp2$pairs, 
#                      levels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "W.M vs CCH.LST")),
#                      labels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "CCH.LST vs W.M")))
# 
# comp2$taxon <- factor(comp2$taxon, 
#                      levels=rev(c("Bacteria", "Fungi", "Arthropoda")),
#                      labels=rev(c("B", "F", "A")))
# 
# g5 <- ggplot(comp2, aes(x=p.adjusted, y=pairs)) + 
#   geom_point(shape = "*", size=6) +
#   #  geom_text(aes(label=taxon), size = 2.5, position = position_jitterdodge(seed=3)) +
#   labs(x = "p-values", y = "") +
#   geom_hline(yintercept=0.05, linetype="dashed") +
#   geom_vline(xintercept = 0.05, linetype="dashed") +
#   #  scale_shape_manual(values = c("B", "F", "A")) +
#   ggtitle("Fungi") +
#   theme_bw() + 
#   theme(
#     text = element_text(size=11),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.text.x = element_text(hjust = 0.5, vjust = 0),
#     axis.text = element_text(size=8),
#     axis.title.x = element_text(vjust = -0.75),
#     legend.title = element_blank(),
#     legend.text = element_text(size=11),
#     legend.position = "bottom",
#     strip.text = element_text(size = 11)) 
# g5
# 
# 
# # Arthropod pairwise permanova comparisons (checking convergence at each stage)
# # only keep comparisons of interest
# comp3 <- comp3[(comp3$pairs=="CCH.E vs W.M" |
#                   comp3$pairs=="CCH.CC vs W.M" |
#                   comp3$pairs=="CCH.EST vs W.M" |
#                   comp3$pairs=="W.M vs CCH.LST"), ]
# 
# # Create factors
# comp3$pairs <- factor(comp3$pairs, 
#                       levels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "W.M vs CCH.LST")),
#                       labels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "CCH.LST vs W.M")))
# 
# comp3$taxon <- factor(comp3$taxon, 
#                       levels=rev(c("Bacteria", "Fungi", "Arthropoda")),
#                       labels=rev(c("B", "F", "A")))
# 
# g6 <- ggplot(comp3, aes(x=p.adjusted, y=pairs)) + 
#   geom_point(shape = "*", size=6) +
#   #  geom_text(aes(label=taxon), size = 2.5, position = position_jitterdodge(seed=3)) +
#   labs(x = "p-values", y = "") +
#   geom_hline(yintercept=0.05, linetype="dashed") +
#   geom_vline(xintercept = 0.05, linetype="dashed") +
#   #  scale_shape_manual(values = c("B", "F", "A")) +
#   ggtitle("Arthropoda") +
#   theme_bw() + 
#   theme(
#     text = element_text(size=11),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     axis.text.x = element_text(hjust = 0.5, vjust = 0),
#     axis.text = element_text(size=8),
#     axis.title.x = element_text(vjust = -0.75),
#     legend.title = element_blank(),
#     legend.text = element_text(size=11),
#     legend.position = "bottom",
#     strip.text = element_text(size = 11)) 
# g6
# 




# OKAY, keep together but use facet wrap to keep taxa separate
# summarize results for pairwise comparisons of beta diversity
comp.all.taxa <- rbind(comp, comp2, comp3)

# only keep comparisons of interest
comp.all.taxa <- comp.all.taxa[(comp.all.taxa$pairs=="CCH.E vs W.M" |
                                  comp.all.taxa$pairs=="CCH.CC vs W.M" |
                                  comp.all.taxa$pairs=="CCH.EST vs W.M" |
                                  comp.all.taxa$pairs=="W.M vs CCH.LST"), ]

# Create factors
comp.all.taxa$pairs <- factor(comp.all.taxa$pairs, 
                              levels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "W.M vs CCH.LST")),
                              labels=rev(c("CCH.E vs W.M", "CCH.CC vs W.M", "CCH.EST vs W.M", "CCH.LST vs W.M")))

comp.all.taxa$taxon <- factor(comp.all.taxa$taxon, levels=c("Bacteria", "Fungi", "Arthropoda"),
                              labels=c("Bacteria", "Fungi", "Arthropoda"))

g7 <- ggplot(comp.all.taxa, aes(x=p.adjusted, y=pairs, shape=taxon)) + 
  geom_point(shape="*", size = 6) +
  #geom_text(aes(label=taxon), size = 2.5, position = position_jitterdodge(seed=3)) +
  labs(x = "p-values", y = "") +
  geom_hline(yintercept=0.05, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  #scale_shape_manual(values = c("B", "F", "A")) +
  ggtitle("") +
  facet_wrap(~taxon) +
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    text = element_text(size=11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(hjust = 0.5, vjust = 0),
    axis.text = element_text(size=8),
    axis.title.x = element_text(vjust = -0.75),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    legend.position = "bottom",
    strip.text = element_text(size = 11, hjust=0))
g7










# put together a plate using richness data
top <- plot_grid(g1, g2, nrow=1, rel_widths = c(1,1), labels=c("(a)","(b)"))
bottom <- plot_grid(l1, l2, nrow=1, rel_widths = c(1,1))
plate <- plot_grid(top, bottom, ncol=1, rel_heights = c(1,0.05))

ggsave("outfiles/Fig2_Richness_Convergence_Recovery.pdf", plate, height = 8, width = 8)

# put together a plate using pairwise-distances, community dissimilarity data
# convergence

partA.plots <- plot_grid(p3, p4, p5, ncol=1)
partA.legend <- plot_grid(l3)
partA <- plot_grid(partA.plots, partA.legend, ncol=1, rel_heights = c(3, 0.1), labels=c("(a)"))


# recovery

partB.plots <- plot_grid(p6, p7, p8, nrow=1, labels = c("(b)"), label_x = 0)
partB.legend <- plot_grid(l4)
partB <- plot_grid(partB.plots, partB.legend, ncol=1, rel_heights=c(1,0.1))

partC <- plot_grid(g7, nrow=1, labels=c("(c)"))


plate2 <- plot_grid(partA, partB, partC, ncol=1, rel_heights = c(3,1,1.4))

ggsave("outfiles/Fig3_BetaDiv_Convergence_Recovery.pdf", plate2, height = 8, width = 8)



# Teresita M. Porter, June 30, 2023
# Rarefaction analyses

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rarecurve
library(purrr) # for map_dfr
library(ggplot2) # ggplot
library(scales) # comma
library(gridExtra) # grid.arrange
library(cowplot) # get_legend
library(RColorBrewer)

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# read in merged taxonomy and metadata from tax_summary_2022-04-27.R
# limited to fire disturbance
a <- read.csv("outfiles/tax_meta1.csv", header=TRUE, stringsAsFactors = FALSE)

# pool results from F230 and BE
a$Amplicon <- gsub("F230", "COI", a$Amplicon)
a$Amplicon <- gsub("BE", "COI", a$Amplicon)
a$GlobalESV <- gsub("F230", "COI", a$GlobalESV)
a$GlobalESV <- gsub("BE", "COI", a$GlobalESV)

# total number of ESVs
length(unique(a$GlobalESV))
# 114776

# total number of reads
# one sample as no sequencing data, ignore it
sum(as.numeric(a$ESVsize), na.rm = TRUE)
# 10,077,185

# Create new column for dcast
a$Sample <- paste(a$SiteCode, a$SiteCodeDM, a$Layer, a$DevStage, a$Replicate, a$Amplicon, sep="_")

# pivot to make sample x ESV matrix (keep each individual sample separate)
A.esv <- reshape2::dcast(a, Sample ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
# [1]    814 114777

# move sample to rownames then delete
rownames(A.esv) <- A.esv$Sample
A.esv$Sample <- NULL

#remove columns with only zeros ****HERE
esv.notnull <- A.esv[, which(colSums(A.esv) !=0)]

#remove rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]

# check distribution of library sizes - may need to toss samples with < 1000 reads
sums <- rowSums(esv.notnull2)
hist(sums)

# one sample with only 555 reads, E90_10_B_Mature_2_ITS
# keep it in because all the curves appear to plateau, even when including the rare ESVs

# #calculate 15th percentile for rrarefy function - may need to analyze each marker separately to avoid such low cutoffs
# esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
# esv.percentile
# # 15% 
# # 6939.75 

# set random seed
set.seed(1234)

# Do rarefaction with pkg 'vegan'
rarecurveout <- rarecurve2(esv.notnull2, 
                           # sample=esv.percentile, 
                           step=500, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rarecurveout, function(x){
  z <- as.data.frame(x)
  z <- data.frame(OTU = z[,1], raw.read = rownames(z))
  z$raw.read <- as.numeric(gsub("N", "",  z$raw.read))
  return(z)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(esv.notnull2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:9]<-c("SiteCode","SiteCode_DM","Layer", "DevStage","Replicate","Amplicon")


# create factors
rare.df$Layer <- factor(rare.df$Layer, levels=c("B", "O", "M"), labels=c("Bryophyte/Litter", "Organic", "Mineral"))
# rare.df$Disturbance <- factor(rare.df$Disturbance, levels=c("Fire", "Harvest", "Salvage"),
#                                labels=c("Wildfire", "Harvest", "Salvage"))
rare.df$DevStage <- factor(rare.df$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                           labels=c("E","CC","EST","LST","M"))
rare.df$Replicate <- factor(rare.df$Replicate, levels=c("1", "2", "3"))
rare.df$Amplicon <- factor(rare.df$Amplicon, levels=c("16S", "ITS", "COI"))

# reds <- brewer.pal(9,"Reds")
# reds <- reds[c(4,6,8)]

# color by amplicon
p1 <- ggplot(data = rare.df) +
  ggtitle("Amplicon") +
  labs(x="Reads", y="SVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Amplicon), size=0.05) +
  scale_x_continuous(label = comma) +
#  scale_color_manual(values = reds) +
#  scale_fill_manual(values = reds) +
  theme_bw() + 
  theme(text = element_text(size=10),
        plot.title = element_text(size=10, hjust=0.01),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        axis.title.x = element_text(vjust = -0.75),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=2)))
p1
ggsave("outfiles/FigS1_rarefaction_amplicon.pdf", p1, width = 4, height = 4)




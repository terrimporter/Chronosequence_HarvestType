# Teresita M. Porter, Apr. 20, 2025
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
library(grid) #viewport

# Read in taxonomic assignments ----

# grab results from MetaWorks1.9.2
BE <- read.csv("infiles/BE_results.csv", header=TRUE, stringsAsFactors = FALSE)
# 4896

F230 <- read.csv("infiles/F230_results.csv", header=TRUE, stringsAsFactors = FALSE)
# 40967

# grab results from MetaWorks1.9.3
ITS <- read.csv("infiles/ITS_results.csv", header=TRUE, stringsAsFactors = FALSE)
# 106141

SSU <- read.csv("infiles/16S_results.csv", header=TRUE, stringsAsFactors = FALSE)
# 340643

# Fix up the fields
# put COI together
COI <- rbind(BE, F230)

# fix up header names
names(COI)[4] <- "seq"
names(ITS)[4] <- "seq"
names(SSU)[4] <- "seq"

# SuperKingdom ~ Domain
names(COI)[9:11] <- c("SuperKingdomDomain", "SuperKingdomDomainRank", "skdBP")
names(SSU)[6:8] <- c("SuperKingdomDomain", "SuperKingdomDomainRank", "skdBP")

# add missing columns
ITS$SuperKingdomDomain <- ""
ITS$SuperKingdomDomainRank <- ""
ITS$skdBP <- ""

SSU$Root <- ""
SSU$RootRank <- ""
SSU$rBP <- ""

SSU$Kingdom <- ""
SSU$KingdomRank <- ""
SSU$kBP <- ""

SSU$Species <- ""
SSU$SpeciesRank <- ""
SSU$sBP <- ""

# Just keep ITS seqs for Fungi
ITS <- ITS[ITS$Kingdom=="k__Fungi",]

# Just keep 16S seqs for Bacteria
SSU <- SSU[SSU$SuperKingdomDomain=="Bacteria",]

# put results together
A <- rbind(COI, ITS, SSU)
# 489,014

# split out amplicon from SampleName
A$Amplicon <- A$SampleName

# grab last part of the SampleName field as the Amplicon
A$Amplicon <- gsub("^.*_", "", A$SampleName)

# Append amplicon to end of global ESV so that ZotuIDs remain unique when datasets are combined
A$GlobalESV <- paste(A$Amplicon, A$GlobalESV, sep="_")

# Remove Amplicon from SampleName
A$SampleName <- unlist(strsplit(A$SampleName, paste("_", A$Amplicon,"$", sep="")))
# remove repetitive part
A$SampleName <- gsub("GRDI_ECO_SM_LV_", "", A$SampleName)
A$SampleName <- gsub("GRDI_ECO_SM_LV", "", A$SampleName)
A$SampleName <- gsub("GRDI_ECO_\\d+_SM_LV_", "", A$SampleName)

# split out IlluminaSample from SampleName
A$IlluminaSample <- A$SampleName

# grab last part of SampleName as the IlluminaSample
A$IlluminaSample <- gsub("^.*_", "", A$IlluminaSample)

# Remove IlluminaSample from SampleName
A$SampleName <- unlist(strsplit(A$SampleName, paste("_", A$IlluminaSample,"$", sep="")))

# Remove Amplicon from SampleName (again)
A$SampleName <- unlist(strsplit(A$SampleName, paste("_", A$Amplicon,"$", sep="")))





# Bring in metadata ----
M <- read.csv("infiles/metadata_2024-11-10_DM.csv", header = TRUE, stringsAsFactors = FALSE)

# Change dashes to underscores to match what's in the taxonomy file
M$SampleIdentifier <- gsub("-","_", M$SampleIdentifier)

# Sanity check, ensure SampleNames from seqdata eq SampleIdentifiers from metadata
A.sorted <- unique(sort(A$SampleName))
# 276
M.sorted <- unique(sort(M$SampleIdentifier))
# 277 one seq file is missing
setdiff(A.sorted, M.sorted) 
setdiff(M.sorted, A.sorted) 
# [1] "CB37_20170705_F_3F", sequences missing for this sample in the metadata







# Create tax_meta.csv ####

# Add metadata to seqfile
# only keep if Fire disturbance for this study
AM <- merge(A, M, by.x="SampleName", by.y="SampleIdentifier", all.y = TRUE)

# remove sample with no reads
AM <- AM[!AM$SampleName=="CB37_20170705_F_3F",]

# check if outfiles dir exists, create if needed
ifelse(!dir.exists(file.path("outfiles")), dir.create(file.path("outfiles")), FALSE)
# save this file for future analyses
write.csv(AM, "outfiles/tax_meta1.csv", quote=FALSE, row.names = FALSE)







# remove ESVs that represent < 0.01% of all reads ----
A.2 <- data.frame(A %>% group_by(GlobalESV) %>% summarize(sum=sum(ESVsize)))

totalESVs <- length(unique(A.2$GlobalESV))
# 114,968

totalReads <- sum(A.2$sum)
# 10,077,185
cutoff <- totalReads*0.0001
# 1007.719

ESVlist <- A.2$GlobalESV[A.2$sum>=cutoff]
# keep 1671

# filter original results to remove the rare tail of ESVs
A.3 <- A[A$GlobalESV %in% ESVlist,]

# Add metadata to seqfile
# only keep if Fire disturbance for this study
AM.2 <- merge(A.3, M, by.x="SampleName", by.y="SampleIdentifier", all.y = TRUE)

# remove sample with no reads
AM.2 <- AM.2[!AM.2$SampleName=="CB37_20170705_F_3F",]

# check if outfiles dir exists, create if needed
ifelse(!dir.exists(file.path("outfiles")), dir.create(file.path("outfiles")), FALSE)
# save this file for future analyses
write.csv(AM.2, "outfiles/tax_meta2.csv", quote=FALSE, row.names = FALSE)

# Figure out remaining SVs and reads
A.4 <- data.frame(AM.2 %>% group_by(GlobalESV) %>% summarize(sum=sum(ESVsize)))

totalESVs <- length(unique(A.4$GlobalESV))
# 1671

totalReads <- sum(A.4$sum)
# 6,642,061




# Table S2 raw sequencing stats ----




AM %>% group_by(Amplicon) %>% dplyr::summarize(sum(ESVsize))
# Amplicon `sum(ESVsize)`
# <chr>             <int>
#   1 16S             3286199
# 2 BE                97266
# 3 F230            2741702
# 4 ITS             3952018

sum(AM$ESVsize)
# 10,077,185


AM %>% group_by(Amplicon) %>% dplyr::summarize(length(unique(GlobalESV)))
# Amplicon `length(unique(GlobalESV))`
# <chr>                          <int>
#   1 16S                            78091
# 2 BE                              1260
# 3 F230                            7733
# 4 ITS                            27884

length(unique(AM$GlobalESV))
# 114,968





# Table S2 sequencing stats more rare removed ----




AM.2 %>% group_by(Amplicon) %>% dplyr::summarize(sum(ESVsize))
# Amplicon `sum(ESVsize)`
# <chr>             <int>
#   1 16S             1612891
# 2 BE                47754
# 3 F230            2100686
# 4 ITS             2880730

sum(AM.2$ESVsize)
# 6,642,061


AM.2 %>% group_by(Amplicon) %>% dplyr::summarize(length(unique(GlobalESV)))
# Amplicon `length(unique(GlobalESV))`
# <chr>                          <int>
#   1 16S                              480
# 2 BE                                17
# 3 F230                             474
# 4 ITS                              700

length(unique(AM.2$GlobalESV))
# 1671



# Teresita M. Porter, April 20, 2025
# Summarize site characteristics for Fig 1 and Table S1

library(dplyr)



# All sites ----

m <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# create factors
m$DevStage <- factor(m$DevStage, levels=c("Establishment",
                                          "CrownClosure","Thinning",
                                          "LateStageThinning","Mature"),
                     labels=c("E","CC","EST","LST","M"))

## Stand age as of 2017 ----
m %>% 
  group_by(DevStage) %>% 
  dplyr::summarize(range=range(StandAge2017))

# DevStage range
# <fct>    <int>
# 1 E           12
# 2 E           15
# 3 CC          20
# 4 CC          25
# 5 EST         37
# 6 EST         41
# 7 LST         50
# 8 LST         66
# 9 M           81
# 10 M           90



# ......................................................................... ----



# Wildfire ####

# Limit to fire sites only
m <- m[m$Disturbance=="Fire",]

# create factors
m$DevStage <- factor(m$DevStage, levels=c("Establishment",
                                          "CrownClosure","Thinning",
                                          "LateStageThinning","Mature"),
                     labels=c("E","CC","EST","LST","M"))

# Fire site characteristics
## Number of sites ----
m %>% group_by(DevStage) %>% dplyr::summarise(num_sites=n_distinct(SiteCode))
# DevStage num_sites
# <fct>        <int>
#   1 E                3
# 2 CC               3
# 3 EST              3
# 4 LST              2
# 5 M                3

## Stand age as of 2017 ----
m %>% group_by(DevStage) %>% dplyr::summarize(range=range(StandAge2017))
# DevStage range
# <fct>    <int>
#   1 E           12
# 2 E           12
# 3 CC          21 <-
# 4 CC          21
# 5 EST         37
# 6 EST         38
# 7 LST         56
# 8 LST         60
# 9 M           81
# 10 M           90

m %>% group_by(DevStage) %>% dplyr::summarize(min=min(StandAge2017), max=max(StandAge2017))

# DevStage   min   max
# <fct>    <int> <int>
#   1 E           12    12
# 2 CC          20    21
# 3 EST         37    38
# 4 LST         56    60
# 5 M           81    90






## Vegetation ----

library(dplyr)

veg <- read.csv("infiles/veg.csv", header=TRUE, stringsAsFactors = FALSE)

# fire disturbance only
veg.fire <- veg[veg$Disturbance=="Fire",]

# list any taxa with 80% + coverage
# establishment
veg.fire.e <- veg.fire[veg.fire$Devstage=="Establishment",]
veg.fire.e.sorted <- head(unique(veg.fire.e[order(-veg.fire.e$Cover_pct),c(1:3)]),15)

# crown closure
veg.fire.cc <- veg.fire[veg.fire$Devstage=="CrownClosure",]
veg.fire.cc.sorted <- head(unique(veg.fire.cc[order(-veg.fire.cc$Cover_pct),c(1:3)]),15)

# thinning
veg.fire.t <- veg.fire[veg.fire$Devstage=="Thinning",]
veg.fire.t.sorted <- head(unique(veg.fire.t[order(-veg.fire.t$Cover_pct),c(1:3)]),15)

# late self-thinning
veg.fire.lst <- veg.fire[veg.fire$Devstage=="LateStageThinning",]
veg.fire.lst.sorted <- head(unique(veg.fire.lst[order(-veg.fire.lst$Cover_pct),c(1:3)]),15)

# mature
veg.fire.m <- veg.fire[veg.fire$Devstage=="Mature",]
veg.fire.m.sorted <- head(unique(veg.fire.m[order(-veg.fire.m$Cover_pct),c(1:3)]),15)




## Site characteristics ----

meta <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# fire disturbance only
meta.fire <- meta[meta$Disturbance=="Fire",]

# create factor
meta.fire$DevStage <- factor(meta.fire$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                             labels=c("E","CC","EST","LST","M"))

### pH ----
#### Layers separate ----
meta.fire %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage Layer `mean(pH)` `sd(pH)`
# <fct>    <chr>        <dbl>      <dbl>
# 1 E        B             5.68     0.120 
# 2 E        M             5.73     0.303 
# 3 E        O             5.88     0.187 
# 4 CC       B             5.19     0.524 
# 5 CC       M             5.43     0.326 
# 6 CC       O             4.96     0.478 
# 7 EST      B             5        0.345 
# 8 EST      M             5.51     0.328 
# 9 EST      O             4.62     0.0866
# 10 LST      B             5.17     0.438 
# 11 LST      M             5.32     0.126 
# 12 LST      O             4.69     0.120 
# 13 M        B             5.12     0.234 
# 14 M        M             5.23     0.156 
# 15 M        O             4.55     0.195 

#### Layers pooled ----
meta.fire %>% group_by(DevStage) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage `mean(pH)` `sd(pH)`
# <fct>           <dbl>      <dbl>
# 1 E                5.76      0.226
# 2 CC               5.19      0.475
# 3 EST              5.04      0.460
# 4 LST              5.06      0.378
# 5 M                4.97      0.357
# 



### TOC ----
#### Layers separate ----
meta.fire %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TOC_gkg), sd(TOC_gkg))
# DevStage Layer `mean(TOC_gkg)` `sd(TOC_gkg)`
# <fct>    <chr>             <dbl>           <dbl>
# 1 E        B                432.             3.57 
# 2 E        M                 17.4            0.747
# 3 E        O                288.            46.9  
# 4 CC       B                500.            15.4  
# 5 CC       M                 12.3            2.73 
# 6 CC       O                485.            22.0  
# 7 EST      B                494.             9.32 
# 8 EST      M                  9.55           0.241
# 9 EST      O                486.             9.31 
# 10 LST      B                502.             6.30 
# 11 LST      M                 16.0            4.47 
# 12 LST      O                491.            30.1  
# 13 M        B                490.            14.4  
# 14 M        M                 10.8            3.29 
# 15 M        O                467.            15.5  


#### Layers pooled ----
meta.fire %>% group_by(DevStage) %>% dplyr::summarize(mean(TOC_gkg), sd(TOC_gkg))

# DevStage `mean(TOC_gkg)` `sd(TOC_gkg)`
# <fct>                <dbl>           <dbl>
# 1 E                     244.            180.
# 2 CC                    333.            231.
# 3 EST                   330.            231.
# 4 LST                   336.            234.
# 5 M                     323.            225.

### TN ----
#### Layers separate ----
meta.fire %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage Layer `mean(TN_gkg)` `sd(TN_gkg)`
# <fct>    <chr>            <dbl>          <dbl>
# 1 E        B               12.9           0.926 
# 2 E        M                0.887         0.0577
# 3 E        O               12.4           1.53  
# 4 CC       B               12.3           4.83  
# 5 CC       M                0.5           0.0300
# 6 CC       O               14.4           4.86  
# 7 EST      B                9.5           1.23  
# 8 EST      M                0.393         0.0589
# 9 EST      O               13.1           1.40  
# 10 LST      B               11.2           1.05  
# 11 LST      M                0.675         0.148 
# 12 LST      O               14.2           0.969 
# 13 M        B                9.90          0.965 
# 14 M        M                0.48          0.177 
# 15 M        O               13.1           1.01  

#### Layers pooled ----
meta.fire %>% group_by(DevStage) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage `mean(TN_gkg)` `sd(TN_gkg)`
# <fct>               <dbl>          <dbl>
# 1 E                    8.58           5.79
# 2 CC                   9.07           7.30
# 3 EST                  7.65           5.54
# 4 LST                  8.67           6.01
# 5 M                    7.81           5.50



# ......................................................................... ----



# Full tree clearcut harvest ----

m <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# Limit to harvest sites only 
m <- m[m$Disturbance=="Harvest",]

# create factors
m$DevStage <- factor(m$DevStage, levels=c("Establishment",
                                          "CrownClosure","Thinning",
                                          "LateStageThinning","Mature"),
                     labels=c("E","CC","EST","LST","M"))

## Number of sites ----
m %>% group_by(DevStage) %>% dplyr::summarise(num_sites=n_distinct(SiteCode))
# DevStage num_sites
# <fct>        <int>
# E                3
# 2 CC               3
# 3 EST              3
# 4 LST              2

## Stand age as of 2017 ----
m %>% group_by(DevStage) %>% dplyr::summarize(range=range(StandAge2017))
# DevStage range
# <fct>    <int>
#   1 E           12
# 2 E           15
# 3 CC          20<-
# 4 CC          25<-
# 5 EST         37<-
# 6 EST         41<-
# 7 LST         50<-
# 8 LST         66

m %>% group_by(DevStage) %>% dplyr::summarize(min=min(StandAge2017), max=max(StandAge2017))

# DevStage   min   max
# <fct>    <int> <int>
#   1 E           12    15
# 2 CC          20<-    25<-
# 3 EST         37<-    41<-
# 4 LST         50<-    66






## Vegetation ----

library(dplyr)

veg <- read.csv("infiles/veg.csv", header=TRUE, stringsAsFactors = FALSE)

# harvest disturbance only
veg.harvest <- veg[veg$Disturbance=="Harvest",]

# list any taxa with 80% + coverage
# establishment
veg.harvest.e <- veg.harvest[veg.harvest$Devstage=="Establishment",]
veg.harvest.e.sorted <- head(unique(veg.harvest.e[order(-veg.harvest.e$Cover_pct),c(1:3)]),15)

# crown closure
veg.harvest.cc <- veg.harvest[veg.harvest$Devstage=="CrownClosure",]
veg.harvest.cc.sorted <- head(unique(veg.harvest.cc[order(-veg.harvest.cc$Cover_pct),c(1:3)]),15)

# thinning
veg.harvest.t <- veg.harvest[veg.harvest$Devstage=="Thinning",]
veg.harvest.t.sorted <- head(unique(veg.harvest.t[order(-veg.harvest.t$Cover_pct),c(1:3)]),15)

# late self-thinning
veg.harvest.lst <- veg.harvest[veg.harvest$Devstage=="LateStageThinning",]
veg.harvest.lst.sorted <- head(unique(veg.harvest.lst[order(-veg.harvest.lst$Cover_pct),c(1:3)]),15)





## Site characteristics ----

meta <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# harvest disturbance only
meta.harvest <- meta[meta$Disturbance=="Harvest",]

# create factor
meta.harvest$DevStage <- factor(meta.harvest$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                             labels=c("E","CC","EST","LST","M"))


### pH ----

#### Layers separate ----
meta.harvest %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage Layer `mean(pH)` `sd(pH)`
# <fct>    <chr>      <dbl>    <dbl>
# 1 E        B           5.34   0.390 
# 2 E        M           5.31   0.212 
# 3 E        O           5.09   0.547 
# 4 CC       B           5.42   0.115 
# 5 CC       M           5.46   0.277 
# 6 CC       O           5.41   0.356 
# 7 EST      B           4.78   0.826 
# 8 EST      M           5.22   0.137 
# 9 EST      O           4.47   0.687 
# 10 LST      B           5.18   0.142 
# 11 LST      M           5.52   0.186 
# 12 LST      O           5.06   0.0493

#### Layers pooled ----
meta.harvest %>% group_by(DevStage) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage `mean(pH)` `sd(pH)`
# <fct>           <dbl>      <dbl>
# 1 E                5.24      0.407
# 2 CC               5.43      0.259
# 3 EST              4.82      0.678
# 4 LST              5.26      0.237




### TOC ----

#### Layers separate ----
meta.harvest %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TOC_gkg), sd(TOC_gkg))

# DevStage Layer `mean(TOC_gkg)` `sd(TOC_gkg)`
# <fct>    <chr>             <dbl>           <dbl>
# 1 E        B                 474.           10.9  
# 2 E        M                  14.1           0.543
# 3 E        O                 428.           53.9  
# 4 CC       B                 499.           14.2  
# 5 CC       M                  13.3           1.58 
# 6 CC       O                 477.           14.3  
# 7 EST      B                 504.           10.8  
# 8 EST      M                  13.8           4.62 
# 9 EST      O                 485.           13.1  
# 10 LST      B                 490.            3.36 
# 11 LST      M                  11.4           5.15 
# 12 LST      O                 478.           27.0  

#### Layers pooled ----
meta.harvest %>% group_by(DevStage) %>% dplyr::summarize(mean(TOC_gkg), sd(TOC_gkg))

# DevStage `mean(TOC_gkg)` `sd(TOC_gkg)`
# <fct>                <dbl>           <dbl>
# 1 E                     306.            213.
# 2 CC                    330.            229.
# 3 EST                   334.            231.
# 4 LST                   326.            230.

### TN ----

#### Layers separate ----
meta.harvest %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage Layer `mean(TN_gkg)` `sd(TN_gkg)`
# <fct>    <chr>          <dbl>        <dbl>
# 1 E        B             10.8         1.46  
# 2 E        M              0.63        0.0606
# 3 E        O             13.6         1.84  
# 4 CC       B             12.1         0.458 
# 5 CC       M              0.613       0.148 
# 6 CC       O             16.3         0.712 
# 7 EST      B             12.6         0.447 
# 8 EST      M              0.553       0.153 
# 9 EST      O             16.5         0.535 
# 10 LST      B             10.9         0.142 
# 11 LST      M              0.415       0.159 
# 12 LST      O             14.2         1.82  

#### Layers pooled ----
meta.harvest %>% group_by(DevStage) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage `mean(TN_gkg)` `sd(TN_gkg)`
# <fct>               <dbl>          <dbl>
# 1 E                    8.35           5.83
# 2 CC                   9.67           6.78
# 3 EST                  9.89           6.94
# 4 LST                  8.52           6.14



# ......................................................................... ----



# Salvage harvest ####

m <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# Limit to salvage sites only 
m <- m[m$Disturbance=="Salvage",]

# create factors
m$DevStage <- factor(m$DevStage, levels=c("Establishment",
                                          "CrownClosure","Thinning",
                                          "LateStageThinning","Mature"),
                     labels=c("E","CC","EST","LST","M"))

## Number of sites ----
m %>% group_by(DevStage) %>% dplyr::summarise(num_sites=n_distinct(SiteCode))
# DevStage num_sites
# <fct>        <int>
#   1 E                3
# 2 CC               3

## Stand age as of 2017 ----
m %>% group_by(DevStage) %>% dplyr::summarize(range=range(StandAge2017))
# DevStage range
# <fct>    <int>
#   1 E           12
# 2 E           13
# 3 CC          20
# 4 CC          20<-

m %>% group_by(DevStage) %>% dplyr::summarize(min=min(StandAge2017), max=max(StandAge2017))

# DevStage   min   max
# <fct>    <int> <int>
# 1 E           12    13
# 2 CC          20    20<-






## Vegetation ----

library(dplyr)

veg <- read.csv("infiles/veg.csv", header=TRUE, stringsAsFactors = FALSE)

# Salvage disturbance only
veg.salvage <- veg[veg$Disturbance=="Salvage",]

# list any taxa with 80% + coverage
# establishment
veg.salvage.e <- veg.salvage[veg.salvage$Devstage=="Establishment",]
veg.salvage.e.sorted <- head(unique(veg.salvage.e[order(-veg.salvage.e$Cover_pct),c(1:3)]),15)

# crown closure
veg.salvage.cc <- veg.salvage[veg.salvage$Devstage=="CrownClosure",]
veg.salvage.cc.sorted <- head(unique(veg.salvage.cc[order(-veg.salvage.cc$Cover_pct),c(1:3)]),15)





## Site characteristics ----

meta <- read.csv("infiles/metadata_2024-11-10_DM.csv", header=TRUE, stringsAsFactors = FALSE)

# salvage disturbance only
meta.salvage <- meta[meta$Disturbance=="Salvage",]

# create factor
meta.salvage$DevStage <- factor(meta.salvage$DevStage, levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                                labels=c("E","CC","EST","LST","M"))

#### pH ----

##### Layers separate ----
meta.salvage %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage Layer `mean(pH)` `sd(pH)`
# <fct>    <chr>        <dbl>      <dbl>
# 1 E        B             5.39     0.135 
# 2 E        M             5.82     0.0854
# 3 E        O             5.50     0.133 
# 4 CC       B             4.89     0.173 
# 5 CC       M             5.46     0.0527
# 6 CC       O             4.68     0.198 

##### Layers pooled ----
meta.salvage %>% group_by(DevStage) %>% dplyr::summarize(mean(pH), sd(pH))

# DevStage `mean(pH_F)` `sd(pH_F)`
# <fct>           <dbl>      <dbl>
# 1 E                5.56      0.215
# 2 CC               5.01      0.366


### TOC ----
##### Layers separate ----
meta.salvage %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TOC_gkg, na.rm = TRUE), sd(TOC_gkg, na.rm = TRUE))

# DevStage Layer `mean(TOC_gkg, na.rm = TRUE)` `sd(TOC_gkg, na.rm = TRUE)`
# <fct>    <chr>                         <dbl>                       <dbl>
# 1 E        B                            470.                         22.6 
# 2 E        M                             15.4                         7.39
# 3 E        O                            412.                         25.5 
# 4 CC       B                            473.                          3.48
# 5 CC       M                              8.30                        1.68
# 6 CC       O                            423.                         43.9 

##### Layers pooled ----
meta.salvage %>% group_by(DevStage) %>% dplyr::summarize(mean(TOC_gkg, na.rm = TRUE), sd(TOC_gkg, na.rm = TRUE))

# DevStage `mean(TOC_gkg, na.rm = TRUE)` `sd(TOC_gkg, na.rm = TRUE)`
# <fct>                            <dbl>                       <dbl>
# 1 E                                 310.<-                        203.<-
# 2 CC                                301.                        214.
# 



### TN ----

##### Layers separate ----
meta.salvage %>% group_by(DevStage, Layer) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage Layer `mean(TN_gkg)` `sd(TN_gkg)`
# <fct>    <chr>            <dbl>          <dbl>
# 1 E        B                9.53          1.89  
# 2 E        M                0.631<-        0.241<- 
# 3 E        O               10.8           1.78  
# 4 CC       B                9.73          1.08  
# 5 CC       M                0.337         0.0557
# 6 CC       O               12.0           1.09  

##### Layers pooled ----
meta.salvage %>% group_by(DevStage) %>% dplyr::summarize(mean(TN_gkg), sd(TN_gkg))

# DevStage `mean(TN_gkg_F)` `sd(TN_gkg_F)`
# <fct>               <dbl>          <dbl>
# 1 E                    7.24<-           4.76<-
# 2 CC                   7.34           5.20







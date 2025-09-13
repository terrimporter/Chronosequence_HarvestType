# Teresita M. Porter, Nov. 3, 2024

library(ggplot2)

# broken due to migration of tiles to new location, install again the usual way, but not from github
library(ggmap)
packageVersion("ggmap")
# Property: ggmap @ https://client.stadiamaps.com/dashboard/#/property/20742/
register_stadiamaps("06ace956-5b20-4ab1-bedf-ee0b159bbd67", 
                    write = TRUE)

# devtools::install_github('oswaldosantos/ggsn')
library(ggsn) # scalebar
library(RColorBrewer)

library(cowplot) #get_legend, plot_grid
library(dplyr)
library(tidyr)

m <- read.csv("infiles/metadata_2024-11-10.csv", header=TRUE, stringsAsFactors = FALSE)
# keep SiteCode, Disturbance, DevStage, Latitude, Longitude
m <- unique(m[,c(3,9,10,12,14,15)])

# Table 1
# create factors
m$DevStage <- factor(m$DevStage, 
                     levels=c("Establishment", "CrownClosure", "Thinning", "LateStageThinning", "Mature"),
                     labels=c("E", "CC", "EST", "LST", "M"))

m$Disturbance <- factor(m$Disturbance, 
                     levels=c("Fire", "Harvest", "Salvage"),
                     labels=c("Wildfire", "Clearcut Harvest", "Salvage Harvest"))

m.table_years_since_disturbance <- m %>%
  dplyr::group_by(DevStage) %>%
  dplyr::summarize(min = min(StandAge2017),
                   max = max(StandAge2017))
m.table_years_since_disturbance 

m.table <- reshape2::dcast(m, DevStage ~ Disturbance)

m.table

bbox <- c(bottom = 47.8, top = max(m$Latitude)+1.2 , 
          right = max(m$Longitude)+1, left = -92.6)

# name of maptype updated
base <- get_stadiamap(bbox = bbox, zoom = 7, maptype = 'stamen_terrain_background') 


# cols for disturbance type
# Star Trek color palette lol (red, blue, amber)
cols <- c("#df0000", "#0099f6", "#f2c300")

map.tmp <- ggmap(base) +
  geom_point(data=m, aes(x=Longitude, y=Latitude, colour = Disturbance, fill = Disturbance, shape=DevStage), cex=2.5) + # plot the points
  # Thunder Bay
  geom_point(aes(x=-89.23477, y=48.3809), colour="black") +
  geom_text(data = NULL, x = -89.23477, y = 48.3809, size = 3, label = "Thunder Bay", hjust=1.2, vjust =1.2 ) +
  # Ignace
  geom_point(aes(x=-91.6587, y=49.4159), fill="black") +
  geom_text(data = NULL, x=-91.6587, y=49.4159, size = 3, label = "Ignace", hjust=1.2, vjust =1.2 ) +
  labs(x="", y="", title="") + # label the axes
  theme_bw() + 
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = c(21:25)) +
  theme(
    text = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key = element_rect(colour = "white"),
    legend.title = element_blank(),
    #      legend.spacing.y = unit(0.01, 'cm'),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-9,-9,-9,-9)
  ) 

map.tmp

l <- get_legend(map.tmp)

map <- ggmap(base) +
  geom_point(data=m, aes(x=Longitude, y=Latitude, fill = Disturbance, shape=DevStage), alpha = 0.5, cex=2.5) + # plot the points
  # Thunder Bay
  geom_point(aes(x=-89.23477, y=48.3809), colour="black") +
  geom_text(data = NULL, x = -89.23477, y = 48.3809, size = 3, label = "Thunder Bay", hjust=1.2, vjust =1.2 ) +
  # Ignace
  geom_point(aes(x=-91.6587, y=49.4159), fill="black") +
  geom_text(data = NULL, x=-91.6587, y=49.4159, size = 3, label = "Ignace", hjust=1.2, vjust =1.2 ) +
  # Lake Nipigon
  # geom_point(aes(x=-91.6587, y=49.4159), fill="black") +
  geom_text(data = NULL, x=-88.484055, y=49.789426, size = 3, label = "Lake\nNipigon", hjust=0.5, vjust =0.5 ) +
  labs(x="", y="", title="") + # label the axes
  theme_bw() + 
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = c(21:25)) +
  theme(
        text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.key = element_rect(colour = "white"),
        legend.title = element_blank(),
  #      legend.spacing.y = unit(0.01, 'cm'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-9,-9,-9,-9)
  ) 

map

map <- map + 
  scalebar(x.min = -89.15, x.max = -87.5,
              y.min = 47.9,  y.max = 51.54807, 
              dist = 50, transform = TRUE, model = 'WGS84', dist_unit = "km", 
              height = 0.03, st.bottom = FALSE, st.dist = 0.04, st.size = 3,
               border.size = 0.7, 
           st.color = "black")
map

############################
# inset map

# bbox <- c(bottom = 48, top = max(m$Latitude)+0.1 , right = -87.76792, left = -92)
bbox2 <- c(bottom = 47.8-7.3, top = max(m$Latitude)+1.2 , right = -87.76792+12.7, left = -92-2.5)

base2 <- get_stadiamap(bbox = bbox2, zoom = 7, maptype = 'stamen_terrain_background') 

map2 <- 
  ggmap(base2, extent = "device") +
  # Ignace
  geom_point(x=-89.4, y=49.36, 
             color="black", fill="transparent", shape=22, cex=10) + # plot the city
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.title= element_blank(),
    legend.title = element_blank(),
  )
map2



# put inset map in the right place
# bbox <- c(bottom = 48, top = max(m$Latitude)+0.1 , right = -87.76792, left = -92)
map3 <- map +
  inset(ggplotGrob(map2), xmin = -92.5, xmax = -90.5, ymin = 48.7, ymax = max(m$Latitude)+2.5)
map3

g <- plot_grid(map3, l, rel_widths = c(1, 0.3))
g

# check if outfiles dir exists, create if needed
ifelse(!dir.exists(file.path("outfiles")), dir.create(file.path("outfiles")), FALSE)
# change to jpg so that gov windows users can see points OMG
ggsave("outfiles/Fig1_map.jpg", g, height = 4, width = 6)
#   Add to legend: © Stadia Maps © OpenMapTiles © OpenStreetMap © Stamen Design

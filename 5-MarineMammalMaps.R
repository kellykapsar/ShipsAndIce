
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(yarrr)
library(ggsn)
library(colorspace)
library(RColorBrewer)
library(scales)
library(viridis)


saveloc <- "../Figures/"

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# all sf == spatial data frame (each row is a unique cell)
allsf <- st_read("../Data_Processed/IceTrafficDataFrame.shp")

# Basemap
# aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% 
  st_transform(AA)
basemap.crop <- st_crop(basemap, st_buffer(allsf, 10000)) %>% 
  st_simplify(dTolerance=1000, preserveTopology = T)

########################### BOWHEAD WHALE ########################### 
# Read in bowhead concentration areas
whale<- st_read("../Data_Raw/Bowhead_RelAbund_Winter/Bowhead_RelAbund_Winter.shp") %>% st_transform(AA)
hicon <- whale %>% filter(area_desc == "High_Conc - Winter")
con <- whale %>% filter(area_desc == "Concentration - Winter")
rang <- whale %>% filter(area_desc == "Range - Winter")

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=rang, fill="#d5ebf1", lwd=0) +
  geom_sf(data=con, fill="#96acf2", lwd=0) +
  geom_sf(data=hicon, fill="#768ce4",lwd=0) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/Bowhead_Winter.png",
       width=10, height=8, units="in", dpi=300)

########################### WALRUS ########################### 
# Read in walrus concentration areas
walrus<- st_read("../Data_Raw/Walrus_WinterSpringRelAbund/OdobenusRosmarus_WinterSpringRelAbund.shp") %>% st_transform(AA)
walcon <- walrus %>% filter(area_desc == "Concentration - Winter & Spring")
walreg <- walrus %>% filter(area_desc == "Regular Use - Winter & Spring")
walrang <- walrus %>% filter(area_desc == "Winter and Spring Range")

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=walrang, fill="#d5ebf1", lwd=0) +
  geom_sf(data=walreg, fill="#a6cee3", lwd=0) +
  geom_sf(data=walcon, fill="#96acf2", lwd=0) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/Walrus_Winter.png",
       width=10, height=8, units="in", dpi=300)

########################### VESSEL TRAFFIC IN BOWHEAD CONCENTRATION AREAS ########################### 

# all sf == spatial data frame (each row is a unique cell)
allsf <- st_read("../Data_Processed/IceTrafficDataFrame.shp")
# studyoutline <- allsf %>% select(id) %>% st_union() %>% st_write("../Data_Processed/StudyOutline.shp", overwrite=T)

alldf <- read.csv("../Data_Processed/IceTrafficDataFrame.csv") %>% dplyr::select(-X)

# ID cells in areas of high concentration of bowheads
hiconcells <- allsf$id[st_intersects(allsf, hicon, sparse=F)]

hicontraff <- alldf[alldf$id %in% hiconcells,]
hicontotals <- hicontraff %>%  filter(month %in% c(11, 12, 1, 2, 3)) %>% group_by(year) %>% summarize(traff=sum(traffic_km))


ggplot(hicontotals, aes(x=year, y=traff, fill=as.factor(year))) +
  geom_bar(stat="identity") +
  scale_fill_manual(labels=2015:2020, values=rev(viridis::viridis(6)), name="Year") +
  scale_x_continuous(breaks= 2015:2020) +
  # scale_y_continuous(breaks=seq(0,300000, by=50000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffi in Bowhead \nWinter High Concentration Areas") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggsave("../Figures/AISInBowheadHiCon_AnnualBar.png", width=9, height=6, units="in")


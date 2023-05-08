########################################################################
# TITLE: Vessel Traffic in Marine Mammal Concentration Areas
#
# DESCRIPTION: This script takes the combined AIS and ice products created
# in script 3 and spatially intersects it with areas of bowhead whale winter
# concentration from the Ecological Atlas of the Bering, Chukchi, and Beaufort 
# Seas produced by Auduobon in 2016. 
#
# CREATED BY: Kelly Kapsar (kelly.kapsar@gmail.com)
# DATE CREATED: 2022-07
# DATE LAST MODIFIED: 2022-08-23
########################################################################

# Load libraries 
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

# Specify projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

###############################################################
# Load data 

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
  # geom_sf(data=hiconnew[4,], fill="red",lwd=0) +
  # geom_sf(data=hiconnew, aes(fill=traff),lwd=0) +
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

# ggsave(filename= "../Figures/Bowhead_Winter.png",
#        width=10, height=8, units="in", dpi=300)

# hicontraffsf <- allsf[st_intersects(allsf, hicon, sparse=F),]
# 
# hicontraffsf$traff.quantile = cut(hicontraffsf$traff, breaks=quantile.interval, include.lowest = TRUE)
# 
# 
# ggplot() +
#   geom_sf(data=st_crop(basemap.crop, st_buffer(hicontraffsf, 100000)), fill="white", color="black", lwd=0.5, alpha = 0.9) +
#   geom_sf(data=hicontraffsf, lwd=0.5, fill=NA) +
#   geom_sf(data=hicon, lwd=0.5) +
#   xlab("") +
#   ylab("") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_bw() +
#   blank()+
#   theme(legend.title = element_text(size = 20),
#         legend.text = element_text(size = 20),
#         legend.position = "left",
#         legend.background = element_rect(fill = "white", color = "black"),
#         axis.ticks = element_blank(),
#         axis.text=element_blank(),
#         # panel.background = element_rect(fill = "lightblue"),
#         panel.border =  element_rect(colour = "black"),
#         panel.grid.major = element_line(colour = "transparent"))

# ggsave(filename= "../Figures/Bowhead_Winter.png",
#        width=10, height=8, units="in", dpi=300)

########################### VESSEL TRAFFIC IN BOWHEAD CONCENTRATION AREAS ########################### 

# all sf == spatial data frame (each row is a unique cell)
allsf <- st_read("../Data_Processed/IceTrafficDataFrame.shp")
# studyoutline <- allsf %>% select(id) %>% st_union() %>% st_write("../Data_Processed/StudyOutline.shp", overwrite=T)

alldf <- read.csv("../Data_Processed/IceTrafficDataFrame.csv") %>% dplyr::select(-X)



alldf$winterid <- ifelse(alldf$year == 2015 & alldf$month < 5, 1, 
                  ifelse(alldf$year == 2015 & alldf$month > 5, 2, 
                  ifelse(alldf$year == 2016 & alldf$month < 5, 2, 
                  ifelse(alldf$year == 2016 & alldf$month > 5, 3, 
                  ifelse(alldf$year == 2017 & alldf$month < 5, 3, 
                  ifelse(alldf$year == 2017 & alldf$month > 5, 4, 
                  ifelse(alldf$year == 2018 & alldf$month < 5, 4, 
                  ifelse(alldf$year == 2018 & alldf$month > 5, 5, 
                  ifelse(alldf$year == 2019 & alldf$month < 5, 5, 
                  ifelse(alldf$year == 2019 & alldf$month > 5, 6, 
                  ifelse(alldf$year == 2020 & alldf$month < 5, 6, 
                  ifelse(alldf$year == 2020 & alldf$month > 5, 7, NA))))))))))))

# ID cells in areas of high concentration of bowheads
hiconnew <- st_cast(hicon, "POLYGON")

st_area(hiconnew[5,])/1000000

hiconnew$traff <- NA
hiconpolydat <- data.frame()

for(i in 1:length(hiconnew$orig_id)){
hiconcells <- allsf$id[st_intersects(allsf, hiconnew[i,], sparse=FALSE)]
hicontraff <- alldf[alldf$id %in% hiconcells,]

hicontotals <- hicontraff %>%  
  filter(month %in% c(11, 12, 1, 2, 3)) %>% 
  group_by(winterid) %>% 
  summarize(traff=sum(traffic_km))

hicontotals$polyid <- i
hiconpolydat <- rbind(hiconpolydat, hicontotals)

hiconnew$traff[i] <- sum(hicontotals$traff)

}


hiconcells <- allsf$id[st_intersects(allsf, hicon, sparse=FALSE)]
hicontraff <- alldf[alldf$id %in% hiconcells,]

hicontotals <- hicontraff %>%  
  filter(month %in% c(11, 12, 1, 2, 3)) %>% 
  group_by(winterid) %>% 
  summarize(traff=sum(traffic_km))

hiconmonth <- hicontraff %>%  
  filter(month %in% c(11, 12, 1, 2, 3)) %>%  
  group_by(year, month, winterid) %>% 
  summarize(traff=sum(traffic_km), iceext = sum(ifelse(icecon >= 15, 1, 0))*622.109950)

hiconmonth$monthorder <- ifelse(hiconmonth$month == 11, 1,
                                ifelse(hiconmonth$month == 12, 2, 
                                       ifelse(hiconmonth$month == 1, 3, 
                                              ifelse (hiconmonth$month == 2, 4, 
                                                      ifelse(hiconmonth$month == 3, 5, NA 
                                )))))

hicon_traffchange <- ((hicontotals$traff[hicontotals$winterid == 6]-hicontotals$traff[hicontotals$winterid == 2])/hicontotals$traff[hicontotals$winterid == 2])*100
hiconarea_km2 <- st_area(hicon)/1e6

# ID cells in areas of  concentration of bowheads
concells <- allsf$id[st_intersects(allsf, con, sparse=F)]
contraff <- alldf[alldf$id %in% concells,]
contotals <- contraff %>%  filter(month %in% c(11, 12, 1, 2, 3)) %>% group_by(winterid) %>% summarize(traff=sum(traffic_km))
con_traffchange <- ((contotals$traff[contotals$winterid == 6]-contotals$traff[contotals$winterid == 2])/contotals$traff[contotals$winterid == 2])*100
conarea_km2 <- st_area(con)/1e6

# ID cells in winter range of bowheads
rangcells <- allsf$id[st_intersects(allsf, rang, sparse=F)]
rangtraff <- alldf[alldf$id %in% rangcells,]
rangtotals <- rangtraff %>%  filter(month %in% c(11, 12, 1, 2, 3)) %>% group_by(winterid) %>% summarize(traff=sum(traffic_km))
rang_traffchange <- ((rangtotals$traff[rangtotals$winterid == 6]-rangtotals$traff[rangtotals$winterid == 2])/rangtotals$traff[rangtotals$winterid == 2])*100
rangarea_km2 <- st_area(rang)/1e6

# Total change in traffic
alltotals <- alldf %>%  filter(month %in% c(11, 12, 1, 2, 3)) %>% group_by(winterid) %>% summarize(traff=sum(traffic_km))

all_traffchange <- ((alltotals$traff[alltotals$winterid == 6]-alltotals$traff[alltotals$winterid == 2])/alltotals$traff[alltotals$winterid == 2])*100

allarea_km2 <- st_area(allsf)/1e6



labs <- c("2014-2015*", "2015-2016", "2016-2017", "2017-2018", "2018-2019", "2019-2020", "2020-2021*")

ggplot(hicontotals, aes(x=winterid, y=traff)) +
  geom_bar(stat="identity") +
  scale_fill_manual(labels=labs,values=rev(viridis::viridis(7)), name="Winter") +
  scale_x_continuous(breaks=1:7, labels=labs) +
  scale_y_continuous(breaks=seq(0,300000, by=50000)) +
  theme_bw(base_size = 20) +
  geom_text(aes(label=c("*", "", "", "","","","*")), position=position_stack(vjust=1), size=5) +
  xlab("") +
  ylab("Total Vessel Traffic in Bowhead \nWinter High Concentration Areas") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggsave("../Figures/AISInBowheadHiCon_km_AnnualBar.png", width=9, height=6, units="in")


ggplot() +
  geom_line(data=hiconmonth, aes(x=monthorder, y =traff, group=winterid, col=as.factor(winterid)),lwd=1.5) +
  scale_color_manual(labels=labs, values=rev(viridis::viridis(7)), name="Winter") + 
  xlab("") +
  scale_x_continuous(breaks=1:5, labels = c("Nov", "Dec", "Jan", "Feb", "Mar") ) +
  scale_y_continuous(name = "Total Vessel Traffic in Bowhead \nWinter High Concentration Areas") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("../Figures/AISinBowhead_MonthlyLines_km.png", width=9, height=6, units="in")

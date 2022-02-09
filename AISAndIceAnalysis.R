
library(raster)
library(stars)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(yarrr)
library(ggsn)


# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Basemap
aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% 
              st_transform(AA)
basemap.crop <- st_crop(basemap, st_buffer(aisbounds, 100000)) %>% 
                  st_simplify(dTolerance=1000, preserveTopology = T)

alldf <- read.csv("../Data_Processed/IceTrafficDataFrame.csv")
# all sf == spatial data frame (each row is a unique cell)
test <- st_read("../Data_Processed/IceTrafficDataFrame.shp") %>% st_intersection(aisbounds)

inice <- alldf[which(alldf$icecon > 0 & alldf$traffic_km > 0),]

hightraffids <- unique(inice$id[which(inice$traffic_km > 10000)])

############################################################################
ggplot(alldf, aes(x=traffic_km)) +
  geom_histogram(bins=50)

ggplot(alldf, aes(x=icecon)) +
  geom_histogram(bins=50)

############################################################################
# Map of study area 
############################################################################ 

p3 <- ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=aisbounds, fill=NA, color="black", lwd=1)+
  geom_sf(data=test, aes(fill=c_2020_01))+
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.05,0.125),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

p3
ggsave(plot=p3, filename= "TestMap.png", 
       width=8, height=8, units="in")

############################################################################
hist(log(alldf$traffic_km))
alldf$log_traffic_km <- log(alldf$traffic_km)

aisinice <- alldf[which(alldf$icecon > 1),]

fit1 <- lm(log_traffic_km ~ icecon, data = aisinice)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(fit1)

ggplot(data=aisinice, aes(x=icecon, y=log(traffic_km))) +
  geom_point() +
  geom_smooth()

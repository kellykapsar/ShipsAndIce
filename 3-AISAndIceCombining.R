########################################################################
# TITLE: Merged AIS and Ice Product
#
# DESCRIPTION: This script takes the sea ice and vessel traffic data sets
# made in scripts 1 & 2 and combines them into a merged product containing
# sea ice concentration, sea ice thickness, vessel traffic (km), and number
# of unique ships for the months of October through April of 2015-2020. Data
# is stored by individual pixel for each pixel of the study area that has a sea 
# ice concentration value of >15% for at least one month during the study period. 
#
# CREATED BY: Kelly Kapsar (kelly.kapsar@gmail.com)
# DATE CREATED: 2022-02
# DATE LAST MODIFIED: 2022-08-23
########################################################################

# Load libraries 
library(raster)
library(stars)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load AIS data 
origmethod <- list.files("../Data_Processed/", pattern="LenAndNShips")
origmethod <- lapply(origmethod, function(x){read.csv(paste0("../Data_Processed/",x))})
origmethod <- do.call(rbind, origmethod)

# Convert distance to km and remove old distance calculation in m
origmethod$traffic_km <- origmethod$length/1000
origmethod$carg_km <- origmethod$CargoDist/1000
origmethod$tank_km <- origmethod$TankDist/1000
origmethod$other_km <- origmethod$OtherDist/1000
origmethod$fish_km <- origmethod$FishDist/1000
origmethod <- origmethod %>% dplyr::select(-length, -CargoDist, -FishDist, -TankDist, -OtherDist) 

allpixels <- origmethod %>% dplyr::select(-X)
allpixels <- allpixels[-which(allpixels$year == -201),]

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Basemap
aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)

# Median sea ice extent 1981-2010
ext <- st_read("../Data_Raw/median_extent_N_03_1981-2010_polyline_v3.0/median_extent_N_03_1981-2010_polyline_v3.0.shp") %>% st_transform(AA)

# Check to make sure no duplicate cell values for each year month combo
qatest <- allpixels %>% group_by(year, month) %>% summarize(ncells=length(id), nuniquecells=length(unique(id)))
which(qatest$ncells != qatest$nuniquecells)


# Create data frames for thickness and concentration 
icesf <- st_read("../Data_Processed/Ice_SMOS.shp")
icesf <- icesf %>% mutate_all(~replace_na(., 0))

condf <-  icesf %>% 
  dplyr::select(c(1:42), id) %>% 
  st_drop_geometry() %>%
  gather(key=yearmon, value=icecon, -id) %>%
  mutate(year = as.numeric(substr(yearmon, start=3, stop=6)),
         month = as.numeric(substr(yearmon, start=8, stop=9)),
         icecon = round(icecon, 2)) %>%
  dplyr::select(-yearmon)

thickdf <-  icesf %>% 
  dplyr::select(c(43:85), id) %>% 
  st_drop_geometry() %>%
  gather(key=yearmon, value=icethick, -id) %>%
  mutate(year = as.numeric(substr(yearmon, start=3, stop=6)),
         month = as.numeric(substr(yearmon, start=8, stop=9)),
         icethick = round(icethick, 2)) %>%
  dplyr::select(-yearmon)

# Join thickness and concentration data frames 
allice <- left_join(condf, thickdf, by=c("year", "month", "id"))

# Assume any ice concentration < 15% is no ice 
allice$icecon <- ifelse(allice$icecon < 15, 0, allice$icecon)
allice$icethick <- ifelse(allice$icecon < 15, 0, allice$icethick)

# Revise sf object to adhere to 15% threshold
iceconsf <- allice %>% 
  mutate(date=as.Date(paste0(year,"-",month,"-1"), format="%Y-%m-%d")) %>% 
  dplyr::select(id, date, icecon) %>%  
  spread(date, icecon)

colnames(iceconsf) <- c(colnames(iceconsf[1]), 
                        paste0("c_", substr(colnames(iceconsf)[2:43], 1,4),"_",substr(colnames(iceconsf[2:43]),6,7)))
allsf <- icesf %>% select(id) %>% left_join(., iceconsf, by=c("id"))

icethicksf <- allice %>% 
  mutate(date=as.Date(paste0(year,"-",month,"-1"), format="%Y-%m-%d")) %>% 
  dplyr::select(id, date, icethick) %>%  
  spread(date, icethick)

colnames(icethicksf) <- c(colnames(icethicksf[1]), 
                          paste0("t_", substr(colnames(icethicksf)[2:43], 1,4),"_",substr(colnames(icethicksf[2:43]),6,7)))
allsf <- allsf %>% left_join(., icethicksf, by=c("id"))



# Remove pixels with no ice for entire study period from data sets
noice <- allice %>% group_by(id) %>% summarize(anyice = ifelse(sum(icecon > 0), TRUE, FALSE)) %>% filter(anyice == FALSE)

allice <- allice[which(!(allice$id %in% noice$id)),]

allsf <- allsf[which(!(allsf$id %in% noice$id)),]

allpixels <- allpixels[which(!(allpixels$id %in% noice$id)),]

# cor.test(allice$icethick, allice$icecon)

# Join ice data with vessel traffic data 
alldf <-left_join(allice, allpixels, by=c("year", "month", "id"))
alldf <- alldf[,c("year", "month", "id", "icecon", "icethick",
                  "traffic_km", "carg_km", "fish_km", "tank_km", "other_km",
                  "nShips", "nCargo", "nFish", "nTank", "nOther")]

# alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$traffic_km <- round(alldf$traffic_km, 2)
alldf$carg_km <- round(alldf$carg_km, 2)
alldf$fish_km <- round(alldf$fish_km, 2)
alldf$tank_km <- round(alldf$tank_km, 2)
alldf$other_km <- round(alldf$other_km, 2)

# write.csv(alldf, "../Data_Processed/TrafficAndIcePixelValues.csv")

# Add traffic onto spatial object for icesf
nShipssf <- allpixels %>% 
  mutate(date=as.Date(paste0(year,"-",month,"-1"), format="%Y-%m-%d")) %>% 
  dplyr::select(id, date, nShips) %>%  
  spread(date, nShips)

colnames(nShipssf) <- c(colnames(nShipssf[1]), 
                        paste0("n_", substr(colnames(nShipssf)[2:73], 1,4),".",substr(colnames(nShipssf[2:73]),6,7)))

allsf <- left_join(allsf, nShipssf, by=c("id"))

# Add nShips onto spatial object 
trafficsf <- allpixels %>% 
  mutate(date=as.Date(paste0(year,"-",month,"-1"), format="%Y-%m-%d")) %>% 
  dplyr::select(id, date, traffic_km) %>%  
  spread(date, traffic_km)

colnames(trafficsf) <- c(colnames(trafficsf[1]), 
                         paste0("v_", substr(colnames(trafficsf)[2:73], 1,4),".",substr(colnames(trafficsf[2:73]),6,7)))

allsf <- left_join(allsf, trafficsf, by=c("id"))

# Replace NAs with zeros 
alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$carg_km[which(is.na(alldf$carg_km))] <- 0
alldf$tank_km[which(is.na(alldf$tank_km))] <- 0
alldf$fish_km[which(is.na(alldf$fish_km))] <- 0
alldf$other_km[which(is.na(alldf$other_km))] <- 0
alldf$nShips[which(is.na(alldf$nShips))] <- 0
alldf$nCargo[which(is.na(alldf$nCargo))] <- 0
alldf$nTank[which(is.na(alldf$nTank))] <- 0
alldf$nFish[which(is.na(alldf$nFish))] <- 0
alldf$nOther[which(is.na(alldf$nOther))] <- 0

# Remove pixels outside the study area
incells <- st_contains(aisbounds, allsf, sparse=FALSE)
allsfin <- allsf[incells,]

# Remove cells outside AIS boundaries from data frame 
alldf <- alldf[which(alldf$id %in% unique(allsfin$id)),]

# alldf == non-spatial data frame (each row is a unique cell, month, year combo)
write.csv(alldf, "../Data_Processed/IceTrafficDataFrame.csv")
# all sf == spatial data frame (each row is a unique cell)
st_write(allsfin, "../Data_Processed/IceTrafficDataFrame.shp")

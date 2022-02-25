
library(raster)
library(stars)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

origmethod <- list.files("../Data_Processed/", pattern="LenAndNShips")
origmethod <- lapply(origmethod, function(x){read.csv(paste0("../Data_Processed/",x))})
origmethod <- do.call(rbind, origmethod)
origmethod$traffic_km <- origmethod$length/1000
origmethod <- origmethod %>% dplyr::select(-length) 

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

# Assume any ice concentration < 10% is no ice 
allice$icecon <- ifelse(allice$icecon < 10, 0, allice$icecon)
allice$icethick <- ifelse(allice$icecon < 10, 0, allice$icethick)

# Revise sf object to adhere to 10% threshold
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

# Joine ice data with vessel traffic data 
alldf <-left_join(allice, allpixels, by=c("year", "month", "id"))
alldf <- alldf[,c("year", "month", "id", "icecon", "icethick","traffic_km", "nShips", "nCargo", "nFish", "nTank", "nOther")]

# alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$traffic_km <- round(alldf$traffic_km, 2)

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
  dplyr::select(-year, -month, -nShips, -nCargo, -nFish, -nTank, -nOther) %>%  
  spread(date, traffic_km)

colnames(trafficsf) <- c(colnames(trafficsf[1]), 
                         paste0("v_", substr(colnames(trafficsf)[2:73], 1,4),".",substr(colnames(trafficsf[2:73]),6,7)))

allsf <- left_join(allsf, trafficsf, by=c("id"))
alldf <- alldf[which(alldf$id %in% unique(allsf$id)),]

alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$nShips[which(is.na(alldf$nShips))] <- 0
alldf$nCargo[which(is.na(alldf$nCargo))] <- 0
alldf$nTank[which(is.na(alldf$nTank))] <- 0
alldf$nFish[which(is.na(alldf$nFish))] <- 0
alldf$nOther[which(is.na(alldf$nOther))] <- 0

# Remove pixels outside the study area
incells <- st_contains(aisbounds, allsf, sparse=FALSE)
allsfin <- allsf[incells,]

# alldf == non-spatial data frame (each row is a unique cell, month, year combo)
write.csv(alldf, "../Data_Processed/IceTrafficDataFrame.csv")
# all sf == spatial data frame (each row is a unique cell)
st_write(allsfin, "../Data_Processed/IceTrafficDataFrame.shp")

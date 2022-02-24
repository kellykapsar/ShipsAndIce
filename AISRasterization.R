
library(raster)
library(sf)
library(dplyr)
library(tidyr)

start <- proc.time()

# Specify year 
year <- 2015

# Load in all vessel shp files 
filelist <- list.files(paste0("D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/", year, "/Vector"), pattern='.shp')

files <- lapply(filelist, 
            function(x){st_read(paste0("D:/AlaskaConservation_AIS_20210225/Data_Processed_HPCC_FINAL/", year,"/Vector/", x))})
# Group all vessel types together for each month 
idx <- paste0(substr(filelist, start=8, stop=11), substr(filelist, start=13, stop=14))
traffic <- lapply(unique(idx), function(x){do.call(rbind, files[which(idx == x)])})

# Load in sea ice data & isolate just cell ids 
icecon <- st_read("../Data_Processed/Ice_SMOS.shp") %>% select(id)


#######################################################################################
# Method 1: For loop
# Takes forever, but uses less memory... 

start <- proc.time()

df <- data.frame(years = c(), months =c(), cells=c(), traffic=c())
for(i in 1:12){
  print(paste0("Processing month: ", i))
  t <- traffic[[i]]
  t$year <- substr(t$AIS_ID, start=11, stop =14)
  t$month <- substr(t$AIS_ID, start=15, stop =16)
  intersections_mat <- sf::st_intersects(t, iceconsf, sparse=F)
  tin <- t[which(rowSums(intersections_mat) > 0),]
  for(j in 1:length(tin$AIS_ID)){
    if(j %% 1000 == 0){
      print(paste0(round(j/length(tin$AIS_ID), 2)*100, "% finished."))
    }
    lenin <- st_intersection(tin[j,], iceconsf) %>% mutate(newlen=st_length(.)) %>% st_drop_geometry() %>% select(id, newlen)
    dftemp <- data.frame(years=rep(t$year[1], length(lenin$id)), 
                         months=rep(t$month[1], length(lenin$id)), 
                         cells= lenin$id,
                         traffic=lenin$newlen)
    df <- rbind(df, dftemp)
  }
}

write.csv(df, paste0("../Data_Processed/TrafficInPixels_",year,"_LoopMethod.csv"))
df2 <- df %>% group_by(years, months, cells) %>% summarize(traffic_km = as.numeric(sum(traffic)/1000))
write.csv(df2, paste0("../Data_Processed/TrafficInPixels_",year,"_LoopMethod_2.csv"))

(proc.time()-start)/60
# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

#################################################################################
# Method 2: Apply 
# Faster, but uses a lot of memory so it doesn't work on 2018-2020

# Crop vessel tracks at pixel boundaries and count total length 
# per pixel in each month of study period 
icepixellength <- function(trafficsf, icesf){
  int <- st_intersection(trafficsf, icesf)
  int$newlen <- st_length(int)
  int$year <- as.numeric(substr(int$AIS_ID, start=11, stop =14))
  int$month <- as.numeric(substr(int$AIS_ID, start=15, stop =16))
  intdf <- int %>%
    st_drop_geometry() %>%
    select(id, newlen, year, month, MMSI_x, AIS_Typ) %>%
    group_by(id, year, month) %>%
    summarize(length=sum(newlen), 
              nShips=length(unique(MMSI_x)), 
              nCargo = length(unique(MMSI_x[which(AIS_Typ == "Cargo")])),
              nFish = length(unique(MMSI_x[which(AIS_Typ == "Fishing")])),
              nTank = length(unique(MMSI_x[which(AIS_Typ == "Tanker")])),
              nOther = length(unique(MMSI_x[which(AIS_Typ == "Other")])))
  return(intdf)
}


# Takes ~15 minutes to run 
trafficcellsa <- lapply(traffic[1:6], function(x){icepixellength(x, iceconsf)})
# Bind together all months of data 
trafficcellsa <- do.call(rbind, trafficcellsa)
# write.csv(trafficcellsa, paste0("../Data_Processed/TrafficInIcePixels_", year, "a.csv"))
rm(trafficcella)

lapply(iceconsf)


traffic7 <- st_intersection(traffic[[7]], iceconsf)
trafficcellsb <- lapply(traffic[7:8], function(x){icepixellength(x, iceconsf)})
# Bind together all months of data 
trafficcellsb <- do.call(rbind, trafficcellsb)
write.csv(trafficcellsb, paste0("../Data_Processed/TrafficInIcePixels_", year, "b.csv"))
rm(trafficcellb)


trafficcellsc <- lapply(traffic[9:12], function(x){icepixellength(x, iceconsf)})
# Bind together all months of data 
trafficcellsc <- do.call(rbind, trafficcellsc)
write.csv(trafficcellsc, paste0("../Data_Processed/TrafficInIcePixels_", year, "c.csv"))


(proc.time()-start)/60





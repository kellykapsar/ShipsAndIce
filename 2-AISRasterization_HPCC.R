
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(stars)
library(foreach)
library(doParallel)

start <- proc.time()

# Specify year 
year <- 2015

# Load in all vessel shp files 
filelist <- list.files(paste0("/mnt/ufs18/home-109/kapsarke/Documents/Data_Processed/Vector/"), pattern='.shp')
filelistyears <- filelist[grep(year, filelist)]

files <- lapply(filelistyears, 
            function(x){st_read(paste0("/mnt/ufs18/home-109/kapsarke/Documents/Data_Processed/Vector/", x))})
# Group all vessel types together for each month 
idx <- paste0(substr(filelistyears, start=8, stop=11), substr(filelist, start=13, stop=14))
traffic <- lapply(unique(idx), function(x){do.call(rbind, files[which(idx == x)])})

# Load in sea ice data & isolate just cell ids 
iceconsf <- st_read("../Data_Raw/Ice_SMOS.shp") %>% select(id)


#######################################################################################
# Method 1: For loop
#######################################################################################
# Takes forever, but uses less memory... 

# registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))
# 
# print(paste0("Traffic list length: ", length(traffic)))
# 
# trafficcells <- foreach(i = 1:12,.packages = c("raster", "sf", "dplyr", "tidyr")) %dopar% {
#   print(paste0("Processing month: ", i))
#   t <- traffic[[i]] %>% st_set_precision(1e6) %>% st_make_valid()
#   df <- data.frame(years = c(), months =c(), cells=c(), traffic=c())
#   t$year <- substr(t$AIS_ID, start=11, stop =14)
#   t$month <- substr(t$AIS_ID, start=15, stop =16)
#   intersections_mat <- sf::st_intersects(t, iceconsf, sparse=F)
#   tin <- t[which(rowSums(intersections_mat) > 0),]
#   for(j in 1:length(tin$AIS_ID)){
#     if(j %% 1000 == 0){
#       print(paste0(round(j/length(tin$AIS_ID), 2)*100, "% finished."))
#     }
#     lenin <- st_intersection(tin[j,], iceconsf) %>% mutate(newlen=st_length(.)) %>% st_drop_geometry() %>% select(id, newlen)
#     dftemp <- data.frame(years=rep(t$year[1], length(lenin$id)),
#                          months=rep(t$month[1], length(lenin$id)),
#                          cells= lenin$id,
#                          traffic=lenin$newlen)
#   df <- rbind(df, dftemp)
#   }
#   return(df)
# }
# 
# trafficcells <- do.call(rbind, trafficcells)
# trafficcellsnew <- trafficcells %>% group_by(years, months, cells) %>% summarize(traffic_km = as.numeric(sum(traffic)/1000))
# 
# write.csv(trafficcellsnew, paste0("../Data_Processed/TrafficInIcePixels_", year, ".csv"))
# 
# 
# 
# (proc.time()-start)/60
# browseURL("https://www.youtube.com/watch?v=K1b8AhIsSYQ")

#################################################################################
# Method 2: Apply 
# Faster, but uses a lot of memory so it doesn't work on 2018-2020

print(paste0("Traffic list length: ", length(traffic)))

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
                      nOther = length(unique(MMSI_x[which(AIS_Typ == "Other")])),
                      CargoDist = sum(newlen[which(AIS_Typ == "Cargo")]),
                      FishDist = sum(newlen[which(AIS_Typ == "Fishing")]),
                      TankDist = sum(newlen[which(AIS_Typ == "Tanker")]),
                      OtherDist = sum(newlen[which(AIS_Typ == "Other")]))
  return(intdf)
}

registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))


trafficcells <- foreach(j = 1:12,.packages = c("raster", "sf", "dplyr", "tidyr")) %dopar% {
  traf <- traffic[[j]] %>% st_set_precision(1) %>% st_make_valid()
  df <- icepixellength(traf, iceconsf)
  return(df)
}
# Takes ~15 minutes to run

# Bind together all months of data
trafficcells <- do.call(rbind, trafficcells)
write.csv(trafficcells, paste0("../Data_Processed/TrafficInIcePixels_LenAndNShips", year, ".csv"))


(proc.time()-start)/60





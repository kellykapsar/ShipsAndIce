########################################################################
# TITLE: AIS Rasterization Script 
#
# DESCRIPTION: This script takes vector data processed as part of the 
# North Pacific and Arctic Vessel Traffic Data Set (2015-2020) and 
# rasterizes it based on the boundaries of the sea ice pixels from the 
# CryoSat2-SMOS merged product from the script "1-IceProcessing.R". 
#
# CREATED BY: Kelly Kapsar (kelly.kapsar@gmail.com)
# DATE CREATED: 2022-02
# DATE LAST MODIFIED: 2022-08-23
#
# NOTE: This script was designed to be run on one year of data at a time
# using Michigan State University's high performance computing system.
# A submission script used to process the data is title "2-AISRasterization.SB".
########################################################################

# Load libraries 
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





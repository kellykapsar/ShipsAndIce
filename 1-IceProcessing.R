
# Load Libraries
library(spatstat)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(sp)
library(dplyr)
library(ncdf4)
library(RNetCDF)
library(stars)

starttot <- proc.time()


wd <- "C:/Users/Kelly Kapsar/OneDrive - Michigan State University/Sync/3-ShipsAndIce/"

# Folder to save processed data 
savedsn <- paste0(wd, "Data_Processed/AIS")

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Study area (bounds for all of AIS data)
aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)

# Basemap
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% st_transform(AA) %>%  st_buffer(0)


#############################################################################
# NSIDC Ice Extent and AIS overlap
#############################################################################
# # Load in all shp files 
# filelist <- list.files(paste0(wd, "Data_Raw/AIS/"), pattern='.tif')
# 
# files <- lapply(filelist, function(x){raster::raster(paste0("../Data_Raw/AIS/", x))})
# 
# # Create indices to join all months
# idx <- rep(1:72, each=4)
# 
# # Create stack and combine all traffic types for each month into one layer
# allstack <- raster::stack(files)
# sumstack <- stackApply(allstack, indices = idx, fun=sum)
# 
# names(sumstack) <- seq(as.Date("2015-01-01"), as.Date("2020-12-01"), by="month")
# 
# raster::crs(sumstack) <- AA
# 
# # Load in sea ice extent data 
# extlist <- as.list(list.files(paste0(wd, "Data_Raw/SeaIceExtent_20220103/"), pattern='.shp'))
# 
# extdates <- lapply(extlist, function(x){as.Date(paste0(substr(x, start=10, stop=15), "01"), format="%Y%m%d")})
# 
# extsf <- lapply(extlist, function(x){st_read(paste0(wd, "/Data_Raw/SeaIceExtent_20220103/", x)) %>% st_transform(AA)})
# 
# 
# vesselinice <- function(trafficraster, extentsf){
#   trafficraster_new <- raster::mask(trafficraster, aisbounds)
#   extentsf_new <- st_crop(extentsf, aisbounds)
#   masktraffic <- raster::mask(trafficraster, extentsf)
#   return(masktraffic)
# }
# 
# shipinice <- lapply(1:72, function(x){vesselinice(sumstack[[x]], extsf[[x]])})
# shipinice <- raster::stack(shipinice)

########################################################################################################
# Load in all SMOS data, reproject, aggregate to monthly average values, and save as rasters 
########################################################################################################

smoslist <- as.list(list.files(paste0(wd, "Data_Raw/SeaIceThickness_SMOS_2015-2020/"), pattern='.nc'))

smosmonths <- lapply(smoslist, function(x){as.numeric(substr(x, start=33, stop=38))})

# Get long and lat matrices (Same for all layers)
ice <- nc_open(paste0(wd,"Data_Raw/SeaIceThickness_SMOS_2015-2020/", smoslist[[1]]))

# Save metadata to a text file
{
  sink(paste0("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", substr(smoslist[[1]], 1, 65),".txt"))
  print(ice)
  sink()
}

# Get projection information 
prj <- ncatt_get(ice, "Lambert_Azimuthal_Grid","proj4_string")$value

# Read lat lon and time for each observation
lon <- ncvar_get(ice, "lon")
lat <- ncvar_get(ice, "lat", verbose = F)

# Flip to correct orientation 
lon2 <- aperm(lon, c(2,1))
lat2 <- aperm(lat, c(2,1))

# https://gis.stackexchange.com/questions/385881/remote-sensing-products-raster-extent-not-agreeing-with-specified-projection-in
llv <- data.frame(long=c(lon2), lat=c(lat2))
llv = st_as_sf(llv, coords=1:2)
st_crs(llv)=4326
llpolar = st_transform(llv,prj)

nc_close()

# Function to process individual netcdfs into raster layers
smosprocess <- function(netcdf, newprj = AA, croparea = aisbounds, llpolar){
  # Open netcdf file 
  
  ice <- nc_open(paste0(wd,"Data_Raw/SeaIceThickness_SMOS_2015-2020/", netcdf))
  
  # Read lat lon and time for each observation
  t <- ncvar_get(ice, "time")
  
  # Read in data from the wind variable and verify the dimensions of the array
  thick.array <- ncvar_get(ice, "analysis_sea_ice_thickness") # 3dim array
  dim(thick.array)
  
  conc.array <- ncvar_get(ice, "sea_ice_concentration") # 3dim array
  dim(conc.array)
  
  # Identify fill value and replace with NA
  fillvalue <- ncatt_get(ice, "sea_ice_concentration","_FillValue")
  fillvalue
  
  thick.array[thick.array == fillvalue$value] <- NA
  conc.array[conc.array == fillvalue$value] <- NA
  
  # Get projection information 
  prj <- ncatt_get(ice, "Lambert_Azimuthal_Grid","proj4_string")$value
  
  # Close netcdf file
  nc_close(ice)
  
  # Flip to realign in correct direction
  conc.array2 <- aperm(conc.array, c(2,1))
  thick.array2 <- aperm(thick.array, c(2,1))
  
  # Make a raster of all ice concentration values 
  icethick <- raster(thick.array2)
  extent(icethick) <- extent(llpolar)
  crs(icethick) <- prj
  
  # icethickbs <- raster::projectRaster(icethick, crs=newprj, res = raster::res(icethick)) %>% 
  #   raster::crop(y=croparea)
  
  
  # Make a raster of all ice thickness values 
  icecon <- raster(conc.array2)
  extent(icecon) <- extent(llpolar)
  crs(icecon) <- prj
  
  # iceconbs <- raster::projectRaster(icecon, crs=newprj, res = raster::res(icecon)) %>% 
  #   raster::crop(y=croparea)
  
  
  # Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
  t2 <- as.POSIXct("1978-01-01 00:00")+as.difftime(t,units="secs")
  t2 <- format(t2, "%Y-%m-%d")
  
  names(icecon) <- t2
  names(icethick) <- t2
  
  return(list(icecon, icethick))
  
}

start <- proc.time()

output <- lapply(1:length(smoslist), function(x){smosprocess(smoslist[[x]], llpolar=llpolar)})

icecon <- sapply(output, "[[", 1) %>% 
  raster::stack() %>% 
  raster::crop(y=st_transform(aisbounds, prj)) %>% 
  raster::projectRaster(crs=AA, res = raster::res(.)) %>% 
  raster::stackApply(indices=smosmonths, fun=mean) %>% 
  stars::st_as_stars() 
icethick <- sapply(output, "[[", 2) %>% 
  raster::stack() %>% 
  raster::crop(y=st_transform(aisbounds, prj)) %>% 
  raster::projectRaster(crs=AA, res = raster::res(.)) %>% 
  raster::stackApply(indices=smosmonths, fun=mean) %>% 
  stars::st_as_stars() 


# Give attribute names to concentration and thickness 
icethick <- setNames(icethick, "thick")
icecon <- setNames(icecon, "con")

# Merge concentration and thickness into one stars object 
icestars <- c(icecon, icethick)

# Get band names (i.e. dates) from stars object 
dates <- st_get_dimension_values(icestars, "band")
# Convert to actual dates
newdates <-  as.Date(paste0(substr(dates, 7,10),"-",substr(dates,11,12),"-01"), format="%Y-%m-%d")
# Set third dimension of stars object to date 
icestars <- st_set_dimensions(icestars, 3, values = newdates, name = "date")


# Convert back to sf object 
icesf <- st_as_sf(icestars) 
# Prep new column names based on dates 
cols <- c(paste0("con_",newdates), paste0("thick_", newdates),"geometry")
# Rename columns (Otherwise it doesn't give the band names as columns for thickness, 
# not sure why - maybe to avoid duplicates?)
colnames(icesf) <- cols
# Create id (matches ID from AISRasterization.R script...)
icesf <- icesf %>% mutate(id=1:nrow(.)) 

write_stars(icestars, "../Data_Processed/Ice_SMOS.tif")
st_write(icesf, "../Data_Processed/Ice_SMOS.shp")

proc.time() - start

########################################################################################################
# Load in all Cryosat data, reproject, aggregate to monthly average values, and save as rasters 
########################################################################################################

cryolist <- as.list(list.files(paste0(wd, "Data_Raw/CryoSat-2_RDEFT4_20220131/"), pattern='.nc'))

cryodates <- lapply(cryolist, function(x){as.Date(substr(x, start=8, stop=15), format="%Y%m%d")})
cryomonths <- lapply(cryolist, function(x){as.numeric(substr(x, start=8, stop=13))})

# Get lat/lon values for all files
# Open netcdf file 

ice <- nc_open(paste0(wd,"Data_Raw/CryoSat-2_RDEFT4_20220131/", cryolist[[1]]))

# Save metadata to a text file
{
  sink(paste0("../Data_Raw/CryoSat-2_RDEFT4_20220131/", cryolist[[1]],".txt"))
  print(ice)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(ice, "lon")
lat <- ncvar_get(ice, "lat", verbose = F)

lon2 <- aperm(lon, c(2,1))
lat2 <- aperm(lat, c(2,1))

# Close netcdf file
nc_close(ice)

# Specify projection (from metadata)
polarstereo <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "

# https://gis.stackexchange.com/questions/385881/remote-sensing-products-raster-extent-not-agreeing-with-specified-projection-in
llv <- data.frame(long=c(lon2), lat=c(lat2))
llv = st_as_sf(llv, coords=1:2)
st_crs(llv)=4326
llpolar = st_transform(llv,polarstereo)


cryoprocess <- function(netcdf, llpolar, lyrname, newprj = AA, croparea = aisbounds){
  # Open netcdf file 
  
  ice <- nc_open(paste0(wd,"Data_Raw/CryoSat-2_RDEFT4_20220131/", netcdf))
  
  # Read in data from the wind variable and verify the dimensions of the array
  cryothick.array <- ncvar_get(ice, "sea_ice_thickness") # 3dim array
  
  cryoconc.array <- ncvar_get(ice, "ice_con") # 3dim array
  
  # Identify fill value and replace with NA
  cryothick.array[cryothick.array == -9999] <- NA
  cryoconc.array[cryoconc.array == -9999] <- NA
  
  # Close netcdf file
  nc_close(ice)
  
  
  # Flip to realign in correct direction
  cryoconc.array2 <- aperm(cryoconc.array, c(2,1))
  cryothick.array2 <- aperm(cryothick.array, c(2,1))
  
  polarstereo <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
  
  # Make a raster of all values 
  icecon <- raster(cryoconc.array2)
  extent(icecon) <- extent(llpolar)
  crs(icecon) <- polarstereo
  
  iceconbs <- raster::projectRaster(icecon, crs=AA) %>% 
    raster::crop(y=aisbounds)
  
  icethick <- raster(cryothick.array2)
  extent(icethick) <- extent(llpolar)
  crs(icethick) <- polarstereo
  
  names(icethick) <- lyrname
  names(icecon) <- lyrname
  
  return(list(icecon, icethick))
}

start <- proc.time()

# cryooutput <- lapply(1:length(cryolist), function(x){cryoprocess(cryolist[[x]], lyrname=cryodates[[x]], llpolar=llpolar)})

cryoicecons <- sapply(cryooutput, "[[", 1) %>% 
  raster::stack() %>% 
  raster::crop(y=st_transform(aisbounds, polarstereo)) %>% 
  raster::projectRaster(crs=AA, res = raster::res(.)) %>% 
  raster::stackApply(indices=cryomonths, fun=mean)


cryoicethicks <- sapply(cryooutput, "[[", 2) %>% 
  raster::stack() %>% 
  raster::crop(y=st_transform(aisbounds, polarstereo)) %>% 
  raster::projectRaster(crs=AA, res = raster::res(.)) %>% 
  raster::stackApply(indices=cryomonths, fun=mean)

writeRaster(cryoicecons, "../Data_Processed/IceConcentration_CryoSat.tif")
saveRDS(cryoicecons, "../Data_Processed/IceConcentration_CryoSat.rds")
writeRaster(cryoicethicks, "../Data_Processed/IceThickness_CryoSat.tif")
saveRDS(cryoicethicks, "../Data_Processed/IceThickness_CryoSat.rds")

proc.time() - start


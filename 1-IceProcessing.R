########################################################################
# TITLE: Ice Processing 
#
# DESCRIPTION: This script takes netcdf files from the CryoSat2-SMOS merged
# product and extracts sea ice thickness and concentration values for the 
# North Pacific study area for which we have AIS data, from 2015-2020. 
# The script then projects the data into Alaska Albers, and converts the 
# files into raster stacks. 
#
# CREATED BY: Kelly Kapsar (kelly.kapsar@gmail.com)
# DATE CREATED: 2022-02
# DATE LAST MODIFIED: 2022-08-23
########################################################################

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

# Folder to save processed data 
savedsn <- paste0("../Data_Processed/AIS")

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Study area (bounds for all of AIS data)
aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)

# Basemap
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% st_transform(AA) %>%  st_buffer(0)


########################################################################################################
# Load in all SMOS data, reproject, aggregate to monthly average values, and save as rasters 
########################################################################################################

smoslist <- as.list(list.files("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", pattern='.nc'))

smosmonths <- lapply(smoslist, function(x){as.numeric(substr(x, start=33, stop=38))})

# Get long and lat matrices (Same for all layers)
ice <- nc_open(paste0("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", smoslist[[1]]))

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

nc_close(ice)

# Function to process individual netcdfs into raster layers
smosprocess <- function(netcdf, newprj = AA, croparea = aisbounds, llpolar){
  # Open netcdf file 
  
  ice <- nc_open(paste0("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", netcdf))
  
  # Read lat lon and time for each observation
  t <- ncvar_get(ice, "time")
  
  # Read in data from the wind variable and verify the dimensions of the array
  thick.array <- ncvar_get(ice, "analysis_sea_ice_thickness") # 3dim array
  dim(thick.array)
  
  conc.array <- ncvar_get(ice, "sea_ice_concentration") # 3dim array
  dim(conc.array)
  
  unc.array <- ncvar_get(ice, "analysis_sea_ice_thickness_unc")
  dim(unc.array)
  
  # Identify fill value and replace with NA
  fillvalue <- ncatt_get(ice, "sea_ice_concentration","_FillValue")
  fillvalue
  
  thick.array[thick.array == fillvalue$value] <- NA
  conc.array[conc.array == fillvalue$value] <- NA
  unc.array[unc.array == fillvalue$value] <- NA
  
  # Get projection information 
  prj <- ncatt_get(ice, "Lambert_Azimuthal_Grid","proj4_string")$value
  
  # Close netcdf file
  nc_close(ice)
  
  # Flip to realign in correct direction
  conc.array2 <- aperm(conc.array, c(2,1))
  thick.array2 <- aperm(thick.array, c(2,1))
  unc.array2 <- aperm(unc.array, c(2,1))
  
  # Make a raster of all ice thickness values 
  icethick <- raster(thick.array2)
  extent(icethick) <- extent(llpolar)
  crs(icethick) <- prj
  
  # icethickbs <- raster::projectRaster(icethick, crs=newprj, res = raster::res(icethick)) %>% 
  #   raster::crop(y=croparea)
  
    # Make a raster of all ice concentration values 
  icecon <- raster(conc.array2)
  extent(icecon) <- extent(llpolar)
  crs(icecon) <- prj
  
  # iceconbs <- raster::projectRaster(icecon, crs=newprj, res = raster::res(icecon)) %>% 
  #   raster::crop(y=croparea)
  
  
  # Make a raster of all ice concentration values 
  iceunc <- raster(unc.array2)
  extent(iceunc) <- extent(llpolar)
  crs(iceunc) <- prj
  
  # iceconbs <- raster::projectRaster(icecon, crs=newprj, res = raster::res(icecon)) %>% 
  #   raster::crop(y=croparea)
  
  
  # Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
  t2 <- as.POSIXct("1978-01-01 00:00")+as.difftime(t,units="secs")
  t2 <- format(t2, "%Y-%m-%d")
  
  names(icecon) <- t2
  names(icethick) <- t2
  names(iceunc) <- t2
  
  return(list(icecon, icethick, iceunc))
  
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
iceunc <- sapply(output, "[[", 3) %>% 
  raster::stack() %>% 
  raster::crop(y=st_transform(aisbounds, prj)) %>% 
  raster::projectRaster(crs=AA, res = raster::res(.)) %>% 
  raster::stackApply(indices=smosmonths, fun=mean) %>% 
  stars::st_as_stars() 


# Give attribute names to concentration and thickness 
icethick <- setNames(icethick, "thick")
icecon <- setNames(icecon, "con")
incunc <- setNames(iceunc, "unc")

# Merge concentration and thickness into one stars object 
icestars <- c(icecon, icethick, iceunc)

# Get band names (i.e. dates) from stars object 
dates <- st_get_dimension_values(icestars, "band")
# Convert to actual dates
newdates <-  as.Date(paste0(substr(dates, 7,10),"-",substr(dates,11,12),"-01"), format="%Y-%m-%d")
# Set third dimension of stars object to date 
icestars <- st_set_dimensions(icestars, 3, values = newdates, name = "date")


# Convert back to sf object 
icesf <- st_as_sf(icestars) 
# Prep new column names based on dates 
cols <- c(paste0("con_",newdates), paste0("thick_", newdates), paste0("unc_", newdates),"geometry")
# Rename columns (Otherwise it doesn't give the band names as columns for thickness, 
# not sure why - maybe to avoid duplicates?)
colnames(icesf) <- cols
# Create id (matches ID from AISRasterization.R script...)
icesf <- icesf %>% mutate(id=1:nrow(.)) 

# Save output
write_stars(icestars, "../Data_Processed/Ice_SMOS.tif")
st_write(icesf, "../Data_Processed/Ice_SMOS_with_unc.shp")

proc.time() - start


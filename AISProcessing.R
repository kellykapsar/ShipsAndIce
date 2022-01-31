
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


# Load in all shp files 
filelist <- list.files(paste0(wd, "Data_Raw/AIS/"), pattern='.tif')

files <- lapply(filelist, function(x){raster::raster(paste0("../Data_Raw/AIS/", x))})

# Create indices to join all months
idx <- rep(1:72, each=4)

# Create stack and combine all traffic types for each month into one layer
allstack <- raster::stack(files)
sumstack <- stackApply(allstack, indices = idx, fun=sum)

names(sumstack) <- seq(as.Date("2015-01-01"), as.Date("2020-12-01"), by="month")

raster::crs(sumstack) <- AA

#############################################################################
# Load in all sea ice extent files
extlist <- as.list(list.files(paste0(wd, "Data_Raw/SeaIceExtent_20220103/"), pattern='.shp'))

extdates <- lapply(extlist, function(x){as.Date(paste0(substr(x, start=10, stop=15), "01"), format="%Y%m%d")})
  
extsf <- lapply(extlist, function(x){st_read(paste0(wd, "/Data_Raw/SeaIceExtent_20220103/", x)) %>% st_transform(AA)})


vesselinice <- function(trafficraster, extentsf){
  trafficraster_new <- raster::mask(trafficraster, aisbounds)
  extentsf_new <- st_crop(extentsf, aisbounds)
  masktraffic <- raster::mask(trafficraster, extentsf)
  return(masktraffic)
}

shipinice <- lapply(1:72, function(x){vesselinice(sumstack[[x]], extsf[[x]])})
shipinice <- raster::stack(shipinice)

# Load in all sea ice thickness files 
thicklist <- as.list(list.files(paste0(wd, "Data_Raw/SeaIceThickness_SMOS_2015-2020/"), pattern='.nc'))

thickdates <- lapply(thicklist, function(x){as.Date(paste0(substr(x, start=33, stop=40), "01"), format="%Y%m%d")})


# Open netcdf file 
ice <- nc_open("c:/Users/Kelly Kapsar/Downloads/downthemall/W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_20150131_20150206_r_v204_01_l4sit.nc")
ice <- open.nc(paste0(wd,"Data_Raw/SeaIceThickness_SMOS_2015-2020/", thicklist[[2]]))

# Save metadata to a text file
{
  sink('../Data_Raw/CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_16269119720374.txt')
  print(wind)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(wind, "lon")
lat <- ncvar_get(wind, "lat", verbose = F)
t <- ncvar_get(wind, "time")

head(lon)

# Read in data from the wind variable and verify the dimensions of the array
wind.array <- ncvar_get(wind, "wind_speed") # 3dim array
dim(wind.array)

# Identify fill value and replace with NA
fillvalue <- ncatt_get(wind, "wind_speed","_FillValue")
fillvalue

wind.array[wind.array == fillvalue$value] <- NA

# Close netcdf file
nc_close(wind)







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

########################################################################################################
# Load in all sea ice thickness files 
thicklist <- as.list(list.files(paste0(wd, "Data_Raw/SeaIceThickness_SMOS_2015-2020/"), pattern='.nc'))

thickdates <- lapply(thicklist, function(x){as.Date(paste0(substr(x, start=33, stop=40), "01"), format="%Y%m%d")})


# Open netcdf file 

ice <- nc_open(paste0(wd,"Data_Raw/SeaIceThickness_SMOS_2015-2020/", thicklist[[2]]))

# Save metadata to a text file
{
  sink(paste0("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", substr(thicklist[[2]], 1, 65),".txt"))
  print(ice)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(ice, "lon")
lat <- ncvar_get(ice, "lat", verbose = F)
t <- ncvar_get(ice, "time")

head(lon)

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


# Make a raster of all values 
icecon <- raster(conc.array2, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +datum=WGS84 +no_defs")) %>% 
  # flip(direction="x") %>%
  flip(direction="y") %>%
  raster::projectRaster(crs=prj) 

test <- raster::projectRaster(icecon, crs=crs(prj))

aisbounds_laea <- st_transform(aisbounds, prj) 

test <- raster::crop(icecon, aisbounds)

test <- raster::projectRaster(from=icecon, to=sumstack[[1]], res=raster::res(icecon))


  
  raster::projectRaster(crs=prj) 

# Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
t2 <- as.POSIXct("1900-01-01 00:00")+as.difftime(t,units="hours")
t2 <- format(t2, "%G-W%V")

######################################################################################
# Load in all sea ice thickness files 
cryolist <- as.list(list.files(paste0(wd, "Data_Raw/CryoSat-2_RDEFT4_20220131/"), pattern='.nc'))

cryodates <- lapply(cryolist, function(x){as.Date(substr(x, start=8, stop=15), format="%Y%m%d")})


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

# Read in data from the wind variable and verify the dimensions of the array
cryothick.array <- ncvar_get(ice, "sea_ice_thickness") # 3dim array
dim(thick.array)

cryoconc.array <- ncvar_get(ice, "ice_con") # 3dim array
dim(conc.array)

# Identify fill value and replace with NA
cryothick.array[cryothick.array == -9999] <- NA
cryoconc.array[cryoconc.array == -9999] <- NA

# Get projection information 
# prj <- ncatt_get(ice, "Lambert_Azimuthal_Grid","proj4_string")$value

# Close netcdf file
nc_close(ice)


# Flip to realign in correct direction
cryoconc.array2 <- aperm(cryoconc.array, c(2,1))
cryothick.array2 <- aperm(cryothick.array, c(2,1))
lon2 <- aperm(lon, c(2,1))
lat2 <- aperm(lat, c(2,1))

polarstereo <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "

# https://gis.stackexchange.com/questions/385881/remote-sensing-products-raster-extent-not-agreeing-with-specified-projection-in
llv <- data.frame(long=c(lon), lat=c(lat))
llv = st_as_sf(llv, coords=1:2)
st_crs(llv)=4326
llpolar = st_transform(llv,polarstereo)
llpolarcoords <- as.data.frame(st_coordinates(llpolar))

# Make a raster of all values 
icecon <- raster(cryoconc.array2)
extent(icecon) <- extent(llpolar)
crs(icecon) <- polarstereo

iceconbs <- raster::projectRaster(icecon, crs=AA) %>% 
  raster::crop(y=aisbounds)

icethick <- raster(cryothick.array2)
extent(icethick) <- extent(llpolar)
crs(icethick) <- polarstereo

icethickbs <- raster::projectRaster(icethick, crs=AA) %>% 
  raster::crop(y=aisbounds)

plot(icethickbs)


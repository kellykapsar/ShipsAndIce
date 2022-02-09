
library(raster)
library(stars)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

origmethod <- list.files("../Data_Processed/", pattern="TrafficInIce")
origmethod <- lapply(origmethod, function(x){read.csv(paste0("../Data_Processed/",x))})
origmethod <- do.call(rbind, origmethod)
origmethod$traffic_km <- origmethod$length/1000
origmethod <- origmethod %>% select(-length) 

loopmethod <- list.files("../Data_Processed/", pattern="LoopMethod")
loopmethod <- lapply(loopmethod, function(x){read.csv(paste0("../Data_Processed/",x))})
loopmethod <- do.call(rbind, loopmethod)
loopmethod <- loopmethod %>% rename(year=years, month=months, id=cells)

allpixels <- rbind(origmethod, loopmethod) %>% select(-X)
allpixels <- allpixels[-which(allpixels$month == 60),] # Manually remove one weird pixel value that idk what happened... 

# Check to make sure no duplicate cell values for each year month combo
qatest <- allpixels %>% group_by(year, month) %>% summarize(ncells=length(id), nuniquecells=length(unique(id)))
which(qatest$ncells != qatest$nuniquecells)


# Create data frames for thickness and concentration 
icesf <- st_read("../Data_Processed/Ice_SMOS.shp")
condf <-  icesf %>% 
            select(c(1:42), id) %>% 
            st_drop_geometry() %>%
            gather(key=yearmon, value=icecon, -id) %>%
            mutate(year = as.numeric(substr(yearmon, start=3, stop=6)),
                   month = as.numeric(substr(yearmon, start=8, stop=9)),
                   icecon = round(icecon, 2)) %>%
            select(-yearmon)

thickdf <-  icesf %>% 
            select(c(43:85), id) %>% 
            st_drop_geometry() %>%
            gather(key=yearmon, value=icethick, -id) %>%
            mutate(year = as.numeric(substr(yearmon, start=3, stop=6)),
                   month = as.numeric(substr(yearmon, start=8, stop=9)),
                   icethick = round(icethick, 2)) %>%
            select(-yearmon)

# Join thickness and concentration data frames 
allice <- left_join(condf, thickdf, by=c("year", "month", "id"))

# cor.test(allice$icethick, allice$icecon)

# Joine ice data with vessel traffic data 
alldf <-left_join(allice, allpixels, by=c("year", "month", "id"))
alldf <- alldf[,c("year", "month", "id", "icecon", "icethick","traffic_km")]

# alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$traffic_km <- round(alldf$traffic_km, 2)

# write.csv(alldf, "../Data_Processed/TrafficAndIcePixelValues.csv")

# Add traffic onto spatial object for icesf
trafficsf <- allpixels %>% 
  mutate(date=as.Date(paste0(year,"-",month,"-1"), format="%Y-%m-%d")) %>% 
  select(-year, -month) %>%  
  spread(date, traffic_km)

colnames(trafficsf) <- c(colnames(trafficsf[1]), 
                         paste0("v_", substr(colnames(trafficsf)[2:73], 1,4),".",substr(colnames(trafficsf[2:73]),6,7)))

allsf <- left_join(icesf, trafficsf, by=c("id"))

# alldf == non-spatial data frame (each row is a unique cell, month, year combo)
write.csv(alldf, "../Data_Processed/IceTrafficDataFrame.csv")
# all sf == spatial data frame (each row is a unique cell)
st_write(allsf, "../Data_Processed/IceTrafficDataFrame.shp")

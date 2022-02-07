
library(raster)
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

# Check to make sure no duplicate cell values for each year month combo
qatest <- allpixels %>% group_by(year, month) %>% summarize(ncells=length(cells), nuniquecells=length(unique(cells)))
which(qatest$ncells != qatest$nuniquecells)

# Add in ice data 
icecon <- readRDS("../Data_Processed/IceConcentration_SMOS.rds")
icecondf <- icecon %>% 
              stars::st_as_stars() %>% 
              st_as_sf() %>% 
              mutate(id=1:nrow(.)) %>% 
              st_drop_geometry() %>% 
              gather(key=yearmon, value=icecon, -id) %>% 
              mutate(year = as.numeric(substr(yearmon, start=7, stop=10)), 
                     month = as.numeric(substr(yearmon, start=11, stop=12)), 
                     icecon = round(icecon, 2)) %>% 
              select(-yearmon)

icethick <- readRDS("../Data_Processed/IceThickness_SMOS.rds")
icethickdf <- icethick %>% 
  stars::st_as_stars() %>% 
  st_as_sf() %>% 
  mutate(id=1:nrow(.)) %>% 
  st_drop_geometry() %>% 
  gather(key=yearmon, value=icethick, -id) %>% 
  mutate(year = as.numeric(substr(yearmon, start=7, stop=10)), 
         month = as.numeric(substr(yearmon, start=11, stop=12)), 
         icethick = round(icethick, 2)) %>% 
  select(-yearmon)


alldf <-left_join(icecondf, allpixels, by=c("year", "month", "id"))
alldf <- alldf[,c("year", "month", "id", "icecon", "traffic_km")]
# alldf$traffic_km[which(is.na(alldf$traffic_km))] <- 0
alldf$traffic_km <- round(alldf$traffic_km, 2)

aisinice <- alldf[which(alldf$icecon > 0 & alldf$traffic_km > 0),]

############################################################################
ggplot(allpixels, aes(x=traffic_km)) +
  geom_histogram(bins=50)

ggplot(icecondf, aes(x=icecon)) +
  geom_histogram(bins=50)

############################################################################
fit1 <- lm(traffic_km ~ icecon, data = aisinice)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(fit1)

ggplot(data=aisinice, aes(x=icecon, y=log(traffic_km))) +
  geom_point() +
  geom_smooth()

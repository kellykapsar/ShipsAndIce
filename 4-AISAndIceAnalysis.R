
library(raster)
library(stars)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(yarrr)
library(ggsn)
library(colorspace)
library(RColorBrewer)
library(scales)
library(viridis)


saveloc <- "../Figures/"

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Median sea ice extent 1981-2010
# ext <- st_read("../Data_Raw/median_extent_N_03_1981-2010_polyline_v3.0/median_extent_N_03_1981-2010_polyline_v3.0.shp") %>% 
#   st_transform(AA) %>% 
#   st_intersection(aisbounds) %>% 
#   st_write("../Data_Processed/MedianMarchIceExtent1981-2010.shp")

#######################################################################
# all sf == spatial data frame (each row is a unique cell)
allsf <- st_read("../Data_Processed/IceTrafficDataFrame.shp")
# studyoutline <- allsf %>% select(id) %>% st_union() %>% st_write("../Data_Processed/StudyOutline.shp", overwrite=T)

alldf <- read.csv("../Data_Processed/IceTrafficDataFrame.csv") %>% dplyr::select(-X)

portcells <- c(1983, 4214, 4168, 2202, 2053, 2054, 2442, 2443, 2516, 2517)

alldf <- alldf[-which(alldf$id %in% portcells),]
allsf <- allsf[-which(allsf$id %in% portcells),]

# Basemap
# aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% 
  st_transform(AA)
basemap.crop <- st_crop(basemap, st_buffer(allsf, 10000)) %>% 
  st_simplify(dTolerance=1000, preserveTopology = T)

#######################################################################
# Call in functions 

# DOESN'T WORK WITH MAP FUNCTION FOR SOME REASON.... 
# Couldn't figure out how to make map function work with calling column names, 
# so just made new functions for each column combo
pizzolato <- function(df, traffcol, icecol){
  traff <- df[{{traffcol}}]
  ice <- df[{{icecol}}]
  if(length(which(traff != 0)) == 0){
    return(NA)
  }
  maxice <- max(ice[which(traff != 0),])
  lotsaice <- which(ice > maxice & traff == 0)
  keeplotsaice <- which(ice[lotsaice,] == min(ice[lotsaice,]))
  lotsaicenew <- lotsaice[-keeplotsaice]
  traffnew <- traff[-lotsaicenew,]
  icenew <- ice[-lotsaicenew,]
  if(length(traffnew) > 0){
    cortest <- cor.test(traffnew, icenew, method="kendall")
    return(cortest)
  }
  if(length(traffnew) <= 0){
    return(NA)
  }
}

# Ice concentration and traffic 
pizzolato_conkm <- function(df){
  if(length(which(df$traffic_km != 0)) == 0){
    return(NA)
  }
  maxice <- max(df$icecon[which(df$traffic_km != 0)])
  lotsaice <- which(df$icecon > maxice & df$traffic_km == 0)
  dfnew <- df
  
  if(length(lotsaice) > 1){
    keeplotsaice <- which(df$icecon[lotsaice] == min(df$icecon[lotsaice]))
    lotsaicenew <- lotsaice[-keeplotsaice]
    dfnew <- df[-lotsaicenew,]
  }
  
  notrafforice <- dfnew[which(dfnew$traffic_km == 0 & dfnew$icecon == 0),]
  
  dfnewest <- dfnew
  
  if(length(notrafforice$year) > 0){
    dfnewest <- dfnew[-which(dfnew$traffic_km == 0 & dfnew$icecon == 0),]
    dfnewest <- rbind(dfnewest, notrafforice[1,])
  }
  
  if(length(dfnewest$traffic_km) >= 20){
    cortest <- cor.test(dfnewest$traffic_km, dfnewest$icecon, method="kendall")
    return(cortest)
  }
  if(length(dfnewest$traffic_km) < 20){
    return(NA)
  }
}

# Ice thickness and traffic
pizzolato_thickkm <- function(df){
  if(length(which(df$traffic_km != 0)) == 0){
    return(NA)
  }
  maxice <- max(df$icethick[which(df$traffic_km != 0)])
  lotsaice <- which(df$icethick > maxice & df$traffic_km == 0)
  keeplotsaice <- which(df$icethick[lotsaice] == min(df$icethick[lotsaice]))
  lotsaicenew <- lotsaice[-keeplotsaice]
  dfnew <- df[-lotsaicenew,]
  if(length(dfnew$traffic_km) > 0){
    cortest <- cor.test(dfnew$traffic_km, dfnew$icethick, method="kendall")
    return(cortest)
  }
  if(length(dfnew$traffic_km) <= 0){
    return(NA)
  }
}

# Ice concentration and number of ships
pizzolato_connship <- function(df){
  if(length(which(df$nShips != 0)) == 0){
    return(NA)
  }
  maxice <- max(df$icecon[which(df$nShips != 0)])
  lotsaice <- which(df$icecon > maxice & df$nShips == 0)
  keeplotsaice <- which(df$icecon[lotsaice] == min(df$icecon[lotsaice]))
  lotsaicenew <- lotsaice[-keeplotsaice]
  dfnew <- df[-lotsaicenew,]
  if(length(dfnew$nShips) > 0){
    cortest <- cor.test(dfnew$nShips, dfnew$icecon, method="kendall")
    return(cortest)
  }
  if(length(dfnew$nShips) <= 0){
    return(NA)
  }
}

# Ice thickness and number of ships
pizzolato_thicknship <- function(df){
  if(length(which(df$nShips != 0)) == 0){
    return(NA)
  }
  maxice <- max(df$icethick[which(df$nShips != 0)])
  lotsaice <- which(df$icethick > maxice & df$nShips == 0)
  keeplotsaice <- which(df$icethick[lotsaice] == min(df$icethick[lotsaice]))
  lotsaicenew <- lotsaice[-keeplotsaice]
  dfnew <- df[-lotsaicenew,]
  if(length(dfnew$nShips) > 0){
    cortest <- cor.test(dfnew$nShips, dfnew$icethick, method="kendall")
    return(cortest)
  }
  if(length(dfnew$nShips) <= 0){
    return(NA)
  }
}

plotKendall <- function(modsf, sigs, savename){
  p3 <- ggplot() +
    geom_sf(data=basemap.crop, fill="#f0f0f0", lwd=0) +
    # geom_sf(data=aisbounds, fill=NA, color="black", lwd=1)+
    geom_sf(data=modsf, aes(fill=estimate)) +
    colorspace::scale_fill_continuous_divergingx("RdBu", name="", rev=T) +
    geom_point(data=sigs, aes(x = sigs[,1], y = sigs[,2]), shape=8, size=0.1) +
    xlab("") +
    ylab("") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    blank()+
    theme(legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = c(0.1,0.2),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.ticks = element_blank(),
          axis.text=element_blank(),
          panel.background = element_rect(fill = "#cfd3d4"),
          panel.border =  element_rect(colour = "black"),
          panel.grid.major = element_line(colour = "transparent"))
  
  # p3
  ggsave(plot=p3, filename= paste0("../Figures/",savename,".png"), 
         width=8, height=8, units="in", dpi = 300)
}

#######################################################################


############################################################################
# Kendall Correlation and Map
############################################################################

alldfnest <- alldf %>% group_by(id) %>% nest()

### Ice con and total vessel traffic 
mods <- alldfnest %>% 
  mutate(cortest = map(data, function(df){pizzolato_conkm(df)}), 
         tidied = map(cortest, broom::tidy)) %>% 
  unnest(tidied)


modelresults <- allsf %>% dplyr::select(id) %>% left_join(., mods[,c("id", "estimate", "statistic", "p.value")])
sigcells <- st_coordinates(st_centroid(modelresults[which(modelresults$p.value < 0.05),])) %>% as.data.frame()
plotKendall(modelresults, sigcells, "KendallCorrelationMap_ConKm_NEW20220412")

# Percent of significant pixels with negative correlation coefficients 
sum(modelresults$estimate[which(modelresults$p.value < 0.05)] <0)/length(modelresults$estimate[which(modelresults$p.value < 0.05)])

# ### Ice con and number of ships
# mods <- alldfnest %>%
#   mutate(cortest = map(data, function(df){pizzolato_connship(df)}),
#          tidied = map(cortest, broom::tidy)) %>%
#   unnest(tidied)
# 
# 
# modelresults <- allsf %>% dplyr::select(id) %>% left_join(., mods[,c("id", "estimate", "statistic", "p.value")])
# sigcells <- st_coordinates(st_centroid(modelresults[which(modelresults$p.value < 0.05),])) %>% as.data.frame()
# plotKendall(modelresults, sigcells, "KendallCorrelationMap_ConnShips")


# ### Ice thickness and number of ships
# mods <- alldfnest %>%
#   mutate(cortest = map(data, function(df){pizzolato_thicknship(df)}),
#          tidied = map(cortest, broom::tidy)) %>%
#   unnest(tidied)
# 
# 
# modelresults <- allsf %>% dplyr::select(id) %>% left_join(., mods[,c("id", "estimate", "statistic", "p.value")])
# sigcells <- st_coordinates(st_centroid(modelresults[which(modelresults$p.value < 0.05),])) %>% as.data.frame()
# plotKendall(modelresults, sigcells, "KendallCorrelationMap_ThicknShips")
# 
# ### Ice thickness and total vessel traffic
# mods <- alldfnest %>%
#   mutate(cortest = map(data, function(df){pizzolato_thickkm(df)}),
#          tidied = map(cortest, broom::tidy)) %>%
#   unnest(tidied)
# 
# 
# modelresults <- allsf %>% dplyr::select(id) %>% left_join(., mods[,c("id", "estimate", "statistic", "p.value")])
# sigcells <- st_coordinates(st_centroid(modelresults[which(modelresults$p.value < 0.05),])) %>% as.data.frame()
# plotKendall(modelresults, sigcells, "KendallCorrelationMap_ThickKm")

############################################################################
# Case study 
############################################################################
# 
alldf %>% filter(id == 1626) %>% group_by(year) %>% summarize(traffic_km = sum(traffic_km))

############################################################################
# Ice and traffic by month line plot 
############################################################################

monthstats <- alldf %>% group_by(year, month) %>% summarize(iceext = sum(ifelse(icecon >= 15, 1, 0))*622.109950, traffic_km = sum(traffic_km), nShips = sum(nShips))

monthstats$timestep <- ifelse(monthstats$month == 10, 1, 
                              ifelse(monthstats$month == 11, 2, 
                                     ifelse(monthstats$month == 12, 3, 
                                            ifelse(monthstats$month == 1, 4, 
                                                   ifelse(monthstats$month == 2, 5,
                                                          ifelse(monthstats$month == 3, 6,
                                                                 ifelse(monthstats$month == 4, 7, NA)))))))

monthstats$traffic_1kkm <- monthstats$traffic_km/1000

# Ice extent (average annual and maximum annual)
annualiceext <- monthstats %>% group_by(year) %>% summarize(maxext = max(iceext))
otherext <- annualiceext$maxext[which(annualiceext$maxext != min(annualiceext$maxext))]
minext <- annualiceext$maxext[which(annualiceext$maxext == min(annualiceext$maxext))]
round(((otherext-minext)/otherext)*100,2) # Percent difference in maximum ice extent between highest and lowest extent years 
studyareasize_km2 <- sum(st_area(allsf))/1000000
minext/studyareasize_km2


coeff <- 2000 # median(monthstats$iceext/monthstats$traffic_1kkm)
ggplot() +
  geom_line(data=monthstats, aes(x=timestep, y =iceext, group=year, col=as.factor(year)),lwd=1.5,lty=2) +
  geom_line(data=monthstats, aes(x=timestep, y = traffic_1kkm*coeff, group=year, col=as.factor(year)),lwd=1.5) +
  scale_color_manual(labels=2015:2020, values=rev(viridis::viridis(6)), name="Year") + 
  xlab("") +
  scale_x_continuous(breaks=1:7, labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ) +
  scale_y_continuous(name = expression("Sea Ice Extent "~(km^2)),
                     sec.axis = sec_axis(~./coeff, name = "Total Vessel Traffic (1000s km)")) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# ggsave("../Figures/IceAndAIS_MonthlyLines.png", width=9, height=6, units="in")


ggplot() +
  geom_line(data=monthstats, aes(x=timestep, y =iceext, group=year, col=as.factor(year)),lwd=1.5,lty=2) +
  geom_line(data=monthstats, aes(x=timestep, y = traffic_1kkm*coeff, group=year), alpha=0,lwd=1.5) +
  scale_color_manual(labels=2015:2020, values=rev(viridis::viridis(6)), name="Year") + 
  xlab("") +
  scale_x_continuous(breaks=1:7, labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ) +
  scale_y_continuous(name = expression("Sea Ice Extent "~(km^2)),
                     sec.axis = sec_axis(~./coeff, name = "Total Vessel Traffic (1000s km)")) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("../Figures/IceAndAIS_MonthlyLines_IceOnly.png", width=9, height=6, units="in")

############################################################################
# Traffic in ice by months - LINE CHART 
############################################################################

inice <- alldf[which(alldf$icecon > 0 & alldf$traffic_km > 0),]

monthinice <- inice %>% group_by(year, month) %>% summarize(traffic_km = sum(traffic_km))

monthinice$timestep <- ifelse(monthinice$month == 10, 1, 
                              ifelse(monthinice$month == 11, 2, 
                                     ifelse(monthinice$month == 12, 3, 
                                            ifelse(monthinice$month == 1, 4, 
                                                   ifelse(monthinice$month == 2, 5,
                                                          ifelse(monthinice$month == 3, 6,
                                                                 ifelse(monthinice$month == 4, 7, NA)))))))


ggplot() +
  geom_line(data=monthinice, aes(x=timestep, y =traffic_km, group=year, col=as.factor(year)),lwd=1.5) +
  scale_color_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") + 
  xlab("") +
  scale_x_continuous(breaks=1:7, labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ) +
  scale_y_continuous(name = "Total Vessel Traffic\nin Ice (km)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# ggsave("../Figures/AISInIce_MonthlyLines.png", width=9, height=6, units="in")

##################################################################################
# Annual traffic totals - BAR CHARTS 
##################################################################################

# Annual traffic totals in ice 
yearinMIZ <- alldf %>% filter(icecon > 0, icecon < 80, traffic_km > 0) %>% 
  group_by(year) %>% 
  summarize(traffic_km = sum(traffic_km),
            totcarg_km = sum(carg_km), 
            totfish_km = sum(fish_km),
            tottank_km = sum(tank_km), 
            totother_km = sum(other_km),
            ncells=n())

ggplot(yearinMIZ, aes(x=year, y=traffic_km, fill=as.factor(year))) +
  geom_bar(stat="identity") +
  scale_fill_manual(labels=2015:2020, values=rev(viridis::viridis(6)), name="Year") +
  scale_x_continuous(breaks= 2015:2020) +
  scale_y_continuous(breaks=seq(0,300000, by=50000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic\nin Marginal Ice Zone (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggsave("../Figures/AISInMIZ_AnnualBar.png", width=9, height=6, units="in")

yearinMIZlong <- yearinMIZlong %>% 
  select(year, totcarg_km, tottank_km, totfish_km, totother_km) %>% 
  gather(-year, key=type, value=distance)

ggplot(yearinMIZlong, aes(x=year, y=distance, fill=type)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=qualitative_hcl(4, palette = "Dynamic"), name="Vessel Type", labels=c("Cargo", "Fishing", "Other", "Tanker")) +
  scale_x_continuous(breaks= 2015:2020) +
  # scale_y_continuous(breaks=seq(0,50000, by=10000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic\nin Marginal Ice Zone (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave("../Figures/AISInMIZIce_Types_AnnualBar.png", width=9, height=6, units="in")

### Total vessel traffic in cells with >50% ice concentration 
yearinlotsaice <- alldf[which(alldf$icecon > 80 & alldf$traffic_km > 0),] %>% 
  group_by(year) %>%
  summarize(traffic_km = sum(traffic_km),
            totcarg_km = sum(carg_km), 
            totfish_km = sum(fish_km),
            tottank_km = sum(tank_km), 
            totother_km = sum(other_km),
            ncells=n())

ggplot(yearinlotsaice, aes(x=year, y=traffic_km, fill=as.factor(year))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=rev(viridis::viridis(6))) +
  scale_x_continuous(breaks= 2015:2020) +
  scale_y_continuous(breaks=seq(0,50000, by=10000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic\nin Pack Ice (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")
ggsave("../Figures/AISIn80Ice_AnnualBar.png", width=9, height=6, units="in")

yearinlotsaicelong <- yearinlotsaice %>% 
  select(year, totcarg_km, tottank_km, totfish_km, totother_km) %>% 
  gather(-year, key=type, value=distance)


ggplot(yearinlotsaicelong, aes(x=year, y=distance, fill=type)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=qualitative_hcl(4, palette = "Dynamic"), name="Vessel Type", labels=c("Cargo", "Fishing", "Other", "Tanker")) +
  scale_x_continuous(breaks= 2015:2020) +
  scale_y_continuous(breaks=seq(0,50000, by=10000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic\nin Pack Ice (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave("../Figures/AISIn80Ice_Types_AnnualBar.png", width=9, height=6, units="in")

##################################################################################
# Growth in vessel traffic by ice con 
##################################################################################

# Change in vessel traffic by ice concentration 
totaltraff <- alldf %>% group_by(year) %>% summarize(traffic_km = sum(traffic_km))
pctchangetotal <- round((totaltraff$traffic_km[totaltraff$year == 2020]-totaltraff$traffic_km[totaltraff$year == 2015])/totaltraff$traffic_km[totaltraff$year == 2020]*100,2)
pctchangeinlotsaice <- round((yearinlotsaice$traffic_km[yearinlotsaice$year == 2020]-yearinlotsaice$traffic_km[yearinlotsaice$year == 2015])/yearinlotsaice$traffic_km[yearinlotsaice$year == 2020]*100,2)
pctchangeinice <- round((yearinice$traffic_km[yearinice$year == 2020]-yearinice$traffic_km[yearinice$year == 2015])/yearinice$traffic_km[yearinice$year == 2020]*100,2)


# Growth in vessel traffic by ice concentration 
growthbyicecon <- alldf %>% 
  filter(year %in% c(2015, 2020)) %>% 
  mutate(icecon = round(icecon, -1)) %>% 
  group_by(icecon, year) %>% 
  summarize(traffic_km = sum(traffic_km)) %>% 
  spread(key=year, value=traffic_km) %>% 
  mutate(pctchange = round((`2020`-`2015`)/`2020`*100,2))

ggplot(growthbyicecon, aes(x=icecon, y=pctchange)) +
  geom_line() +
  scale_x_continuous(breaks= seq(0,100,10)) +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Change in vessel activity\n(2015-2020)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
# ggsave("../Figures/AISIbyIceCon_PctChange.png", width=9, height=6, units="in")

### Total vessel traffic in cells with 100% ice concentration is zero!  
yearinallice <- alldf[which(alldf$icecon >= 100  & alldf$traffic_km > 0),] %>% 
  group_by(year) %>% 
  summarize(traffic_km = sum(traffic_km))

##################################################################################
# Vessel traffic in occupied cells x ice concentration - BOXPLOT
##################################################################################

trafficbyicecon <- alldf[which(alldf$icecon > 0  & alldf$traffic_km > 0),] %>% 
  filter(icecon > 0) %>% 
  mutate(icecon = round(icecon, -1)) 

ggplot(trafficbyicecon, aes(x=icecon, y=traffic_km, group=icecon)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Vessel Traffic\nin Occupied Cells (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
# ggsave("../Figures/VesselTrafficbyIceCon_Boxplot.png", width=9, height=6, units="in")

##################################################################################
# Proportion of cells occupied x ice concentration - LINE PLOT 
##################################################################################

propocc <- alldf %>%  
  mutate(icecon = round(icecon, 0)) %>% 
  group_by(icecon) %>% 
  summarize(ncells=n(), nocc=length(which(nShips > 0))) %>% 
  mutate(propocc = round(nocc/ncells*100, 2))

propocc$propocc[which(propocc$icecon == 0)]

propocc <- propocc %>% filter(icecon > 0)

cor.test(propocc$propocc, propocc$icecon, method= "spearman")

ggplot(propocc, aes(x=icecon, y=propocc)) +
  geom_line() +
  scale_x_continuous(breaks= seq(0,100,10), expand = c(0, 0), limits = c(0,NA)) +
  scale_y_continuous(breaks=seq(0,60,10), expand = c(0, 0), limits = c(0,60)) +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Occupied Cells (%)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        axis.line = element_line())
# ggsave("../Figures/PropOccbyIceCon_Lineplot.png", width=9, height=6, units="in")

##################################################################################
# Total traffic by ice concentration - LINE PLOT
##################################################################################

# Line plot: Mean vessel traffic in occupied cells x ice concentration 
trafficbyiceconsumms <- alldf[which(alldf$icecon > 0  & alldf$traffic_km > 0),] %>% 
  filter(icecon > 0) %>% 
  mutate(icecon = round(icecon, 0)) %>% 
  group_by(icecon) %>% 
  summarize(traffic_km = median(traffic_km), nShips= median(nShips), ncells=n())

ggplot(trafficbyiceconsumms, aes(x=icecon, y=nShips)) +
  geom_line() +
  scale_x_continuous(breaks= seq(0,100,10)) +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Median Vessel Traffic\nin Occupied Cells (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
# ggsave("../Figures/nShipsbyIceCon.png", width=9, height=6, units="in")

test <- alldf %>% 
  filter(icecon > 0) %>% 
  mutate(icecon = round(icecon, -1)) %>% 
  group_by(icecon) %>% 
  summarize(traf_km = mean(traffic_km), 
            tottraff_km = sum(traffic_km),
            totcarg_km = sum(carg_km), 
            totfish_km = sum(fish_km),
            tottank_km = sum(tank_km), 
            totother_km = sum(other_km),
            ncells=n())

testlong <- test %>% select(icecon, totcarg_km, tottank_km, totfish_km, totother_km) %>% gather(-icecon, key=type, value=distance)
testlong %>% group_by(type) %>% summarize(distance=sum(distance))

ggplot(testlong, aes(x=icecon, y=distance, group=type, color=type)) +
  geom_line( lwd=2) +
  scale_color_manual(values=qualitative_hcl(4, palette = "Dynamic"), name="Vessel Type", labels=c("Cargo", "Fishing", "Other", "Tanker")) +
  scale_x_continuous(breaks= seq(0,100,10), expand = c(0, 0), limits = c(0,NA)) +
  # scale_y_continuous(breaks=seq(0,250000,50000), expand=c(0,0),limits=c(0, 200000)) +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Vessel Traffic\nin Sea Ice (km)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave("../Figures/VesselTrafficbyIceCon_Types_Lineplot.png", width=9, height=6, units="in")

ggplot(test, aes(x=icecon, y=tottraff_km)) +
  geom_line() +
  scale_x_continuous(breaks= seq(0,100,10), expand = c(0, 0), limits = c(0,NA)) +
  scale_y_continuous(breaks=seq(0,50000,10000), expand = c(0, 0), limits=c(0, NA)) + 
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Vessel Traffic\nin Occupied Cells (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        panel.border = element_blank(), 
        axis.line = element_line())
ggsave("../Figures/VesselTrafficbyIceCon_Lineplot.png", width=9, height=6, units="in")

##################################################################################
# Total traffic by ice concentration - MAPS 
##################################################################################

#### Marginal ice zone 
iceMIZpixels <- alldf %>% filter(icecon > 15, icecon < 80, traffic_km > 0) %>% 
  group_by(id) %>% 
  summarize(traffic_km = sum(traffic_km))
iceMIZsf <- allsf %>% filter(id %in% iceMIZpixels$id) %>% mutate(iceMIZtraff = iceMIZpixels$traffic_km)

min = min(iceMIZsf$iceMIZtraff)
max = max(iceMIZsf$iceMIZtraff)
diff <- max - min
std = sd(iceMIZsf$iceMIZtraff)

equal.interval = round(seq(min, max, by = diff/6), 0)
quantile.interval = round(quantile(iceMIZsf$iceMIZtraff, probs=seq(0, 1, by = 1/6)), 0)
std.interval = round(c(seq(min, max, by=std), max), 0)
natural.interval = round(classInt::classIntervals(iceMIZsf$iceMIZtraff, n = 6, style = 'jenks')$brks,0)

iceMIZsf$iceMIZtraff.equal = cut(iceMIZsf$iceMIZtraff, breaks=equal.interval, include.lowest = TRUE)
iceMIZsf$iceMIZtraff.quantile = cut(iceMIZsf$iceMIZtraff, breaks=quantile.interval, include.lowest = TRUE)
iceMIZsf$iceMIZtraff.std = cut(iceMIZsf$iceMIZtraff, breaks=std.interval, include.lowest = TRUE)
iceMIZsf$iceMIZtraff.natural = cut(iceMIZsf$iceMIZtraff, breaks=natural.interval, include.lowest = TRUE)

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=st_buffer(allsf$geometry, 10000), fill=NA) +
  geom_sf(data=iceMIZsf,aes(fill=iceMIZtraff.quantile)) +
  scale_fill_manual(values=brewer.pal(7,"YlOrRd"), name="Total Traffic (km)") +
  # geom_sf(data=ice50sf$geometry[ice50sf$ice50traff > 1000], fill="red") +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/IniceMIZ_Quantile.png",
       width=10, height=8, units="in", dpi=300)

## Pack ice 
ice80pixels <- alldf %>% filter(icecon > 80, traffic_km > 0) %>% group_by(id) %>% summarize(traffic_km = sum(traffic_km))
ice80sf <- allsf %>% filter(id %in% ice80pixels$id) %>% mutate(ice80traff = ice80pixels$traffic_km)

min = min(ice80sf$ice80traff)
max = max(ice80sf$ice80traff)
diff <- max - min
std = sd(ice80sf$ice80traff)

equal.interval = round(seq(min, max, by = diff/6), 0)
quantile.interval = round(quantile(ice80sf$ice80traff, probs=seq(0, 1, by = 1/6)), 0)
std.interval = round(c(seq(min, max, by=std), max), 0)
natural.interval = round(classInt::classIntervals(ice80sf$ice80traff, n = 6, style = 'jenks')$brks,0)

ice80sf$ice80traff.equal = cut(ice80sf$ice80traff, breaks=equal.interval, include.lowest = TRUE)
ice80sf$ice80traff.quantile = cut(ice80sf$ice80traff, breaks=quantile.interval, include.lowest = TRUE)
ice80sf$ice80traff.std = cut(ice80sf$ice80traff, breaks=std.interval, include.lowest = TRUE)
ice80sf$ice80traff.natural = cut(ice80sf$ice80traff, breaks=natural.interval, include.lowest = TRUE)

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  # geom_sf(data=allsf$geometry, fill=NA) +
  geom_sf(data=ice80sf,aes(fill=ice80traff.quantile)) +
  scale_fill_manual(values=brewer.pal(7,"YlOrRd"), name="Total Traffic (km)") +
  # geom_sf(data=ice50sf$geometry[ice50sf$ice50traff > 1000], fill="red") +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/Inice80_Quantile.png",
       width=10, height=8, units="in", dpi=300)

##################################################################################
# Total traffic in ice - MAPS 
##################################################################################

#### Marginal ice zone 
icepixels <- alldf %>% filter(icecon > 15, traffic_km > 0) %>% 
  group_by(id) %>% 
  summarize(traffic_km = sum(traffic_km))
icesf <- allsf %>% filter(id %in% icepixels$id) %>% mutate(icetraff = icepixels$traffic_km)

min = min(icesf$icetraff)
max = max(icesf$icetraff)
diff <- max - min
std = sd(icesf$icetraff)

equal.interval = round(seq(min, max, by = diff/6), 0)
quantile.interval = round(quantile(icesf$icetraff, probs=seq(0, 1, by = 1/6)), 0)
std.interval = round(c(seq(min, max, by=std), max), 0)
natural.interval = round(classInt::classIntervals(icesf$icetraff, n = 6, style = 'jenks')$brks,0)

icesf$icetraff.equal = cut(icesf$icetraff, breaks=equal.interval, include.lowest = TRUE)
icesf$icetraff.quantile = cut(icesf$icetraff, breaks=quantile.interval, include.lowest = TRUE)
icesf$icetraff.std = cut(icesf$icetraff, breaks=std.interval, include.lowest = TRUE)
icesf$icetraff.natural = cut(icesf$icetraff, breaks=natural.interval, include.lowest = TRUE)

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=icesf,aes(fill=icetraff.quantile)) +
  scale_fill_manual(values=brewer.pal(7,"YlOrRd"), name="Total Traffic (km)") +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/Inice_Quantile.png",
       width=10, height=8, units="in", dpi=300)

# Read in bowhead concentration areas
whale<- st_read("../Data_Raw/Bowhead_RelAbund_Winter/Bowhead_RelAbund_Winter.shp")
hicon <- whale %>% filter(area_desc == "High_Conc - Winter")
con <- whale %>% filter(area_desc == "Concentration - Winter")
rang <- whale %>% filter(area_desc == "Range - Winter")

ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=hicon, fill="darkblue", alpha = 0.5, lwd=0) +
  geom_sf(data=con, fill="blue", alpha = 0.5, lwd=0) +
  geom_sf(data=rang, fill="lightblue", alpha = 0.5, lwd=0) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "left",
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

ggsave(filename= "../Figures/Bowhead_Winter.png",
       width=10, height=8, units="in", dpi=300)

############################################################################
# Histograms and distribution fitting 
############################################################################

library(fitdistrplus)
# Traffic 
ggplot(inice, aes(x=traffic_km)) +
  geom_histogram(bins=50)

fittraff <- fitdist(inice$traffic_km, distr="gamma", method="mle")
summary(fittraff)
plot(fittraff)

# Number of ships 
ggplot(inice, aes(x=nShips)) +
  geom_histogram(bins=50)

fitships <- fitdist(inice$nShips, distr="gamma", method="mle")
summary(fitships)
plot(fitships)

# Ice concentration 
ggplot(inice, aes(x=icecon)) +
  geom_histogram(bins=50)

fiticecon <- fitdist(inice$icecon, distr="gamma", method="mle")
summary(fiticecon)
plot(fiticecon)

# Ice thickness
ggplot(inice, aes(x=icethick)) +
  geom_histogram(bins=50)

thick <- inice[!is.na(inice$icethick),]

fitthick <- fitdist(thick$icethick, distr="gamma", method="mle")
summary(fitthick)
plot(fitthick)

############################################################################
# Map of study area 
############################################################################ 

p3 <- ggplot() +
  # geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  geom_sf(data=aisbounds, fill=NA, color="black", lwd=1)+
  geom_sf(data=allsf[which(allsf$c_2020_03 > 1),], fill="green") +
  geom_sf(data=allsf[which(allsf$id %in% highshipids),], fill="red")+
  geom_sf(data=allsf[which(allsf$id %in% hightraffids),], fill="blue")+
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.05,0.125),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

p3
ggsave(plot=p3, filename= "TestMap.png", 
       width=8, height=8, units="in")

############################################################################
# Pixel timeseries 
############################################################################
alldf$date <- as.Date(paste0(alldf$year, "-", alldf$month, "-1"), format="%Y-%m-%d")

pixeltimeseries <- function(longdata, var, varlab, savedata=TRUE, plotdata=TRUE, name, saveloc){
  # Create custom labels for plot
  df.labels <- longdata[longdata$id == 1, ] %>% mutate( label = ifelse( month == 1, 
                                                                        format( date, "%Y"), 
                                                                        "" ) ) 
  
  p1 <- ggplot(longdata, aes_string(x = "date", group="id", y=var)) + 
    scale_x_date(breaks = df.labels$date,  labels =  df.labels$label, expand=c(0,1)) +
    scale_y_continuous(expand= c(0, 0)) +
    xlab(label = "") + 
    ylab(label = varlab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          text = element_text(size = 25), 
          plot.margin = margin(10, 5, 5, 5),
          panel.grid = element_blank(),
          axis.ticks.length=unit(.15, "in")) +
    geom_line(size=0.2, alpha=0.4) 
  
  ifelse(plotdata == TRUE, return(p1), NA)
  
  ifelse(savedata == TRUE, ggsave(plot=p1, filename= paste0(saveloc,"PixelLinePlot_",name,".png"), 
                                  width=12, height=8, units="in"), NA)
}

pixeltimeseries(alldf, plotdata = F, savedata=T, saveloc=saveloc, 
                name="Traffic_km", var= "traffic_km", varlab = "Total Vessel Traffic (km)")

pixeltimeseries(alldf, plotdata = F, savedata=T, saveloc=saveloc, 
                name="nShips", var= "nShips", varlab = "Number of Ships")


############################################################################
# Non-parametric correlations across all data 
############################################################################

cor.test(inice$traffic_km, inice$icecon, method="kendall")
cor.test(inice$traffic_km, inice$icethick, method="kendall")
cor.test(inice$nShips, inice$icecon, method="kendall")
cor.test(inice$nShips, inice$icethick, method="kendall")

hist(log(inice$traffic_km))
inice$log_traffic_km <- log(inice$traffic_km)
inice$sqrt_icecon <- sqrt(inice$icecon)

fit1 <- lm(log_traffic_km ~ sqrt_icecon, data = inice)

assumps <- augment(fit1)

qqnorm(assumps$.std.resid)
qqline(assumps$.std.resid)


library(nortest)
ad.test(assumps$.std.resid)
ggplot(data = assumps, aes(y = .std.resid, x = .fitted)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -2) +
  geom_hline(yintercept = 2) +
  ggtitle("Standardized Residuals vs Fits")



test <- glm(traffic_km ~ icecon, data=inice, family=Gamma)
t <- augment(test)
summary(test)

library(segmented)
test <- segmented(fit1, npsi=1)

# breakpoint
test$psi

# Slope
slope(test)

fitteddata <- fitted(test)
my.model <- data.frame(sqrticecon = inice$sqrt_icecon, logtraff=fitteddata)





ggplot() + 
  geom_point(data=inice, aes(x=sqrt_icecon, y=log_traffic_km)) +
  geom_line(data=my.model, aes(x= sqrticecon, y=logtraff), colour="red", lwd=2)




# https://rpubs.com/mpfoley73/495822#:~:text=In%20general%2C%20transforming%20the%20response,predictor%20variables%20corrects%20non%2Dlinearity.
fit2 <- glm(log_traffic_km ~ icethick, data = inice)
fit3 <- lm(log_traffic_km ~ icecon + icethick, data = inice)



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
ggplotRegression(fit2)


ggplot(data=inice, aes(x=icecon, y=log(traffic_km))) +
  geom_point() +
  geom_smooth()

ggplot(data=inice, aes(x=icecon, y=traffic_km)) +
  geom_point() +
  geom_smooth()

ggplot(data=inice, aes(x=icethick, y=log(traffic_km))) +
  geom_point() +
  geom_smooth()

ggplot(data=inice, aes(x=icethick, y=nShips)) +
  geom_point() +
  geom_smooth()

ggplot(data=inice, aes(x=icecon, y=nShips)) +
  geom_point() +
  geom_smooth()


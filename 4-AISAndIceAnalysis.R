
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


saveloc <- "../Figures/"

# Projection (Alaska Albers)
AA <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Basemap
aisbounds <- st_read("../Data_Raw/ais_reshape.shp") %>% st_transform(AA)
basemap <- read_sf("../Data_Raw/AK_CAN_RUS/AK_CAN_RUS.shp") %>% 
  st_transform(AA)
basemap.crop <- st_crop(basemap, st_buffer(aisbounds, 100000)) %>% 
  st_simplify(dTolerance=1000, preserveTopology = T)

# Median sea ice extent 1981-2010
# ext <- st_read("../Data_Raw/median_extent_N_03_1981-2010_polyline_v3.0/median_extent_N_03_1981-2010_polyline_v3.0.shp") %>% 
#   st_transform(AA) %>% 
#   st_intersection(aisbounds) %>% 
#   st_write("../Data_Processed/MedianMarchIceExtent1981-2010.shp")

# all sf == spatial data frame (each row is a unique cell)
allsf <- st_read("../Data_Processed/IceTrafficDataFrame.shp")





# plot(basemap.crop$geometry)
# plot(test$geometry, add=T, col="red")
# plot(allsf$geometry, add=T, border="green", col=NA)
# 
# 
# linecells <- unlist(st_intersects(st_union(ext), allsf))
# plot(allsf$geometry[linecells], add=T, border="purple", col=NA)
# 
# online <- as.data.frame(st_coordinates(st_centroid(allsf$geometry[linecells])))
# 
# allcoords <- as.data.frame(st_coordinates(st_centroid(allsf)))
# 
# aboveline <- function(linecentroids, cellcentroid){
#   liney <- linecentroids[which(linecentroids$X == cellcentroid$X),]
# }
# 
# 
# plot(allsf$geometry[1], col="red", add=T)
# length(unique(t$X))
# length(unique(t$Y))
# 
# 
# length(which(round(allcoords$X,0) %in% round(online$X, 0)))
# 
# 
# plot(allsf$geometry[which(round(allcoords$X,0) %in% round(online$X, 0))], add=T, border="red")


alldf <- read.csv("../Data_Processed/IceTrafficDataFrame.csv") %>% dplyr::select(-X)


fulldf <- expand.grid(year=2015:2020, month=1:12, id=unique(alldf$id))
fulldf <- left_join(fulldf, alldf, by=c("year", "month", "id"))

fulldf$traffic_km[which(is.na(fulldf$traffic_km))] <- 0
fulldf$nShips[which(is.na(fulldf$nShips))] <- 0
fulldf$nCargo[which(is.na(fulldf$nCargo))] <- 0
fulldf$nTank[which(is.na(fulldf$nTank))] <- 0
fulldf$nFish[which(is.na(fulldf$nFish))] <- 0
fulldf$nOther[which(is.na(fulldf$nOther))] <- 0


# DECIDED NOT TO RUN TREND TEST BECAUSE WOULDN'T EXPECT CLIMATE-RELATED TRENDS OVER SIX YEARS... 
# NOT SURE IF THAT'S A COMPLETELY VALID REASON.
# Calculate de-trended values for sea ice concentration/thickness and vessel traffic 
# test <- fulldf[-which(fulldf$month %in% c(5:9)),]
# 
# test$date <- as.Date(paste0(test$year,"-",test$month,"-01"), format="%Y-%m-%d")
# test <- test[order(test$date),]
# 
# df <- data.frame(date= seq(as.Date("2015-01-01"),as.Date("2020-12-01"), by=c("month")), timestep=1:72)
# test <- left_join(test, df, by=c("date"))
# 
# temp <- EnvStats::kendallSeasonalTrendTest(traffic_km ~ month + year, 
#                                            data = subset(test, id == 2279))
# 
# library(zyp)
# testone <- test[which(test$id == 2279),]
# t <- zyp::zyp.trend.vector(testone$traffic_km, testone$timestep, method=c("zhang"))
# 
# preds <- (t["trend"]*testone$timestep + t["intercept"])
# predslinear <- (t["linear"]*testone$timestep + t["intercept"])
# 
# resids <- testone$traffic_km-preds
# 
# plot(testone$timestep, testone$traffic_km)
# lines(testone$timestep, preds, type="b", col="red")
# lines(testone$timestep, predslinear, type="b", col="blue")
# legend("topleft", legend=c("Zhang method", "Linear regression"),
#        col=c("red", "blue"), lty = 1:2, cex=0.8)

start <- proc.time()
fulldfnest <- fulldf %>% group_by(id) %>% nest()

pizzolato <- function(df){
  if(length(which(df$traffic_km != 0)) == 0){
    return(NA)
  }
  maxice <- max(df$icecon[which(df$traffic_km != 0)])
  lotsaice <- which(df$icecon > maxice & df$traffic_km == 0)
  keeplotsaice <- which(df$icecon[lotsaice] == min(df$icecon[lotsaice]))
  lotsaicenew <- lotsaice[-keeplotsaice]
  dfnew <- df[-lotsaicenew,]
  if(length(dfnew$traffic_km) > 25){
    cortest <- cor.test(dfnew$traffic_km, dfnew$icecon, method="kendall")
    return(cortest)
  }
  if(length(dfnew$traffic_km) <= 25){
    return(NA)
  }
}

mods <- fulldfnest %>% mutate(cortest = map(data, function(df){pizzolato(df)}), 
                              tidied = map(cortest, broom::tidy)) %>% 
  unnest(tidied)
proc.time() - start

modelresults <- allsf %>% dplyr::select(id) %>% left_join(., mods[,c("id", "estimate", "statistic", "p.value")])
sigcells <- st_coordinates(st_centroid(modelresults[which(modelresults$p.value < 0.05),])) %>% as.data.frame()

# studyoutline <- allsf %>% select(id) %>% st_union() %>% st_write("../Data_Processed/StudyOutline.shp")

p3 <- ggplot() +
  geom_sf(data=basemap.crop, fill="white", color="black", lwd=0.5, alpha = 0.9) +
  # geom_sf(data=aisbounds, fill=NA, color="black", lwd=1)+
  geom_sf(data=modelresults, aes(fill=estimate)) +
  colorspace::scale_fill_continuous_divergingx("RdBu", name="", rev=T) +
  geom_point(data=sigcells, aes(x = sigcells[,1], y = sigcells[,2]), shape=8, size=0.1) +
  xlab("") +
  ylab("") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  blank()+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.1,0.2),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        # panel.background = element_rect(fill = "lightblue"),
        panel.border =  element_rect(colour = "black"),
        panel.grid.major = element_line(colour = "transparent"))

# p3
ggsave(plot=p3, filename= "../Figures/KendallCorrelationMap_25obs.png", 
       width=8, height=8, units="in", dpi = 300)



############################################################################
# Ice and traffic line plot 
############################################################################

monthstats <- alldf %>% group_by(year, month) %>% summarize(iceext = sum(ifelse(icecon >= 10, 1, 0))*622.109950, traffic_km = sum(traffic_km))

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


coeff <- 3000 # median(monthstats$iceext/monthstats$traffic_1kkm)
library(RColorBrewer)
ggplot() +
  geom_line(data=monthstats, aes(x=timestep, y =iceext, group=year, col=as.factor(year)),lwd=1.5,lty=2) +
  geom_line(data=monthstats, aes(x=timestep, y = traffic_1kkm*coeff, group=year, col=as.factor(year)),lwd=1.5) +
  scale_color_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") + 
  xlab("") +
  scale_x_continuous(breaks=1:7, labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ) +
  scale_y_continuous(name = expression("Sea Ice Extent "~(km^2)),
                     sec.axis = sec_axis(~./coeff, name = "Total Shipping Traffic (1000s km)")) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave("../Figures/IceAndAIS_MonthlyLines.png", width=9, height=6, units="in")

############################################################################
# Traffic in ice plot 
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


library(RColorBrewer)
ggplot() +
  geom_line(data=monthinice, aes(x=timestep, y =traffic_km, group=year, col=as.factor(year)),lwd=1.5) +
  scale_color_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") + 
  xlab("") +
  scale_x_continuous(breaks=1:7, labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ) +
  scale_y_continuous(name = "Total Vessel Traffic\nin Ice (km)") +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("../Figures/AISInIce_MonthlyLines.png", width=9, height=6, units="in")

# Annual traffic totals in ice 
yearinice <- monthinice %>% group_by(year) %>% summarize(traffic_km = sum(traffic_km))

ggplot(yearinice, aes(x=year, y=traffic_km, fill=as.factor(year))) +
  geom_bar(stat="identity") +
  scale_fill_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") +
  scale_x_continuous(breaks= 2015:2020) +
  scale_y_continuous(breaks=seq(0,300000, by=50000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic in Ice (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggsave("../Figures/AISInIce_AnnualBar.png", width=9, height=6, units="in")

### Total vessel traffic in cells with >50% ice concentration 
yearinlotsaice <- alldf[which(alldf$icecon > 50 & alldf$traffic_km > 0),] %>% group_by(year) %>% summarize(traffic_km = sum(traffic_km))

ggplot(yearinlotsaice, aes(x=year, y=traffic_km, fill=as.factor(year))) +
  geom_bar(stat="identity") +
  scale_fill_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") +
  scale_x_continuous(breaks= 2015:2020) +
  scale_y_continuous(breaks=seq(0,50000, by=10000)) +
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("Total Vessel Traffic\nin >50% Ice (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggsave("../Figures/AISInIce_AnnualBar.png", width=9, height=6, units="in")

### Total vessel traffic in cells with 100% ice concentration is zero!  
yearinallice <- alldf[which(alldf$icecon >= 99  & alldf$traffic_km > 0),] %>% group_by(year) %>% summarize(traffic_km = sum(traffic_km))


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

# Box plot: vessel traffic in occupied cells x ice concentration 
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

# Proportion of cells occupied x ice concentration 
propocc <- alldf %>% 
  mutate(icecon = round(icecon, 0)) %>% 
  group_by(icecon) %>% 
  summarize(ncells=n(), nocc=length(which(nShips > 0))) %>% 
  mutate(propocc = round(nocc/ncells*100, 2))

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

# Total traffic by ice concentration 
test <- alldf %>% filter(icecon > 0) %>% mutate(icecon = round(icecon, 0)) %>% group_by(icecon) %>% summarize(traffic_km = mean(traffic_km), ncells=n())

ggplot(test, aes(x=icecon, y=traffic_km)) +
  geom_line() +
  # scale_fill_manual(labels=2015:2020, values=brewer.pal(7,"YlGnBu")[2:7], name="Year") +
  # scale_x_continuous(breaks= 2015:2020) +
  # scale_y_continuous(breaks=seq(0,50000, by=10000)) +
  theme_bw(base_size = 20) +
  xlab("Sea Ice Concentration (%)") +
  ylab("Mean Vessel Traffic\nin Occupied Cells (km)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")

############################




hightraffids <- unique(inice$id[which(inice$traffic_km > 1000)])
highshipids <- unique(inice$id[which(inice$nShips > 100)])
hightraff <- inice[which(inice$traffic_km > 1000),]

library(fitdistrplus)

############################################################################
# Exploring data distributions

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
############################################################################
# Pixel timeseries 
############################################################################
fulldf$date <- as.Date(paste0(fulldf$year, "-", fulldf$month, "-1"), format="%Y-%m-%d")

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

pixeltimeseries(fulldf, plotdata = F, savedata=T, saveloc=saveloc, 
                name="Traffic_km", var= "traffic_km", varlab = "Total Vessel Traffic (km)")

pixeltimeseries(fulldf, plotdata = F, savedata=T, saveloc=saveloc, 
                name="nShips", var= "nShips", varlab = "Number of Ships")


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


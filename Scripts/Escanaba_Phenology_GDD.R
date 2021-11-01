#Use GAMM to model walleye spawning phenology as a function of degree days following Pederson et al.


library(tidyverse)
library(climwin)
library(lubridate)
library(mgcv)
library(geosphere)
library(ggridges)

#Read in walleye spawning phenology data and summarize to daily catches
###Walleye data
dat <- read.csv("AllEscanabaWAE_1946_2019.csv")
dat
dat$Sex <- ifelse(dat$Sex %in% c(1, "Male","M"), "Male", ifelse(dat$Sex %in% c(2, "Female","F"), "Female",
                                                                ifelse(dat$Sex %in% c(3, "Unknown","U"), "Unknown", NA)))

dat <- mutate_at(dat, .vars="Date", .funs=as.Date, format="%m/%d/%Y")
dat <- as_tibble(dat)
dat <- uncount(dat, Count)
dat

#Fix lengths
dat$Length <- ifelse(is.na(dat$Length), dat$TLbin, dat$Length)
#Test test
#Summarize catch by date
wae <- 
  dat %>%
  filter(CorrGear==2, Month %in% c(3:5), Year != 1948) %>%
  group_by(Year, DOY, Date) %>%
  summarize(Total=n(), Females=sum(Sex=="Female", na.rm=T), PropFemale=Females/Total) %>%
  group_by(Year) %>%
  mutate(cumTotal = cumsum(Total), cumFemale=cumsum(Females), propcumTotal=cumTotal/sum(Total), propcumFemale=cumFemale/sum(Females),
         PropFemaleTotal = Females/sum(Total), scDOY = scale(DOY, center=T, scale=F), fempres = ifelse(Females > 0, 1, 0))
wae <- wae %>% group_by(Year) %>% mutate(scFem = scale(Females))


ggplot(data=wae[abs(wae$scDOY)<=15,], aes(x=scDOY, y=Females, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=filter(wae, abs(scDOY)<=15), aes(x=DOY, y=Year, group=Year, height=scFem+2,fill=stat(x))) + 
  geom_density_ridges_gradient(stat="identity", scale=6) + 
  scale_fill_viridis_c(name="DOY", option="D") + 
  scale_y_continuous(breaks=seq(1944,2040,4)) + theme_bw() + scale_x_continuous(breaks=seq(50,150,5))

ggplot(data=wae[abs(wae$scDOY)<=15,], aes(x=scDOY, y=scFem, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

#Environmental data
#Bring in environmental data
#Now Environmental data - start with precipitation
#Precipitation
precip <- read_csv("NOAA_NCDC_NorthernWIWeather_1940_2020.csv", 
                   col_types=cols(
                     STATION = col_character(),
                     NAME = col_character(),
                     LATITUDE = col_double(),
                     LONGITUDE = col_double(),
                     ELEVATION = col_double(),
                     DATE = col_date(format="%m/%d/%Y"),
                     AWND = col_double(),
                     MDWM = col_double(),
                     PRCP = col_double(),
                     SNOW = col_double(),
                     SNWD = col_double(),
                     TAVG = col_double(),
                     TMAX = col_double(),
                     TMIN = col_double(),
                     TSUN = col_double(),
                     WDMV = col_double()))
meanprecip <-
  precip %>%
  filter(NAME %in% c("EAGLE RIVER, WI US", "MINOCQUA, WI US", "PHELPS, WI US",
                     "RAINBOW RSVR LK TOMAHAWK, WI US", "REST LAKE, WI US")) %>%
  group_by(DATE) %>%
  summarize(meanPrecip=mean(PRCP, na.rm=T), meanSnow=mean(SNOW, na.rm=T), 
            meanMaxTemp=mean(TMAX, na.rm=T), meanMinTemp=mean(TMIN, na.rm=T), meanTemp = mean(c(meanMaxTemp, meanMinTemp)), N=sum(!is.na(PRCP)), meanWind=mean(WDMV, na.rm=T))
meanprecip

#Photoperiod

Photo <- data.frame(Date = seq(min(precip$DATE), as.Date("12/31/2020", format="%m/%d/%Y"), by=1), Year = NA, DOY = NA, DayLen = NA)
Photo$Year <- year(Photo$Date)
Photo$DOY <- as.numeric(strftime(Photo$Date, format="%j"))
Photo$DayLen <- daylength(Photo$Date, lat = 46.06413190)
Photo

#Water temperature (under ice modeled from Sparkling lake under ice buoy)
temps <- read_csv(file="ModeledEscanabaWaterTemps_1956.2020.csv",
                  col_types=cols(
                    Date=col_date(format="%m/%d/%Y"),
                    Notes = col_character()))
temps

#Ice data
#Ice on and off dates and water temps during them
ice <- read_csv("Escanaba_IceDates_1956_2020.csv",
                col_types=cols(
                  IceOn=col_date(format="%m/%d/%Y"),
                  IceOff=col_date(format="%m/%d/%Y")))
ice$IceOnDOY <- as.numeric(strftime(ice$IceOn, format="%j"))
ice$IceOffDOY <- as.numeric(strftime(ice$IceOff, format="%j"))
ice

#Function to calculate degree days from start date to end date by year
#x=temperature file with date and temperature
#y=fish dates
#start=starting date
start.month=1
start.day=1
gdd <- function(x, y, start.month, start.day) {  
  gdd <- matrix(NA, ncol=4, nrow=nrow(y))
  
  for (i in 1:nrow(y)) {
    if(start.month > 5){
    start = as.Date(paste0(start.month,"/",start.day,"/",y$Year[i]-1),format="%m/%d/%Y")
    } else {
      start = as.Date(paste0(start.month,"/",start.day,"/",y$Year[i]),format="%m/%d/%Y")
    }
    t = filter(temps, Date %in% seq(start,y$Date[i],1))
    gdd[i,1] = month(y$Date[i])
    gdd[i,2] = day(y$Date[i])
    gdd[i,3] = year(y$Date[i])
    gdd[i,4] = ifelse(length(t$PredWaterTemp)>0,sum(t$PredWaterTemp),NA)
  }
  colnames(gdd) <- c("Month","Day","Year","GDD")
  gdd <- as_tibble(gdd)
  gdd$Date <- with(gdd,as.Date(paste0(Month,"/",Day,"/",Year), format="%m/%d/%Y"))
  return(gdd)
}

Jan1GDD <- gdd(x=temps, y=sel.wae, start.month=1, start.day=1)
Oct1GDD <- gdd(x=temps, y=sel.wae, start.month=10, start.day=1)
names(Jan1GDD)[4] <- "Jan1GDD"
names(Oct1GDD)[4] <- "Oct1GDD"
Jan1GDD

Mar15GDD <- gdd(x=temps, y=sel.wae, start.month=3, start.day=15)

sel.wae <- left_join(sel.wae, select(Jan1GDD, Date, Jan1GDD)) %>% left_join(select(Oct1GDD, Date, Oct1GDD))
sel.wae <- left_join(sel.wae, select(Mar15GDD, Date, GDD))

ggplot(data=sel.wae, aes(x=Jan1GDD, y=scFem, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=sel.wae, aes(x=Oct1GDD, y=scFem, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=sel.wae, aes(x=GDD, y=scFem, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

#Write random effects gam to predict scaled female catch as function of DOY and year, with DOY nested in year?
#Pederson et al paper
hist(wae$scFem)
filter(sel.wae, Year==1985)

plot(PredWaterTemp ~ Date, data=temps[temps$Year==1985,])
lines(MeanWaterTemp ~ Date, data=temps[temps$Year==1985,])
lines(log(Females) ~ Date, data=sel.wae[sel.wae$Year==1985,], col="red", type="l")
lines(MeanAirTemp ~ Date, data=temps[temps$Year==1985,], col="blue")

##Model GI
wae$fYear <- as.factor(wae$Year)
sel.wae <- subset(wae, abs(scDOY) <= 15)

modGI.7.gaus <- gam(scFem ~ s(Jan1GDD, bs="tp", m=2) + s(Jan1GDD, by=fYear, k=7, bs="tp", m=2) + s(fYear, bs="re"), 
                    data=sel.wae, method="REML", family="gaussian")
summary(modGI.7.gaus)
plot(sel.wae$scFem[!is.na(sel.wae$Jan1GDD)] ~ predict(modGI.7.gaus))

plot(modGI.7.gaus)




pred.dat <- tibble(Jan1GDD = sel.wae$Jan1GDD[!(is.na(sel.wae$Jan1GDD))], fYear=sel.wae$fYear[!(is.na(sel.wae$Jan1GDD))])
pred.dat$pred <- predict(modGI.7.gaus, newdata=pred.dat)

ggplot(data=sel.wae, aes(x=Jan1GDD, y=scFem, group=Year)) + 
  geom_point() + 
  geom_line(data=pred.dat, aes(x=Jan1GDD, y=pred, group=fYear), inherit.aes=F) + 
  facet_wrap(~fYear, scales="free")

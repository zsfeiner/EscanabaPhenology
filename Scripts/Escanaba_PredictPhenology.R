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

#Write random effects gam to predict scaled female catch as function of DOY and year, with DOY nested in year?
#Pederson et al paper
hist(wae$scFem)

##Model GI
wae$fYear <- as.factor(wae$Year)
sel.wae <- subset(wae, abs(scDOY) <= 15)

modGI.7.gaus <- gam(scFem ~ s(scDOY, bs="tp", m=2) + s(scDOY, by=fYear, k=7, bs="tp", m=1) + s(fYear, bs="re"), 
               data=sel.wae, method="REML", family="gaussian")
summary(modGI.7.gaus)
plot(sel.wae$scFem ~ predict(modGI.7.gaus))

pred.dat.all <- data.frame(scDOY = rep(seq(-10,10,0.1),61), fYear=as.factor(rep(unique(wae$Year),each=201)))
pred.dat.all$pred <- predict(modGI.7.gaus, newdata=pred.dat.all)
pred.dat.all

unscale <- wae %>%
  group_by(fYear) %>%
  summarize(meanDOY = mean(DOY), sdDOY = sd(DOY), meanFem=mean(Females), sdFem=sd(Females))
unscale

#Unscale DOY and total female catches
all.pred <- left_join(pred.dat.all, unscale) %>%
  mutate(DOY = scDOY + meanDOY, PredFem = pred * sdFem + meanFem, PredFem=ifelse(PredFem<0,0,PredFem))
all.pred

ggplot(data=sel.wae, aes(x=DOY, y=scFem, group=fYear)) + 
  geom_point() + geom_line() +
  geom_vline(data=maxDOY, aes(xintercept=maxDOY), col="red", lwd=1) + 
  geom_line(data=unsc.pred, aes(x=DOY, y=pred, group=fYear), inherit.aes=F, col="blue", lwd=1) +
  facet_wrap(~fYear, scales="free") + theme_bw()


ggplot(data=sel.wae, aes(x=DOY, y=Females, group=fYear)) + 
  geom_point() + geom_line() +
  #geom_vline(data=maxDOY, aes(xintercept=maxDOY), col="red", lwd=1) + 
  geom_line(data=all.pred, aes(x=DOY, y=PredFem, group=fYear), inherit.aes=F, col="blue", lwd=1) +
  facet_wrap(~fYear, scales="free") + theme_bw()

maxDOY <- all.pred %>%
  group_by(fYear) %>%
  filter(PredFem == max(PredFem)) %>%
  mutate(Year = as.numeric(as.character(fYear))) %>%
  select(Year, DOY, PredFem, everything())
maxDOY

plot(DOY ~ Year, maxDOY)


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

#Get recruitment
YOYPE <- read_csv("EscanabaStockRecruitPE.csv")

#Get recruitment sampling dates

YOYdates <- read_csv("WAE_YOY_CompositeThrough2019.csv")
YOYdates <- YOYdates %>%
  filter(MWBC == 2339900) %>% mutate(Date = as.Date(DATE, format="%m/%d/%Y")) %>%
  group_by(YEAR) %>%
  summarize(Date=mean(Date))
YOYdates


dat <- maxDOY %>%
  ungroup(.) %>%
  mutate(FemDOY = DOY, rFemDOY = round(FemDOY)) %>%
  select(Year, FemDOY, rFemDOY, PredFem) %>%
  left_join(ice) %>%
  left_join(select(Photo, Year, DOY, DayLen), by=c("Year"="Year","rFemDOY"="DOY")) %>%
  left_join(select(temps, Year, DOY, PredWaterTemp), by=c("Year"="Year","rFemDOY"="DOY")) %>%
  left_join(select(YOYPE, Year, AdultPE, Age0PE))
dat

#Changes in spawn timing, water temp at spawn, and light at spawn
plot(FemDOY ~ Year, dat, type="b"); abline(reg=lm(FemDOY ~ Year, dat), col="blue")
plot(PredWaterTemp ~ Year, dat); abline(reg=lm(PredWaterTemp ~ Year, dat), col="blue")
plot(DayLen ~ Year, dat); abline(reg=lm(DayLen ~ Year, dat), col="blue")

#Correlation to ice-on, ice-off, ice-duration
plot(FemDOY ~ IceOnDOY, dat); abline(reg=lm(FemDOY ~ IceOnDOY, dat), col="blue")
plot(FemDOY ~ IceOffDOY, dat); abline(reg=lm(FemDOY ~ IceOffDOY, dat), col="blue")
plot(FemDOY ~ IceDuration, dat); abline(reg=lm(FemDOY ~ IceDuration, dat), col="blue")

#Recruitment vs spawn time, winter, and temp at spawn
plot(Age0PE ~ FemDOY, dat); abline(reg=lm(Age0PE ~ FemDOY, dat), col="blue") 
plot(Age0PE ~ DayLen, dat); abline(reg=lm(Age0PE ~ DayLen, dat), col="blue") 
plot(Age0PE ~ PredWaterTemp, dat); abline(reg=lm(Age0PE ~ PredWaterTemp, dat), col="blue") 
plot(Age0PE ~ IceOnDOY, dat); abline(reg=lm(Age0PE ~ IceOnDOY, dat), col="blue") 
plot(Age0PE ~ IceOffDOY, dat); abline(reg=lm(Age0PE ~ IceOffDOY, dat), col="blue") 
plot(Age0PE ~ IceDuration, dat); abline(reg=lm(Age0PE ~ IceDuration, dat), col="blue") 


#Climate windows analysis

#Create GDD0
temps$GDD0 <- ifelse(temps$PredWaterTemp > 0, temps$PredWaterTemp, 0)

#Find latest spawning date
max(dat$rFemDOY)
as.Date(max(dat$rFemDOY)-1, origin="2013-01-01")
as.Date(median(dat$rFemDOY)-1, origin="2013-01-01")

xvar <- inner_join(select(temps, GDD0, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))

xvar

dat$FemDate <- as.Date(dat$rFemDOY-1, origin=paste0(dat$Year,"-01-01"))
dat$FemDate

biol <- filter(dat, Year > 1955)

baseline <- lm(rFemDOY ~ 1, data=biol)
baseline

fem.win.GDD0 <- slidingwin(exclude = c(10,30),
                    xvar=list(GDD0=xvar$GDD0),# Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                    cdate=xvar$Date,
                    bdate=biol$FemDate,
                    baseline=baseline,
                    type="absolute",
                    refday=c(15,5),
                    stat=c("sum","CV","mean"),
                    func=c("lin"),
                    range=c(365, 0),
                    cinterval="day",
                    cmissing="method1",
                    k=5)
fem.win.GDD0$combos

fem.win.precip <- slidingwin(exclude = c(10,30),
                           xvar=list(Precip=xvar$meanPrecip),# DayLen=xvar$Photoperiod),
                           cdate=xvar$Date,
                           bdate=biol$FemDate,
                           baseline=baseline,
                           type="absolute",
                           refday=c(15,5),
                           stat=c("sum","mean"),
                           func=c("lin"),
                           range=c(365, 0),
                           cinterval="day",
                           cmissing="method1",
                           k=5)
fem.win.precip$combos

fem.win.daylen <- slidingwin(exclude = c(10,30),
                             xvar=list(DayLen=xvar$Photoperiod),
                             cdate=xvar$Date,
                             bdate=biol$FemDate,
                             baseline=baseline,
                             type="absolute",
                             refday=c(15,5),
                             stat=c("sum","mean"),
                             func=c("lin"),
                             range=c(365, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5)
fem.win.daylen$combos




#Randomize sum of GDD0
fem.win.randGDD0 <- randwin(exclude = c(10,30),
                      xvar=list(GDD0=xvar$GDD0),
                      cdate=xvar$Date,
                      bdate=biol$FemDate,
                      baseline=baseline,
                      type="absolute",
                      refday=c(15,5),
                      stat=c("sum"),
                      func=c("lin"),
                      range=c(365, 0),
                      cinterval="day",
                      cmissing="method1",
                      k=0, repeats=100)

#Randomize mean precip
fem.win.randPrecip <- randwin(exclude = c(10,30),
                            xvar=list(Precip=xvar$meanPrecip),
                            cdate=xvar$Date,
                            bdate=biol$FemDate,
                            baseline=baseline,
                            type="absolute",
                            refday=c(15,5),
                            stat=c("mean"),
                            func=c("lin"),
                            range=c(365, 0),
                            cinterval="day",
                            cmissing="method1",
                            k=0, repeats=100)


fem.win.precip$combos

fem.win.precip[[2]]
pvalue(dataset=fem.win.precip[[2]]$Dataset, datasetrand=fem.win.randPrecip[[1]],
  metric='AIC', sample.size=57)
plothist(dataset=fem.win.precip[[2]]$Dataset, datasetrand=fem.win.randPrecip[[1]])
plotdelta(dataset=fem.win.precip[[2]]$Dataset)
plotweights(dataset=fem.win.precip[[2]]$Dataset)
plotbetas(dataset=fem.win.precip[[2]]$Dataset)
plotwin(dataset=fem.win.precip[[2]]$Dataset)
plotbest(dataset=fem.win.precip[[2]]$Dataset,
         bestmodel=fem.win.precip[[2]]$BestModel,
         bestmodeldata=fem.win.precip[[2]]$BestModelData)


fem.win.GDD0$combos

fem.win.GDD0[[1]]
pvalue(dataset=fem.win.GDD0[[1]]$Dataset, datasetrand=fem.win.randGDD0[[1]],
       metric='AIC', sample.size=57)
plothist(dataset=fem.win.GDD0[[1]]$Dataset, datasetrand=fem.win.randGDD0[[1]])
plotdelta(dataset=fem.win.GDD0[[1]]$Dataset)
plotweights(dataset=fem.win.GDD0[[1]]$Dataset)
plotbetas(dataset=fem.win.GDD0[[1]]$Dataset)
plotwin(dataset=fem.win.GDD0[[1]]$Dataset)
plotbest(dataset=fem.win.GDD0[[1]]$Dataset,
         bestmodel=fem.win.GDD0[[1]]$BestModel,
         bestmodeldata=fem.win.GDD0[[1]]$BestModelData)



#Try on recruitment using relative window before and after spawning
xvar <- inner_join(select(temps, GDD0, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))
xvar

dat$FemDate <- as.Date(dat$rFemDOY-1, origin=paste0(dat$Year,"-01-01"))
dat$FemDate

r.biol <- filter(dat, Year > 1955, !is.na(Age0PE))
r.biol
r.baseline <- lm(Age0PE ~ 1, data=r.biol)
r.baseline

CV <- function(x) {sd(x)/mean(x)}

#Thermal habitat before spawning
rec.win.before.GDD0 <- slidingwin(exclude = c(10,1),
                      xvar=list(GDD0=xvar$GDD0),#, Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                      cdate=xvar$Date,
                      bdate=r.biol$FemDate,
                      baseline=r.baseline,
                      type="relative",
                      refday=c(15,5),
                      stat=c("sum","CV","mean"),#,"slope","mean"),
                      func=c("lin"),
                      range=c(300, 0),
                      cinterval="day",
                      cmissing="method1",
                      k=5)
rec.win.before.GDD0$combos

#Thermal habitat after spawning
rec.win.after.GDD0 <- slidingwin(exclude = c(10,-230),
                                 xvar=list(GDD0=xvar$GDD0),#, Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                                 cdate=xvar$Date,
                                 bdate=r.biol$FemDate,
                                 baseline=r.baseline,
                                 type="relative",
                                 refday=c(15,5),
                                 stat=c("sum", "CV","mean"),#,"slope","mean"),
                                 func=c("lin"),
                                 range=c(0, -230),
                                 cinterval="day",
                                 cmissing="method1",
                                 k=5)
rec.win.after.GDD0$combos

#Precipitation after spawning
rec.win.after.precip <- slidingwin(exclude = c(10,-230),
                                   xvar=list(Precip=xvar$meanPrecip),
                                   cdate=xvar$Date,
                                   bdate=r.biol$FemDate,
                                   baseline=r.baseline,
                                   type="relative",
                                   refday=c(15,5),
                                   stat=c("sum", "mean"),#,"slope","mean"),
                                   func=c("lin"),
                                   range=c(0, -230),
                                   cinterval="day",
                                   cmissing="method1",
                                   k=5)

#Precipitation before spawning
rec.win.before.precip <- slidingwin(exclude = c(10,1),
                                    xvar=list(Precip=xvar$meanPrecip),
                                    cdate=xvar$Date,
                                    bdate=r.biol$FemDate,
                                    baseline=r.baseline,
                                    type="relative",
                                    refday=c(15,5),
                                    stat=c("sum", "mean"),#,"slope","mean"),
                                    func=c("lin"),
                                    range=c(300, 0),
                                    cinterval="day",
                                    cmissing="method1",
                                    k=5)

rec.win.before.precip$combos
rec.win.after.precip$combos



rec.win.before.GDD0_rand <- randwin(exclude = c(10,1),
                                      xvar=list(GDD0=xvar$GDD0),
                                      cdate=xvar$Date,
                                      bdate=r.biol$FemDate,
                                      baseline=r.baseline,
                                      type="relative",
                                      refday=c(15,5),
                                      stat=c("CV", 'sum'),
                                      func=c("lin"),
                                      range=c(300, 0),
                                      cinterval="day",
                                      cmissing="method1",
                                      k=0, repeats=100)

rec.win.after.GDD0_rand <- randwin(exclude = c(10,-230),
                                     xvar=list(GDD0=xvar$GDD0),
                                     cdate=xvar$Date,
                                     bdate=r.biol$FemDate,
                                     baseline=r.baseline,
                                     type="relative",
                                     refday=c(15,5),
                                     stat=c("CV", "sum"),
                                     func=c("lin"),
                                     range=c(0, -230),
                                     cinterval="day",
                                     cmissing="method1",
                                     k=0, repeats=100)

rec.win.before.GDD0$combos
rec.win.before.GDD0_rand$combos

###GDD Before spawning
#Check GDD CV
pvalue(dataset=rec.win.before.GDD0[[2]]$Dataset, datasetrand=rec.win.before.GDD0_rand[[1]],
       metric='AIC', sample.size=nrow(r.biol))

#Check GDD sum
pvalue(dataset=rec.win.before.GDD0[[1]]$Dataset, datasetrand=rec.win.before.GDD0_rand[[2]],
       metric='AIC', sample.size=nrow(r.biol))
plothist(dataset=rec.win.before.GDD0[[1]]$Dataset, datasetrand=rec.win.before.GDD0_rand[[2]])
plotdelta(dataset=rec.win.before.GDD0[[1]]$Dataset)
plotweights(dataset=rec.win.before.GDD0[[1]]$Dataset)
plotbetas(dataset=rec.win.before.GDD0[[1]]$Dataset)
plotwin(dataset=rec.win.before.GDD0[[1]]$Dataset)
plotbest(dataset=rec.win.before.GDD0[[1]]$Dataset,
         bestmodel=rec.win.before.GDD0[[1]]$BestModel,
         bestmodeldata=rec.win.before.GDD0[[1]]$BestModelData)


###GDD after spawning
rec.win.after.GDD0$combos
rec.win.after.GDD0_rand$combos

#Check GDD CV
pvalue(dataset=rec.win.after.GDD0[[2]]$Dataset, datasetrand=rec.win.after.GDD0_rand[[1]],
       metric='AIC', sample.size=nrow(r.biol))

#Check GDD sum
pvalue(dataset=rec.win.after.GDD0[[1]]$Dataset, datasetrand=rec.win.after.GDD0_rand[[2]],
       metric='AIC', sample.size=nrow(r.biol))

plothist(dataset=rec.win.after.GDD0[[1]]$Dataset, datasetrand=rec.win.after.GDD0_rand[[2]])
plotdelta(dataset=rec.win.after.GDD0[[1]]$Dataset)
plotweights(dataset=rec.win.after.GDD0[[1]]$Dataset)
plotbetas(dataset=rec.win.after.GDD0[[1]]$Dataset)
plotwin(dataset=rec.win.after.GDD0[[1]]$Dataset)
plotbest(dataset=rec.win.after.GDD0[[1]]$Dataset,
         bestmodel=rec.win.after.GDD0[[1]]$BestModel,
         bestmodeldata=rec.win.after.GDD0[[1]]$BestModelData)





rec.win.after.precip <- slidingwin(exclude = c(10,-230),
                                 xvar=list(Precip=xvar$meanPrecip),
                                 cdate=xvar$Date,
                                 bdate=r.biol$FemDate,
                                 baseline=r.baseline,
                                 type="relative",
                                 refday=c(15,5),
                                 stat=c("sum", "mean"),#,"slope","mean"),
                                 func=c("lin"),
                                 range=c(0, -230),
                                 cinterval="day",
                                 cmissing="method1",
                                 k=5)

rec.win.before.precip <- slidingwin(exclude = c(10,1),
                                   xvar=list(Precip=xvar$meanPrecip),
                                   cdate=xvar$Date,
                                   bdate=r.biol$FemDate,
                                   baseline=r.baseline,
                                   type="relative",
                                   refday=c(15,5),
                                   stat=c("sum", "mean"),#,"slope","mean"),
                                   func=c("lin"),
                                   range=c(300, 0),
                                   cinterval="day",
                                   cmissing="method1",
                                   k=5)

rec.win.before.precip$combos
rec.win.after.precip$combos
rec.win.before.precip_rand <- randwin(exclude = c(10,1),
                                      xvar=list(Precip=xvar$meanPrecip),
                                      cdate=xvar$Date,
                                      bdate=r.biol$FemDate,
                                      baseline=r.baseline,
                                      type="relative",
                                      refday=c(15,5),
                                      stat=c("sum","mean"),
                                      func=c("lin"),
                                      range=c(300, 0),
                                      cinterval="day",
                                      cmissing="method1",
                                      k=0, repeats=100)

rec.win.after.precip_rand <- randwin(exclude = c(10,-230),
                                     xvar=list(Precip=xvar$meanPrecip),
                                     cdate=xvar$Date,
                                     bdate=r.biol$FemDate,
                                     baseline=r.baseline,
                                     type="relative",
                                     refday=c(15,5),
                                     stat=c("sum","mean"),
                                     func=c("lin"),
                                     range=c(0, -230),
                                     cinterval="day",
                                     cmissing="method1",
                                     k=0, repeats=100)
#Precipitation before spawning
rec.win.before.precip$combos
rec.win.before.precip_rand$combos

pvalue(dataset=rec.win.before.precip[[1]]$Dataset, datasetrand=rec.win.before.precip_rand[[1]],
       metric='AIC', sample.size=55)
pvalue(dataset=rec.win.before.precip[[2]]$Dataset, datasetrand=rec.win.before.precip_rand[[2]],
       metric='AIC', sample.size=55)

plothist(dataset=rec.win.before.precip[[1]]$Dataset, datasetrand=rec.win.before.precip_rand[[1]])
plotdelta(dataset=rec.win.before.precip[[1]]$Dataset)
plotweights(dataset=rec.win.before.precip[[1]]$Dataset)
plotbetas(dataset=rec.win.before.precip[[1]]$Dataset)
plotwin(dataset=rec.win.before.precip[[1]]$Dataset)
plotbest(dataset=rec.win.before.precip[[1]]$Dataset,
         bestmodel=rec.win.before.precip[[1]]$BestModel,
         bestmodeldata=rec.win.before.precip[[1]]$BestModelData)

#Precipitation after spawning
rec.win.after.precip$combos
rec.win.after.precip_rand$combos

pvalue(dataset=rec.win.after.precip[[1]]$Dataset, datasetrand=rec.win.after.precip_rand[[1]],
       metric='AIC', sample.size=55)
pvalue(dataset=rec.win.after.precip[[2]]$Dataset, datasetrand=rec.win.after.precip_rand[[2]],
       metric='AIC', sample.size=55)

plothist(dataset=rec.win.after.precip[[1]]$Dataset, datasetrand=rec.win.after.precip_rand[[1]])
plotdelta(dataset=rec.win.after.precip[[1]]$Dataset)
plotweights(dataset=rec.win.after.precip[[1]]$Dataset)
plotbetas(dataset=rec.win.after.precip[[1]]$Dataset)
plotwin(dataset=rec.win.after.precip[[1]]$Dataset)
plotbest(dataset=rec.win.after.precip[[1]]$Dataset,
         bestmodel=rec.win.after.precip[[1]]$BestModel,
         bestmodeldata=rec.win.after.precip[[1]]$BestModelData)
# 
# #Split dataset into thirds?
# #Try on recruitment using relative window before and after spawning
# #Nothing
# #Do last third of recruitment years
# nrow(r.biol)/3
# 
# plot(Age0PE ~ Year, data=r.biol)
# lines(Age0PE ~ Year, data=r.biol[1:18,], col=1)
# lines(Age0PE ~ Year, data=r.biol[19:37,], col=2)
# lines(Age0PE ~ Year, data=r.biol[38:55,], col=3)
# lines(x=c(r.biol$Year[c(1,18)]), y=rep(mean(r.biol$Age0PE[1:18]),2), col=1, lty=2)
# lines(x=c(r.biol$Year[c(19,37)]), y=rep(mean(r.biol$Age0PE[19:37]),2), col=2, lty=2)
# lines(x=c(r.biol$Year[c(38,55)]), y=rep(mean(r.biol$Age0PE[38:55]),2), col=3, lty=2)
# 
# r.biol.late <- r.biol[38:55, ]
# r.biol.late
# r.late.baseline <- lm(Age0PE ~ 1, data=r.biol.late)
# r.late.baseline
# 
# #Thermal habitat before spawning
# rec.late.before <- slidingwin(exclude = c(10,1),
#                                   xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                                   cdate=xvar$Date,
#                                   bdate=r.biol.late$FemDate,
#                                   baseline=r.late.baseline,
#                                   type="relative",
#                                   stat=c("sum","CV"),
#                                   func=c("lin"),
#                                   range=c(300, 0),
#                                   cinterval="day",
#                                   cmissing="method1",
#                                   k=5)
# 
# rec.late.after <- slidingwin(exclude = c(10,-230),
#                               xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                               cdate=xvar$Date,
#                               bdate=r.biol.late$FemDate,
#                               baseline=r.late.baseline,
#                               type="relative",
#                               stat=c("sum","CV"),
#                               func=c("lin"),
#                               range=c(0, -230),
#                               cinterval="day",
#                               cmissing="method1",
#                               k=5)
# 
# rec.late.before$combos
# rec.late.after$combos
# 
# plotwin(dataset=rec.late.before[[2]]$Dataset)
# plotwin(dataset=rec.late.after[[3]]$Dataset)
# 
# ###Do middle third of recruitment years
# #Do last third of recruitment years
# nrow(r.biol)/3
# 
# plot(Age0PE ~ Year, data=r.biol)
# lines(Age0PE ~ Year, data=r.biol[1:18,], col=1)
# lines(Age0PE ~ Year, data=r.biol[19:37,], col=2)
# lines(Age0PE ~ Year, data=r.biol[38:55,], col=3)
# lines(x=c(r.biol$Year[c(1,18)]), y=rep(mean(r.biol$Age0PE[1:18]),2), col=1, lty=2)
# lines(x=c(r.biol$Year[c(19,37)]), y=rep(mean(r.biol$Age0PE[19:37]),2), col=2, lty=2)
# lines(x=c(r.biol$Year[c(38,55)]), y=rep(mean(r.biol$Age0PE[38:55]),2), col=3, lty=2)
# 
# r.biol.mid <- r.biol[19:37, ]
# r.biol.mid
# r.mid.baseline <- lm(Age0PE ~ 1, data=r.biol.mid)
# r.mid.baseline
# 
# #Thermal habitat before spawning
# rec.mid.before <- slidingwin(exclude = c(10,1),
#                               xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                               cdate=xvar$Date,
#                               bdate=r.biol.mid$FemDate,
#                               baseline=r.mid.baseline,
#                               type="relative",
#                               stat=c("sum","CV"),
#                               func=c("lin"),
#                               range=c(300, 0),
#                               cinterval="day",
#                               cmissing="method1",
#                               k=5)
# 
# rec.mid.after <- slidingwin(exclude = c(10,-230),
#                              xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                              cdate=xvar$Date,
#                              bdate=r.biol.mid$FemDate,
#                              baseline=r.mid.baseline,
#                              type="relative",
#                              stat=c("sum","CV"),
#                              func=c("lin"),
#                              range=c(0, -230),
#                              cinterval="day",
#                              cmissing="method1",
#                              k=5)
# 
# rec.mid.before$combos
# rec.mid.after$combos
# 
# plotwin(dataset=rec.mid.before[[2]]$Dataset)
# plotwin(dataset=rec.mid.after[[3]]$Dataset)
# 
# #Quickly look at absolute day, May 1
# #Thermal habitat before spawning
# ##Nothing
# rec.abs.before <- slidingwin(exclude = c(10,1),
#                               xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                               cdate=xvar$Date,
#                               bdate=r.biol$FemDate,
#                               baseline=r.baseline,
#                               type="absolute",
#                               refday=c(1,5),
#                               stat=c("sum","CV"),
#                               func=c("lin"),
#                               range=c(300, 0),
#                               cinterval="day",
#                               cmissing="method1",
#                               k=5)
# 
# rec.abs.after <- slidingwin(exclude = c(10,-230),
#                              xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
#                              cdate=xvar$Date,
#                              bdate=r.biol$FemDate,
#                              baseline=r.baseline,
#                              type="absolute",
#                              refday=c(1,5),
#                              stat=c("sum","CV"),
#                              func=c("lin"),
#                              range=c(0, -230),
#                              cinterval="day",
#                              cmissing="method1",
#                              k=5)
# 
# rec.abs.before$combos
# rec.abs.after$combos

r.biol <- left_join(r.biol, YOYdates, by=c("Year"="YEAR"))
r.biol
r.biol2 <- filter(r.biol, !is.na(Date))
r.baseline2 <- lm(Age0PE ~ 1, r.biol2)

rec.abs.beforefall2 <- slidingwin(exclude = c(10,1),
                             xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                             cdate=xvar$Date,
                             bdate=r.biol2$Date,
                             baseline=r.baseline2,
                             type="absolute",
                             refday=c(15,10),
                             stat=c("mean"),
                             func=c("lin"),
                             range=c(230, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=0)
rec.abs.beforefall$combos

plotdelta(dataset=rec.abs.beforefall[[3]]$Dataset)
plotweights(dataset=rec.abs.beforefall[[3]]$Dataset)
plotbetas(dataset=rec.abs.beforefall[[3]]$Dataset)
plotwin(dataset=rec.abs.beforefall[[3]]$Dataset)
plotbest(dataset=rec.abs.beforefall[[3]]$Dataset,
         bestmodel=rec.abs.beforefall[[3]]$BestModel,
         bestmodeldata=rec.abs.beforefall[[3]]$BestModelData)

####Log scale recruitment
#####
lnr.baseline <- lm(log(Age0PE) ~ 1, data=r.biol)
lnr.baseline

CV <- function(x) {sd(x)/mean(x)}

#Thermal habitat before spawning
rec.win.before.GDD0.lnr <- slidingwin(exclude = c(10,1),
                                  xvar=list(GDD0=xvar$GDD0),#, Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                                  cdate=xvar$Date,
                                  bdate=r.biol$FemDate,
                                  baseline=lnr.baseline,
                                  type="relative",
                                  refday=c(15,5),
                                  stat=c("sum","CV","mean"),#,"slope","mean"),
                                  func=c("lin"),
                                  range=c(300, 0),
                                  cinterval="day",
                                  cmissing="method1",
                                  k=0)
rec.win.before.GDD0$combos
rec.win.before.GDD0.lnr$combos


#Thermal habitat after spawning
rec.win.after.GDD0.lnr <- slidingwin(exclude = c(10,-230),
                                 xvar=list(GDD0=xvar$GDD0),#, Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                                 cdate=xvar$Date,
                                 bdate=r.biol$FemDate,
                                 baseline=lnr.baseline,
                                 type="relative",
                                 refday=c(15,5),
                                 stat=c("sum", "CV","mean"),#,"slope","mean"),
                                 func=c("lin"),
                                 range=c(0, -230),
                                 cinterval="day",
                                 cmissing="method1",
                                 k=0)
rec.win.after.GDD0.lnr$combos
rec.win.after.GDD0$combos


plotdelta(dataset=rec.win.after.GDD0.lnr[[2]]$Dataset)
plotweights(dataset=rec.win.after.GDD0.lnr[[2]]$Dataset)
plotbetas(dataset=rec.win.after.GDD0.lnr[[2]]$Dataset)
plotwin(dataset=rec.win.after.GDD0.lnr[[2]]$Dataset)
plotbest(dataset=rec.win.after.GDD0[[2]]$Dataset,
         bestmodel=rec.win.after.GDD0[[2]]$BestModel,
         bestmodeldata=rec.win.after.GDD0[[2]]$BestModelData)


##Quick look at lengths and relationship with phenology
yoylen <- read_csv("C:/Users/feinezs/Documents/WAE Spawning Phenology/YOYWAEData_Escanaba/Escanaba_FallShockingLengths_1990-2019.csv")
yoylen$Date <- as.Date(yoylen$Date, format="%m/%d/%Y")
yoylen <- filter(yoylen, Age==0, Species=="WAE")
plot(Length.mm ~ Year, data=yoylen)

yoylen <- yoylen %>%
  left_join(select(dat, Year, FemDOY, IceOnDOY, IceOffDOY, FemDate))
plot(Length.mm ~ FemDOY, data=yoylen)
yoylen

yoylen <- yoylen %>%
  mutate(yoyDOY = as.numeric(strftime(Date, format="%j"))) %>%
  mutate(growDays = yoyDOY - FemDOY) %>%
  mutate(growthrate = Length.mm/growDays) %>% rowwise() %>%
  mutate(GDD = ifelse(Year<2020, sum(temps$GDD0[temps$Date %in% seq(FemDate, Date, 1)]), NA))

yoylen$GDD

temps$GDD0[temps$Date %in% (seq(yoylen$FemDate[1],yoylen$Date[1], 1))]

plot(Length.mm ~ growDays, data=yoylen)
plot(Length.mm ~ IceOffDOY, data=yoylen)
plot(growthrate ~ FemDOY, data=yoylen)
plot(growthrate ~ IceOffDOY, data=yoylen)
plot(growthrate ~ GDD, data=yoylen)
plot(Length.mm ~ GDD, data=yoylen)
plot(growthrate ~ Year, data=yoylen)

plot(Length.mm ~ yoyDOY, data=yoylen)

yoylen
yoylen$yoyDOY

meanyoylen <- yoylen %>%
  group_by(Date) %>%
  summarize(Length = mean(Length.mm, na.rm=T), FemDOY = mean(FemDOY), FemDate=mean(FemDate), growDays=mean(growDays), growthrate=mean(growthrate), GDD=mean(GDD))
meanyoylen


plot(Length ~ growDays, data=meanyoylen)
plot(Length ~ IceOffDOY, data=meanyoylen)
plot(Length ~ year(Date), data=meanyoylen)
plot(growthrate ~ FemDOY, data=meanyoylen)
plot(growthrate ~ IceOffDOY, data=meanyoylen)
plot(growthrate ~ GDD, data=meanyoylen)
plot(Length ~ GDD, data=meanyoylen)
plot(growthrate ~ Year, data=meanyoylen)

meanyoylen <- meanyoylen %>%
  mutate(Year=year(Date)) %>%
  left_join(select(dat, AdultPE, Age0PE, Year))
meanyoylen

####Climwin for YOY length post spawn relative window

#Try on recruitment using relative window before and after spawning
xvar <- inner_join(select(temps, GDD0, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))
xvar

meanyoylen <- filter(meanyoylen, year(Date) < 2020)
len.baseline <- lm(Length ~ 1, data=meanyoylen)

#Thermal habitat before spawning
yoylen.win <- slidingwin(exclude = c(10,1),
                         xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                         cdate=xvar$Date,
                         bdate=meanyoylen$Date,
                         baseline=len.baseline,
                         type="relative",
                         #refday=c(15,5),
                         stat=c("sum", "CV", "mean"),
                         func=c("lin"),
                         range=c(180, 0),
                         cinterval="day",
                         cmissing="method1",
                         k=5)
yoylen.win$combos

plotdelta(dataset=yoylen.win[[1]]$Dataset)
plotweights(dataset=yoylen.win[[1]]$Dataset)
plotbetas(dataset=yoylen.win[[1]]$Dataset)
plotwin(dataset=yoylen.win[[1]]$Dataset)
plotbest(dataset=yoylen.win[[1]]$Dataset,
         bestmodel=yoylen.win[[1]]$BestModel,
         bestmodeldata=yoylen.win[[1]]$BestModelData)

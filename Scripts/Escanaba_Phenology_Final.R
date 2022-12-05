library(tidyverse)
library(climwin)
library(lubridate)
library(mgcv)
library(geosphere)
library(ggridges)

#Read in walleye spawning phenology data and summarize to daily catches
###Walleye data
dat <- read.csv("AllEscanabaWAE_1946_2021.csv")
dat
dat$Sex <- ifelse(dat$Sex %in% c(1, "Male","M"), "Male", ifelse(dat$Sex %in% c(2, "Female","F"), "Female",
                                                                ifelse(dat$Sex %in% c(3, "Unknown","U"), "Unknown", NA)))

dat <- mutate_at(dat, .vars="Date", .funs=as.Date, format="%m/%d/%Y")
dat <- as_tibble(dat)
dat <- uncount(dat, Count)
dat

#Fix lengths
dat$Length <- ifelse(is.na(dat$Length), dat$TLbin, dat$Length)


#Summarize catch by date
wae <- 
  dat %>%
  filter(CorrGear==2, Month %in% c(3:5), Year != 1948) %>%
  group_by(Year, DOY, Date) %>%
  summarize(Total=n(), Females=sum(Sex=="Female", na.rm=T), PropFemale=Females/Total)

ggplot(data=wae, aes(x=DOY, y=Females, group=Year)) + 
  geom_point() + 
  facet_wrap(~Year, scales="free")

ggplot(data=wae, aes(x=DOY, y=Year, group=Year, height=Females,fill=stat(x))) + 
  geom_density_ridges_gradient(stat="identity", scale=6) + 
  scale_fill_viridis_c(name="DOY", option="D") + 
  scale_y_continuous(breaks=seq(1944,2040,4)) + theme_bw() + scale_x_continuous(breaks=seq(50,150,5))


#Write random effects gam to predict female catch as function of scaled DOY and year
#Pederson et al paper
hist(wae$Females)

#Filter out observations more than 15 days from the center and after 1955 (when ice records start)
wae$fYear <- as.factor(wae$Year)
sel.wae <- wae %>%
  filter(Year > 1955) %>%
  group_by(Year) %>%
  filter(abs(mean(DOY) - DOY) < 15) %>%
  mutate(scDOY = scale(DOY, center=T, scale=F))
sel.wae

hist(sel.wae$Females)
hist(sel.wae$DOY)
summary(sel.wae$DOY)

#Fit poisson models to check model fit on raw DOY
# modG_pois <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))
# modGS_pois <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="cc"),m=2) +
#                     s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))
# modGI_pois <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, by=fYear, k=6,  bs="cc",m=2) +
#                     s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))
# 
# modG_pois$aic
# modGS_pois$aic
# modGI_pois$aic


#Fit quasi-poisson models, better fit
modG <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(-15,15)))

modGS <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="cc"),m=2) +
               s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

modGI <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, by=fYear, k=6,  bs="cc",m=2) +
                    s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

summary(modG)
summary(modGS)
summary(modGI)

library(gratia)

draw(modGS)
draw(modG)
draw(modGI)

par(mfrow=c(2,2))
gam.check(modG)
gam.check(modGS)
gam.check(modGI)

AIC(modG)


lim.pred.data <- sel.wae %>%
  group_by(fYear) %>%
  summarize(first=min(DOY)-20, last=max(DOY)+20) %>%
  group_by(fYear) %>%
  group_modify(~ tibble(DOY = seq(.$first, .$last))) %>%
  ungroup()




#pred.data <- tibble(DOY=rep(c(80:150), 61), fYear=as.factor(rep(unique(sel.wae$fYear), each=length(80:150))))
#pred.data
modG.pred <- as_tibble(predict(modG, newdata=lim.pred.data, se.fit=T, type="response"))
modG.pred

modGS.pred <- as_tibble(predict(modGS, newdata=lim.pred.data, type="response", se.fit=T))
modGS.pred

modGI.pred <- as_tibble(predict(modGI, newdata=lim.pred.data, type="response", se.fit=T))
modGI.pred

lim.pred.data$modGfit <- modG.pred$fit
lim.pred.data$modGse <- modG.pred$se.fit
lim.pred.data$modGSfit <- modGS.pred$fit
lim.pred.data$modGSse <- modGS.pred$se.fit
lim.pred.data$modGIfit <- modGI.pred$fit
lim.pred.data$modGIse <- modGI.pred$se.fit


ggplot(data=lim.pred.data, aes(x=DOY, y=modGfit, group=fYear)) + 
  geom_ribbon(aes(x=DOY, ymin=modGfit-modGse, ymax=modGfit+modGse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")


ggplot(data=lim.pred.data, aes(x=DOY, y=modGSfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGSfit-modGSse, ymax=modGSfit+modGSse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")

ggplot(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGIfit-modGIse, ymax=modGIfit+modGIse, group=fYear), alpha=0.4) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + #ylim(c(-2,3)) +
  theme_classic() + facet_wrap(~fYear, scales="free")

#plot all predictions
ggplot(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + 
  geom_line(aes(x=DOY, y=modGfit, group=fYear), color="red") + 
  geom_line(aes(x=DOY, y=modGSfit, group=fYear), color="green") +
  theme_classic() + facet_wrap(~fYear, scales="free")


#Quasi-poisson model GI is the best fit and provides reasonable predictions
#Extract predictions of day with highest predicted number of spawning females
#Generate 
lim.pred.data

esc.phen <- lim.pred.data %>%
  group_by(fYear) %>%
  summarize(SpawnDOY = DOY[which.max(modGIfit)])

library(HDInterval)

hdi.year <- tibble()
for (i in unique(lim.pred.data$fYear)) {
  set <- filter(lim.pred.data, fYear==i)
  reps <- lapply(set, rep, set$modGIfit)
  hdi.year <- bind_rows(hdi.year, HDInterval::hdi(reps$DOY,credMass=0.80))
}
hdi.year

esc.phen <- bind_cols(esc.phen, hdi.year) %>%
  mutate(Duration = upper-lower)

esc.phen

plot(Duration ~ fYear, data=esc.phen)
plot(Duration ~ SpawnDOY, data=esc.phen)

ggplot(sel.wae, aes(x=DOY, y=Females, group=fYear)) + 
  geom_point() + 
  geom_vline(data=esc.phen, aes(xintercept=lower, group=fYear)) + 
  geom_vline(data=esc.phen, aes(xintercept=upper, group=fYear)) + 
  facet_wrap(~fYear, scales='free') + 
  geom_line(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear))
par(mfrow=c(1,1))
plot(SpawnDOY ~ as.numeric(fYear), data=esc.phen, type="b")

########Compare predicted water temps against USGS model - does pretty well, ice season is more iffy#######
########Test water temps from my model against USGS lake model

#Predicted water temperatures for Escanaba

predtemps <- read_csv("ModeledEscanabaWaterTemps_1956.2020.csv")
predtemps$Date <- as.Date(predtemps$Date, format="%m/%d/%Y")

usgstemps <- read_csv("Escanaba_USGS_Temperatures.csv")


combtemps <- left_join(usgstemps, predtemps, by=c("date"="Date"))
combtemps

select(combtemps, date, temp_0, PredWaterTemp)
cor(combtemps$temp_0, combtemps$PredWaterTemp)


plot(PredWaterTemp ~ temp_0, data=filter(combtemps, month(date) %in% c(12,1,2,3)))



######
#######Use climate windows package to predict model-predicted female catches based on environmental variables
######

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



#Climate windows analysis

#Create GDD0
temps$GDD0 <- ifelse(temps$PredWaterTemp > 0, temps$PredWaterTemp, 0)
temps$GDD5 <- ifelse(temps$PredWaterTemp >= 5, temps$PredWaterTemp, 0)


##Model predicted female walleye abundance data, scaled and centered within years
pred.wae <- select(lim.pred.data, fYear, DOY, modGIfit) %>%
  mutate(Date = as.Date(DOY-1, origin=paste0(fYear, "-01-01"))) %>%
  group_by(fYear) %>%
  mutate(PredWaeCount = round(modGIfit, 0), scPredFem = scale(modGIfit, center=F, scale=T), propFem = PredWaeCount/sum(PredWaeCount)) %>%
  filter(fYear != 2021)

ggplot(pred.wae, aes(x=DOY, y=scPredFem, group=fYear)) + 
  geom_line()

#Set up environmental variables
xvar <- inner_join(select(temps, GDD0, GDD5, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))
xvar

baseline <- glm(scPredFem ~ 1, data=pred.wae, family="Gamma")
baseline

#Use mean, CV of water temperature, sum of GDD5, and sum of precip
fem.win.temp <- slidingwin(exclude = c(10,30),
                           xvar=list(Temp=xvar$PredWaterTemp),
                           cdate=xvar$Date,
                           bdate=pred.wae$Date,
                           baseline=baseline,
                           type="relative",
                           stat=c("mean"),
                           func=c("quad"),
                           range=c(270, 0),
                           cinterval="day",
                           cmissing="method1",
                           k=5)

##Assess temp
fem.win.temp$combos
fem.win.temp[[1]]

#pvalue(dataset=fem.win.temp[[1]]$Dataset, datasetrand=fem.win.rand.temp[[1]],
#       metric='C', sample.size=57)
#plothist(dataset=fem.win.temp[[1]]$Dataset, datasetrand=fem.win.rand.temp[[1]])
plotdelta(dataset=fem.win.temp[[1]]$Dataset)
plotweights(dataset=fem.win.temp[[1]]$Dataset)
plotbetas(dataset=fem.win.temp[[1]]$Dataset)
plotwin(dataset=fem.win.temp[[1]]$Dataset)
plotbest(dataset=fem.win.temp[[1]]$Dataset,
         bestmodel=fem.win.temp[[1]]$BestModel,
         bestmodeldata=fem.win.temp[[1]]$BestModelData)


fem.win.GDD5 <- slidingwin(exclude = c(10,30),
                           xvar=list(Temp=xvar$GDD5),# Precip=xvar$meanPrecip, DayLen=xvar$Photoperiod),
                           cdate=xvar$Date,
                           bdate=pred.wae$Date,
                           baseline=baseline,
                           type="relative",
                           stat=c("sum"),
                           func=c("quad"),
                           range=c(270, 0),
                           cinterval="day",
                           cmissing="method1",
                           k=5)
##Assess GDD
fem.win.GDD5$combos
fem.win.GDD5[[1]]
#pvalue(dataset=fem.win.GDD5[[1]]$Dataset, datasetrand=fem.win.rand.GDD5[[1]],
#       metric='C', sample.size=57)
#plothist(dataset=fem.win.GDD5[[1]]$Dataset, datasetrand=fem.win.rand.GDD5[[1]])
plotdelta(dataset=fem.win.GDD5[[1]]$Dataset)
plotweights(dataset=fem.win.GDD5[[1]]$Dataset)
plotbetas(dataset=fem.win.GDD5[[1]]$Dataset)
plotwin(dataset=fem.win.GDD5[[1]]$Dataset)
plotbest(dataset=fem.win.GDD5[[1]]$Dataset,
         bestmodel=fem.win.GDD5[[1]]$BestModel,
         bestmodeldata=fem.win.GDD5[[1]]$BestModelData)


fem.win.precip <- slidingwin(exclude = c(10,30),
                             xvar=list(Precip=xvar$meanPrecip),
                             cdate=xvar$Date,
                             bdate=pred.wae$Date,
                             baseline=baseline,
                             type="relative",
                             stat=c("sum"),
                             func=c("lin"),
                             range=c(270, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5)

##Assess GDD
fem.win.precip$combos
fem.win.precip[[1]]
#pvalue(dataset=fem.win.precip[[1]]$Dataset, datasetrand=fem.win.rand.precip[[1]],
#       metric='C', sample.size=57)
#plothist(dataset=fem.win.precip[[1]]$Dataset, datasetrand=fem.win.rand.precip[[1]])
plotdelta(dataset=fem.win.precip[[1]]$Dataset)
plotweights(dataset=fem.win.precip[[1]]$Dataset)
plotbetas(dataset=fem.win.precip[[1]]$Dataset)
plotwin(dataset=fem.win.precip[[1]]$Dataset)
plotbest(dataset=fem.win.precip[[1]]$Dataset,
         bestmodel=fem.win.precip[[1]]$BestModel,
         bestmodeldata=fem.win.precip[[1]]$BestModelData)

fem.win.daylen <- slidingwin(exclude = c(10,30),
                             xvar=list(Photoperiod=xvar$Photoperiod),
                             cdate=xvar$Date,
                             bdate=pred.wae$Date,
                             baseline=baseline,
                             type="relative",
                             stat=c("mean"),
                             func=c("quad"),
                             range=c(270, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5)

##Assess GDD
fem.win.daylen$combos
fem.win.daylen[[1]]
#pvalue(dataset=fem.win.daylen[[1]]$Dataset, datasetrand=fem.win.rand.daylen[[1]],
#       metric='C', sample.size=57)
#plothist(dataset=fem.win.daylen[[1]]$Dataset, datasetrand=fem.win.rand.daylen[[1]])
plotdelta(dataset=fem.win.daylen[[1]]$Dataset)
plotweights(dataset=fem.win.daylen[[1]]$Dataset)
plotbetas(dataset=fem.win.daylen[[1]]$Dataset)
plotwin(dataset=fem.win.daylen[[1]]$Dataset)
plotbest(dataset=fem.win.daylen[[1]]$Dataset,
         bestmodel=fem.win.daylen[[1]]$BestModel,
         bestmodeldata=fem.win.daylen[[1]]$BestModelData)


#Randomize relationships to test for significance
fem.win.rand.temp <- randwin(exclude = c(10,30),
                             xvar=list(Temp=xvar$PredWaterTemp),
                             cdate=xvar$Date,
                             bdate=pred.wae$Date,
                             baseline=baseline,
                             type="relative",
                             stat=c("mean"),
                             func=c("quad"),
                             range=c(270, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5, repeats=10)

fem.win.rand.GDD5 <- randwin(exclude = c(10,30),
                             xvar=list(Temp=xvar$GDD5),
                             cdate=xvar$Date,
                             bdate=pred.wae$Date,
                             baseline=baseline,
                             type="relative",
                             stat=c("sum"),
                             func=c("quad"),
                             range=c(270, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5, repeats=10)

fem.win.rand.precip <- randwin(exclude = c(10,30),
                             xvar=list(Temp=xvar$meanPrecip),
                             cdate=xvar$Date,
                             bdate=pred.wae$Date,
                             baseline=baseline,
                             type="relative",
                             stat=c("sum"),
                             func=c("lin"),
                             range=c(270, 0),
                             cinterval="day",
                             cmissing="method1",
                             k=5, repeats=10)

fem.win.rand.daylen <- randwin(exclude = c(10,30),
                               xvar=list(Temp=xvar$Photoperiod),
                               cdate=xvar$Date,
                               bdate=pred.wae$Date,
                               baseline=baseline,
                               type="relative",
                               stat=c("sum"),
                               func=c("quad"),
                               range=c(270, 0),
                               cinterval="day",
                               cmissing="method1",
                               k=5, repeats=10)


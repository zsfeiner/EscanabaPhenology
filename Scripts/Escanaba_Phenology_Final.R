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

wae <- wae %>%
  group_by(Year)%>%
  mutate(diff=abs(DOY-mean(DOY)))

ggplot(data=filter(wae, diff < 15), aes(x=DOY, y=Year, group=Year, height=Females,fill=after_stat(x))) + 
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
#modG <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(-15,15)))

#modGS <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="cc"),m=2) +
              # s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

modGI <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, by=fYear, k=6,  bs="cc",m=2) +
                    s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

#summary(modG)
#summary(modGS)
summary(modGI)

library(gratia)

#draw(modGS)
#draw(modG)
draw(modGI)

#par(mfrow=c(2,2))
#gam.check(modG)
#gam.check(modGS)
gam.check(modGI)

#AIC(modG)


lim.pred.data <- sel.wae %>%
  group_by(fYear) %>%
  summarize(first=min(DOY)-20, last=max(DOY)+20) %>%
  group_by(fYear) %>%
  group_modify(~ tibble(DOY = seq(.$first, .$last))) %>%
  ungroup()




#pred.data <- tibble(DOY=rep(c(80:150), 61), fYear=as.factor(rep(unique(sel.wae$fYear), each=length(80:150))))
#pred.data
#modG.pred <- as_tibble(predict(modG, newdata=lim.pred.data, se.fit=T, type="response"))
#modG.pred

#modGS.pred <- as_tibble(predict(modGS, newdata=lim.pred.data, type="response", se.fit=T))
#modGS.pred

modGI.pred <- as_tibble(predict(modGI, newdata=lim.pred.data, type="response", se.fit=T))
modGI.pred

#lim.pred.data$modGfit <- modG.pred$fit
#lim.pred.data$modGse <- modG.pred$se.fit
#lim.pred.data$modGSfit <- modGS.pred$fit
#lim.pred.data$modGSse <- modGS.pred$se.fit
lim.pred.data$modGIfit <- modGI.pred$fit
lim.pred.data$modGIse <- modGI.pred$se.fit


# ggplot(data=lim.pred.data, aes(x=DOY, y=modGfit, group=fYear)) + 
#   geom_ribbon(aes(x=DOY, ymin=modGfit-modGse, ymax=modGfit+modGse, group=fYear), alpha=0.4) +
#   geom_line(color="blue") + #ylim(c(-2,3)) +
#   geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
#   theme_classic() + facet_wrap(~fYear, scales="free")
# 
# 
# ggplot(data=lim.pred.data, aes(x=DOY, y=modGSfit, group=fYear)) + 
#   #geom_ribbon(aes(x=DOY, ymin=modGSfit-modGSse, ymax=modGSfit+modGSse, group=fYear), alpha=0.4) +
#   geom_line(color="blue") + #ylim(c(-2,3)) +
#   geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
#   theme_classic() + facet_wrap(~fYear, scales="free")

ggplot(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGIfit-modGIse, ymax=modGIfit+modGIse, group=fYear), alpha=0.4) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + #ylim(c(-2,3)) +
  theme_classic() + facet_wrap(~fYear, scales="free")

# #plot all predictions
# ggplot(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
#   geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
#   geom_line(color="blue") + 
#   geom_line(aes(x=DOY, y=modGfit, group=fYear), color="red") + 
#   geom_line(aes(x=DOY, y=modGSfit, group=fYear), color="green") +
#   theme_classic() + facet_wrap(~fYear, scales="free")


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

ggplot(filter(temps, Year %in% c(2010:2013), DOY %in% c(80:150)), 
       aes(x=as.Date(Date), y=PredWaterTemp, group=Year)) + 
  geom_line(lwd=1, color="blue") +  
  facet_wrap(~Year, nrow=6, scales="free_x") + theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), panel.border=element_rect(fill=NA),
        axis.text.x=element_text(angle=90)) + 
  ylab("Water temperature") + xlab("Date")

format(temps$Date[1:10], "%m/%d")

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


##Make a gamm to see if can detect effects on catches of fish reasonably well
xvar
sel.wae

env.wae <- left_join(sel.wae, xvar) %>% filter(!is.na(PredWaterTemp))
env.wae

plot(Females ~ PredWaterTemp, data=env.wae)
plot(Females ~ meanPrecip, data=env.wae)
plot(Females ~ Photoperiod, data=env.wae)

mod.fem <- gam(Females ~ s(PredWaterTemp, bs="cs", k=6) + s(Photoperiod, bs="cs") +
                 s(meanPrecip, bs='cs') + ti(Photoperiod, PredWaterTemp, bs='cs') +
               s(fYear, bs="re"), data=env.wae, family=quasipoisson(), method="REML")
mod.fem.int <- gam(Females ~ s(meanPrecip, bs='cs') + te(Photoperiod, PredWaterTemp, bs='cs') +
                     s(fYear, bs="re"), data=env.wae, family=quasipoisson(), method="REML")

mod.fem0 <- gam(Females ~ s(PredWaterTemp, bs="cs", k=6) + s(Photoperiod, bs="cs") +
                 s(meanPrecip, bs='cs') +
                 s(fYear, bs="re"), data=env.wae, family=quasipoisson(), method="REML")
summary(mod.fem)
hist(resid(mod.fem))
gam.check(mod.fem)
plot(mod.fem)
BIC(mod.fem, mod.fem0)
anova(mod.fem, mod.fem.int, mod.fem0)

summary(mod.fem)
summary(mod.fem.int)
summary(mod.fem0)
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, theta=45, phi=25, color="topo", ticktype="detailed", xlab="Water temperature", zlab="Response")
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, theta=225, phi=25, color="topo", ticktype="detailed", xlab="Water temperature", zlab="Response")
plot(mod.fem, scheme=1, residuals=T, select=1, cex=3, ylab="Partial effect", xlab="Water temperature")
plot(mod.fem, scheme=1, residuals=T, select=2, cex=3, ylab="Partial effect")


plot(Photoperiod ~ DOY, data=env.wae)
plot(Females ~ DOY, data=env.wae)
plot(Females ~ Photoperiod, data=env.wae)

vis.gam(mod.fem.int, theta=120)
plot(mod.fem.int, scheme=1, too.far=0.1)
plot(mgcViz::getViz(mod.fem))

gampredplot <- ggplot(env.wae, aes(x=DOY, y=Females, group=fYear)) + 
  geom_point() + 
  #geom_vline(data=esc.phen, aes(xintercept=lower, group=fYear)) + 
  #geom_vline(data=esc.phen, aes(xintercept=upper, group=fYear)) + 
  facet_wrap(~fYear, scales='free', ncol=10) + 
  #geom_line(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  geom_line(data=env.wae, aes(x=DOY, y=predict(mod.fem, type="response")), color="blue", inherit.aes=F, lwd=1) + 
  geom_ribbon(data=env.wae, aes(x=DOY, ymin=predict(mod.fem, type="response") - predict(mod.fem, type="response", se=T)$se.fit*1.96,
                                ymax=predict(mod.fem, type="response") + predict(mod.fem, type="response", se=T)$se.fit*1.96),
              color=NA, fill="blue", alpha=0.3) +
  theme_classic()

ggsave("Gam_Pred_plot.png", gampredplot, units="in", dpi=300, height=13, width=20, scale=0.8)

ggplot(env.wae, aes(x=DOY, y=Females)) + 
  geom_point() + 
  facet_wrap(~Year, scales='free') + 
  geom_line(aes(x=DOY, y=predict(mod.fem, type="response"))) + 
  geom_line(aes(x=DOY, y=predict(mod.fem0, type="response")), color="blue")

#Set up environmental variables
xvar <- inner_join(select(temps, GDD0, GDD5, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))
xvar

#Predictions for recruitment
YOYPE
YOYPE$lnRS <- log(YOYPE$Age0PE/YOYPE$AdultPE)
plot(lnRS ~ AdultPE, data=YOYPE)

esc.phen$Year <- as.numeric(as.character(esc.phen$fYear))
rec.dat <- left_join(YOYPE, esc.phen, by=c("Year"="Year"))
rec.dat

rec.dat <- rec.dat %>%
  filter(!is.na(SpawnDOY), Year < 2021) %>%
  mutate(SpawnDate = as.Date(SpawnDOY, origin=paste0((Year-1),"-12-31")))
rec.dat


xvar
r.baseline <- lm(lnRS ~ 1, data=rec.dat)

#Thermal habitat before spawning
rec.win.all <- slidingwin(exclude = c(10,-150),
                                  xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                                  cdate=xvar$Date,
                                  bdate=rec.dat$SpawnDate,
                                  baseline=r.baseline,
                                  type="relative",
                                  stat=c("mean", "slope"),
                                  func=c("lin", "quad"),
                                  range=c(120, -150),
                                  cinterval="day",
                                  cmissing="method1")
rec.win.all$combos

rec.win.all[[1]]
plotdelta(dataset=rec.win.all[[1]]$Dataset)
plotweights(dataset=rec.win.all[[1]]$Dataset)
plotbetas(dataset=rec.win.all[[1]]$Dataset)
plotwin(dataset=rec.win.all[[1]]$Dataset)
plotbest(dataset=rec.win.all[[1]]$Dataset,
         bestmodel=rec.win.all[[1]]$BestModel,
         bestmodeldata=rec.win.all[[1]]$BestModelData)


#Randomize all data
rec.win.all.rand <- randwin(exclude = c(10,-150),
                             xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                             cdate=xvar$Date,
                             bdate=rec.dat$SpawnDate,
                             baseline=r.baseline,
                             type="relative",
                             stat=c("mean","slope"),
                             func=c("lin", "quad"),
                             range=c(120, -150),
                             cinterval="day",
                             cmissing="method1", repeats=100)

pvalue(dataset=rec.win.all[[2]]$Dataset, datasetrand=rec.win.all_rand[[1]],
       metric='AIC', sample.size=55)
plothist(dataset=rec.win.all[[1]]$Dataset, datasetrand=rec.win.all_rand[[1]])

#Do three year groups in 17-19 yr clusters
rec.dat <- rec.dat %>%
  mutate(period = ifelse(Year %in% c(1958:1983), "past",ifelse(Year %in% c(1984:2002), "mid", "present")))

ggplot(rec.dat, aes(x=Year, y=lnRS, color=period)) + 
  geom_point()


pastbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="past"))
rec.win.past <- slidingwin(exclude = c(10,-150),
                            xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                            cdate=xvar$Date,
                            bdate=filter(rec.dat, period=="past")$SpawnDate,
                            baseline=pastbaseline,
                            type="relative",
                            stat=c("mean", "slope"),
                            func=c("lin", "quad"),
                            range=c(120, -150),
                            cinterval="day",
                            cmissing="method1")

midbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="mid"))
rec.win.mid <- slidingwin(exclude = c(10,-150),
                           xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                           cdate=xvar$Date,
                           bdate=filter(rec.dat, period=="mid")$SpawnDate,
                           baseline=midbaseline,
                           type="relative",
                           stat=c("mean","slope"),
                           func=c("lin", "quad"),
                           range=c(120, -150),
                           cinterval="day",
                           cmissing="method1")

presentbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="present"))
rec.win.present <- slidingwin(exclude = c(10,-150),
                           xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                           cdate=xvar$Date,
                           bdate=filter(rec.dat, period=="present")$SpawnDate,
                           baseline=presentbaseline,
                           type="relative",
                           stat=c("mean","slope"),
                           func=c("lin", "quad"),
                           range=c(120, -150),
                           cinterval="day",
                           cmissing="method1")

rec.win.past$combos
rec.win.mid$combos
rec.win.present$combos

plotdelta(dataset=rec.win.past[[1]]$Dataset)
plotdelta(dataset=rec.win.mid[[1]]$Dataset)
plotdelta(dataset=rec.win.present[[1]]$Dataset)

plotdelta(dataset=rec.win.past[[2]]$Dataset)
plotdelta(dataset=rec.win.mid[[2]]$Dataset)
plotdelta(dataset=rec.win.present[[2]]$Dataset)

plotdelta(dataset=rec.win.past[[3]]$Dataset)
plotdelta(dataset=rec.win.mid[[3]]$Dataset)
plotdelta(dataset=rec.win.present[[3]]$Dataset)

plotdelta(dataset=rec.win.past[[4]]$Dataset)
plotdelta(dataset=rec.win.mid[[4]]$Dataset)
plotdelta(dataset=rec.win.present[[4]]$Dataset)

plotbest(dataset=rec.win.past[[4]]$Dataset,
         bestmodel=rec.win.past[[4]]$BestModel,
         bestmodeldata=rec.win.past[[4]]$BestModelData)
plotbest(dataset=rec.win.mid[[4]]$Dataset,
         bestmodel=rec.win.mid[[4]]$BestModel,
         bestmodeldata=rec.win.mid[[4]]$BestModelData)
plotbest(dataset=rec.win.present[[4]]$Dataset,
         bestmodel=rec.win.present[[4]]$BestModel,
         bestmodeldata=rec.win.present[[4]]$BestModelData)
rec.win.past[[1]]
rec.win.mid[[1]]
rec.win.present[[1]]

pastbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="past"))
rec.win.past.rand <- randwin(exclude = c(10,-150),
                           xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                           cdate=xvar$Date,
                           bdate=filter(rec.dat, period=="past")$SpawnDate,
                           baseline=pastbaseline,
                           type="relative",
                           stat=c("mean","slope"),
                           func=c("lin", "quad"),
                           range=c(120, -150),
                           cinterval="day",
                           cmissing="method1", repeats=100)

midbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="mid"))
rec.win.mid.rand <- randwin(exclude = c(10,-150),
                          xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                          cdate=xvar$Date,
                          bdate=filter(rec.dat, period=="mid")$SpawnDate,
                          baseline=midbaseline,
                          type="relative",
                          stat=c("mean","slope"),
                          func=c("lin", "quad"),
                          range=c(120, -150),
                          cinterval="day",
                          cmissing="method1", repeats=100)

presentbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="present"))
rec.win.present.rand <- randwin(exclude = c(10,-150),
                              xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                              cdate=xvar$Date,
                              bdate=filter(rec.dat, period=="present")$SpawnDate,
                              baseline=presentbaseline,
                              type="relative",
                              stat=c("mean","slope"),
                              func=c("lin", "quad"),
                              range=c(120, -150),
                              cinterval="day",
                              cmissing="method1", repeats=100)

#Past
rec.win.past
pvalue(dataset=rec.win.past[[1]]$Dataset, datasetrand=rec.win.past.rand[[1]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.past[[2]]$Dataset, datasetrand=rec.win.past.rand[[2]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.past[[3]]$Dataset, datasetrand=rec.win.past.rand[[3]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.past[[4]]$Dataset, datasetrand=rec.win.past.rand[[4]],
       metric='AIC', sample.size=19)
plothist(dataset=rec.win.past[[1]]$Dataset, datasetrand=rec.win.past.rand[[1]])
plothist(dataset=rec.win.past[[4]]$Dataset, datasetrand=rec.win.past.rand[[4]])

#mid
pvalue(dataset=rec.win.mid[[1]]$Dataset, datasetrand=rec.win.mid.rand[[1]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.mid[[2]]$Dataset, datasetrand=rec.win.mid.rand[[2]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.mid[[3]]$Dataset, datasetrand=rec.win.mid.rand[[3]],
       metric='AIC', sample.size=19)
pvalue(dataset=rec.win.mid[[4]]$Dataset, datasetrand=rec.win.mid.rand[[4]],
       metric='AIC', sample.size=19)
plothist(dataset=rec.win.mid[[1]]$Dataset, datasetrand=rec.win.mid.rand[[1]])

#present
pvalue(dataset=rec.win.present[[1]]$Dataset, datasetrand=rec.win.present.rand[[1]],
       metric='AIC', sample.size=17)
pvalue(dataset=rec.win.present[[2]]$Dataset, datasetrand=rec.win.present.rand[[2]],
       metric='AIC', sample.size=17)
pvalue(dataset=rec.win.present[[3]]$Dataset, datasetrand=rec.win.present.rand[[3]],
       metric='AIC', sample.size=17)
pvalue(dataset=rec.win.present[[4]]$Dataset, datasetrand=rec.win.present.rand[[4]],
       metric='AIC', sample.size=17)
plothist(dataset=rec.win.present[[3]]$Dataset, datasetrand=rec.win.present.rand[[3]])

plotbest(dataset=rec.win.past[[4]]$Dataset,
         bestmodel=rec.win.past[[4]]$BestModel,
         bestmodeldata=rec.win.past[[4]]$BestModelData)
plotbest(dataset=rec.win.mid[[4]]$Dataset,
         bestmodel=rec.win.mid[[4]]$BestModel,
         bestmodeldata=rec.win.mid[[4]]$BestModelData)
plotbest(dataset=rec.win.present[[4]]$Dataset,
         bestmodel=rec.win.present[[4]]$BestModel,
         bestmodeldata=rec.win.present[[4]]$BestModelData)

plotdelta(rec.win.present[[4]]$Dataset)
plotwin(rec.win.present[[4]]$Dataset)
rec.win.past$combos
rec.win.mid$combos
rec.win.present$combos




plotwin(rec.win.past[[1]]$Dataset)
plotwin(rec.win.past[[2]]$Dataset)
plotwin(rec.win.past[[3]]$Dataset)
plotwin(rec.win.past[[4]]$Dataset)

plotwin(rec.win.mid[[1]]$Dataset)
plotwin(rec.win.mid[[2]]$Dataset)
plotwin(rec.win.mid[[3]]$Dataset)
plotwin(rec.win.mid[[4]]$Dataset)

plotwin(rec.win.present[[1]]$Dataset)
plotwin(rec.win.present[[2]]$Dataset)
plotwin(rec.win.present[[3]]$Dataset)
plotwin(rec.win.present[[4]]$Dataset)

xvar
plot(meanPrecip ~ Date, xvar)

rec.win.before$combos
plotwin(filter(rec.win.before[[1]]$Dataset, deltaAICc <= -2), cw=1)
ggplot(data=filter(rec.win.before[[1]]$Dataset, deltaAICc <= -2) %>%
         select(WindowOpen, WindowClose) %>%
         pivot_longer(cols=c(WindowOpen,WindowClose)), aes(x=value, y=name)) + 
  geom_boxplot



#Examine what happens if zoom into each window individually, do 20-80 days post spawning
#Doesn't matter, still not significant
rec.win.2080 <- slidingwin(exclude = c(5,-50),
                             xvar=list(GDD0=xvar$meanPrecip),
                             cdate=xvar$Date,
                             bdate=rec.dat$SpawnDate,
                             baseline=r.baseline,
                             type="relative",
                             stat=c("mean","slope"),
                             func=c("lin"),
                             range=c(-20, -80),
                             cinterval="day",
                             cmissing="method1")
rec.win.2080$combos

rec.win.2080[[1]]
pvalue(dataset=rec.win.2080[[1]]$Dataset, datasetrand=rec.win.2080.rand[[1]],
       metric='AIC', sample.size=55)
plothist(dataset=rec.win.2080[[1]]$Dataset, datasetrand=rec.win.2080.rand[[1]])
plotdelta(dataset=rec.win.2080[[2]]$Dataset)
plotweights(dataset=rec.win.2080[[1]]$Dataset)
plotbetas(dataset=rec.win.2080[[1]]$Dataset)
plotwin(dataset=rec.win.2080[[1]]$Dataset)
plotbest(dataset=rec.win.2080[[2]]$Dataset,
         bestmodel=rec.win.2080[[2]]$BestModel,
         bestmodeldata=rec.win.2080[[2]]$BestModelData)

#Randomize all data
rec.win.2080.rand <- randwin(exclude = c(5,-50),
                               xvar=list(GDD0=xvar$GDD0),
                               cdate=xvar$Date,
                               bdate=rec.dat$SpawnDate,
                               baseline=r.baseline,
                               type="relative",
                               stat=c("mean"),
                               func=c("lin"),
                               range=c(-20, -80),
                               cinterval="day",
                               cmissing="method1", repeats=100)

#Perform on wind on shorter timeseries - not a lot there, hold back unless someone asks for it
wind <- read_csv("ntl17_1_v19.csv")
wind
colSums(is.na(wind))
filter(wind, is.na(avg_wind_speed))
plot(is.na(avg_wind_speed) ~ sampledate, wind)

wind.xvar <- left_join(xvar, select(wind, sampledate, avg_wind_speed), by=c("Date"="sampledate")) %>%
  filter(!is.na(avg_wind_speed))
wind.xvar

wind.base <- lm(lnRS ~ 1, data=filter(rec.dat, Year > 1989))

rec.wind <- slidingwin(exclude = c(10,-150),
                              xvar=list(wind=wind.xvar$avg_wind_speed),
                              cdate=wind.xvar$Date,
                              bdate=filter(rec.dat, Year > 1989)$SpawnDate,
                              baseline=wind.base,
                              type="relative",
                              stat=c("mean","slope"),
                              func=c("lin", "quad"),
                              range=c(120, -150),
                              cinterval="day",
                              cmissing="method2")
rec.wind$combos
pvalue(dataset=rec.wind[[1]]$Dataset, datasetrand=rec.wind.rand[[1]],
       metric='AIC', sample.size=55)
plothist(dataset=rec.wind[[1]]$Dataset, datasetrand=rec.wind.rand[[1]])
plotdelta(dataset=rec.wind[[4]]$Dataset)
plotweights(dataset=rec.wind[[1]]$Dataset)
plotbetas(dataset=rec.wind[[1]]$Dataset)
plotwin(dataset=rec.wind[[1]]$Dataset)
plotbest(dataset=rec.wind[[1]]$Dataset,
         bestmodel=rec.wind[[1]]$BestModel,
         bestmodeldata=rec.wind[[1]]$BestModelData)

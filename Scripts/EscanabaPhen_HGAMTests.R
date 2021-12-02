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


ggplot(data=wae, aes(x=DOY, y=Females, group=Year)) + 
  geom_point() + 
  facet_wrap(~Year, scales="free")

ggplot(data=wae, aes(x=DOY, y=Females, group=Year, color=Year)) + 
  geom_line(lwd=1) + geom_point(color="black") + 
  theme_classic() + scale_color_viridis_c() + facet_wrap(~Year, scales="free")


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


modG <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

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


lim.pred.data <- sel.wae %>%
  group_by(fYear) %>%
  summarize(first=min(DOY)-10, last=max(DOY)+10) %>%
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



#Try with GDD as predictor
#Water temperature (under ice modeled from Sparkling lake under ice buoy)
temps <- read_csv(file="ModeledEscanabaWaterTemps_1956.2020.csv",
                  col_types=cols(
                    Date=col_date(format="%m/%d/%Y"),
                    Notes = col_character()))
temps

#Function to calculate degree days from start date to end date by year
#x=temperature file with date and temperature
#y=fish dates
#start=starting date

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
names(Mar15GDD)[4] <- "Mar15GDD"

sel.wae <- left_join(sel.wae, select(Jan1GDD, Date, Jan1GDD)) %>% left_join(select(Oct1GDD, Date, Oct1GDD))
sel.wae <- left_join(sel.wae, select(Mar15GDD, Date, Mar15GDD))

ggplot(data=sel.wae, aes(x=Jan1GDD, y=Females, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=sel.wae, aes(x=Oct1GDD, y=Females, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=sel.wae, aes(x=Mar15GDD, y=Females, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")


modG.GDD <- gam(Females ~ s(Jan1GDD, bs="tp", k=7, m=2) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

modGS.GDD <- gam(Females ~ s(Jan1GDD, bs="tp", k=6) + s(Jan1GDD, fYear, k=6, bs="fs", xt=list(bs="tp"),m=2) +
               s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

modGI.GDD <- gam(Females ~ s(Jan1GDD, bs="tp", k=6) + s(Jan1GDD, by=fYear, k=6,  bs="tp",m=2) +
               s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

summary(modG.GDD)
summary(modGS.GDD)
summary(modGI.GDD)

library(gratia)

draw(modGS.GDD)
draw(modG.GDD)
draw(modGI.GDD)

par(mfrow=c(2,2))
gam.check(modG.GDD)
gam.check(modGS.GDD)
gam.check(modGI.GDD)


#Predictions
GDD.pred.data <- sel.wae %>%
  filter(!is.na(Jan1GDD)) %>%
  group_by(fYear) %>%
  summarize(first=min(Jan1GDD, na.rm=T)-100, last=max(Jan1GDD, na.rm=T)+300) %>%
  group_by(fYear) %>%
  group_modify(~ tibble(Jan1GDD = seq(.$first, .$last))) %>%
  ungroup()




#pred.data <- tibble(DOY=rep(c(80:150), 61), fYear=as.factor(rep(unique(sel.wae$fYear), each=length(80:150))))
#pred.data
modG.GDD.pred <- as_tibble(predict(modG.GDD, newdata=GDD.pred.data, se.fit=T, type="response"))
modG.GDD.pred

modGS.GDD.pred <- as_tibble(predict(modGS.GDD, newdata=GDD.pred.data, type="response", se.fit=T))
modGS.GDD.pred

modGI.GDD.pred <- as_tibble(predict(modGI.GDD, newdata=GDD.pred.data, type="response", se.fit=T))
modGI.GDD.pred

GDD.pred.data$modGfit <- modG.GDD.pred$fit
GDD.pred.data$modGse <- modG.GDD.pred$se.fit
GDD.pred.data$modGSfit <- modGS.GDD.pred$fit
GDD.pred.data$modGSse <- modGS.GDD.pred$se.fit
GDD.pred.data$modGIfit <- modGI.GDD.pred$fit
GDD.pred.data$modGIse <- modGI.GDD.pred$se.fit


ggplot(data=GDD.pred.data, aes(x=Jan1GDD, y=modGfit, group=fYear)) + 
  geom_ribbon(aes(x=Jan1GDD, ymin=modGfit-modGse, ymax=modGfit+modGse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=Jan1GDD, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")


ggplot(data=GDD.pred.data, aes(x=Jan1GDD, y=modGSfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGSfit-modGSse, ymax=modGSfit+modGSse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=Jan1GDD, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")

ggplot(data=GDD.pred.data, aes(x=Jan1GDD, y=modGIfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGIfit-modGIse, ymax=modGIfit+modGIse, group=fYear), alpha=0.4) +
  geom_point(data=sel.wae, aes(x=Jan1GDD, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + #ylim(c(-2,3)) +
  theme_classic() + facet_wrap(~fYear, scales="free")

#plot all predictions
ggplot(data=GDD.pred.data, aes(x=Jan1GDD, y=modGIfit, group=fYear)) + 
  geom_point(data=sel.wae, aes(x=Jan1GDD, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + 
  geom_line(aes(x=Jan1GDD, y=modGfit, group=fYear), color="red") + 
  geom_line(aes(x=Jan1GDD, y=modGSfit, group=fYear), color="green") +
  theme_classic() + facet_wrap(~fYear, scales="free")

#Try temperature

sel.wae
temp
temps

sel.wae <- left_join(sel.wae, select(temps, Date, PredWaterTemp), by="Date")
sel.wae

ggplot(sel.wae, aes(x=PredWaterTemp, y=Females, group=fYear)) + geom_point() + facet_wrap(~fYear, scales="free")

modG.temp <- gam(Females ~ s(PredWaterTemp, bs="tp", k=7, m=2) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

modGS.temp <- gam(Females ~ s(PredWaterTemp, bs="tp", k=6) + s(PredWaterTemp, fYear, k=6, bs="fs", xt=list(bs="tp"),m=2) +
                   s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

modGI.temp <- gam(Females ~ s(PredWaterTemp, bs="tp", k=6) + s(PredWaterTemp, by=fYear, k=6,  bs="tp",m=2) +
                   s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML")

summary(modG.temp)
summary(modGS.temp)
summary(modGI.temp)

library(gratia)

draw(modGS.temp)
draw(modG.temp)
draw(modGI.temp)

par(mfrow=c(2,2))
gam.check(modG.temp)
gam.check(modGS.temp)
gam.check(modGI.temp)


#Predictions
temp.pred.data <- sel.wae %>%
  filter(!is.na(PredWaterTemp)) %>%
  group_by(fYear) %>%
  summarize(first=min(PredWaterTemp, na.rm=T)-5, last=max(PredWaterTemp, na.rm=T)+5) %>%
  group_by(fYear) %>%
  group_modify(~ tibble(PredWaterTemp = seq(.$first, .$last))) %>%
  ungroup()




#pred.data <- tibble(DOY=rep(c(80:150), 61), fYear=as.factor(rep(unique(sel.wae$fYear), each=length(80:150))))
#pred.data
modG.temp.pred <- as_tibble(predict(modG.temp, newdata=temp.pred.data, se.fit=T, type="response"))
modG.temp.pred

modGS.temp.pred <- as_tibble(predict(modGS.temp, newdata=temp.pred.data, type="response", se.fit=T))
modGS.temp.pred

modGI.temp.pred <- as_tibble(predict(modGI.temp, newdata=temp.pred.data, type="response", se.fit=T))
modGI.temp.pred

temp.pred.data$modGfit <- modG.temp.pred$fit
temp.pred.data$modGse <- modG.temp.pred$se.fit
temp.pred.data$modGSfit <- modGS.temp.pred$fit
temp.pred.data$modGSse <- modGS.temp.pred$se.fit
temp.pred.data$modGIfit <- modGI.temp.pred$fit
temp.pred.data$modGIse <- modGI.temp.pred$se.fit


ggplot(data=temp.pred.data, aes(x=PredWaterTemp, y=modGfit, group=fYear)) + 
  geom_ribbon(aes(x=PredWaterTemp, ymin=modGfit-modGse, ymax=modGfit+modGse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=PredWaterTemp, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")


ggplot(data=temp.pred.data, aes(x=PredWaterTemp, y=modGSfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGSfit-modGSse, ymax=modGSfit+modGSse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=PredWaterTemp, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")

ggplot(data=temp.pred.data, aes(x=PredWaterTemp, y=modGIfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGIfit-modGIse, ymax=modGIfit+modGIse, group=fYear), alpha=0.4) +
  geom_point(data=sel.wae, aes(x=PredWaterTemp, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + #ylim(c(-2,3)) +
  theme_classic() + facet_wrap(~fYear, scales="free")

#plot all predictions
ggplot(data=temp.pred.data, aes(x=PredWaterTemp, y=modGIfit, group=fYear)) + 
  geom_point(data=sel.wae, aes(x=PredWaterTemp, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue") + 
  geom_line(aes(x=PredWaterTemp, y=modGfit, group=fYear), color="red") + 
  geom_line(aes(x=PredWaterTemp, y=modGSfit, group=fYear), color="green") +
  theme_classic() + facet_wrap(~fYear, scales="free")

modGI.all <- gam(Females ~ s(PredWaterTemp, bs="tp", k=6) + s(DOY, bs="cc", k=6) + #s(PredWaterTemp, by=fYear, k=6,  bs="tp",m=2) +
                  s(DOY, by=fYear, k=6, bs="cc", m=2)  + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

ggplot(sel.wae, aes(x=DOY, y=Females, group=fYear)) + 
  geom_point() + ylim(c(0,175)) +
  geom_smooth(method="lm", se=F, color="gray80") + theme_classic()





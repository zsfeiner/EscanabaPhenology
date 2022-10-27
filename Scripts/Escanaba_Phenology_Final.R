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
  summarize(Total=n(), Females=sum(Sex=="Female", na.rm=T), PropFemale=Females/Total) %>%
  group_by(Year) %>%
  mutate(cumTotal = cumsum(Total), cumFemale=cumsum(Females), propcumTotal=cumTotal/sum(Total), propcumFemale=cumFemale/sum(Females),
         PropFemaleTotal = Females/sum(Total), scDOY = scale(DOY, center=T, scale=F), fempres = ifelse(Females > 0, 1, 0))
wae <- wae %>% group_by(Year) %>% mutate(scFem = scale(Females))


ggplot(data=wae, aes(x=DOY, y=Females, group=Year)) + 
  geom_point() + 
  facet_wrap(~Year, scales="free")

ggplot(data=wae[abs(wae$scDOY)<=15,], aes(x=scDOY, y=Females, group=Year)) + 
  geom_point() + geom_smooth(se=F) +
  facet_wrap(~Year, scales="free")

ggplot(data=filter(wae, abs(scDOY)<=15), aes(x=DOY, y=Year, group=Year, height=scFem+2,fill=stat(x))) + 
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
  mutate(scDOY = scale(DOY))

hist(sel.wae$Females)
hist(sel.wae$DOY)
summary(sel.wae$DOY)

#Fit poisson models to check model fit on raw DOY
modG_pois <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))
modGS_pois <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="cc"),m=2) +
                    s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))
modGI_pois <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, by=fYear, k=6,  bs="cc",m=2) +
                    s(fYear, bs="re"), data=sel.wae, family=poisson(), method="REML", knots=list(DOY=c(80,150)))

modG_pois$aic
modGS_pois$aic
modGI_pois$aic


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
gam.check(modG_pois)
gam.check(modGS_pois)
gam.check(modGI_pois)

AIC(modG)


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

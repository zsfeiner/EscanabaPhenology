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


modG <- gam(Females ~ s(DOY, bs="tp", k=7) + s(fYear, bs="re"), data=sel.wae, family=poisson, method="REML")

modGS <- gam(Females ~ s(DOY, bs="tp", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="tp"),m=2) +
               s(fYear, bs="re"), data=sel.wae, family=poisson, method="REML")

modGI.gamma <- gam(Females+0.1 ~ s(DOY, bs="tp", k=7) + s(DOY, by=fYear, k=7,  bs="tp",m=2) +
               s(fYear, bs="re"), data=sel.wae, family=Gamma(link="log"), method="REML")
summary(modGI)
sel.wae$predGI <- predict(modGI, type="response")

ggplot(data=sel.wae, aes(x=DOY, y=Females, group=Year, color=Year)) + 
  geom_point(color="black") + geom_line(aes(x=DOY, y=predGI)) + 
  theme_classic() + facet_wrap(~Year, scales="free")


gam.check(modGS)
AIC(modGS, modG)

library(gratia)

draw(modGS)
draw(modG)

pred.data <- tibble(DOY=rep(c(90:130), 61), fYear=as.factor(rep(unique(sel.wae$fYear), each=length(90:130))))
pred.data
modG.pred <- as_tibble(predict(modG, newdata=pred.data, se.fit=T, type="response"))
modG.pred

modGS.pred <- as_tibble(predict(modGS, newdata=pred.data, type="response", se.fit=T))
modGS.pred

modGI.pred <- as_tibble(predict(modGI, newdata=pred.data, type="response", se.fit=T))
modGI.pred

pred.data$modGfit <- modG.pred$fit
pred.data$modGse <- modG.pred$se.fit
pred.data$modGSfit <- modGS.pred$fit
pred.data$modGSse <- modGS.pred$se.fit
pred.data$modGIfit <- modGI.pred$fit
pred.data$modGIse <- modGI.pred$se.fit


summary(modG)

ggplot(data=pred.data, aes(x=DOY, y=modGfit, group=fYear)) + 
  geom_ribbon(aes(x=DOY, ymin=modGfit-modGse, ymax=modGfit+modGse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")


ggplot(data=pred.data, aes(x=DOY, y=modGSfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGSfit-modGSse, ymax=modGSfit+modGSse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")

ggplot(data=pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  #geom_ribbon(aes(x=DOY, ymin=modGIfit-modGIse, ymax=modGIfit+modGIse, group=fYear), alpha=0.4) +
  geom_line(color="blue") + #ylim(c(-2,3)) +
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  theme_classic() + facet_wrap(~fYear, scales="free")


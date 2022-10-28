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


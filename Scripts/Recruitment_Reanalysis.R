#Revised recruitment windows analysis using straight YOY PE

#load libraries
library(tidyverse)
library(climwin)
library(lubridate)
library(mgcv)
library(geosphere)
library(ggridges)

#Read in walleye spawning phenology data and summarize to daily catches
###Walleye data - available upon request, see .gitignore
dat <- read.csv("./Data/AllEscanabaWAE_1946_2021.csv")
dat

#Recategorize sex
dat$Sex <- ifelse(dat$Sex %in% c(1, "Male","M"), "Male", ifelse(dat$Sex %in% c(2, "Female","F"), "Female",
                                                                ifelse(dat$Sex %in% c(3, "Unknown","U"), "Unknown", NA)))
#Clean up data and convert counts to individual fish
dat <- mutate_at(dat, .vars="Date", .funs=as.Date, format="%m/%d/%Y")
dat <- as_tibble(dat)
dat <- uncount(dat, Count)
dat

#Fix lengths
dat$Length <- ifelse(is.na(dat$Length), dat$TLbin, dat$Length)

ggplot(data=filter(dat, CorrGear==2, Month %in% c(3:5), Year != 1948, Sex=="Female"), aes(y=Length, x=Year)) + 
  geom_point() + geom_smooth(method="loess")


ggplot(data=filter(dat, CorrGear==2, Month %in% c(3:5), Year != 1948, Sex=="Female"),
       aes(x=Length)) + geom_histogram() + facet_wrap(~Year, scales="free_y")


#Summarize catch by date
wae <- 
  dat %>%
  filter(CorrGear==2, Month %in% c(3:5), Year != 1948) %>%
  group_by(Year, DOY, Date) %>%
  summarize(Total=n(), Females=sum(Sex=="Female", na.rm=T), PropFemale=Females/Total, length=mean(Length[Sex=="Female"], na.rm=T))

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


#Fit quasi-poisson models, better fit - take a long time, just fit modGI for now
#modG <- gam(Females ~ s(DOY, bs="cc", k=7) + s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(-15,15)))

#modGS <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, fYear, k=6, bs="fs", xt=list(bs="cc"),m=2) +
# s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

modGI <- gam(Females ~ s(DOY, bs="cc", k=6) + s(DOY, by=fYear, k=6,  bs="cc",m=2) +
               s(fYear, bs="re"), data=sel.wae, family=quasipoisson(), method="REML", knots=list(DOY=c(80,150)))

#Model comparison
anova(modG, modGS, modGI)

#Model explanatory power
summary(modG) #0.567, 60.7%
summary(modGS) #0.76, 82.4%
summary(modGI) #(0.805, 85.2%)

plot(modGI)

#par(mfrow=c(2,2))
#gam.check(modG)
#gam.check(modGS)
gam.check(modGI)


#Create dataset for predictions
lim.pred.data <- sel.wae %>%
  group_by(fYear) %>%
  summarize(first=min(DOY)-20, last=max(DOY)+20) %>%
  group_by(fYear) %>%
  group_modify(~ tibble(DOY = seq(.$first, .$last))) %>%
  ungroup()

#plot observed vs predicted and histogram of modGI residuals
par(mfrow=c(2,1))
plot(sel.wae$Females , predict(modGI, type="response"), xlab="Observed female catch", ylab="Predicted female catch")
abline(a=0, b=1)
hist(resid(modGI), xlab="Residual", main="")

#Generate model predictions
modGI.pred <- as_tibble(predict(modGI, newdata=lim.pred.data, type="response", se.fit=T))
modGI.pred

lim.pred.data$modGIfit <- modGI.pred$fit
lim.pred.data$modGIse <- modGI.pred$se.fit

#Create predicted timeseries plot for supplemental material
modgipredplot <- ggplot(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  geom_point(data=sel.wae, aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue", lwd=1) + #ylim(c(-2,3)) +
  theme(base_size=20) + theme_bw() + facet_wrap(~fYear, scales="free_y", ncol=4) + 
  xlab("DOY") + ylab("Catch (# female walleye)") 

#Create predicted timeseries plot for supplemental material
ggplot(data=filter(lim.pred.data, fYear %in% c(1957,1958,2019)), aes(x=DOY, y=modGIfit, group=fYear)) + 
  geom_point(data=filter(sel.wae, fYear %in% c(1957,1958,2019)), aes(x=DOY, y=Females, group=fYear), color="black", inherit.aes=F) + 
  geom_line(color="blue", lwd=1) +  theme_classic(base_size=20) +
  theme(panel.background = element_rect(fill=NA, color="black")) + facet_wrap(~fYear, scales="free_y", ncol=1) + 
  xlab("DOY") + ylab("Catch (# female walleye)") 


modgipredplot
ggsave("./Manuscript/CJFAS Submission/Revisions/ModGIPredPlot_Supp.png", width=14, height=20, units="in", dpi=400, scale=0.75)  

######
#######Use climate windows package to predict model-predicted female catches based on environmental variables
######

#########Bring in environmental data#############
###Precipitation###
precip <- read_csv("./Data/NOAA_NCDC_NorthernWIWeather_1940_2020.csv", 
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

#Generate table for sites used in analysis
precip.sites <- precip %>%
  filter(NAME %in% c("EAGLE RIVER, WI US", "MINOCQUA, WI US", "PHELPS, WI US",
                     "RAINBOW RSVR LK TOMAHAWK, WI US", "REST LAKE, WI US")) %>%
  filter(year(DATE) > 1946) %>%
  group_by(STATION, NAME, LATITUDE, LONGITUDE) %>%
  summarize(coverage = sum(!is.na(PRCP))/n())
print.data.frame(precip.sites)


#Summarize precipitation as daily mean across all available sites
meanprecip <-
  precip %>%
  filter(NAME %in% c("EAGLE RIVER, WI US", "MINOCQUA, WI US", "PHELPS, WI US",
                     "RAINBOW RSVR LK TOMAHAWK, WI US", "REST LAKE, WI US")) %>%
  group_by(DATE) %>%
  summarize(meanPrecip=mean(PRCP, na.rm=T), meanSnow=mean(SNOW, na.rm=T), 
            meanMaxTemp=mean(TMAX, na.rm=T), meanMinTemp=mean(TMIN, na.rm=T), meanTemp = mean(c(meanMaxTemp, meanMinTemp)), N=sum(!is.na(PRCP)), meanWind=mean(WDMV, na.rm=T))
meanprecip

newprecip <- meanprecip %>% mutate(DOY = as.numeric(strftime(DATE, format="%j"))) %>%
  mutate(fYear=year(DATE))

ggplot(filter(newprecip, fYear %in% c(1957,1958,2019)), aes(x=DOY, y=meanPrecip)) +
  geom_point() + facet_wrap(~fYear, ncol=1, scales="free") + theme_classic(base_size=20) + labs(y="Mean precipitation (mm)") + 
  theme(panel.border=element_rect(fill=0, color="black"))
    


####Photoperiod####
#Bring in photoperiod using geosphere package
Photo <- data.frame(Date = seq(min(precip$DATE), as.Date("12/31/2020", format="%m/%d/%Y"), by=1), Year = NA, DOY = NA, DayLen = NA)
Photo$Year <- year(Photo$Date)
Photo$DOY <- as.numeric(strftime(Photo$Date, format="%j"))
Photo$DayLen <- daylength(Photo$Date, lat = 46.06413190)
Photo


######Water temperature (under ice modeled from Sparkling lake under ice buoy, see PredictWaterTemps.R)#####
temps <- read_csv(file="./Data/ModeledEscanabaWaterTemps_1956.2020.csv",
                  col_types=cols(
                    Date=col_date(format="%m/%d/%Y"),
                    Notes = col_character()))
temps

newtemps <- temps %>% mutate(DOY = as.numeric(strftime(Date, format="%j")))
newtemps <- newtemps %>% mutate(model = ifelse(temp.iced=="Ice", "Predicted","Observed"))

ggplot(filter(newtemps, Year %in% c(1957,1958,2019)), aes(x=DOY, y=PredWaterTemp, color=model)) +
  geom_point() + facet_wrap(~Year, ncol=1, scales="free") + theme_classic(base_size=20) + labs(y="Mean water temperature") + 
  scale_color_manual(values=c("red","blue")) + theme(panel.background=element_rect(fill=NA, color="black"), legend.position="none")

#####Ice data##### - Not needed
#Ice on and off dates and water temps during them
#ice <- read_csv("./Data/Escanaba_IceDates_1956_2020.csv",
#                col_types=cols(
#                  IceOn=col_date(format="%m/%d/%Y"),
#                  IceOff=col_date(format="%m/%d/%Y")))
#ice$IceOnDOY <- as.numeric(strftime(ice$IceOn, format="%j"))
#ice$IceOffDOY <- as.numeric(strftime(ice$IceOff, format="%j"))
#ice

#Get walleye recruitment data, available on request
YOYPE <- read_csv("./Data/EscanabaStockRecruitPE.csv")

#Get recruitment sampling dates, available on request
YOYdates <- read_csv("./Data/WAE_YOY_CompositeThrough2019.csv")
YOYdates <- YOYdates %>%
  filter(MWBC == 2339900) %>% mutate(Date = as.Date(DATE, format="%m/%d/%Y")) %>%
  group_by(YEAR) %>%
  summarize(Date=mean(Date))
YOYdates


#######################################################################
############Walleye spawning phenology analysis########################
#######################################################################

#Create GDD0
temps$GDD0 <- ifelse(temps$PredWaterTemp > 0, temps$PredWaterTemp, 0)
temps$GDD5 <- ifelse(temps$PredWaterTemp >= 5, temps$PredWaterTemp, 0)

##Make a gamm to see if can detect effects on catches of fish reasonably well
#Set up environmental variables
xvar <- inner_join(select(temps, GDD0, GDD5, PredWaterTemp, Photoperiod, Date), select(meanprecip, DATE, meanPrecip), by=c("Date"="DATE"))
xvar

sel.wae

plot(PredWaterTemp ~ Date, xvar)

#Join environmental variables to walleye catches and drop missing water temperature data
env.wae <- left_join(sel.wae, xvar) %>% filter(!is.na(PredWaterTemp))
env.wae


#GAMM for effects of water temperature, photoperiod, and precipitation on female catches
mod.fem <- gam(Females ~ s(PredWaterTemp, bs="cs", k=6) + s(Photoperiod, bs="cs") +
                 s(meanPrecip, bs='cs') + ti(Photoperiod, PredWaterTemp, bs='cs') +
                 s(fYear, bs="re"), data=env.wae, family=nb(link="log"), method="REML")

summary(mod.fem)
hist(resid(mod.fem))
gam.check(mod.fem)
plot(mod.fem)

#Check correlation of water temperature and photoperiod
cor.test(env.wae$PredWaterTemp, env.wae$Photoperiod)


#View interaction from a few different angles
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, theta=45, phi=25, color="topo", ticktype="detailed", xlab="Water temperature", zlab="Response")
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, theta=225, phi=25, color="topo", ticktype="detailed", xlab="Water temperature", zlab="Response")
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, xlab="Water temperature", plot.type="contour", main="", contour.col="gray20",
        color="topo")

#Create manuscript figure 1
png(file="./Manuscript/CJFAS Submission/Revisions/Fig1_GAMMResult.png", res=500, width=7, height=6, units="in", pointsize=12)
par(mfrow=c(2,2), mai = c(0.55, 0.75, 0.1, 0.1), mgp=c(2,0.75,0))
plot(mod.fem, scheme=1, residuals=T, select=1, cex=3, ylab="Partial effect", xlab=expression("Water temperature ("*degree*"C)"))
mtext(side=3, text="a)", adj=0.015, padj=2)
plot(mod.fem, scheme=1, residuals=T, select=2, cex=3, ylab="Partial effect")
mtext(side=3, text="b)", adj=0.015, padj=2)
plot(mod.fem, scheme=1, residuals=T, select=3, cex=3, ylab="Partial effect", xlab="Mean precipitation (mm)")
mtext(side=3, text="c)", adj=0.015, padj=2)
vis.gam(mod.fem, view=c("PredWaterTemp","Photoperiod"), type="link", too.far=0.1, xlab=expression("Water temperature ("*degree*"C)"), plot.type="contour", main="", contour.col="gray20", ylab="Photoperiod (hr)", color="topo")
mtext(side=3, text="d)", adj=0.015, padj=2)
dev.off()

#Create manuscript figure 2
gampredplot <- ggplot(env.wae, aes(x=DOY, y=Females, group=fYear)) + 
  geom_point() +
  #geom_vline(data=esc.phen, aes(xintercept=lower, group=fYear)) + 
  #geom_vline(data=esc.phen, aes(xintercept=upper, group=fYear)) + 
  facet_wrap(~fYear, scales='free', ncol=6) + 
  #geom_line(data=lim.pred.data, aes(x=DOY, y=modGIfit, group=fYear)) + 
  geom_line(data=env.wae, aes(x=DOY, y=predict(mod.fem, type="response")), color="blue", inherit.aes=F, lwd=1) + 
  geom_ribbon(data=env.wae, aes(x=DOY, ymin=predict(mod.fem, type="response") - predict(mod.fem, type="response", se=T)$se.fit*1.96,
                                ymax=predict(mod.fem, type="response") + predict(mod.fem, type="response", se=T)$se.fit*1.96),
              color=NA, fill="blue", alpha=0.3) +
  theme_classic(base_size=14) + scale_x_continuous(breaks=seq(50,200,2)) + 
  theme(axis.text.x=element_text(size=9))
gampredplot
ggsave("./Manuscript/CJFAS Submission/Revisions/Fig2_Gam_Pred_plot.png", gampredplot, units="in", dpi=300, height=20, width=20, scale=0.8)


##########################################################################
#################CLIMATE WINDOWS ANALYSIS OF RECRUITMENT##################
##########################################################################

##Bring in age0 PEs and CVs, available on request - calculate them here
newPE_raw <- read_csv(file = "Data/EscanabaAge0PE_raw.csv")
newPE <- 
newPE_raw %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  summarize(CM = sum(C.M, na.rm=T),
            R = sum(R, na.rm=T)) %>%
  mutate(PE = CM/(R+1), invPE = (R+1)/CM, invVar=R/(CM^2)) %>%
  mutate(CV = sqrt(invVar)/invPE*100) %>%
  mutate(PE.sd = (CV/100) * PE)

print(newPE, n=Inf)

#Plot timeseries of PEs
ggplot(newPE, aes(x=Year, y=PE)) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(x=Year, ymax=PE+PE.sd, ymin=PE-PE.sd))

#Drop 2009 - bad PE, no recaps
newPE <- filter(newPE, Year != 2009)
newPE

plot(newPE$CV)

#Predictions for recruitment - use straight PEs
newPE

esc.phen <- lim.pred.data %>%
  group_by(fYear) %>%
  summarize(SpawnDOY = DOY[which.max(modGIfit)])

#Create recruitment dataset
esc.phen$Year <- as.numeric(as.character(esc.phen$fYear))
rec.dat <- left_join(newPE, esc.phen, by=c("Year"="Year"))
rec.dat

print(rec.dat, n=Inf)

#Set reference date as spawn date
rec.dat <- rec.dat %>%
  filter(!is.na(SpawnDOY), Year < 2021) %>%
  mutate(SpawnDate = as.Date(SpawnDOY, origin=paste0((Year-1),"-12-31")))
rec.dat
rec.dat <- rename(rec.dat, "Age0PE" = "PE")

#Weight by CV
rec.dat$weights <- 1/rec.dat$CV
xvar

#Set climwin baseline as intercept only model
r.baseline <- lm(Age0PE ~ 1, data=rec.dat, weights=weights)
r.baseline

ggplot(newPE, aes(x=Year, y=PE)) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(x=Year, ymax=PE+PE.sd, ymin=PE-PE.sd)) + 
  geom_hline(yintercept=coef(r.baseline), lty=2) + theme_classic(base_size=20) + 
  labs(y="Age-0 abundance")

#Perform climwin analysis using mean and slope of temperature and precipitation
#slidingwin automatically does all combinations of functions and stats, but
#we will ignore the slope of precipitation for interpretation
#Include all windows 120 days before spawning and 150 days after, up to 10 days in duration
a <- Sys.time()

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
#Examine combos
rec.win.all$combos

nrow(rec.win.all[[1]]$Dataset)


#Randomize all data for calculation of PdAICc
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

#Print p values and plot histograms
for (i in c(1,2,3,5,6,7)) {
  print(pvalue(dataset=rec.win.all[[i]]$Dataset, datasetrand=rec.win.all.rand[[i]],
               metric='AIC', sample.size=53))
  plothist(dataset=rec.win.all[[i]]$Dataset, datasetrand=rec.win.all.rand[[i]])  
}


#For supplemental diagram, just fit a window 20-40 days post spawning
#Just for figure making, not used
# diag.mean.lin <- singlewin(xvar=list(GDD0=xvar$GDD0),
#                       cdate=xvar$Date,
#                       bdate=rec.dat$SpawnDate,
#                       baseline=r.baseline,
#                       type="relative",
#                       stat=c("mean"),
#                       func=c("lin"),
#                       range=c(40, 20),
#                       cinterval="day",
#                       cmissing="method1")
# 
# diag.mean.quad <- singlewin(xvar=list(GDD0=xvar$GDD0),
#                            cdate=xvar$Date,
#                            bdate=rec.dat$SpawnDate,
#                            baseline=r.baseline,
#                            type="relative",
#                            stat=c("mean"),
#                            func=c("quad"),
#                            range=c(40, 20),
#                            cinterval="day",
#                            cmissing="method1")
# 
# diag.slope.lin <- singlewin(xvar=list(GDD0=xvar$GDD0),
#                            cdate=xvar$Date,
#                            bdate=rec.dat$SpawnDate,
#                            baseline=r.baseline,
#                            type="relative",
#                            stat=c("slope"),
#                            func=c("lin"),
#                            range=c(40, 20),
#                            cinterval="day",
#                            cmissing="method1")
# 
# diag.slope.quad <- singlewin(xvar=list(GDD0=xvar$GDD0),
#                            cdate=xvar$Date,
#                            bdate=rec.dat$SpawnDate,
#                            baseline=r.baseline,
#                            type="relative",
#                            stat=c("slope"),
#                            func=c("quad"),
#                            range=c(40, 20),
#                            cinterval="day",
#                            cmissing="method1")
# 
# diag.precip.mean.lin <- singlewin(xvar=list(Precip=xvar$meanPrecip),
#                            cdate=xvar$Date,
#                            bdate=rec.dat$SpawnDate,
#                            baseline=r.baseline,
#                            type="relative",
#                            stat=c("mean"),
#                            func=c("lin"),
#                            range=c(40, 20),
#                            cinterval="day",
#                            cmissing="method1")
# 
# 
# diag.precip.mean.quad <- singlewin(xvar=list(Precip=xvar$meanPrecip),
#                             cdate=xvar$Date,
#                             bdate=rec.dat$SpawnDate,
#                             baseline=r.baseline,
#                             type="relative",
#                             stat=c("mean"),
#                             func=c("quad"),
#                             range=c(40, 20),
#                             cinterval="day",
#                             cmissing="method1")
# 
# 
# 
# plotbest(diag.mean.lin$Dataset, diag.mean.lin$BestModel, diag.mean.lin$BestModelData)
# plotbest(diag.mean.quad$Dataset, diag.mean.quad$BestModel, diag.mean.quad$BestModelData)
# plotbest(diag.slope.lin$Dataset, diag.slope.lin$BestModel, diag.slope.lin$BestModelData)
# plotbest(diag.slope.quad$Dataset, diag.slope.quad$BestModel, diag.slope.quad$BestModelData)
# plotbest(diag.precip.mean.lin$Dataset, diag.precip.mean.lin$BestModel, diag.precip.mean.lin$BestModelData)
# plotbest(diag.precip.mean.quad$Dataset, diag.precip.mean.quad$BestModel, diag.precip.mean.quad$BestModelData)
# 
# diag.mean.lin
# diag.mean.quad
# diag.slope.lin
# diag.slope.quad
# diag.precip.mean.lin
# diag.precip.mean.quad

b <- Sys.time() - a

#Make giant figure for all outputs of full timeseries model
library(ggpubr)
alltime_suppl_plot <- ggpubr::ggarrange(
  plotwin(rec.win.all[[1]]$Dataset) + xlab(expression("Linear mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[1]]$Dataset),
  plotweights(rec.win.all[[1]]$Dataset),
  plothist(dataset=rec.win.all[[1]]$Dataset, datasetrand=rec.win.all.rand[[1]]),
  plotwin(rec.win.all[[2]]$Dataset) + xlab(expression("Linear mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[2]]$Dataset),
  plotweights(rec.win.all[[2]]$Dataset),
  plothist(dataset=rec.win.all[[2]]$Dataset, datasetrand=rec.win.all.rand[[2]]),
  plotwin(rec.win.all[[3]]$Dataset) + xlab(expression("Linear slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[3]]$Dataset),
  plotweights(rec.win.all[[3]]$Dataset),
  plothist(dataset=rec.win.all[[3]]$Dataset, datasetrand=rec.win.all.rand[[3]]),
  plotwin(rec.win.all[[5]]$Dataset) + xlab(expression("Quadratic mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[5]]$Dataset),
  plotweights(rec.win.all[[5]]$Dataset),
  plothist(dataset=rec.win.all[[5]]$Dataset, datasetrand=rec.win.all.rand[[5]]),
  plotwin(rec.win.all[[6]]$Dataset) + xlab(expression("Quadratic mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[6]]$Dataset),
  plotweights(rec.win.all[[6]]$Dataset),
  plothist(dataset=rec.win.all[[6]]$Dataset, datasetrand=rec.win.all.rand[[6]]),
  plotwin(rec.win.all[[7]]$Dataset) + xlab(expression("Quadratic slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.all[[7]]$Dataset),
  plotweights(rec.win.all[[7]]$Dataset),
  plothist(dataset=rec.win.all[[7]]$Dataset, datasetrand=rec.win.all.rand[[7]]),
  nrow=6, ncol=4)
rec.win.all$combos
alltime_suppl_plot

ggsave("./Manuscript/CJFAS Submission/Revisions/AllTimeResults_SupplPlot_revised.png", alltime_suppl_plot,
       height=20,width=20, units="in", dpi=300, scale=1)

###Make diagram plot
ggpubr::ggarrange(
 
  plotdelta(rec.win.all[[1]]$Dataset),

  plothist(dataset=rec.win.all[[1]]$Dataset, datasetrand=rec.win.all.rand[[1]]),
 
  plotdelta(rec.win.all[[2]]$Dataset),
 
  plothist(dataset=rec.win.all[[2]]$Dataset, datasetrand=rec.win.all.rand[[2]]),

  plotdelta(rec.win.all[[3]]$Dataset),
 
  plothist(dataset=rec.win.all[[3]]$Dataset, datasetrand=rec.win.all.rand[[3]]),

  plotdelta(rec.win.all[[5]]$Dataset),
 
  plothist(dataset=rec.win.all[[5]]$Dataset, datasetrand=rec.win.all.rand[[5]]),
 
  plotdelta(rec.win.all[[6]]$Dataset),
 
  plothist(dataset=rec.win.all[[6]]$Dataset, datasetrand=rec.win.all.rand[[6]]),

  plotdelta(rec.win.all[[7]]$Dataset),
 
  plothist(dataset=rec.win.all[[7]]$Dataset, datasetrand=rec.win.all.rand[[7]]),
  nrow=2, ncol=6)



#Do three year groups in 17-19 yr clusters to test for non-stationarity in climate windows
rec.dat <- rec.dat %>%
  mutate(period = ifelse(Year %in% c(1958:1983), "past",ifelse(Year %in% c(1984:2002), "mid", "present")))

#Show different periods
recplot <- ggplot(rec.dat, aes(x=Year, y=Age0PE, color=period)) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(x=Year, ymax=Age0PE+PE.sd, ymin=Age0PE-PE.sd)) + 
  theme_classic() + ylab("Age-0 walleye") + labs(color="Period")

ggsave("./Manuscript/CJFAS Submission/Revisions/Recruitmentperiods_SupplPlot_revised.png", recplot,
       height=8,width=10, units="in", dpi=400, scale=0.7)

#Perform climwin on past data following other methods from above
pastbaseline <- lm(Age0PE ~ 1, weights=weights, data=filter(rec.dat, period=="past"))
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

#Perform climwin on mid data following other methods from above
midbaseline <- lm(Age0PE ~ 1, weights=weights,  data=filter(rec.dat, period=="mid"))
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

#Perform climwin on late (present) data following other methods from above
presentbaseline <- lm(Age0PE ~ 1, weights=weights, data=filter(rec.dat, period=="present"))
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

#Examine model combinations and outputs
rec.win.past$combos
rec.win.mid$combos
rec.win.present$combos

plotweights(dataset=rec.win.past[[1]]$Dataset)
plotweights(dataset=rec.win.past[[2]]$Dataset)
plotweights(dataset=rec.win.past[[3]]$Dataset)
plotweights(dataset=rec.win.past[[5]]$Dataset)
plotweights(dataset=rec.win.past[[6]]$Dataset)
plotweights(dataset=rec.win.past[[7]]$Dataset)

plotweights(dataset=rec.win.mid[[1]]$Dataset)
plotweights(dataset=rec.win.mid[[2]]$Dataset)
plotweights(dataset=rec.win.mid[[3]]$Dataset)
plotweights(dataset=rec.win.mid[[5]]$Dataset)
plotweights(dataset=rec.win.mid[[6]]$Dataset)
plotweights(dataset=rec.win.mid[[7]]$Dataset)

plotweights(dataset=rec.win.present[[1]]$Dataset) #temp mean lin
plotweights(dataset=rec.win.present[[2]]$Dataset) #Precip mean lin
plotweights(dataset=rec.win.present[[3]]$Dataset) #temp slope lin
plotweights(dataset=rec.win.present[[5]]$Dataset) #temp mean quad
plotweights(dataset=rec.win.present[[6]]$Dataset) #precip mean quad
plotweights(dataset=rec.win.present[[7]]$Dataset) #temp slope quad


#Computationally the randomization process takes 10s of hours
#Just perform randomization on best performing window, statistic, and function for each variables
#Only consider mean precipitation, but consider mean and slope of GDD0

rec.win.past$combos

#Use rec.win.past[[6]] which is quadratic mean precip
#Use rec.win.past[[7]] which is quadratic slope GDD0

pastbaseline <- lm(Age0PE ~ 1, weights=weights, data=filter(rec.dat, period=="past"))

#Assess quadratic mean precip - not significant
plotdelta(rec.win.past[[6]]$Dataset)
plotweights(rec.win.past[[6]]$Dataset)
plotwin(rec.win.past[[6]]$Dataset)
plotbest(dataset=rec.win.past[[6]]$Dataset,
         bestmodel=rec.win.past[[6]]$BestModel,
         bestmodeldata=rec.win.past[[6]]$BestModelData)
rec.win.past.rand_Precip.mean.quad <- randwin(exclude = c(10,-150),
                                              xvar=list(Precip=xvar$meanPrecip),
                                              cdate=xvar$Date,
                                              bdate=filter(rec.dat, period=="past")$SpawnDate,
                                              baseline=pastbaseline,
                                              type="relative",
                                              stat=c("mean"),
                                              func=c("quad"),
                                              range=c(120, -150),
                                              cinterval="day",
                                              cmissing="method1", repeats=100)
pvalue(dataset=rec.win.past[[6]]$Dataset, datasetrand=rec.win.past.rand_Precip.mean.quad[[1]],
       metric='AIC', sample.size=19)


#Assess quadratic slope GDD0 - not significant
rec.win.past$combos[7,]
plotdelta(rec.win.past[[7]]$Dataset)
plotbest(dataset=rec.win.past[[7]]$Dataset,
         bestmodel=rec.win.past[[7]]$BestModel,
         bestmodeldata=rec.win.past[[7]]$BestModelData)
rec.win.past.rand_Temp.slope.quad <- randwin(exclude = c(10,-150),
                                             xvar=list(GDD0=xvar$GDD0),
                                             cdate=xvar$Date,
                                             bdate=filter(rec.dat, period=="past")$SpawnDate,
                                             baseline=pastbaseline,
                                             type="relative",
                                             stat=c("slope"),
                                             func=c("quad"),
                                             range=c(120, -150),
                                             cinterval="day",
                                             cmissing="method1", repeats=100)
pvalue(dataset=rec.win.past[[7]]$Dataset, datasetrand=rec.win.past.rand_Temp.slope.quad[[1]],
       metric='AIC', sample.size=19)


#Repeat for mid
rec.win.mid$combos[6:7,]
#Use rec.win.mid[[6]] which is quad mean precip
#Use rec.win.mid[[7]] which is quad slope GDD0

##Assess quadratic mean precip - not significant
midbaseline <- lm(Age0PE ~ 1, weights=weights, data=filter(rec.dat, period=="mid"))
plotdelta(rec.win.mid[[6]]$Dataset)
plotbest(dataset=rec.win.mid[[6]]$Dataset,
         bestmodel=rec.win.mid[[6]]$BestModel,
         bestmodeldata=rec.win.mid[[6]]$BestModelData)
rec.win.mid.rand_Precip.mean.quad <- randwin(exclude = c(10,-150),
                                            xvar=list(Precip=xvar$meanPrecip),
                                            cdate=xvar$Date,
                                            bdate=filter(rec.dat, period=="mid")$SpawnDate,
                                            baseline=midbaseline,
                                            type="relative",
                                            stat=c("mean"),
                                            func=c("quad"),
                                            range=c(120, -150),
                                            cinterval="day",
                                            cmissing="method1", repeats=100)
pvalue(dataset=rec.win.mid[[6]]$Dataset, datasetrand=rec.win.mid.rand_Precip.mean.quad[[1]],
       metric='AIC', sample.size=18)

#Assess quadratic slope GDD0 - not significant
plotdelta(rec.win.mid[[7]]$Dataset)
plotbest(dataset=rec.win.mid[[7]]$Dataset,
         bestmodel=rec.win.mid[[7]]$BestModel,
         bestmodeldata=rec.win.mid[[7]]$BestModelData)
rec.win.mid.rand_Temp.slope.quad <- randwin(exclude = c(10,-150),
                                           xvar=list(GDD0=xvar$GDD0),
                                           cdate=xvar$Date,
                                           bdate=filter(rec.dat, period=="mid")$SpawnDate,
                                           baseline=midbaseline,
                                           type="relative",
                                           stat=c("slope"),
                                           func=c("quad"),
                                           range=c(120, -150),
                                           cinterval="day",
                                           cmissing="method1", repeats=100)
pvalue(dataset=rec.win.mid[[7]]$Dataset, datasetrand=rec.win.mid.rand_Temp.slope.quad[[1]],
       metric='AIC', sample.size=18)

#Repeat for present
rec.win.present$combos[6:7,]
#Use rec.win.present[[6]] which is quadratic mean precip
#Use rec.win.present[[7]] which is quadratic slope GDD0

#Assess quadratic mean precip - not significant
presentbaseline <- lm(Age0PE ~ 1, weights=weights, data=filter(rec.dat, period=="present"))
plotdelta(rec.win.present[[6]]$Dataset)
plotwin(rec.win.present[[6]]$Dataset)
plotweights(rec.win.present[[6]]$Dataset)
plotbest(dataset=rec.win.present[[6]]$Dataset,
         bestmodel=rec.win.present[[6]]$BestModel,
         bestmodeldata=rec.win.present[[6]]$BestModelData)
rec.win.present.rand_Precip.mean.quad <- randwin(exclude = c(10,-150),
                                                 xvar=list(Precip=xvar$meanPrecip),
                                                 cdate=xvar$Date,
                                                 bdate=filter(rec.dat, period=="present")$SpawnDate,
                                                 baseline=presentbaseline,
                                                 type="relative",
                                                 stat=c("mean"),
                                                 func=c("quad"),
                                                 range=c(120, -150),
                                                 cinterval="day",
                                                 cmissing="method1", repeats=100)
pvalue(dataset=rec.win.present[[6]]$Dataset, datasetrand=rec.win.present.rand_Precip.mean.quad[[1]],
       metric='AIC', sample.size=16)

#Assess quadratic slope GDD0 - not significant
plotdelta(rec.win.present[[7]]$Dataset)
plotbest(dataset=rec.win.mid[[7]]$Dataset,
         bestmodel=rec.win.mid[[7]]$BestModel,
         bestmodeldata=rec.win.mid[[7]]$BestModelData)
rec.win.present.rand_Temp.slope.quad <- randwin(exclude = c(10,-150),
                                                xvar=list(GDD0=xvar$GDD0),
                                                cdate=xvar$Date,
                                                bdate=filter(rec.dat, period=="present")$SpawnDate,
                                                baseline=presentbaseline,
                                                type="relative",
                                                stat=c("slope"),
                                                func=c("quad"),
                                                range=c(120, -150),
                                                cinterval="day",
                                                cmissing="method1", repeats=100)
pvalue(dataset=rec.win.present[[7]]$Dataset, datasetrand=rec.win.present.rand_Temp.slope.quad[[1]],
       metric='AIC', sample.size=16)

# #Further evaluate rec.win.present[[6]] effect of mean precip in late time period
# ##NO LONGER SIGNIFICANT - COMMENTED OUT
# plotweights(rec.win.present[[6]]$Dataset)
# 
# #Pull out model fits
# presdat <- rec.win.present[[6]]$Dataset
# 
# #12 models comprise 95% AICc weight set
# cumsum(presdat$ModWeight)
# presdatset <- presdat[1:12,]
# presdatset
# 
# #Model averaged coefficients and median climate window open and close
# weighted.mean(presdatset$ModelBeta, w=presdatset$ModWeight)
# weighted.mean(presdatset$ModelBetaQ, w=presdatset$ModWeight)
# weighted.mean(presdatset$ModelInt, w=presdatset$ModWeight)
# median(presdatset$WindowOpen)
# median(presdatset$WindowClose)
# 
# #Create manuscript Figure 3 - results of rec.win.present[[6]]
# deltaplot <- plotdelta(rec.win.present[[6]]$Dataset)
# windowplot <- plotwin(dataset=rec.win.present[[6]]$Dataset) + 
#   ggtitle(paste("Climate window range for", (0.95 * 100), "% confidence set"))
# histplot <- plothist(dataset=rec.win.present[[6]]$Dataset, datasetrand=rec.win.present.rand_Precip.mean.quad[[1]])
# bestplot <- plotbest(dataset=rec.win.present[[6]]$Dataset,
#                      bestmodel=rec.win.present[[6]]$BestModel,
#                      bestmodeldata=rec.win.present[[6]]$BestModelData) + ylab("Recruitment response") + xlab("Mean precipitation")
# 
# bestdat <- rec.win.present[[6]]$BestModelData
# bestdat$pred <- predict(rec.win.present[[6]]$BestModel)
# bestdat$UCI <- predict(rec.win.present[[6]]$BestModel, se=T)$se * 1.96 + bestdat$pred
# bestdat$LCI <- bestdat$pred - predict(rec.win.present[[6]]$BestModel, se=T)$se * 1.96
# 
# predat <- tibble(climate=seq(min(bestdat$climate), max(bestdat$climate), 0.1)) %>%
#   mutate("I(climate^2)" = climate^2)
# 
# preds <- predict(rec.win.present[[6]]$BestModel, newdata=predat, se=T)
# 
# predat <- predat %>%
#   mutate(pred=preds$fit,
#          UCI=preds$fit + preds$se.fit*1.96,
#          LCI=preds$fit - preds$se.fit*1.96)
# predat
# 
# bestmodplot <- ggplot(bestdat, aes(x=climate, y=yvar)) + 
#   geom_point(size=3, alpha=0.5, shape=21, fill="dark grey") + 
#   geom_line(data=predat, aes(x=climate, y=pred), inherit.aes=F, lwd=1) + 
#   geom_ribbon(data=predat, aes(x=climate, ymax=UCI, ymin=LCI), inherit.aes=F, alpha=0.3) + 
#   labs(title = "Output of best model", x = "Mean precipitation (mm)", 
#        y = "ln(R/S)") +
#   theme_climwin() 
# 
# library(ggpubr)
# 
# fig3 <- ggarrange(deltaplot, windowplot, histplot, bestmodplot, labels=c("a)", "b)","c)","d)"), nrow=2, ncol=2)
# fig3
# ggsave("./Manuscript/Figures/Fig3_RecentPrecipResult.png", fig3, dpi=500, units="in", width=7, height=6, scale=1.5)



##########################CREATE SUPPLEMENTAL FIGURES#########################
past_suppl_plot <- ggpubr::ggarrange(
  plotwin(rec.win.past[[1]]$Dataset) + xlab(expression("Linear mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[1]]$Dataset),
  plotweights(rec.win.past[[1]]$Dataset),
  #plothist(dataset=rec.win.past[[1]]$Dataset, datasetrand=rec.win.past.rand[[1]]),
  plotwin(rec.win.past[[2]]$Dataset) + xlab(expression("Linear mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[2]]$Dataset),
  plotweights(rec.win.past[[2]]$Dataset),
  #plothist(dataset=rec.win.past[[2]]$Dataset, datasetrand=rec.win.past.rand[[2]]),
  plotwin(rec.win.past[[3]]$Dataset) + xlab(expression("Linear slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[3]]$Dataset),
  plotweights(rec.win.past[[3]]$Dataset),
  #plothist(dataset=rec.win.past[[3]]$Dataset, datasetrand=rec.win.past.rand[[3]]),
  plotwin(rec.win.past[[5]]$Dataset) + xlab(expression("Quadratic mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[5]]$Dataset),
  plotweights(rec.win.past[[5]]$Dataset),
  #plothist(dataset=rec.win.past[[5]]$Dataset, datasetrand=rec.win.past.rand[[5]]),
  plotwin(rec.win.past[[6]]$Dataset) + xlab(expression("Quadratic mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[6]]$Dataset),
  plotweights(rec.win.past[[6]]$Dataset),
  #plothist(dataset=rec.win.past[[6]]$Dataset, datasetrand=rec.win.past.rand[[6]]),
  plotwin(rec.win.past[[7]]$Dataset) + xlab(expression("Quadratic slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.past[[7]]$Dataset),
  plotweights(rec.win.past[[7]]$Dataset),
  #plothist(dataset=rec.win.past[[7]]$Dataset, datasetrand=rec.win.past.rand[[7]]),
  nrow=6, ncol=3)
past_suppl_plot
mid_suppl_plot <- ggpubr::ggarrange(
  plotwin(rec.win.mid[[1]]$Dataset) + xlab(expression("Linear mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[1]]$Dataset),
  plotweights(rec.win.mid[[1]]$Dataset),
  #plothist(dataset=rec.win.mid[[1]]$Dataset, datasetrand=rec.win.mid.rand[[1]]),
  plotwin(rec.win.mid[[2]]$Dataset) + xlab(expression("Linear mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[2]]$Dataset),
  plotweights(rec.win.mid[[2]]$Dataset),
  #plothist(dataset=rec.win.mid[[2]]$Dataset, datasetrand=rec.win.mid.rand[[2]]),
  plotwin(rec.win.mid[[3]]$Dataset) + xlab(expression("Linear slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[3]]$Dataset),
  plotweights(rec.win.mid[[3]]$Dataset),
  #plothist(dataset=rec.win.mid[[3]]$Dataset, datasetrand=rec.win.mid.rand[[3]]),
  plotwin(rec.win.mid[[5]]$Dataset) + xlab(expression("Quadratic mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[5]]$Dataset),
  plotweights(rec.win.mid[[5]]$Dataset),
  #plothist(dataset=rec.win.mid[[5]]$Dataset, datasetrand=rec.win.mid.rand[[5]]),
  plotwin(rec.win.mid[[6]]$Dataset) + xlab(expression("Quadratic mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[6]]$Dataset),
  plotweights(rec.win.mid[[6]]$Dataset),
  #plothist(dataset=rec.win.mid[[6]]$Dataset, datasetrand=rec.win.mid.rand[[6]]),
  plotwin(rec.win.mid[[7]]$Dataset) + xlab(expression("Quadratic slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.mid[[7]]$Dataset),
  plotweights(rec.win.mid[[7]]$Dataset),
  #plothist(dataset=rec.win.mid[[7]]$Dataset, datasetrand=rec.win.mid.rand[[7]]),
  nrow=6, ncol=3)
mid_suppl_plot

present_suppl_plot <- ggpubr::ggarrange(
  plotwin(rec.win.present[[1]]$Dataset) + xlab(expression("Linear mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[1]]$Dataset),
  plotweights(rec.win.present[[1]]$Dataset),
  #plothist(dataset=rec.win.present[[1]]$Dataset, datasetrand=rec.win.present.rand[[1]]),
  plotwin(rec.win.present[[2]]$Dataset) + xlab(expression("Linear mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[2]]$Dataset),
  plotweights(rec.win.present[[2]]$Dataset),
  #plothist(dataset=rec.win.present[[2]]$Dataset, datasetrand=rec.win.present.rand[[2]]),
  plotwin(rec.win.present[[3]]$Dataset) + xlab(expression("Linear slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[3]]$Dataset),
  plotweights(rec.win.present[[3]]$Dataset),
  #plothist(dataset=rec.win.present[[3]]$Dataset, datasetrand=rec.win.present.rand[[3]]),
  plotwin(rec.win.present[[5]]$Dataset) + xlab(expression("Quadratic mean GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[5]]$Dataset),
  plotweights(rec.win.present[[5]]$Dataset),
  #plothist(dataset=rec.win.present[[5]]$Dataset, datasetrand=rec.win.present.rand[[5]]),
  plotwin(rec.win.present[[6]]$Dataset) + xlab(expression("Quadratic mean precipitation")) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[6]]$Dataset),
  plotweights(rec.win.present[[6]]$Dataset),
  #plothist(dataset=rec.win.present[[6]]$Dataset, datasetrand=rec.win.present.rand[[6]]),
  plotwin(rec.win.present[[7]]$Dataset) + xlab(expression("Quadratic slope GDD"[0])) + theme(axis.title.y = element_text(size=17)),
  plotdelta(rec.win.present[[7]]$Dataset),
  plotweights(rec.win.present[[7]]$Dataset),
  #plothist(dataset=rec.win.present[[7]]$Dataset, datasetrand=rec.win.present.rand[[7]]),
  nrow=6, ncol=3)
present_suppl_plot

ggsave("./Manuscript/CJFAS Submission/Revisions/PastResults_SupplPlot.png", past_suppl_plot,
       height=20,width=20, units="in", dpi=300, scale=1)
ggsave("./Manuscript/CJFAS Submission/Revisions/MidResults_SupplPlot.png", mid_suppl_plot,
       height=20,width=20, units="in", dpi=300, scale=1)
ggsave("./Manuscript/CJFAS Submission/Revisions/PresentResults_SupplPlot.png", present_suppl_plot,
       height=20,width=20, units="in", dpi=300, scale=1)


#P hist plots for tested relationships - Table 3
rec.win.past$combos

phistplots <- ggarrange(
  plothist(dataset=rec.win.past[[7]]$Dataset, datasetrand=rec.win.past.rand_Temp.slope.quad[[1]]) + ylab("Past (1958-1983)") + theme(axis.title.y = element_text(size=17)) + 
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","   GDD"[0]~"slope")))), hjust=0.05, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  plothist(dataset=rec.win.past[[6]]$Dataset, datasetrand=rec.win.past.rand_Precip.mean.quad[[1]]) +
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","              Mean precipitation")))), hjust=0.27, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  
  plothist(dataset=rec.win.mid[[7]]$Dataset, datasetrand=rec.win.mid.rand_Temp.slope.quad[[1]]) + ylab("Mid (1984-2002)") + theme(axis.title.y = element_text(size=17)) +
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","   GDD"[0]~"slope")))), hjust=0.05, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  plothist(dataset=rec.win.mid[[6]]$Dataset, datasetrand=rec.win.mid.rand_Precip.mean.quad[[1]]) + 
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","              Mean precipitation")))), hjust=0.27, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  
  plothist(dataset=rec.win.present[[7]]$Dataset, datasetrand=rec.win.present.rand_Temp.slope.quad[[1]]) + ylab("Present (2003-2019)") + theme(axis.title.y = element_text(size=17)) +
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","   GDD"[0]~"slope")))), hjust=0.05, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  plothist(dataset=rec.win.present[[6]]$Dataset, datasetrand=rec.win.present.rand_Precip.mean.quad[[1]]) + 
    geom_text(aes(x=-Inf,y=Inf,label=as.character(expression(atop("Quadratic","              Mean precipitation")))), hjust=0.27, vjust=1.2, fontface="plain", size=8, inherit.aes=F, parse=T),
  
  ncol=2, nrow=3
)

phistplots

ggsave("./Manuscript/CJFAS Submission/Revisions/Phist_SupplPlot.png", phistplots,
       height=15,width=15, units="in", dpi=300, scale=1)

##END





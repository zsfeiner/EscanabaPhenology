library(tidyverse)
library(lubridate)
library(viridis)
library(zoo)
library(gamm4)


#####################################################################
###Winter temperature breakpoint model from underice observations in Sparkling Lake, WI
SL.IceDays <- read_delim("./Data/SparklingLakeIceDates.txt", delim="\t")
SL.IceDays <- mutate_at(SL.IceDays, .vars=vars(ends_with("Date")), .funs=as.Date, format="%m/%d/%Y")
SL.IceDays$iceyear <- SL.IceDays$iceoff_year


wint <- read_delim("./Data/SparklingLakeWinterTemps.txt", delim="\t")  
wint$sampledate <- as.Date(wint$sampledate, format="%m/%d/%Y")  
wint <- left_join(wint, dplyr::select(SL.IceDays, IceOnDate, IceOffDate, iceyear), by=c("IceYear"="iceyear"))
wint$Iced <- with(wint, ifelse(sampledate >= (IceOnDate-1) & sampledate <= (IceOffDate+1), "ice","open"))
wint

wint <- filter(wint, Iced=="ice")
wint$WaterTemp[wint$WaterTemp > 3.7 & wint$IceYear==2006] <- NA
wint$WaterTemp <- na.approx(wint$WaterTemp)
wint$avg_air_temp <- na.approx(wint$avg_air_temp)

ggplot(wint, aes(x=sampledate, y=WaterTemp, group=IceYear)) + 
  geom_line() + geom_point() + facet_wrap(~IceYear, scales="free")
ggplot(wint, aes(x=sampledate, y=avg_air_temp, group=IceYear)) + 
  geom_line() + geom_point() + facet_wrap(~IceYear, scales="free")
ggplot(wint, aes(x=sampledate, y=avg_wind_speed_1m, group=IceYear)) + 
  geom_line() + geom_point() + facet_wrap(~IceYear, scales="free")

DaysFrozen <- c(1:sum(wint$IceYear==2000), 
                1:sum(wint$IceYear==2001),
                1:sum(wint$IceYear==2002),
                1:sum(wint$IceYear==2003),
                1:sum(wint$IceYear==2004),
                1:sum(wint$IceYear==2005),
                1:sum(wint$IceYear==2006))
wint$DaysFrozen <- DaysFrozen

fromThaw <- c(sum(wint$IceYear==2000):1, 
              sum(wint$IceYear==2001):1,
              sum(wint$IceYear==2002):1,
              sum(wint$IceYear==2003):1,
              sum(wint$IceYear==2004):1,
              sum(wint$IceYear==2005):1,
              sum(wint$IceYear==2006):1)
wint$fromThaw <- fromThaw


ggplot(wint, aes(x=fromThaw, y=WaterTemp, group=IceYear)) + 
  geom_line(aes(color=IceYear))

names(wint)
wint$f.year <- as.factor(wint$IceYear)
wint$AboveFreezing <- ifelse(wint$avg_air_temp > 0, wint$avg_air_temp, 0)
wint <- 
  wint %>%  group_by(IceYear) %>%
  mutate(cumAboveFreezing = cumsum(AboveFreezing), day1temp=WaterTemp[which(DaysFrozen==1)], cumTemp = cumsum(avg_air_temp),
         tendaytemp=rollapply(avg_air_temp, width=21, FUN=sum, align="right", partial=T),
         tendaywarm=rollapply(AboveFreezing, width=21, FUN=sum, align="right", partial=T))
par(mfrow=c(1,1))
plot(WaterTemp ~ cumAboveFreezing, data=wint, type="l")
plot(WaterTemp ~ cumTemp, data=wint)
plot(WaterTemp ~ tendaytemp, data=wint)
plot(WaterTemp ~ tendaywarm, data=wint)


watertemppreds <- ggplot(dat=wint, aes(x=DaysFrozen, y=WaterTemp, group=IceYear)) + 
  facet_wrap(~IceYear)+
  geom_point() + 
  geom_line(aes(y=predict(WT.mod$gam)), color="blue") + 
  geom_ribbon(aes(ymax=predict(WT.mod$gam)+predict(WT.mod$gam, se=T)$se.fit*1.96,
                  ymin=predict(WT.mod$gam)-predict(WT.mod$gam, se=T)$se.fit*1.96),
              color=NA, alpha=0.3, fill="blue") + theme_bw() + theme(panel.border=element_rect(color="black", fill=NA)) + 
  xlab("Days frozen") + 
  ylab("Water temperature")

ggsave("./Manuscript/SparklingWaterTempPreds.png", watertemppreds, dpi=300, width=10, height=10, units="in")
unique(wint$IceYear)
test <- gamm4(WaterTemp ~ s(avg_air_temp) + s(cumAboveFreezing) + s(fromThaw) + s(tendaytemp) + s(tendaywarm) + day1temp, data=wint)
test
summary(test$gam)
par(mfrow=c(3,2)); plot(test$gam, pages=1, scheme=1, residuals=T)
hist(resid(test))
plot(WaterTemp ~ fromThaw, data=wint, col=wint$IceYear)
points(predict(test$gam) ~ wint$fromThaw, col=wint$IceYear, pch="-")
plot(predict(test$gam) ~ wint$WaterTemp)
predict(test$gam)

WT.mod <- gamm4(WaterTemp ~ s(avg_air_temp) + s(cumAboveFreezing) + s(fromThaw) + s(tendaytemp) + s(tendaywarm) + day1temp, data=wint)
WT.mod
summary(WT.mod$gam)
plot(WT.mod$gam, pages=1, scheme=1, nrow=3)
plot(WaterTemp ~ fromThaw, data=wint, col=wint$IceYear)
points(predict(WT.mod$gam) ~ wint$fromThaw, col=wint$IceYear, pch="-")
plot(predict(WT.mod$gam) ~ wint$WaterTemp)

Sparkling_obsvpred <- ggplot(dat=wint, aes(x=WaterTemp, y=predict(WT.mod$gam))) + 
  geom_point() + theme_bw() + xlab("Observed water temperature") + ylab("Predicted water temperature")
hist(resid(WT.mod$gam))
ggsave("./Manuscript/SparklingObsvPred.png", Sparkling_obsvpred, dpi=300, width=10, height=10, units="in")


#####################################
#Daily water temperature in Escanaba
library(geosphere)
temps <- read_csv("./Data/EscanabaTemps_1956_2020.csv",
                  col_types = cols(
                    Date = col_date(format="%m/%d/%Y"),
                    Notes = col_character()))
temps

#Fill in missing air temperatures from weather station data
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
meanprecip <-
  precip %>%
  filter(NAME %in% c("EAGLE RIVER, WI US", "MINOCQUA, WI US", "PHELPS, WI US",
                     "RAINBOW RSVR LK TOMAHAWK, WI US", "REST LAKE, WI US")) %>%
  group_by(DATE) %>%
  summarize(meanPrecip=mean(PRCP, na.rm=T), meanSnow=mean(SNOW, na.rm=T), 
            meanMaxTemp=mean(TMAX, na.rm=T), meanMinTemp=mean(TMIN, na.rm=T), meanTemp = mean(c(meanMaxTemp, meanMinTemp)), N=sum(!is.na(PRCP)), meanWind=mean(WDMV, na.rm=T))
meanprecip

temps$LowAirTempC[is.na(temps$LowAirTempC)] <- meanprecip$meanMinTemp[meanprecip$DATE %in% temps$Date[is.na(temps$LowAirTempC)]]
temps$HighAirTempC[is.na(temps$HighAirTempC)] <- meanprecip$meanMaxTemp[meanprecip$DATE %in% temps$Date[is.na(temps$HighAirTempC)]]

temps[,-c(1:4,13)] <- apply(temps[,-c(1:4,13)], 2, FUN=na.approx, na.rm=F, maxgap=2)
temps$Photoperiod <- daylength(lat=46.06413190, doy=temps$Date)

temps

plot(Photoperiod ~ Date, filter(temps, Year==2017))

#Add day of year
temps$DOY <- as.numeric(strftime(temps$Date, format="%j"))
temps$DOY

#Predict overwinter water temperatures based on patterns observed in sparkling lake
temps
temps$iceyear <- temps$Year
temps$iceyear[temps$Month %in% c("Jul","Aug","Sep","Oct","Nov","Dec")] <- temps$Year[temps$Month %in% c("Jul","Aug","Sep","Oct","Nov","Dec")]+1

#Means of temps
temps$MeanWaterTemp <- rowMeans(select(temps, WaterTempC.AM, WaterTempC.PM), na.rm=T)
temps$MeanAirTemp <- rowMeans(select(temps, LowAirTempC, HighAirTempC), na.rm=T)
sum(is.na(temps$MeanWaterTemp))
sum(is.na(temps$MeanAirTemp))
temps[is.na(temps$MeanWaterTemp),]


pred.temps <- 
  temps %>%
  #filter(iceyear != 1956) %>%
  group_by(iceyear) %>%
  mutate(temp.iced = ifelse(MeanWaterTemp > 0, "Open","Ice"),
         AboveFreezing = ifelse(MeanAirTemp>0, MeanAirTemp, 0),
         threeweektemp=rollapply(MeanAirTemp, width=21, FUN=sum, align="right", partial=T),
         threeweekwarm=rollapply(AboveFreezing, width=21, FUN=sum, align="right", partial=T)) %>%
  group_by(iceyear, temp.iced) %>%
  mutate(fromThaw = ifelse(temp.iced=="Open", 0, rev(row_number())), cumAboveFreezing=cumsum(AboveFreezing)) %>%
  group_by(iceyear) %>%
  mutate(day1temp = ifelse(iceyear>1956, MeanWaterTemp[which.max(fromThaw)-1], NA),
         PredWaterTemp = ifelse(MeanWaterTemp > 0, MeanWaterTemp,
                                predict(WT.mod$gam, 
                                        newdata=data.frame(day1temp=day1temp,
                                                           avg_air_temp=MeanAirTemp,
                                                           cumAboveFreezing=cumAboveFreezing,
                                                           fromThaw=fromThaw,
                                                           tendaytemp=threeweektemp,
                                                           tendaywarm=threeweekwarm))))

Escanaba_WTmod <- ggplot(filter(pred.temps, MeanWaterTemp>0), aes(x=DOY, y=PredWaterTemp, group=iceyear, color=Year)) +
  geom_line() + facet_wrap(~Year) + theme_minimal() + ylab("Predicted water temperature") + 
  geom_point(data=filter(pred.temps, MeanWaterTemp <= 0), aes(x=DOY, y=PredWaterTemp), color="red", size=0.5, inherit.aes=F)
ggsave("./Manuscript/Escanaba_modeledWT.png", Escanaba_WTmod, dpi=300, width=10, height=7, units="in", scale=1)


print(pred.temps[is.na(pred.temps$MeanWaterTemp),], width=Inf)

select(pred.temps, PredWaterTemp)
pred.temps
write.csv(pred.temps, "ModeledEscanabaWaterTemps_1956.2020.csv", quote=F, row.names=F)


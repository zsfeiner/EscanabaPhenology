library(patchwork)
library(tidyverse)
library(brms)
library(ggplot2)

set.seed(666)

# generate some skewed data
data <- rskew_normal(1e4, mu = 0, sigma = 1, alpha = -5)

# fitting a brms model with a Gaussian likelihood
model_normal <- brm(data ~ 1, family = gaussian(), data = data)

# fitting a brms model with a skew normal likelihood
model_skew <- brm(data ~ 1, family = skew_normal(), data = data)

# posterior predictive checking
pp_check(model_skew, nsamples = 1e2)

model_skew



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
fem <- filter(dat, Sex=="Female", !is.na(DOY), CorrGear==2, Month %in% c(3:5))
hist(fem$DOY)
hist(data)

ggplot(data=fem, aes(x=DOY, group=Year)) + 
  geom_density() + facet_wrap(~Year, scales="free_y")


fem$Year <- as.factor(fem$Year)
seldata <- filter(fem, Year %in% c(2003:2006))
hist(seldata$DOY)

formula <- bf(DOY ~ 0+(1|Year), sigma ~ 0+(1|Year))#, alpha ~ 0+(1|Year))

#myprior <-  get_prior(DOY ~ 0+intercept+(1|Year), alpha ~ 0+intercept+(1|Year), family=skew_normal(), data=seldata)

myprior <- get_prior(formula, family=gaussian(), data=fem)

code <- make_stancode(formula, prior=myprior, family=skew_normal(), data=seldata)
data <- make_standata(formula, prior=myprior, family=skew_normal(), data=seldata)

test_skew <- brm(formula, family=gaussian(), data=fem, prior=myprior,
                 cores=3, chains=3, iter=1000, thin=3, inits=0)

summary(test_skew)
test_skew$fit
pp_check(test_skew, ndraws=100)
plot(test_skew)
test_skew$model



plot(density(rskew_normal(n=1000, mu=116.18, sigma=1.5, alpha=-88.18)))
lines(density(filter(seldata, Year==2003)$DOY))

plot(density(rskew_normal(n=1000, mu=89.30, sigma=0.98, alpha=0.01)), xlim=c(80,100))
lines(density(filter(seldata, Year==2012)$DOY))

plot(density(rskew_normal(n=1000, mu=130.19, sigma=0.07, alpha=0.01)))
lines(density(filter(seldata, Year==2013)$DOY))

plot(density(rskew_normal(n=200, mu=89.30, sigma=0.98, alpha=0.01)), xlim=c(80,100))
lines(density(filter(seldata, Year==2012)$DOY))

plot(table(round(rskew_normal(n=100, mu=116.18, sigma=1.5, alpha=40),0)), xlim=c(80,140))
lines(x=unique(filter(seldata, Year==2003)$DOY), y=table(filter(seldata, Year==2003)$DOY))

plot(table(round(rskew_normal(n=200, mu=130.28, sigma=0.89, alpha=300),0)), xlim=c(80,140))
lines(x=unique(filter(seldata, Year==2013)$DOY), y=table(filter(seldata, Year==2013)$DOY))


      
library(sn)
params <- cp2dp(c(116, 1.87, -0.09), "SN")
params

sims <- rsn(n=1000, dp=params)
sims

hist(sims, freq=F)
hist(filter(seldata, Year==2003)$DOY, add=T, freq=F, col="blue")

plot(density(sims), lty=2, col="red")
lines(density(filter(seldata, Year==2003)$DOY))

params2015 <- cp2dp(c(108, 1.58, -0.19), family="SN")
sims2015 <- rsn(n=1000, dp=params2015)

plot(density(sims), lty=2, col="red", xlim=c(100,125))
lines(density(sims2015), lty=2, col="blue")
lines(density(filter(seldata, Year==2003)$DOY), col="red")
lines(density(filter(seldata, Year==2015)$DOY), col="blue")

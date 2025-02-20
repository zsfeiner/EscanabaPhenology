#Perform climwin on late (present) data following other methods from above
presentbaseline <- lm(lnRS ~ 1, data=filter(rec.dat, period=="present"))
rec.win.present_null <- slidingwin(exclude = c(10,0),
                              xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                              cdate=xvar$Date,
                              bdate=filter(rec.dat, period=="present")$SpawnDate,
                              baseline=presentbaseline,
                              type="relative",
                              stat=c("mean","slope"),
                              func=c("lin", "quad"),
                              range=c(60, 3),
                              cinterval="day",
                              cmissing="method1")

rec.win.present$combos

rec.dat$lnAdultPE <- log(rec.dat$AdultPE)
ricker <- lm(lnRS ~ lnAdultPE, data=filter(rec.dat, period=="mid"))
summary(ricker)

ggplot(rec.dat, aes(x=lnAdultPE, y=lnRS, color=period)) + 
  geom_point() + geom_smooth(method="lm")


presentbaseline <- lm(lnRS ~ lnAdultPE, data=filter(rec.dat, period=="past"))
summary(presentbaseline)
rec.win.present.ricker <- slidingwin(exclude = c(10,0),
                              xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                              cdate=xvar$Date,
                              bdate=filter(rec.dat, period=="mid")$SpawnDate,
                              baseline=lm(lnRS ~ lnAdultPE, data=filter(rec.dat, period=="mid")),
                              type="relative",
                              stat=c("mean","slope"),
                              func=c("lin", "quad"),
                              range=c(60, 3),
                              cinterval="day",
                              cmissing="method1")
rec.win.present.ricker$combos

ggplot(rec.dat, aes(x=AdultPE, y=Age0PE, color=Year)) + 
  geom_point() + geom_line()

plotbest(rec.win.present.ricker[[4]]$Dataset, rec.win.present.ricker[[4]]$BestModel, rec.win.present.ricker[[4]]$BestModelData)

presentbaseline <- lm(log(Age0PE) ~ 1, data=filter(rec.dat, period=="present"))

rec.win.present.R <- slidingwin(exclude = c(10,0),
                                     xvar=list(GDD0=xvar$GDD0, Precip=xvar$meanPrecip),
                                     cdate=xvar$Date,
                                     bdate=filter(rec.dat, period=="present")$SpawnDate,
                                     baseline=presentbaseline,
                                     type="relative",
                                     stat=c("mean","slope"),
                                     func=c("lin", "quad"),
                                     range=c(60, 3),
                                     cinterval="day",
                                     cmissing="method1")
rec.win.present.R$combos



# load in the included demonstration dataset data("demo.siber.data")
siber.example <- createSiberObject(demo.siber.data) 

# The first ellipse is referenced using a character string representation 
# where in "x.y", "x" is the community, and "y" is the group within that 
# community.
ellipse1 <- "1.2" 

# Ellipse two is similarly defined: community 1, group3 
ellipse2 <- "1.3"

# the overlap betweeen the corresponding 95% prediction ellipses is given by: 
maxLikOverlap(ellipse2, ellipse1, siber.example,
                                   p.interval = 0.95, n = 100, do.plot=T)
ellipse95.overlap

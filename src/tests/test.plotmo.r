# test.plotmo.R: regression tests for plotmo
# Many of thewse  tests are culled from man page examples and modified to try to confuse plotmo.
# Many of the plots are plotted twice so you can visually check by comparing
# plots in the same window, they should be substantially the same.
# Stephen Milborrow, Petaluma Jan 2007

Trace = FALSE

dopar <- function(nrows, ncols, caption = "")
{
    cat("                             ", caption, "\n")
    earth:::make.space.for.sub.caption(caption)
    par(mfrow=c(nrows, ncols))
    par(mar = c(3, 3, 1.7, 0.5))
    par(mgp = c(1.6, 0.6, 0))
    par(cex = 0.7)
}
library(earth)
data(ozone1)

caption <- "basic earth test of plotmo"
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, degree1=2, degree2=4, sub.caption=caption, trace=Trace)

caption <- "test 5 x 5 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=51, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", sub.caption=caption, trace=Trace)

caption <- "test 4 x 4 layout with ylab"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=30, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="ozone", sub.caption=caption, trace=Trace)

caption <- "test 3 x 3 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=16, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", sub.caption=caption, trace=Trace)

caption <- "test 2 x 2 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=9, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", sub.caption=caption, trace=Trace)

caption <- "test 1 x 1 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=4, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", sub.caption=caption, trace=Trace)

caption <- "test plotmo basic params"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(3,2,caption)
set.seed(1) # needed for reproducibility because of sample for rug in plotmo
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, sub.caption=caption,
        main="test main", xlab="test xlab", ylab="test ylab", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=F, degree2=4, grid.func=mean, col.persp="white", ngrid=10, phi=40, trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, lty.degree1=2, col.degree1=2, nrug=50, degree2=F, main="nrug=50", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, main="nrug=-1", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, nrug=20, ndegree1=50, degree2=F, main="ndegree1=50 nrug=20", trace=Trace)

caption <- "test plotmo ylim"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(3,3,caption)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, sub.caption=caption, xlab="ylim=default", trace=Trace)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=NA, xlab="ylim=NA", trace=Trace)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=c(0,20), xlab="ylim=c(0,20)", trace=Trace)

# term.plot calls predict.earth with an se parameter, even with termplot(se=FALSE)

caption <- "basic earth test against termplot"
dopar(4,4,caption)
earth:::make.space.for.sub.caption("test caption1")
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, do.par=FALSE, ylim=NA, sub.caption=caption, degree2=FALSE, trace=Trace)
cat("Ignore two warnings: predict.earth ignored argument \"se\"\n")
termplot(a)

caption <- "test change order of earth predictors"
dopar(4,4,caption)
a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
plotmo(a, do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)
termplot(a)

a.global <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
func1 <- function()
{
    caption <- "test call plotmo from within a function, global dataframe"
    dopar(4,4,caption)
    a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
    plotmo(a,        do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)
    plotmo(a.global, do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)

    caption <- "test call plotmo from within a function, local dataframe"
    dopar(4,4,caption)
    ozone1.local <- ozone1[,c(1,2,3,4,5,6,7,8,10)]  # drop vis
    a <- earth(doy ~ humidity + temp + wind, data=ozone1.local, degree=2)
    plotmo(a,        do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)
    plotmo(a.global, do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)

    caption <- "test call plotmo from within a function, local dataframe using x,y interface"
    dopar(4,4,caption)
    x <- ozone1.local[,c(4,5,3)]    # humidty temp wind
    y <- ozone1.local[,9]           # doy
    a <- earth(x, y, degree=2)
    plotmo(a,        do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)
    plotmo(a.global, do.par=FALSE, ylim=NA, sub.caption=caption, degree2=c(1,3), trace=Trace)
}
func1()

caption <- "test earth formula versus to x,y model"
dopar(4,4,caption)
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, do.par=FALSE, sub.caption=caption, trace=Trace)
a <- earth(ozone1[, -1], ozone1[,1], degree=2)
plotmo(a, do.par=FALSE, trace=Trace)

# single predictor
caption <- "test earth(O3~wind, data=ozone1, degree=2), single predictor"
dopar(2,2,caption)
a <- earth(O3~wind, data=ozone1, degree=2)
plotmo(a, trace=Trace)

caption <- "test lm(log(doy) ~ vh+wind+humidity+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, col.response=3, pch.response=20, trace=Trace)
termplot(a)

caption <- "test lm(log(doy) ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, col.resp=3, pch.response=20, clip=FALSE, trace=Trace)
termplot(a)

caption <- "test lm(doy ~ (vh+wind+humidity)^2, data=ozone1)"
dopar(4,3,caption)
a <- lm(doy ~ (vh+wind+humidity)^2, data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NULL, trace=Trace)
# termplot(a) # termplot fails with Error in `[.data.frame`(mf, , i): undefined columns selected

caption <- "test lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,3,caption)
a <- lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NULL, trace=Trace)
termplot(a) # termplot draws a funky second wind plot

caption <- "test lm with data=ozone versus attach(ozone)"
dopar(4,3,caption)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, degree1=c(1,2,4,5), trace=Trace)
attach(ozone1)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh))
plotmo(a, do.par=FALSE, degree1=c(1,2,4,5), trace=Trace)
detach(ozone1)

# commented out because "$" in names is not yet supported
# a <- lm(log(ozone1$doy) ~ I(ozone1$vh*ozone1$wind) + log(ozone1$ibh))
# plotmo(a, trace=Trace)

set.seed(1)
caption <- "test glm, glm(cbind(damage, 6-damage) ~ temp, family=binomial, data=orings)"
dopar(2,2,caption)
library(faraway)
data(orings)
a <- lm(damage/6 ~ temp, data=orings)
plotmo(a, do.par=FALSE, sub.caption=caption, col.response="pink", clip=FALSE, nrug=-1, ylim=c(0,1),
    main="lm(damage/6 ~ temp, data=orings)", trace=Trace)
a <- glm(cbind(damage, 6-damage) ~ temp, family=binomial, data=orings)
plotmo(a, do.par=FALSE, sub.caption=caption, col.response="pink", clip=FALSE, nrug=-1, ylim=c(0,1),
    main="glm(cbind(damage, 6-damage) ~ temp, family=binomial, data=orings)", trace=Trace)
termplot(a)

set.seed(1)
caption <- "test glm(lot2~log(u),data=clotting,family=Gamma)"
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
plotmo(a, do.par=FALSE, sub.caption=caption, col.response=3, clip=FALSE, nrug=-1, trace=Trace)
termplot(a)

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test plot.gam, with mgcv::gam(y ~ s(x) + s(x,z)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z <- runif(n);
y <- test1(cbind(x,z)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z))
plotmo(a, do.par=FALSE, type2="contour", sub.caption=caption, col.response=3, func=test1, col.func=4, pch.func=1, trace=Trace)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", ylim=NA, trace=Trace)
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)

n<-400
sig<-2
set.seed(1)
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f <- f0(x0) + f1(x1) + f2(x2)
e <- rnorm(n, 0, sig)
y <- f + e
test.func <- function(x) f0(x[,1]) + f1(x[,2]) + f2(x[,3])
library(mgcv)
caption <- "test mgcv::gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20)) (and extra persp plot)"
dopar(3,3,caption)
a <- gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20))
plot(a, select=2)
plot(a, select=3)
plot(a, select=1)
plotmo(a, do.par=FALSE, type2="contour", sub.caption=caption, xlab=NULL, main="", func=test.func, trace=Trace)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, theta=-35, trace=Trace)

set.seed(1)
caption <- "test plot.gam, with mgcv::gam(doy~s(wind)+s(humidity,wind)+s(vh)+temp,data=ozone1)"
dopar(3,3,caption)
a <- gam(doy ~ s(wind) + s(humidity,wind) + s(vh) + temp, data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, type2="contour", degree1=c(1,3), xlab=NULL, main="", clip=FALSE, trace=Trace)
plot(a, select=1)
plot(a, select=3)
plot(a, select=2)
plot(a, select=4)

detach("package:mgcv")
library(gam)
caption <- "test gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
set.seed(1)
dopar(3,2,caption)
data(airquality)
airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, col.response=3, trace=Trace)
# $$ termplot gives fishy looking wind plot, plotmo looks ok
termplot(a)
detach("package:gam")

library(mda)
caption <- "test mars and earth (expect not a close match)"
dopar(6,3,caption)
a <- mars( ozone1[, -1], ozone1[,1], degree=2)
b <- earth(ozone1[, -1], ozone1[,1], degree=2)
plotmo(a, do.par=FALSE, sub.caption=caption, trace=Trace)
plotmo(b, do.par=FALSE, trace=Trace)

caption <- "test mars and mars.to.earth(mars) (expect no degree2 for mars)"
dopar(6,3,caption)
a <- mars(ozone1[, -1], ozone1[,1], degree=2)
b <- mars.to.earth(a)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, trace=Trace)
plotmo(b, do.par=FALSE, ylim=NA, trace=Trace)

# test inverse.func and func

caption <- "test inverse.func=exp"
a <- lm(log(Volume) ~ Girth + Height + I(Girth*Height), data=trees)
my.func <- function(x) -60 + 5 * x[,1] + x[,2] / 3
plotmo(a, sub.caption=caption, inverse.func = exp, col.response = "pink", func=my.func, col.func="grey", ndegree1=1000, type2="p", trace=Trace)

# se testing

caption = "se=2, lm(doy~., data=ozone1) versus termplot"
dopar(6,3,caption)
a <- lm(doy~., data=ozone1)
plotmo(a, se=2, do.par=FALSE, trace=Trace, sub.caption=caption)
termplot(a, se=2)

caption <- "test different se options, se=2, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(4,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, se=2, trace=Trace)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, se=2, col.shade.se=0, col.se=1, trace=Trace)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, se=2, col.se=1, trace=Trace)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NULL, se=2, col.se=1, trace=Trace)

caption = "se=2, earth(doy~humidity+temp+wind, data=ozone1) versus termplot (expect no se lines)"
dopar(3,2,caption)
a <- earth(doy~humidity + temp + wind, data=ozone1, degree=2)
cat("Ignore warning: predict.earth ignored argument \"se\"\n")
termplot(a)
cat("Ignore two warnings: predict.earth ignored argument \"se.fit\"\n")
plotmo(a, se=2, do.par=FALSE, ylim=NA, degree2=c(1:2), clip=FALSE, sub.caption=caption, trace=Trace)

caption <- "test se=2, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(2,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, se=2, trace=Trace)
termplot(a, se=TRUE)

caption <- "test se=2 and inverse.func, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(3,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, se=2, trace=Trace)
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NULL, se=2, inverse.func=exp, trace=Trace)
termplot(a, se=TRUE)

caption <- "test se=3, glm(lot2~log(u),data=clotting,family=Gamma)"
set.seed(1)
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
plotmo(a, do.par=FALSE, sub.caption=caption, col.response=4, pch.response=7, clip=FALSE, nrug=-1, se=3, trace=Trace)
termplot(a, se=TRUE)

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test se=2, plot.gam, with mgcv::gam(y ~ s(x) + s(x,z)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z <- runif(n);
y <- test1(cbind(x,z)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z))
plotmo(a, do.par=FALSE, type2="contour", sub.caption=caption, col.response=3, func=test1, col.func=4, se=2, trace=Trace)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", col.image=topo.colors(10), 
        ylim=NA, se=2, trace=Trace, main="topo.colors")
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)

detach("package:mgcv")
library(gam)
set.seed(1)
caption <- "test se=2, gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
dopar(3,2,caption)
data(airquality)
airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
cat("Ignore three warnings: No standard errors (currently) for gam predictions with newdata\n")
plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, col.response=3, se=2, trace=Trace)
termplot(a)
detach("package:gam")

# test factors by changing wind to a factor

ozone1[,"wind"] <- factor(ozone1[,"wind"])

# commented out because factors are not yet supported by plotmo.earth
# caption <- "test wind=factor, earth(O3 ~ ., data=ozone1)"
# a <- earth(doy ~ ., data=ozone1)
# set.seed(1)
# dopar(4,3,caption)
# plotmo(a, col.response="gray", se=2, nrug=-1, do.par=FALSE, sub.caption=caption, trace=Trace)
# termplot(a)

caption <- "test wind=factor, lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone1)"
a <- lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone1)
set.seed(1)
dopar(4,3,caption)
plotmo(a, col.response="gray", se=2, nrug=-1, do.par=FALSE, sub.caption=caption, trace=Trace)
termplot(a, se=TRUE)

caption <- "test test se options like col.se"
dopar(2,2,caption)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, se=2, sub.caption=caption, trace=Trace)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, se=2, lty.se=1, col.se=2, trace=Trace)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, se=2, lty.se=1, col.shade=0, trace=Trace)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, se=2, lty.se=3, col.shade="gray", trace=Trace)

caption <- "test wind=factor, glm(y ~ i + j, family=poisson())"
y <- c(18,17,15,20,10,20,25,13,12)
i <- gl(3,1,9)
j <- gl(3,3)
a <- glm(y ~ i + j, family=poisson())
set.seed(1)
dopar(2,2,caption)
plotmo(a, do.par=F, se=2, nrug=-1, sub.caption=caption, trace=Trace)
termplot(a, se=1, rug=T)

if(length(grep("package:gam", search())))
   detach("package:gam")
caption <- "test wind=factor, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone1)"
library(mgcv)
a <- gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone1)
plotmo(a, se=1, sub.caption=caption, trace=Trace)
caption <- "test wind=factor, clip=TRUE, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone1)"
plotmo(a, se=1, sub.caption=caption, clip=FALSE, trace=Trace)
# termplot doesn't work here so code commented out
# dopar(3,3,caption)
# plotmo(a, do.par=FALSE, trace=Trace)
# termplot(a)

# commented out because se is not supported for gam::predict.gam
# detach("package:mgcv")
# library(gam)
# caption <- "test wind=factor, gam:gam(doy ~ vh + wind + lo(humidity) + lo(vh) + temp, data=ozone1)"
# a <- gam(doy ~ vh + wind + lo(humidity) + lo(vh) + temp, data=ozone1)
# set.seed(1)
# dopar(4,3,caption)
# plotmo(a, do.par=FALSE, sub.caption=caption, ylim=NA, col.response=3, se=2, subcaption=caption, trace=Trace)
# termplot(a, se=1)
# detach("package:gam")

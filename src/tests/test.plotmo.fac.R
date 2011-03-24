# test.plotmo.fac.R: test factor plotting in plotmo. This also tests swapxy, xflip, and yflip
# Stephen Milborrow, Berea Mar 2011

library(earth)
library(rpart)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur
if(!interactive())
    postscript(paper="letter")

cat("==test plotmo with factors==\n")
test.fac.with.rpart <- function(ngrid2=20)
{
    et <- etitanic

    col.response <- as.numeric(et$sex)+2
    et$pclass.fac <- et$pclass
    et$parch.num <- et$parch
    parch.fac <- et$parch
    parch.fac[parch.fac > 3] <- 3
    et$parch.fac <- factor(parch.fac, labels=c( "lev0", "lev1", "lev2", "lev3"))
    et$pclass.num <- as.numeric(et$pclass)
    et$pclass <- et$sex <- et$age <- et$sibsp <- et$parch <- NULL
    cat("names(et):", names(et), "\n") # survived pclass.fac parch.num parch.fac pclass.num

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(4,5))
    par(mar = c(2, 2, 3, 0.5), cex=.6)

    # numeric x numeric
    a2 <- rpart(survived ~ pclass.num+parch.num, data=et)
    set.seed(145)
    plotmo(a2, do.par=F, type2="im",
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a2, do.par=F, type2="con", degree1=NA,
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a2, do.par=F, type2="persp", degree1=NA,
           ngrid2=100, theta=NA, ticktype="d", border=NA, cex.lab=.8, ntick=2)

    # factor x numeric
    a3 <- rpart(survived ~ pclass.fac+parch.num, data=et)
    set.seed(145)
    plotmo(a3, do.par=F, type2="im",
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a3, do.par=F, type2="con", degree1=NA, ylim=c(0,1),
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a3, do.par=F, type2="persp", degree1=NA,
           ngrid2=100, theta=NA, ticktype="d", border=NA, cex.lab=.8, ntick=2)

    # numeric x factor
    a4 <- rpart(survived ~ pclass.num+parch.fac, data=et)
    set.seed(145)
    plotmo(a4, do.par=F, type2="im",
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a4, do.par=F, type2="con", degree1=NA,
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a4, do.par=F, type2="persp", degree1=NA,
           ngrid2=100, theta=NA, ticktype="d", border=NA, cex.lab=.8, ntick=2)

    # factor x factor
    a5 <- rpart(survived ~ pclass.fac+parch.fac, data=et)
    set.seed(145)
    plotmo(a5, do.par=F, type2="im",
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a5, do.par=F, type2="con", degree1=NA,
           col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
    set.seed(145)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           ngrid2=100, theta=NA, ticktype="d", border=NA, cex.lab=.8, ntick=2)
    par(mfrow=c(1,1))
}
test.fac.with.rpart()
cat("==test plotmo swapxy with factors==\n")
test.swapxy.with.rpart <- function(ngrid2=20)
{
    et <- etitanic[c(1:50,300:350,600:650),]

    col.response <- as.numeric(et$sex)+2
    et$pclass.fac <- et$pclass
    et$parch.num <- et$parch
    parch.fac <- et$parch
    parch.fac[parch.fac > 2] <- 2
    et$parch.fac <- factor(parch.fac, labels=c( "lev0", "lev1", "lev2"))
    et$pclass.num <- as.numeric(et$pclass)
    et$pclass <- et$sex <- et$age <- et$sibsp <- et$parch <- NULL
    cat("names(et):", names(et), "\n") # survived pclass.fac parch.num parch.fac pclass.num

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(4,4))
    par(mar = c(2, 2, 5, 0.5), cex=.6)

    # factor x factor
    a5 <- rpart(survived ~ pclass.fac+parch.fac, data=et)
    for(swapxy in c(F,T)) {
        for(xflip in c(F,T))
            for(yflip in c(F,T)) {
                set.seed(145)
                plotmo(a5, do.par=F, type2="im", degree1=NA,
                       swapxy=swapxy, xflip=xflip, yflip=yflip,
                       main=paste("swapxy=", swapxy, "\nxflip=", xflip, "\nyflip=", yflip),
                       col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
                set.seed(145)
                plotmo(a5, do.par=F, type2="con", degree1=NA,
                       swapxy=swapxy, xflip=xflip, yflip=yflip,
                       main=paste("swapxy=", swapxy, "\nxflip=", xflip, "\nyflip=", yflip),
                       col.response=col.response, cex.response=.3, jitter.response=1, pch.response=20)
            }
    }
    par(mfrow=c(2,2))
    set.seed(146)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           swapxy=FALSE, main=paste("swapxy=", FALSE),
           ngrid2=100, theta=NA, ticktype="d", cex.lab=.8, ntick=5)
    set.seed(146)
    plotmo(a5, do.par=F, type2="persp", degree1=NA,
           swapxy=TRUE, main=paste("swapxy=", TRUE),
           ngrid2=100, theta=NA, ticktype="d", cex.lab=.8, ntick=5)
}
test.swapxy.with.rpart()

aflip <- earth(O3~vh + wind + humidity + temp, data=ozone1, degree=2)
col.response<- ifelse(ozone1$O3 == 38, "red", "pink")

# test xflip arg, degree1 plots
par(mfrow=c(2,2))
set.seed(102)
plotmo(aflip, degree1=1:2, degree2=0, do.par=F, col.response=col.response, pch.response=20, nrug=-1, ylab="O3")
plotmo(aflip, degree1=1:2, degree2=F, do.par=F, col.response=col.response, pch.response=20, nrug=-1, ylab="O3", xflip=T, main="xflip=TRUE, degree1 plots")

col.response<- ifelse(ozone1$O3 == 1, "green", "pink")

# test flip args, type2=persp
par(mfrow=c(2,2))
plotmo(aflip, degree1=0, degree2=2, do.par=F, ticktype="d")
plotmo(aflip, degree1=0, degree2=2, do.par=F, tickt="d", swapxy=T, main="swapxy=TRUE")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

# test swapxy args, type2=image
par(mfrow=c(3,3))

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, main="test swapxy on image plots\nreference plot")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, swapxy=T, main="swapxy=T")
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, xflip=T, main="xflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, yflip=T, main="yflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, xflip=T, yflip=T, main="xflip=T, yflip=T")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, swapxy=T, xflip=T, main="swapxy=T, xflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, swapxy=T, yflip=T, main="swapxy=T, yflip=T")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="im", col.response=col.response, pch.response=20, swapxy=T, xflip=T, yflip=T, main="swapxy=T, xflip=T, yflip=T")

# test flip args, type2=contour
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, main="test flip on contour plots\nreference plot")
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, swapxy=T)
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, xflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, yflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, xflip=T, yflip=T)

plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, swapxy=T, xflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, swapxy=T, yflip=T)
plotmo(aflip, degree1=0, degree2=2, do.par=F, type2="con", col.response=col.response, pch.response=20, swapxy=T, xflip=T, yflip=T)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

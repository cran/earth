# test.earth.full.R: test earth

library(earth)
library(mda)
data(ozone1)
data(trees)

print.time <- FALSE         # FALSE for no time results (for diff against reference)
plot.it.default <- TRUE     # TRUE to do plots too, FALSE for speed
options.old <- options()
options(digits = 2)

#--- test examples from man pages ------------------------------------------------------------

cat("--- earth.Rd ----------------------\n")
a <- earth(Volume ~ ., data = trees)
print(summary(a, digits = 2))

set.seed(1)
train.subset <- rep(FALSE, nrow(trees))
train.subset[sample(1:nrow(trees), .8 * nrow(trees))] <- TRUE
test.subset <- !train.subset
a <- earth(Volume~., data=trees[train.subset, ])
yhat <- predict(a, newdata=trees[test.subset, ])
y <- trees[test.subset, "Volume"]
rsq <- 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
print(rsq, digits=2)

cat("--- format.earth.Rd ----------------------\n")
as.func <- function( # convert expression string to func
               object, digits = 8, use.names = TRUE, ...)
  eval(parse(text=paste(
    "function(x)\n",
    "{\n",
    "if(is.vector(x))\n",
    "  x <- matrix(x, nrow = 1, ncol = length(x))\n",
    "with(as.data.frame(x),\n",
    format(object, digits = digits, use.names = use.names, ...),
    ")\n",
    "}\n", sep = "")))
a <- earth(Volume ~ ., data = trees)
my.func <- as.func(a, use.names = FALSE)
print(my.func(c(10,80)))     # yields 17.76888
print(predict(a, c(10,80)))  # yields 17.76888, but is slower
a <- earth(Volume ~ ., data = trees)
cat(format(a))
cat(format(a, use.names = FALSE, add.labels = TRUE))
cat("--- get.nterms.per.degree.Rd ----------------------\n")
a <- earth(O3 ~ ., data = ozone1, degree = 2)
print(get.nterms.per.degree(a))
cat("--- get.nused.preds.per.subset.Rd ----------------------\n")
a <- earth(O3 ~ ., data = ozone1, degree = 2)
print(get.nused.preds.per.subset(a$dirs, a$selected.terms))
print(get.nused.preds.per.subset(a$dirs, a$prune.terms))
cat("--- mars.to.earth.Rd ----------------------\n")
a <- mars(trees[,-3], trees[,3])
a <- mars.to.earth(a)
print(summary(a, digits = 2))
cat("--- model.matrix.Rd ----------------------\n")
a <- earth(Volume ~ ., data = trees)
print(summary(a, decomp = "none"))
bx <- model.matrix(a)
a.lm <- lm(trees$Volume ~ bx[,-1])
print(summary(a.lm))
cat("--- plot.earth.models.Rd ----------------------\n")
a1 <- earth(O3 ~ ., data = ozone1, degree = 2)
a2 <- earth(O3 ~ .-wind, data = ozone1, degree = 2, nk = 31)
a3 <- earth(O3 ~ .-humidity, data = ozone1, degree = 2, nk = 31)
plot.earth.models(list(a1,a2,a3), rlim=c(.6,.8), cum.grid="n")
cat("--- plot.earth.Rd ----------------------\n")
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(a, col.npreds = 1, cum.grid = "none")
cat("--- predict.earth.Rd ----------------------\n")
a <- earth(Volume ~ ., data = trees)
predict(a)
predict(a, c(10,80))
cat("--- reorder.earth.Rd ----------------------\n")
a <- earth(O3 ~ ., data = ozone1, degree = 2)
print(reorder(a, decomp = "none"))
print(a$selected.terms[reorder(a)])
cat("--- update.earth.Rd ----------------------\n")
a <- earth(Volume~ ., data = trees)
print(summary(a, digits = 2))
a <- earth(O3 ~ ., data = ozone1, degree = 2)
print(summary(update(a, formula = O3 ~ . - temp)))
print(summary(update(a, nprune = 8)))

#--- test model building capabilities --------------------------------------------------------

cat("--- test model building capabilities ----------------------\n")
itest <- 0
N <- 100
set.seed(1)
x1 <- runif(N, -1, 1)
x2 <- runif(N, -1, 1)
x3 <- runif(N, -1, 1)
x4 <- runif(N, -1, 1)
x5 <- runif(N, -1, 1)
x6 <- runif(N, -1, 1)
x7 <- runif(N, -1, 1)
x8 <- runif(N, -1, 1)
x9 <- runif(N, -1, 1)
x10 <- runif(N, -1, 1)

make.func <- function(
    obj      = stop("no 'obj' arg"),
    digits   = 14,
    use.names = TRUE,   # use predictor names, else "x[,1]" etc
    ...)                # extra args passed onto format, eg add.labels=TRUE
{
    s <- paste(
        "function(x)\n",
        "{\n",
        "if(is.vector(x))\n",
        "  x <- matrix(x, nrow=1, ncol=length(x))\n",
        "with(as.data.frame(x),\n",
        format(obj, digits=digits, use.names=use.names, ...),
        ")\n",
        "}\n", sep="")

    eval.parent(parse(text=s))
}

# this cross checks that RSq and GRSq claimed by 
# the model versus an independent calc of RSq and GRSq

test.model.rsq <- function(object, x, y, MarsFunc, nCases, nUsedTerms, penalty, RefFunc=NULL, ...)
{
    y1 <- RefFunc(x, ...)
    rss <- sum((y1 - MarsFunc(x))^2)
    rss.null <- sum((y - mean(y))^2)
    gcv.null <- earth:::get.gcv(rss.null, 1, 0, penalty, nCases)
    gcv <- earth:::get.gcv(rss, nUsedTerms, (nUsedTerms-1)/2, penalty, nCases)
    # $$ these happen, look into it
    if(is.finite(object$rsq))
        if(!isTRUE(all.equal(object$rsq, 1 - rss/rss.null)))
            cat("\nWarning: RSq mismatch object$rsq", object$rsq, "calculated RSq", 1 - rss/rss.null)
        else if(!isTRUE(all.equal(object$grsq, 1 - gcv/gcv.null)))
            cat("\nWarning GRSq mismatch object$grsq", object$grsq, "calculated GRSq", 1 - gcv/gcv.null)
}

# this uses the global x and d
test.earth <- function(itest, func, degree=2, nk=51, plotit=plot.it.default, test.rsq=TRUE, trace=0)
{
    cat(sprintf("%-3d", itest), sprintf("%-32s", deparse(substitute(func))), 
        "degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk))
    gc()
    earthTime <- system.time(fite <- earth(d[,-1, drop=FALSE], d[,1], degree=degree,
                                           trace=trace, nk=nk, pmethod="b", fast.k=-1))
    funca <- make.func(fite)
    nCases <- nrow(x)
    penalty <- ifelse(degree>1,3,2)
    nUsedTerms <- sum(fite$selected.terms!=0)
    cat(" nTerms",  sprintf("%-2d", nUsedTerms), "of", sprintf("%-3d ", ncol(fite$bx)))
    if(print.time)
        cat(" time", earthTime[1])
    cat("GRSq", sprintf("%4.2g", fite$grsq))
    sub.caption <- paste(itest, ": ", deparse(substitute(func)), " degree=", degree, " nk=", nk, sep="")
    if(test.rsq)
    test.model.rsq(fite, x=x, y=d[,1], MarsFunc=funca, nCases=nCases,
                   nUsedTerms=nUsedTerms, penalty=penalty, RefFunc=func)
    if(plotit) {
        plotmo(fite, func=func, sub.caption=sub.caption)
        plot(fite, nresiduals=500, sub.caption=sub.caption)
    }
    cat("\n")
    fite
}

ozone.test <- function(itest, sModel, x, y, degree=2, nk=51,
                    plotit=plot.it.default, trace=0, col.loess="lightblue")
{
    fite <- earth(x, y, degree=degree, nk=nk, trace=trace)
    fitm <- mars(x, y, degree=degree, nk=nk)

    cat(sprintf("%-3d", itest),
        sprintf("%-32s", sModel), 
        "degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk),
        "nTerms",  sprintf("%-2d", sum(fite$selected.terms != 0)), 
        "of", sprintf("%-3d", sum(fite$selected.terms)),
        "GRSq", sprintf("%4.2g", fite$grsq),
        "GRSq ratio", fite$grsq/mars.to.earth(fitm)$grsq,
        "\n")
    sub.caption <- paste(itest, ": ", sModel, " degree=", degree, " nk=", nk, sep="")
    if(plotit) {
        fitme <- mars.to.earth(fitm)
        plotmo(fite, sub.caption=paste("EARTH", sub.caption))
        plotmo(fitme, sub.caption=paste("MARS", sub.caption))
        plot(fite, nresiduals=500, col.loess=col.loess, sub.caption=paste("EARTH", sub.caption))
        plot(fitme, sub.caption=paste("MARS", sub.caption))
        fitme <- update(fitme)  # generate model selection data
        plot.earth.models(list(fite, fitme), sub.caption=paste(itest, ": Compare earth to mars ", sModel, sep=""))
    }
    fite
}

funcNoise <- function(x)    # noise
{
    rnorm(1)
}
x <- cbind(              x1)
d <- cbind(funcNoise(x), x1)
# plotit=FALSE because there is only an intercept
itest <- itest+1; test.earth(itest, funcNoise, nk=5,  degree=1, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=5,  degree=2, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=51, degree=1, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=51, degree=2, plotit=FALSE, test.rsq=FALSE)

func1 <- function(x)
{
    sin(3 * x[,1]) + x[,2]
}
x <- cbind(          x1, x2)
d <- cbind(func1(x), x1, x2)
itest <- itest+1; test.earth(itest, func1, nk=5,   degree=1)
itest <- itest+1; test.earth(itest, func1, nk=5,   degree=2)
itest <- itest+1; test.earth(itest, func1, nk=51, degree=1)
itest <- itest+1; test.earth(itest, func1, nk=51, degree=2)

func7 <- function(x)    # just one predictor
{
    sin(5 * x[,1])
}
x <- cbind(x1)
d <- cbind(func7(x), x1)
itest <- itest+1; a<-test.earth(itest, func7, nk=5, degree=1)
itest <- itest+1; test.earth(itest, func7, nk=5,   degree=2)
itest <- itest+1; test.earth(itest, func7, nk=51, degree=1)
itest <- itest+1; test.earth(itest, func7, nk=51, degree=2)

func8 <- function(x)
{
    ret <- 0
    for (i in 1:5)
        ret <- ret + sin(2 * x[,i])
    ret + x[,1]*cos(4 * x[,2]) + (x[,3]-2)* x[,4]
}

func8noise <- function(x)
{
    func8(x) + rnorm(nrow(x),0,1)
}

x <- cbind(          x1,  x2,  x3,  x4,  x5)
d <- cbind(func8(x), x1,  x2,  x3,  x4,  x5)
itest <- itest+1; test.earth(itest, func8, nk=11, degree=1)
itest <- itest+1; test.earth(itest, func8, nk=11, degree=2)
itest <- itest+1; test.earth(itest, func8, nk=11, degree=10)
itest <- itest+1; test.earth(itest, func8, nk=51, degree=1)
itest <- itest+1; test.earth(itest, func8, nk=51, degree=2)
itest <- itest+1; test.earth(itest, func8, nk=51, degree=10)
itest <- itest+1; test.earth(itest, func8noise, nk=11, degree=1,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, func8noise, nk=11, degree=2,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, func8noise, nk=11, degree=10, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, func8noise, nk=51, degree=1,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, func8noise, nk=51, degree=2,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, func8noise, nk=51, degree=10, test.rsq=FALSE)

eqn56 <- function(x) # Friedman MARS paper equation 56
{
    0.1 * exp(4*x[,1]) +
    4 / (1 + exp(-20*(x[,2]-0.5))) +
    3 * x[,3] +
    2 * x[,4] +
    x[,5]
}

eqn56noise <- function(x)
{
    eqn56(x) + rnorm(nrow(x),0,1)
}

x <- cbind(          x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 )
d <- cbind(eqn56(x), x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 )
itest <- itest+1; test.earth(itest, eqn56, nk=11, degree=1)
itest <- itest+1; test.earth(itest, eqn56, nk=11, degree=2)
itest <- itest+1; test.earth(itest, eqn56, nk=11, degree=10)
itest <- itest+1; test.earth(itest, eqn56, nk=51, degree=1)
itest <- itest+1; test.earth(itest, eqn56, nk=51, degree=2)
itest <- itest+1; test.earth(itest, eqn56, nk=51, degree=10)
itest <- itest+1; test.earth(itest, eqn56noise, nk=11, degree=1,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, eqn56noise, nk=11, degree=2,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, eqn56noise, nk=11, degree=10, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, eqn56noise, nk=51, degree=1,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, eqn56noise, nk=51, degree=2,  test.rsq=FALSE)
itest <- itest+1; test.earth(itest, eqn56noise, nk=51, degree=10, test.rsq=FALSE)

robotArm <- function(x) # Friedman Fast MARS paper
{
    l1     <- x[,1]
    l2     <- x[,2]
    theta1 <- x[,3]
    theta2 <- x[,4]
    phi    <- x[,5]

    x1 <- l1 * cos(theta1) - l2 * cos(theta1 + theta2) * cos(phi)
    y <-  l1 * sin(theta1) - l2 * sin(theta1 + theta2) * cos(phi)
    z <-  l2 *  sin(theta2) * sin(phi)

    sqrt(x1^2 + y^2 + z^2)
}
N1 <- 100
set.seed(1)
x1. <- runif(N1, -1, 1)
x2. <- runif(N1, -1, 1)
x3. <- runif(N1, -1, 1)
x4. <- runif(N1, -1, 1)
x5. <- runif(N1, -1, 1)

x <- cbind(             (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
d <- cbind(robotArm(x), (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
colnames(x) <- c("l1", "l2", "theta1", "theta2", "phi")
colnames(d) <- c("arm", "l1", "l2", "theta1", "theta2", "phi")
itest <- itest+1; test.earth(itest, robotArm, nk=51, degree=1)
itest <- itest+1; test.earth(itest, robotArm, nk=51, degree=10)
itest <- itest+1; test.earth(itest, robotArm, nk=201, degree=1)
itest <- itest+1; test.earth(itest, robotArm, nk=201, degree=10)

# tests with ozone data

data(ozone1)
attach(ozone1)

x <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- itest+1; ozone.test(itest, "doy ~ wind+humidity+temp+vis", x, y, degree=1, nk=21)

x <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- itest+1; ozone.test(itest, "doy ~ wind+humidity+temp+vis", x, y, degree=2, nk=21)

# this is a basic test of RegressAndFix (because this generates lin dep bx cols)

cat("--Expect warning from mda::mars: NAs introduced by coercion\n") # why do we get a warning?
x <- cbind(wind, exp(humidity))
y <- doy
# col.loess is 0 else get loess errors
itest <- itest+1; ozone.test(itest, "doy ~ wind+exp(humidity)", x, y, degree=1, nk=21, col.loess=0)

x <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; ozone.test(itest, "O3~.", x, y, degree=2, nk=21)

x <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; a<-ozone.test(itest, "O3~., nk=51", x, y, degree=2, nk=51)

detach(ozone1)

options(options.old)

#--- subset ----------------------------------------------------------------------------------

set.seed(9)
train.subset <- sample(1:nrow(ozone1), .8 * nrow(ozone1))
test.subset <- (1:nrow(ozone1))[-train.subset]

# all the following models should be identical
a <- earth(ozone1[,-1], ozone1[,1], subset=train.subset, nprune=7, degree=2)
print(a)
a <- earth(ozone1[train.subset,-1], ozone1[train.subset,1], nprune=7, degree=2)
print(a)
a <- earth(O3 ~ ., data=ozone1, subset=train.subset, nprune=7, degree=2)
print(a)

y <- ozone1[test.subset, 1]
yhat <- predict(a, newdata = ozone1[test.subset, -1])
print(1 - sum((y - yhat)^2)/sum((y - mean(y))^2)) # print RSquared

#--- ppenalty and update----------------------------------------------------------------------

a <- earth(O3 ~ ., data=ozone1, degree=2)
print(update(a, ppenalty = -1))
print(update(a, ppenalty = 10))
a <- earth(O3 ~ ., data=ozone1, nk=31, pmethod="n", degree=2)
a.none <- print(update(a, nprune=10, pmethod="n"))
print(update(a.none, pmethod="b"))
print(update(a.none, nprune=4, pmethod="e"))
a.updated <- update(a.none, nprune=10, pmethod="b")
print(a.updated)
a.backwards <- update(a, nprune=10, pmethod="b")
print(a.backwards)
stopifnot(all.equal(a.updated$bx, a.backwards$bx))
a <- earth(O3 ~ ., data=ozone1, nk=31, nprune=10, pmethod="b", degree=2)
print(a)
stopifnot(all.equal(a$bx, a.backwards$bx))

#--- fast mars -------------------------------------------------------------------------------

print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = -1, fast.beta = 1))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = -1, fast.beta = 0))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 1))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 0))

# --- plot.earth and plot.earth.models--------------------------------------------------------

a <- earth(O3 ~ ., data=ozone1) # formula interface

plot(a, sub.caption="plot.earth test 1", col.rsq=3, col.loess=4, col.qq="pink", 
    col.vline=1, col.npreds=0, nresiduals=100, cum.grid="grid")

set.seed(1)
plot(a, sub.caption="plot.earth test 2", which=c(3,4,1), rlim=c(.2,.9), 
    jitter=.01, id.n=5, legend.pos=c(10,.4), pch=20, lty.vline=1)

plot(a, sub.caption="plot.earth test 3", which=c(2), main="test main")

a1 <- earth(ozone1[,c(2:4,10)], ozone1[,1])     # x,y interface

plot(a, sub.caption="plot.earth test 4", id.n=1)

set.seed(1)
plot.earth.models(a, which=1, rlim=c(.4,.8), jitter=.01)

plot.earth.models(a1)

plot.earth.models(list(a, a1), col.cum=c(3,4),  col.grsq=c(1,2), col.rsq=c(3,4),
    col.npreds=1, col.vline=1, lty.vline=3,
    col.grid="pink", cum.grid="grid", legend.pos=c(5,.4), 
    legend.text=c("a", "b", "c"))

# --- extractAIC.earth -----------------------------------------------------------------------

a <-earth(O3 ~ ., data=ozone1, degree=2)
cat("Ignore 10 warnings: extractAIC.earth: using GCV instead of AIC\n")
print(drop1(a))

#--- tests/test.earth.R ----------------------------------------------------------------------

cat("--- ../../tests/test.earth.R -------------------------\n")
source("../../tests/test.earth.R")


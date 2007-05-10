# test.earth.full.R: test earth

print(R.version.string)

library(earth)
library(mda)
data(ozone1)
data(trees)

print.time <- FALSE         # FALSE for no time results (for diff against reference)
plot.it.default <- TRUE     # TRUE to do plots too, FALSE for speed
options.old <- options()
options(digits = 3)

#--- test examples from man pages ------------------------------------------------------------

cat("--- earth.Rd -----------------------------\n")
example(earth)

a <- earth(mpg ~ ., data = mtcars, pmethod = "none", trace = 4)

set.seed(1)
train.subset <- sample(1:nrow(trees), .8 * nrow(trees))
test.subset <- (1:nrow(trees))[-train.subset]
a <- earth(Volume ~ ., data = trees[train.subset, ])
yhat <- predict(a, newdata = trees[test.subset, ])
y <- trees$Volume[test.subset]
print(1 - sum((y - yhat)^2) / sum((y - mean(y))^2)) # print R-Squared
cat("--- print.default of earth object---------\n")
print.default(a, digits=3)
cat("--- done print.default of earth object----\n")
plot(a)
library(mda)
(a <- fda(Species~., data=iris, method=earth, keepxy=TRUE))
plot(a)
print(summary(a$fit))
plot(a$fit)
plotmo(a$fit, ycolumn=1, ylim=c(-1.5,1.5), clip=FALSE)
plotmo(a$fit, ycolumn=2, ylim=c(-1.5,1.5), clip=FALSE)
a <- update(a, nk=3)	# not on man page
print(a)
print(summary(a$fit))

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
example(format.earth)
cat("--- get.nterms.per.degree.Rd ----------------------\n")
example(get.nterms.per.degree)
cat("--- get.nused.preds.per.subset.Rd ----------------------\n")
example(get.nused.preds.per.subset)
cat("--- mars.to.earth.Rd ----------------------\n")
example(mars.to.earth)
cat("--- plot.earth.models.Rd ----------------------\n")
example(plot.earth.models)
cat("--- plot.earth.Rd ----------------------\n")
example(plot.earth)
cat("--- predict.earth.Rd ----------------------\n")
example(predict.earth)
cat("--- reorder.earth.Rd ----------------------\n")
example(reorder.earth)
cat("--- update.earth.Rd ----------------------\n")
example(update.earth)

cat("--- test predict.earth -------------------\n")
a <- earth(Volume ~ ., data = trees)
predict(a, c(10,80))
predict(a)
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
predict(a, xpredict)
colnames(xpredict) <- c("Girth", "Height")
predict(a, xpredict)
predict(a, as.data.frame(xpredict))
# reverse dataframe columns (and their names), predict should deal with it correctly
xpredict <- as.data.frame(cbind(xpredict[,2], xpredict[,1]))
colnames(xpredict) <- c("Height", "Girth")
predict(a, xpredict)
# repeat but with x,y (not formula) call to earth
x1 <- cbind(trees$Girth, trees$Height)
colnames(x1) <- c("Girth", "Height")
a <- earth(x1, trees$Volume)
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
predict(a, xpredict)
colnames(xpredict) <- c("Girth", "Height")
predict(a, xpredict)
predict(a, as.data.frame(xpredict))
cat("--Expect warning from predict.earth: the variable names in 'data' do not match those in 'object'\n")
xpredict <- as.data.frame(cbind(xpredict[,2], xpredict[,1]))
colnames(xpredict) <- c("Height", "Girth")
predict(a, xpredict)

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
    gcv.null <- earth:::get.gcv(rss.null, 1, penalty, nCases)
    gcv <- earth:::get.gcv(rss, nUsedTerms, penalty, nCases)
    if(is.finite(object$rsq))
        if(!isTRUE(all.equal(object$rsq, 1 - rss/rss.null)))
            cat("\nWarning: RSq mismatch object$rsq", object$rsq, "calculated RSq", 1 - rss/rss.null)
        else if(!isTRUE(all.equal(object$grsq, 1 - gcv/gcv.null)))
            cat("\nWarning GRSq mismatch object$grsq", object$grsq, "calculated GRSq", 1 - gcv/gcv.null)
}

# this uses the global matrix data.global (data.global[,1] is the response)

test.earth <- function(itest, func, degree=2, nk=51, plotit=plot.it.default, test.rsq=TRUE, trace=0)
{
    cat("itest", sprintf("%-3d", itest), sprintf("%-32s", deparse(substitute(func))),
        "degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk))
    if(trace)
        cat("\n")
    gc()
    earthTime <- system.time(fite <- earth(data.global[,-1], data.global[,1],
                                        degree=degree, trace=trace, nk=nk, pmethod="b", fast.k=-1))
    funca <- make.func(fite)
    nCases <- nrow(data.global)
    penalty <- ifelse(degree>1,3,2)
    nUsedTerms <- sum(fite$selected.terms!=0)
    cat(" nTerms",  sprintf("%-2d", nUsedTerms), "of", sprintf("%-3d ", nrow(fite$dirs)))
    if(print.time)
        cat(" time", earthTime[1])
    cat("GRSq", sprintf("%4.2g", fite$grsq))
    caption <- paste("itest ", itest, ": ", deparse(substitute(func)),
                        " degree=", degree, " nk=", nk, sep="")
    if(test.rsq)
        test.model.rsq(fite, x=data.global[,-1, drop=FALSE], y=data.global[,1], MarsFunc=funca,
            nCases=nCases, nUsedTerms=nUsedTerms, penalty=penalty, RefFunc=func)
    if(plotit) {
        #$$ plotmo needs column names on x matrix wich are dropped if data.global[,-1] is a vector
        if(!is.vector(data.global[,-1]))
            plotmo(fite, func=func, caption=caption)
        plot(fite, nresiduals=500, caption=caption)
    }
    cat("\n")
    fite
}

ozone.test <- function(itest, sModel, x, y, degree=2, nk=51,
                    plotit=plot.it.default, trace=0, col.loess="lightblue")
{
    fite <- earth(x, y, degree=degree, nk=nk, trace=trace)
    fitm <- mars(x, y, degree=degree, nk=nk)

    cat("itest",
        sprintf("%-3d", itest),
        sprintf("%-32s", sModel),
        "degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk),
        "nTerms",  sprintf("%-2d", sum(fite$selected.terms != 0)),
        "of", sprintf("%-3d", nrow(fite$dirs)),
        "GRSq", sprintf("%4.2g", fite$grsq),
        "GRSq ratio", fite$grsq/mars.to.earth(fitm)$grsq,
        "\n")
    caption <- paste("itest ", itest, ": ", sModel, " degree=", degree, " nk=", nk, sep="")
    if(plotit) {
        fitme <- mars.to.earth(fitm)
        plotmo(fite, caption=paste("EARTH", caption))
        plotmo(fitme, caption=paste("MARS", caption))
        plot(fite, nresiduals=500, col.loess=col.loess, caption=paste("EARTH", caption))
        plot(fitme, caption=paste("MARS", caption))
        fitme <- update(fitme)  # generate model selection data
        plot.earth.models(list(fite, fitme), caption=paste(itest, ": Compare earth to mars ", sModel, sep=""))
    }
    fite
}

funcNoise <- function(x)    # noise
{
    rnorm(1)
}
x <- cbind(x1)
data.global <- cbind(funcNoise(x), x1)
# plotit=FALSE because there is only an intercept
itest <- itest+1; test.earth(itest, funcNoise, nk=5,  degree=1, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=5,  degree=2, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=51, degree=1, plotit=FALSE, test.rsq=FALSE)
itest <- itest+1; test.earth(itest, funcNoise, nk=51, degree=2, plotit=FALSE, test.rsq=FALSE)

func1 <- function(x)
{
    sin(3 * x[,1]) + x[,2]
}
x.global <- cbind(                    x1, x2)
data.global <- cbind(func1(x.global), x1, x2)
itest <- itest+1; test.earth(itest, func1, nk=5,   degree=1)
itest <- itest+1; test.earth(itest, func1, nk=5,   degree=2)
itest <- itest+1; test.earth(itest, func1, nk=51, degree=1)
itest <- itest+1; test.earth(itest, func1, nk=51, degree=2)

func7 <- function(x)    # just one predictor
{
    sin(5 * x[,1])
}
x.global <- cbind(                    x1)
data.global <- cbind(func7(x.global), x1)
itest <- itest+1; test.earth(itest, func7, nk=5,  degree=1)
itest <- itest+1; test.earth(itest, func7, nk=5,  degree=2)
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

x.global <- cbind(                    x1,  x2,  x3,  x4,  x5)
data.global <- cbind(func8(x.global), x1,  x2,  x3,  x4,  x5)
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

neg.eqn56noise <- function(x)
{
-eqn56noise(x)
}

x.global <- cbind(                    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 )
data.global <- cbind(eqn56(x.global), x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 )
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

x.global <- cbind(                       (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
data.global <- cbind(robotArm(x.global), (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
colnames(x.global) <- c("l1", "l2", "theta1", "theta2", "phi")
colnames(data.global) <- c("arm", "l1", "l2", "theta1", "theta2", "phi")
itest <- itest+1; test.earth(itest, robotArm, nk=51, degree=1)
itest <- itest+1; test.earth(itest, robotArm, nk=51, degree=10)
itest <- itest+1; test.earth(itest, robotArm, nk=201, degree=1)
itest <- itest+1; test.earth(itest, robotArm, nk=201, degree=10)

cat("--- tests with ozone data ----------------------\n")

data(ozone1)
attach(ozone1)

x.global <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- itest+1; ozone.test(itest, "doy ~ wind+humidity+temp+vis", x.global, y, degree=1, nk=21)

x.global <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- itest+1; ozone.test(itest, "doy ~ wind+humidity+temp+vis", x.global, y, degree=2, nk=21)

# this is a basic test of RegressAndFix (because this generates lin dep bx cols)

cat("--Expect warning from mda::mars: NAs introduced by coercion\n") # why do we get a warning?
x.global <- cbind(wind, exp(humidity))
y <- doy
# col.loess is 0 else get loess errors
itest <- itest+1; ozone.test(itest, "doy ~ wind+exp(humidity)", x.global, y, degree=1, nk=21, col.loess=0)

x.global <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; ozone.test(itest, "O3~.", x.global, y, degree=2, nk=21)

x.global <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; ozone.test(itest, "O3~., nk=51", x.global, y, degree=2, nk=51)

detach(ozone1)
options(options.old)

cat("--- fast mars -----------------------------------\n")

print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = -1, fast.beta = 1))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = -1, fast.beta = 0))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 1))
print(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 0))

cat("--- plot.earth and plot.earth.models ------------\n")

a <- earth(O3 ~ ., data=ozone1) # formula interface

plot(a, caption="plot.earth test 1", col.rsq=3, col.loess=4, col.qq="pink",
    col.vline=1, col.npreds=0, nresiduals=100, cum.grid="grid")

set.seed(1)
plot(a, caption="plot.earth test 2", which=c(3,4,1), rlim=c(.2,.9),
    jitter=.01, id.n=5, legend.pos=c(10,.4), pch=20, lty.vline=1)

plot(a, caption="plot.earth test 3", which=c(2), main="test main")

a1 <- earth(ozone1[,c(2:4,10)], ozone1[,1])     # x,y interface

plot(a, caption="plot.earth test 4", id.n=1)

set.seed(1)
plot.earth.models(a, which=1, rlim=c(.4,.8), jitter=.01)

plot.earth.models(a1)

plot.earth.models(list(a, a1), col.cum=c(3,4),  col.grsq=c(1,2), col.rsq=c(3,4),
    col.npreds=1, col.vline=1, lty.vline=3,
    legend.pos=c(5,.4), legend.text=c("a", "b", "c"))

cat("--- test multiple responses ---------------------\n")

# this uses the global matrix data.global (data.global[,1:2] is the response)

test.earth.two.responses <- function(itest, func1, func2,
    degree=2, nk=51, plotit=plot.it.default, test.rsq=TRUE, trace=0)
{
    if(typeof(func1) == "character")
        funcnames <- paste("multiple responses", func1, func2)
    else
        funcnames <- paste("multiple responses", deparse(substitute(func1)), deparse(substitute(func2)))
    cat("itest", sprintf("%-3d", itest), funcnames,
        " degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk), "\n\n")
    gc()
    fite <- earth(data.global[,c(-1,-2), drop=FALSE], data.global[,1:2],
                degree=degree, trace=trace, nk=nk, pmethod="b", fast.k=-1)
    print(fite)
    caption <- paste("itest ", itest, ": ", funcnames, " degree=", degree, " nk=", nk, sep="")
    if(plotit) {
        if(typeof(func1) == "character") {
            plotmo(fite, caption=caption)
            plotmo(fite, ycolumn=2)
        } else {
            plotmo(fite, func=func1, caption=caption)
            plotmo(fite, func=func2, ycolumn=2)
        }
        plot(fite, caption=caption)
        plot(fite, ycolumn=2)
    }
    cat("\n")
    fite
}

x.global <- cbind(                                     x1, x2)
data.global <- cbind(func1(x.global), func7(x.global), x1, x2)
itest <- itest+1; a <- test.earth.two.responses(itest, func1, func7, nk=51, degree=1)
print(summary(a))
plotmo(a, ycolumn=1)     # test generation of caption based on response name
plotmo(a, ycolumn=2)
plot(a, ycolumn=1)
plot(a, ycolumn=2)

x.global <- cbind(                                     x1, x2)
data.global <- cbind(func1(x.global), func7(x.global), x1, x2)
itest <- itest+1; a <- test.earth.two.responses(itest, func1, func7, nk=51, degree=3)
print(summary(a))

x.global <- cbind(                                           x1, x2, x3, x4, x5)
data.global <- cbind(eqn56(x.global), neg.eqn56noise(x.global), x1, x2, x3, x4, x5)
itest <- itest+1; a <- test.earth.two.responses(itest, eqn56, neg.eqn56noise, nk=51, degree=1)

x.global <- cbind(                                           x1, x2, x3, x4, x5)
data.global <- cbind(eqn56(x.global), neg.eqn56noise(x.global), x1, x2, x3, x4, x5)
itest <- itest+1; a <- test.earth.two.responses(itest, eqn56, neg.eqn56noise, nk=51, degree=2)

N1 <- 100
set.seed(1)
x1. <- runif(N1, -1, 1)
x2. <- runif(N1, -1, 1)
x3. <- runif(N1, -1, 1)
x4. <- runif(N1, -1, 1)
x5. <- runif(N1, -1, 1)

x.global <- cbind(                                        (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
data.global <- cbind(robotArm(x.global), eqn56(x.global), (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
colnames(x.global)    <- c(                "l1", "l2", "theta1", "theta2", "phi")
colnames(data.global) <- c("arm", "eqn56", "l1", "l2", "theta1", "theta2", "phi")
itest <- itest+1; test.earth.two.responses(itest, robotArm, eqn56, nk=51, degree=1)
itest <- itest+1; test.earth.two.responses(itest, robotArm, eqn56, nk=51, degree=10)
itest <- itest+1; test.earth.two.responses(itest, robotArm, eqn56, nk=201, degree=1)
itest <- itest+1; test.earth.two.responses(itest, robotArm, eqn56, nk=201, degree=2)
itest <- itest+1; test.earth.two.responses(itest, robotArm, eqn56, nk=201, degree=10)

attach(ozone1)
x.global <- cbind(                wind, humidity, temp, ibh, dpg, ibt, vis)
data.global <- cbind(O3, doy, vh, wind, humidity, temp, ibh, dpg, ibt, vis)
itest <- itest+1; test.earth.two.responses(itest, "O3", "doy", nk=51, degree=2)
detach(ozone1)

cat("--- formula based multiple response -------------\n")

a2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=2)
plotmo(a2)                  # test generation of caption based on response name
plotmo(a2, ycolumn=1)
plotmo(a2, ycolumn=2)
plot(a2)
plot(a2, ycolumn=1)
plot(a2, ycolumn=2)

cat("--- test plot.earth.models with multiple responses ---\n")

a <- earth(O3 ~ ., data=ozone1, degree=2)
a2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=2)
b2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=1)
plot.earth.models(list(a, a2), caption="plot.earth.models with multiple responses, list(a,a2)")
plot.earth.models(list(a2, a), caption="plot.earth.models with multiple responses, list(a2,a)")
plot.earth.models(list(a2, b2), caption="plot.earth.models with multiple responses, list(a2,b2)")

cat("--- subset --------------------------------------\n")

set.seed(9)
train.subset <- sample(1:nrow(ozone1), .8 * nrow(ozone1))
test.subset <- (1:nrow(ozone1))[-train.subset]

# all the following models should be identical
a <- earth(ozone1[,-1], ozone1[,1], subset=train.subset, nprune=7, degree=2)
print(a)
plotmo(a, caption="test subset: earth(ozone1[,-1], ozone1[,1], subset=train.subset)")

a <- earth(ozone1[train.subset,-1], ozone1[train.subset,1], nprune=7, degree=2)
print(a)
plotmo(a, caption="test subset: earth(ozone1[train.subset,-1], ozone1[train.subset,1]")

a <- earth(O3 ~ ., data=ozone1, subset=train.subset, nprune=7, degree=2)
print(a)
plotmo(a, caption="test subset: earth(O3 ~ ., data=ozone1, subset=train.subset")

y <- ozone1[test.subset, 1]
yhat <- predict(a, newdata = ozone1[test.subset, -1])
print(1 - sum((y - yhat)^2)/sum((y - mean(y))^2)) # print RSquared

cat("--- ppenalty and update -------------------------\n")

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
print(all.equal(a.updated$bx, a.backwards$bx))
a <- earth(O3 ~ ., data=ozone1, nk=31, nprune=10, pmethod="b", degree=2)
print(a)
print(all.equal(a$bx, a.backwards$bx))

cat("--- Force.xtx.prune -----------------------------\n")

a <- earth(Volume ~ ., data = trees)
a1 <- earth(Volume ~ ., data = trees, Force.xtx.prune=TRUE)
print(all.equal(a$bx, a1$bx))

a <- earth(O3 ~ ., data = ozone1, nk=51)
a1 <- earth(O3 ~ ., data = ozone1, nk=51, Force.xtx.prune=TRUE)
print(all.equal(a$bx, a1$bx))

cat("--- extractAIC.earth ----------------------------\n")

a <-earth(O3 ~ ., data=ozone1, degree=2)
cat("Ignore 10 warnings: extractAIC.earth: using GCV instead of AIC\n")
print(drop1(a))

cat("--- fda and mda with earth -----------------------------------\n")

am <- fda(Species ~ ., data=iris, method=mars, degree=1, keepxy=TRUE)
print(am)
a <- fda(Species ~ ., data=iris, method=earth, degree=1, keepxy=TRUE)
print(a)
print(confusion(a))
par(mar=c(3, 3, 2, .5))  # small margins and text to pack figs in
par(mgp=c(1.6, 0.6, 0))  # flatten axis elements
par(oma=c(0,0,4,0))      # make space for caption
layout(rbind(c(1,1,0,0), c(2,3,4,5), c(6,7,8,9)), heights=c(2,1,1))
plot(a)
plotmo(a$fit, ycolum=1, ylim=c(-1.5,1.5), clip=FALSE, do.par=FALSE)
plotmo(a$fit, ycolum=2, ylim=c(-1.5,1.5), clip=FALSE, do.par=FALSE)
mtext("fda test", outer=TRUE, font=2, line=1.5, cex=1)

data(glass)
set.seed(123)
samp <- sample(c(1:214), size=100, replace=FALSE)
glass.train <- glass[samp,]
glass.test <- glass[-samp,]
am <- mda(Type ~ ., data=glass.train, method=mars,  keepxy=TRUE, degree=2)
a <-  mda(Type ~ ., data=glass.train, method=earth, keepxy=TRUE, degree=2)
print(am)
print(a)
cat("mda with mars  ", attr(confusion(am), "error"), "\n")
cat("mda with earth ", attr(confusion(a),  "error"), "\n")
plotmo(a$fit, ycolum=9, clip=FALSE)

cat("--- ../../tests/test.earth.R -------------------------\n")

source("../../tests/test.earth.R")

if(!interactive())
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher)

# test.mods.R: test earth's ability to build various models

library(earth)

options(digits=3)
options(warn=2) # treat warnings as errors

PLOT <- FALSE
PRINT <- FALSE
PRINT.DATA <- FALSE
SUMMARY <- FALSE
TIME <- FALSE
TRACE <- 0
FORCE.WEIGHTS <- FALSE

itest <- 0

printf <- function(format, ...) cat(sprintf(format, ...)) # like c printf

ss <- function(x) # sum of squares
{
    sum(as.vector(x^2))
}
test.earth <- function(func, x, xtest, npreds, nk=NULL, degree=2, ...)
{
    itest <<- itest + 1
    set.seed(itest)
    x <- x[, 1:npreds, drop=FALSE]
    y <- func(x)
    nk <- if(is.null(nk)) min(200, max(20, 2 * ncol(x))) + 1 else nk
    printf("TEST %-2g n %g %-28s nk %-3g degree %-2g ",
        itest, nrow(x), deparse(substitute(func)), nk, degree)
    gc()
    earth.time <- system.time(mod <- earth(x, y, nk=nk, degree=degree,
                              trace=TRACE, Force.weights=FORCE.WEIGHTS, ...))

    ytest <- func(xtest)
    fitted <- predict(mod, xtest)
    stopifnot(length(fitted) == nrow(xtest))
    test.rsq <- 1 - ss(ytest - fitted) / ss(ytest - mean(ytest))

    time.string <- if(TIME) sprintf(" [time %5.3f]", earth.time[3]) else ""
    printf(" nterms %-3g%s grsq % 4.2f test.rsq % 4.2f\n",
        ncol(mod$bx), time.string, mod$grsq, test.rsq)
    if(PRINT) {
        print(mod)
        printf("\n")
    }
    if(SUMMARY) {
        print(summary(mod))
        printf("\n")
    }
    if(PRINT.DATA) {
        print(cbind(y, x))
        printf("\n")
    }
    if(PLOT) {
        nterms2 <- 0
        nterms <- earth:::get.nterms.per.degree(mod)
        if(length(nterms) >= 3)
            nterms2 <- nterms[3] # number of degree2 terms
        nrow <- 3
        par(mfrow=c(nrow, nrow), mar = c(1.5, 1.5, 3, .5), mgp = c(1.5, .5, 0), cex.main=1)
        oma <- par("oma") # space for caption
        oma[3] <- 3
        par(oma=oma)
        col.response <- "red" # ifelse(mod$leverages<.25, "green", "red")
        plotmo(mod, trace=-1, do.par=FALSE, col.response=col.response, npoints=500, cex.response=.8)
        caption <- sprintf("test %g n %d %s nk %g degree %g grsq %.2f test.rsq %.2f",
                    itest, nrow(x), deparse(substitute(func)), nk, degree, mod$grsq, test.rsq)
        mtext(caption, outer=TRUE, font=2, line=1, cex=1)
        plotmo(mod, trace=-1, do.par=FALSE, col.response=col.response, cex.response=.8,
               degree1=0, type2="image", npoints=500)
        col.response <- ifelse(mod$leverages<.25, "green", "red")
        plot(mod, versus=4, caption="", col.points=col.response, npoints=1000) # leverage plot
    }
    mod
}
testn <- function(n)
{
    itest <<- 0
    max.ncol <- 10
    set.seed(n)
    x <- matrix(runif(max.ncol * n, -1, 1), ncol=max.ncol)
    x <- x[order(x[,1]), , drop=FALSE] # sort first column for convenience
    colnames(x) <- paste("x", 1:ncol(x), sep="")

    xtest <- matrix(runif(max.ncol * 10000, -1, 1), ncol=max.ncol)
    xtest <- xtest[order(xtest[,1]), , drop=FALSE]
    colnames(xtest) <- c(paste("x", 1:max.ncol, sep=""))

    univariate <- function(x)
    {
        x[,1] + .3 * rnorm(nrow(x))
    }
    test.earth(univariate, x, xtest, 1, degree=1)

    bivariate <- function(x)
    {
        x[,1] + x[,2] + .3 * rnorm(nrow(x))
    }
    test.earth(bivariate, x, xtest, 2)

    bivariate.with.interaction <- function(x)
    {
        x[,1] + x[,2] + (x[,1] * x[,2]) + .3 * rnorm(nrow(x))
    }
    test.earth(bivariate.with.interaction, x, xtest, 2)

    trivariate <- function(x)
    {
        x[,1] + x[,2] + x[,3] + .3 * rnorm(nrow(x))
    }
    test.earth(trivariate, x, xtest, 3)

    sq <- function(x) x * x

    trivariate.with.interaction <- function(x)
    {
        x[,1] + x[,2] + sq(x[,3]) + (x[,1] * x[,2]) + .2 * rnorm(nrow(x))
    }
    test.earth(trivariate.with.interaction, x, xtest, 3)

    trivariate.two.interactions <- function(x)
    {
        x[,1] + x[,2] + sq(x[,3]) + (x[,1] * x[,2]) + sq(x[,1] * sq(x[,3])) + .1 * rnorm(nrow(x))
    }
    test.earth(trivariate.two.interactions, x, xtest, 3)
    printf("\n")

    sin.5.times.x1 <- function(x)
    {
        sin(5 * x[,1])
    }
    test.earth(sin.5.times.x1, x, xtest, 1, nk=5, degree=1)
    test.earth(sin.5.times.x1, x, xtest, 1, nk=5)
    test.earth(sin.5.times.x1, x, xtest, 1, degree=1)
    test.earth(sin.5.times.x1, x, xtest, 1)
    printf("\n")

    sin.3x1.times.x2 <- function(x)
    {
        sin(3 * x[,1]) + x[,2]
    }
    test.earth(sin.3x1.times.x2, x, xtest, 2, nk=5, degree=1)
    test.earth(sin.3x1.times.x2, x, xtest, 2, nk=5)
    test.earth(sin.3x1.times.x2, x, xtest, 2, degree=1)
    test.earth(sin.3x1.times.x2, x, xtest, 2)
    printf("\n")

    if(n < 100)
        return(invisible())

    pure.noise <- function(x)
    {
        rnorm(nrow(x))
    }
    test.earth(pure.noise, x, xtest, 1, nk=5,  degree=1)
    test.earth(pure.noise, x, xtest, 1, nk=5)
    test.earth(pure.noise, x, xtest, 1, nk=51, degree=1)
    test.earth(pure.noise, x, xtest, 1, nk=51)
    test.earth(pure.noise, x, xtest, 5, nk=5,  degree=1)
    test.earth(pure.noise, x, xtest, 5, nk=5)
    test.earth(pure.noise, x, xtest, 5, nk=51, degree=1)
    # TODO n=100, with weights test.rsq -1.34, without weights test.rsq 0.01
    test.earth(pure.noise, x, xtest, 5, nk=51)
    printf("\n")

    five.predictors <- function(x)
    {
        y <- 0
        for (i in 1:5)
            y <- y + sin(2 * x[,i])
        y + x[,1] * cos(4 * x[,2]) + (x[,3]-2)* x[,4]
    }
    test.earth(five.predictors, x, xtest, 5, degree=1)
    test.earth(five.predictors, x, xtest, 5)
    # TODO n=100, with weights test.rsq 0.69, without weights test.rsq 0.91
    test.earth(five.predictors, x, xtest, 5, degree=10)
    test.earth(five.predictors, x, xtest, 5, nk=51, degree=1)
    test.earth(five.predictors, x, xtest, 5, nk=51)
    test.earth(five.predictors, x, xtest, 5, nk=51, degree=3)
    printf("\n")

    # gives negative test.rsq
    #
    # five.predictors.noise <- function(x)
    # {
    #     five.predictors(x) + .01 * nrow(x)
    # }
    # test.earth(five.predictors.noise, x, xtest, 5, degree=1)
    # test.earth(five.predictors.noise, x, xtest, 5)
    # test.earth(five.predictors.noise, x, xtest, 5, degree=3)
    # test.earth(five.predictors.noise, x, xtest, 5, nk=51, degree=1)
    # test.earth(five.predictors.noise, x, xtest, 5, nk=51)
    # test.earth(five.predictors.noise, x, xtest, 5, nk=51, degree=3)
    # printf("\n")

    eqn56 <- function(x) # Friedman MARS paper equation 56
    {
        0.1 * exp(4 * x[,1]) +
        4 / (1 + exp(-20 * (x[,2] - 0.5))) +
        3 * x[,3] +
        2 * x[,4] +
        x[,5]
    }
    test.earth(eqn56, x, xtest, 5, degree=1)
    test.earth(eqn56, x, xtest, 5)
    test.earth(eqn56, x, xtest, 5, degree=3)
    test.earth(eqn56, x, xtest, 5, nk=51, degree=1)
    test.earth(eqn56, x, xtest, 5, nk=51)
    test.earth(eqn56, x, xtest, 5, nk=51, degree=3)
    printf("\n")

    eqn56.noise <- function(x)
    {
        eqn56(x) + 1 * rnorm(nrow(x))
    }
    test.earth(eqn56.noise, x, xtest, 5, degree=1)
    test.earth(eqn56.noise, x, xtest, 5)
    test.earth(eqn56.noise, x, xtest, 5, degree=3)
    test.earth(eqn56.noise, x, xtest, 5, nk=51, degree=1)
    test.earth(eqn56.noise, x, xtest, 5, nk=51)
    test.earth(eqn56.noise, x, xtest, 5, nk=51, degree=3)
    printf("\n")

    eqn56.noise.plus.extra.preds <- function(x)
    {
        eqn56(x) + 1 * rnorm(nrow(x))
    }
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10, degree=1)
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10)
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10, degree=3)
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10, nk=51, degree=1)
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10, nk=51)
    test.earth(eqn56.noise.plus.extra.preds, x, xtest, 10, nk=51, degree=3)
    printf("\n")

    # this has linear preds in both 1 and 2 degree terms
    test.earth(eqn56, x, xtest, 5, linpreds=c("^x1$","x3","5"))
    test.earth(eqn56, x, xtest, 5, linpreds=c(3,5))

    # check symmetry by using negative of eqn56 (may not be completely
    # symmetric because of fencepost odd/even alignment issues)
    neg.eqn56 <- function(x)
    {
        -eqn56(x)
    }
    test.earth(neg.eqn56, x, xtest, 5, linpreds=c(3,5))
    printf("\n")

    robot.arm <- function(x) # Friedman Fast MARS paper
    {
        l1     <- x[,1]
        l2     <- x[,2]
        theta1 <- x[,3]
        theta2 <- x[,4]
        phi    <- x[,5]

        x1 <- l1 * cos(theta1) - l2 * cos(theta1 + theta2) * cos(phi)
        y <-  l1 * sin(theta1) - l2 * sin(theta1 + theta2) * cos(phi)
        z <-  l2 * sin(theta2) * sin(phi)

        sqrt(x1^2 + y^2 + z^2)
    }

    x[,1] <- (x[,1] + 1) / 2    # l1  0..1
    x[,2] <- (x[,2] + 1) / 2    # l2  0..1
    x[,3] <- pi * (x[,3] + 1)   # theta1
    x[,4] <- pi * (x[,4] + 1)   # theta2
    x[,5] <- pi * x[,5] / 2     # phi
    colnames(x) <- c("l1", "l2", "theta1", "theta2", "phi", paste("x", 6:ncol(x), sep=""))

    xtest[,1] <- (xtest[,1] + 1) / 2    # l1  0..1
    xtest[,2] <- (xtest[,2] + 1) / 2    # l2  0..1
    xtest[,3] <- pi * (xtest[,3] + 1)   # theta1
    xtest[,4] <- pi * (xtest[,4] + 1)   # theta2
    xtest[,5] <- pi * xtest[,5] / 2     # phi
    colnames(xtest) <- c("l1", "l2", "theta1", "theta2", "phi", paste("x", 6:ncol(x), sep=""))

    test.earth(robot.arm, x, xtest, 5, nk=51,  degree=3)
    test.earth(robot.arm, x, xtest, 5, nk=201, degree=3)
    test.earth(robot.arm, x, xtest, 5, nk=201, degree=5)

    invisible()
}
if(PLOT) {
    printf("opening test.mods.pdf\n")
    pdf("test.mods.pdf")
}
start.time <- proc.time()
testn(30)
testn(100)
testn(500)
printf("[total time %.3f]\n", (proc.time() - start.time)[3])
if(PLOT)
    dev.off()
if(!interactive())
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs

# earth.varmod.R: build variance models for estimating prediction intervals
#
# TODO Extend the coverage table (print.inconf.tab) to show percentages in
# lower and upper intervals, so the user can check for asymmetry of the
# residuals.
#
# TODO Add QQ plot for prediction intervals, a "PIQ plot"
#
# TODO Consider making the code automatically detect non-monotonicity and
# issuing a warning.  Probably only possible for univariate models.
#
# TODO Could maybe prevent "Error in numericDeriv" by internally passing
# an explicit derivative in the call to nls.

VARMOD.METHODS <- c("const", "power", "power0",
                    "lm", "rlm", "earth", "gam",
                    "x.lm", "x.rlm", "x.earth", "x.gam")

TRACE.VARMOD         <- .3
TRACE.VARMOD.DETAILS <- .31 # will also cause plotting

# varmod returns a "varmod" object which can be used to
# estimate prediction intervals.
#
# y is the observed response (it is a n x 1 matrix)

varmod <- function(parent,
                   method, exponent, conv, clamp, minspan,
                   trace, x, y, meanfit, model.var, ...)
{
    UseMethod("varmod")
}
varmod.earth <- function(parent,
                      method, exponent, conv, clamp, minspan,
                      trace, x, y, meanfit, model.var, ...)
{
    check.classname(parent, deparse(substitute(parent)), "earth")
    varmod.internal(parent,
                    method, exponent, conv, clamp, minspan,
                    trace, x, y, meanfit, model.var, ...)
}
varmod.default <- function(parent,
                           method, exponent, conv, clamp, minspan,
                           trace, x, y, meanfit, model.var, ...)
{
    warning0("varmod.default: varmods are not supported for \"",
             class(parent)[1], "\" objects",
             "\nForging on anyway")

    varmod.internal(parent,
                    method, exponent, conv, clamp, minspan,
                    trace, x, y, meanfit, model.var, ...)
}
varmod.internal <- function(parent,
                            method, exponent=1, conv=1, clamp=.1, minspan=-5,
                            trace=0, parent.x=NULL, parent.y=NULL,
                            meanfit, model.var,
                            ...)
{
    # The following constant was an argument to earth but I removed it
    # and hardcoded it here for simplicity in the earth interface.
    # We use lambda to transform the squared residuals as follows:
    #     transformed.resids = squared.resids ^ (lambda / 2)
    # So with lambda=1, we transform to absolute residuals, and if
    # lambda=2, then there is no transform.  We call the transformed
    # residuals the abs.resids in the code (which is actually the correct
    # nomenclature only when lambda is 1).  See also get.resids.name.
    lambda <- 1

    trace <- check.trace.arg(trace)
    check.lambda.arg(lambda)
    check.exponent.arg(exponent, method)
    check.conv.arg(conv)
    check.clamp.arg(clamp)
    stopifnot(!is.null(parent.x))
    stopifnot(is.matrix(parent.x))
    stopifnot(!is.null(parent.y))
    stopifnot(is.matrix(parent.y))
    if(NCOL(parent.y) != 1)
        stop0("varmod.method cannot be used on multiple response models ",
              "(implementation restriction)")

    if(trace >= TRACE.VARMOD) {
        printf("\nvarmod method=\"%s\" exponent=%g lambda=%g conv=%g clamp=%g minspan=%g:\n",
            method, exponent, lambda, conv, clamp, minspan)
        if(trace == TRACE.VARMOD.DETAILS) {
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            par(mfrow=c(2, 3), mar=c(3, 3, 3, 1), mgp=c(1.5, 0.5, 0))
        }
    }
    n <- nrow(parent.x)
    df <- length(parent$selected.terms)
    squared.resids <- n / (n - df) * (parent.y - meanfit)^2 + model.var
    abs.resids <- squared.resids ^ (lambda / 2) # by default lambda=1, so take sqrt here

    residmod <- iterate.residmod(parent, abs.resids,
                                 method, exponent, lambda, conv, clamp, minspan,
                                 trace, parent.x, parent.y, ...)
    if(trace >= TRACE.VARMOD)
        printf("\n")

    varmod             <- NULL
    varmod$call        <- make.call.generic(match.call(expand.dots=TRUE), "varmod")
    varmod$parent      <- parent
    varmod$method      <- method
    varmod$exponent    <- exponent
    varmod$lambda      <- lambda
    varmod$residmod    <- residmod
    varmod$min.sd      <- get.min.sd(residmod, lambda, clamp)
    varmod$model.var   <- model.var
    varmod$meanfit     <- meanfit
    varmod$abs.resids  <- abs.resids # transformed residuals (actually on abs when lambda is 1)
    varmod$parent.x    <- parent.x
    varmod$parent.y    <- parent.y
    class(varmod)      <- "varmod"
    varmod$iter.rsq    <- get.varmod.weighted.rsq(varmod)
    varmod$iter.stderr <- get.varmod.stderr(varmod)
    varmod
}
get.varmod.weighted.rsq <- function(object, y) # return NULL if can't get rsq
{
    check.classname(object, deparse(substitute(object)), "varmod")
    if(object$method == "const")
        return(NULL)
    residmod  <- object$residmod
    residuals <- residmod$residuals
    fitted    <- residmod$fitted.values
    weights   <- weights(residmod)
    if(is.null(residuals) || is.null(fitted) || is.null(weights))
        return(NULL)
    mss <- sum(weights * (fitted - weighted.mean(fitted, weights))^2)
    rss <- sum(weights * residuals^2)
    mss / (mss + rss)
}
get.varmod.stderr <- function(object)
{
    check.classname(object, deparse(substitute(object)), "varmod")
    residmod <- object$residmod
    if(class(residmod)[1] %in% c("lm", "rlm", "nls")) {
        coef <- summary(residmod)$coefficients
        stopifnot(!is.null(coef[,"Std. Error"]))
        coef[,"Std. Error"]
    } else if(class(residmod)[1] == "gam" && "package:gam" %in% search()) {
        coef <- coefficients(summary.glm(residmod))
        stopifnot(!is.null(coef[,"Std. Error"]))
        coef[,"Std. Error"]
    } else if(class(residmod)[1] == "earth") {
        coef <- summary(lm(residmod$y ~ residmod$bx))$coefficients
        coef[,"Std. Error"]
    } else
        NULL
}
iterate.residmod <- function(parent, abs.resids,
                            method, exponent, lambda, conv, clamp, minspan,
                            trace, parent.x, parent.y, ...)
{
    varmod    <- NULL
    max.iter  <- 30
    weights   <- rep(1, nrow(parent.y))
    residmod  <- NULL
    trace.tab <- NULL # to trace, build a dataframe and print it when we're done

    for(iter in 1:max.iter) {
        residmod <- get.residmod(method, exponent, minspan, parent.x, parent.y,
                                 abs.resids, weights, trace, iter, parent, residmod, ...)

        coef.change <- get.coef.change(method, iter, residmod)

        # fill in enough of varmod for predict.varmod in get.residmod.weights
        varmod$parent   <- parent
        varmod$method   <- method
        varmod$exponent <- exponent
        varmod$lambda   <- lambda
        varmod$residmod <- residmod
        varmod$min.sd   <- get.min.sd(residmod, lambda, clamp)
        class(varmod)   <- "varmod"
        weights         <- get.residmod.weights(varmod, iter, trace)

        trace.tab <- trace.residmod(trace.tab, trace, iter, max.iter, residmod, varmod,
                                    coef, coef.change, parent.y, weights, exponent)

        if(possibly.break(iter, max.iter, coef.change, method, conv))
            break

    }
    if(iter == max.iter)
        warning0(sprintf(
            "varmod did not converge after %d iters (final coefchange%% %.2f)",
            iter, mean(abs(coef.change))))
    if(trace >= TRACE.VARMOD) {
        if(trace == TRACE.VARMOD.DETAILS && inherits(residmod, "nls"))
            printf("\n")
        print(trace.tab[1:iter,], row.names=FALSE, digits=2)
    }
    residmod
}
TRACE.NCOEF <- 0

blank.plot <- function(main=NULL)
{
    plot(0, 0, col=0, bty="n", xlab="", ylab="",
         xaxt="n", yaxt="n", main=main)
}
trace.residmod <- function(trace.tab, trace, iter, max.iter, residmod, varmod,
                           coef, coef.change, parent.y, weights, exponent)
{
    if(trace >= TRACE.VARMOD) {
        if(trace == TRACE.VARMOD.DETAILS) {
            if(inherits(residmod, "nls"))
                printf("\n")
            plotmo::plotmo(varmod, type="abs.residual",
                main="residmod first predictor",
                do.par=FALSE, degree1=1, degree2=0, trace=-1,
                col.response=1, col.degree1="lightblue", lwd.degree1=3)
        }
        coef <- coef(residmod)
        if(iter == 1) {
            unlockBinding("TRACE.NCOEF", asNamespace("earth"))
            TRACE.NCOEF <<- length(coef) # needed because nbr of earth coefs can change
            lockBinding("TRACE.NCOEF", asNamespace("earth"))
            trace.tab <- as.data.frame(matrix(NA, nrow=max.iter, ncol=3+TRACE.NCOEF))
            colnames(trace.tab) <- c("    iter", "weight.ratio", "coefchange%",
                     fix.coef.names(names(coef), colnames(parent.y), exponent))
        }
        trace.tab[iter,] <-
            c(iter, max(weights) / min(weights),
              mean(abs(coef.change)), c(coef, repl(NA, TRACE.NCOEF))[1:TRACE.NCOEF])
    }
    trace.tab
}
get.residmod <- function(method, exponent, minspan, parent.x, parent.y,
                         abs.resids, weights, trace, iter, parent, prev.residmod, ...)
{
    switch(method,
        const   = residmod.const(parent.x, parent.y, abs.resids,
                                 weights, trace),

        power   = residmod.power(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, parent, prev.residmod, iter, ...),

        power0  = residmod.power0(exponent, parent.x, parent.y, abs.resids,
                                  weights, trace, parent, prev.residmod, iter, ...),

        lm      = residmod.lm(exponent, parent.x, parent.y, abs.resids,
                              weights, trace, parent, ...),

        rlm     = residmod.rlm(exponent, parent.x, parent.y, abs.resids,
                               weights, trace, parent, ...),

        earth   = residmod.earth(exponent, minspan, parent.x, parent.y, abs.resids,
                                 weights, trace, parent, ...),

        gam     = residmod.gam(exponent, parent.x, parent.y, abs.resids,
                               weights, trace, parent, ...),

        x.lm    = residmod.x.lm(exponent, parent.x, parent.y, abs.resids,
                                weights, trace, ...),

        x.rlm   = residmod.x.rlm(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, ...),

        x.earth = residmod.x.earth(exponent, minspan, parent.x, parent.y, abs.resids,
                                   weights, trace, ...),

        x.gam   = residmod.x.gam(exponent, parent.x, parent.y, abs.resids,
                                 weights, trace, ...),

        stop0("illegal varmod.method \"", method, "\""))
}
PREV.COEF <- NULL

get.coef.change <- function(method, iter, residmod) # returns percentages for each coef
{
    # TODO use of abs means we think coef unchanged if sign but not abs value changes
    coef <- abs(coef(residmod))

    unlockBinding("PREV.COEF", asNamespace("earth"))
    if(iter == 1)
        PREV.COEF <<- coef  # note <<- not <-

    if(method %in% c("earth", "x.earth")) # see comments in possibly.break
        if(length(PREV.COEF) != length(coef))
            return(9999)

    denom <- PREV.COEF

    # Prevent divide by near zero.  The code below downweights extremely
    # small coefs and prevents them from dominating if rest are large.
    min.coef <- .01 * max(coef)
    denom[denom < min.coef] <- min.coef

    coef.change <- 100 * (coef - PREV.COEF) / denom

    PREV.COEF <<- coef      # note <<- not <-
    lockBinding("PREV.COEF", asNamespace("earth"))
    coef.change # a percentage for each coef
}
possibly.break <- function(iter, max.iter, coef.change, method, conv)
{
    if(conv < 0)
        iter == -conv
    else {
        method == "const" ||

        # TODO Since the earth basis funcs can change, looking at changes
        #      in the coefs can't be used to determine convergence.
        #      So for now, always do 1 iter.

        (method %in% c("earth", "x.earth")) ||

        # TODO following will sometimes create an intercept only model.
        # (method %in% c("earth", "x.earth") && iter >= 2) ||

        iter > 1 && mean(abs(coef.change)) < conv
    }
}
plot.residmod.weights <- function(w, main, min=NA, max=NA, median=NA) # for debugging
{
    plot(w, type="l", main=main, ylim=c(0, max(w, if(is.na(max)) 0 else max)))
    if(!is.na(min)) {
        abline(h=min, col=2)
        abline(h=max, col=2)
        abline(h=median, col=2)
    }
    else
        legend("topright", sprintf("max/min %.0f", max(w) / min(w)))
    lines(w) # replot over other annotations
}
# clamp to prevent extreme weights after squaring and inverse in get.residmod.weights

clamp.se <- function(se, iter, trace)
{
    median <- median(se)
    min <- median / 5   # 5 seemed ok with some limited simulation studies
    max <- 5 * median

    if(trace == TRACE.VARMOD.DETAILS)
        plot.residmod.weights(se, sprintf("iter %d: se", iter), min, max, median)

    se[se < min] <- min
    se[se > max] <- max

    se
}
# The variance for a regression on absolute residuals is proportional to
# the square of the regression model predicted value (Carrol and Ruppert
# book Section 3.3.3 and Table 3.3).

get.residmod.weights <- function(object, iter=0, trace=0)
{
    check.classname(object, deparse(substitute(object)), "varmod")

    # square to convert se to variance, inverse to convert variance to weight
    weights <- 1 / clamp.se(predict.varmod(object, type="se"), iter, trace)^2

    weights <- weights / mean(weights) # normalization not strictly necessary, may help numerics

    if(trace == TRACE.VARMOD.DETAILS)
        plot.residmod.weights(weights, "weights")

    weights
}
# we calculate LAMBDA.FACTOR only when necessary because the calculation
# can be slow,  hence we need the following global variables

LAMBDA <- LAMBDA.FACTOR <- -999

update.LAMBDA.FACTOR <- function(lambda, trace)
{
    approx.equal <- function(x, y)
    {
        # allow for limited precision in doubles, also allows .33 for 1/3
        abs(x - y) < 1e-2
    }
    #--- update.LAMBDA.FACTOR starts here ---
    if(lambda != LAMBDA) {
        unlockBinding("LAMBDA", asNamespace("earth"))
        unlockBinding("LAMBDA.FACTOR", asNamespace("earth"))

        LAMBDA <<- lambda   # note <<- not <-

        # some values have been precalculated
        if(approx.equal(lambda, 2))
            LAMBDA.FACTOR <<- 1
        # sqrt(pi / 2) = 1.2533, ratio mean dev to stddev, Geary 1935
        else if(approx.equal(lambda, 1))
            LAMBDA.FACTOR <<- sqrt(pi / 2)
        # (residuals^2)^(1/3) is approx normal by the Wilson-Hilferty
        # transform, although the left tail will still be short
        else if(approx.equal(lambda, 2/3))
            LAMBDA.FACTOR <<- 1.2464
        else
        {
            rnorm(1) # seems to be necessary to make .Random.seed available
            old.seed <- .Random.seed
            set.seed(1) # for reproducibility, 1e6 could be bigger but then slow
            LAMBDA.FACTOR <<- 1 / mean(rnorm(1e6)^2 ^ (lambda/2))
            set.seed(old.seed)
        }
        lockBinding("LAMBDA", asNamespace("earth"))
        lockBinding("LAMBDA.FACTOR", asNamespace("earth"))

        if(trace >= TRACE.VARMOD)
            printf("lambda %g LAMBDA.FACTOR %g\n", lambda, LAMBDA.FACTOR)
    }
}
# scale a prediction by the residmod back to a standard deviation

to.sd <- function(abs.resids, lambda, trace=0)
{
    update.LAMBDA.FACTOR(lambda, trace)
    # pmax is necessary to prevent e.g. sqrt of neg prediction from residmod
    (LAMBDA.FACTOR * pmax(abs.resids, 0)) ^ (1 / lambda)
}
get.min.sd <- function(residmod, lambda, clamp=.1)
{
    predict <- predict(residmod)
    predict <- predict[predict > 0]
    stopifnot(length(predict) > 0)
    stopifnot(clamp >= 0 && clamp <= 1)
    clamp * mean(to.sd(predict, lambda, 0))
}
check.lambda.arg <- function(lambda)
{
    check.numeric.scalar(lambda)
    if(lambda < 0.25 || lambda > 2)
        stop0("varmod.lambda must be between 0.25 and 2 ",
              "(you have lambda=", lambda, ")")
}
# TRUE if estimation of variance depends only on the fitted response (not on x)
method.uses.fitted.response <- function(method)
{
    method %in% c("power", "power0", "lm", "rlm", "earth", "gam")
}
check.exponent.arg <- function(exponent, method)
{
    check.numeric.scalar(exponent)
    # TODO following restriction could be lifted but currently only partially implemented
    if(exponent != 1 && !method.uses.fitted.response(method))
        stop0("varmod.exponent argument is not allowed with method=\"", method, "\"\n",
              "(varmod.exponent is only allowed for varmod.methods that depend only ",
              "on the fitted response)")
    if(exponent < .1 || exponent > 5)
        stop0("varmod.exponent must be between .1 and 5 ",
              "(you have varmod.exponent=", exponent, ")")
}
check.conv.arg <- function(conv)
{
    err <- function(conv)
        stop0("varmod.conv must be a negative integer or a percent between 0 and 100\n",
              "       (you have varmod.conv=", conv, ")")

    check.numeric.scalar(conv)
    if(conv < 0) {
        if(floor(conv) != conv) # conv is negative, check that it is an integer
            err(conv)
    } else if(conv == 0 || conv > 100)
        err(conv)
}
check.clamp.arg <- function(clamp)
{
    check.numeric.scalar(clamp)
    if(clamp < 0 || clamp > 1)
        stop0("varmod.clamp must be between 0 and 1, you have varmod.clamp=", clamp)
}
residmod.const <- function(parent.x, parent.y, abs.resids, weights, trace)
{
    # Predictions can be handled in a simple consistent way in
    # residmod.predict if instead of calculating the variance directly
    # here, we achieve the same result by building an intercept-only model
    # which always predicts mean(abs.resids).
    #
    # The conversion to a dataframe is necessary if the user later calls
    # plot(parent$varmod$residmod) or plotmo(parent$varmod$residmod).
    # Note that plotmo will call predict.varmod via predict.earth.

    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    lm(abs.resids~1, data=data, weights=weights, y=TRUE)
}
apply.exponent <- function(yhat, exponent, xname=colnames(yhat))
{
    stopifnot(!is.null(xname))
    check.vec(yhat, "yhat")
    # exponents of neg numbers are allowed only for integer exponents
    if(floor(exponent) != exponent) {
        check.that.most.are.positive(
            yhat, xname, sprintf("exponent=%g", exponent), "nonpositive")
        yhat[yhat < 0] <- 0 # don't want to take say sqrt of a neg number
    }
    yhat ^ exponent
}
nls.wrapper <- function(formula, data, start, weights, abs.resids, trace)
{
    # We use algorithm="port" below because the default algorithm more often causes
    # "Error in numericDeriv: Missing value or an infinity produced"
    # Also, on test data we sometimes need more iterations than the default 50

    mod <- nls(formula,
               data=data, start=start, weights=weights,
               trace=(trace == TRACE.VARMOD.DETAILS),
               algorithm="port", control=list(maxiter=100))

    # make model data available for plotmo and plotmor
    mod$y             <- abs.resids
    mod$fitted.values <- predict(mod)

    mod
}
estimate.power.start.values <- function(prev.residmod, abs.resids, data, weights, trace, iter)
{
    if(is.null(prev.residmod)) { # first iteration in iterate.residmod?
        # use a linear model to estimate the start values
        lm <- lm(abs.resids~., data=data, weights=weights)
        coefs <- coef(lm)
        if(trace == TRACE.VARMOD.DETAILS) {
            plotmo::plotmo(lm, col.response=2, do.par=F, trace=-1,
                main=sprintf("iter 1: lm for start vals\nvarmod.method=\"power\""))
            plot(lm, which=1)
            blank.plot()
        }
        start <- list(coefs[1], coefs[2], exponent=1)
        if(trace >= TRACE.VARMOD)
            printf(
                "\n     start: (Intercept)=%.2g    coef=%.2g    exponent=%.2g\n\n",
                start[[1]], start[[2]], start[[3]])
    } else { # not first iteration
        # use previous model values as starting values
        coefs <- coef(prev.residmod)
        stopifnot(length(coefs) == 3)
        start <- list(coefs[1], coefs[2], coefs[3])
    }
    names(start) <- c("(Intercept)", "coef", "exponent")
    if(trace == TRACE.VARMOD.DETAILS)
        cat(sprintf("iter %d  RSS:   ", iter), names(start), "\n")
    start
}
residmod.power <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, parent, prev.residmod, iter, ...)
{
    if(exponent != 1) # TODO allow this?
        stop0("the exponent argument cannot be used with varmod.method=\"power\"")
    parent.fit <- predict(parent)
    check.that.most.are.positive(
        parent.fit, "predict(parent)", "varmod.method=\"power\"", "nonpositive")
    parent.fit[parent.fit < 0] <- 0 # force negative values to zero
    data <- data.frame(abs.resids, apply.exponent(parent.fit, exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    start <- estimate.power.start.values(prev.residmod, abs.resids, data, weights, trace, iter)
    nls.wrapper(abs.resids~`(Intercept)` + coef * FITTED^exponent,
                data, start, weights, abs.resids, trace)
}
estimate.power0.start.values <- function(prev.residmod, abs.resids, data, weights, trace, iter)
{
    if(is.null(prev.residmod)) { # first iteration in iterate.residmod?
        # use a linear model to estimate the start values
        lm <- lm(abs.resids~.-1, data=data, weights=weights)
        coefs <- coef(lm)
        if(trace == TRACE.VARMOD.DETAILS) {
            plotmo::plotmo(lm, col.response=2, do.par=F, trace=-1,
                main=sprintf("iter 1: lm for start vals\nvarmod.method=\"power0\""))
            plot(lm, which=1)
            blank.plot()
        }
        start <- list(coefs[1], exponent=1)
        if(trace >= TRACE.VARMOD)
            printf(
                "\n     start: coef=%.2g    exponent=%.2g\n\n",
                start[[1]], start[[2]])
    } else { # not first iteration
        # use previous model values as starting values
        coefs <- coef(prev.residmod)
        stopifnot(length(coefs) == 2)
        start <- list(coefs[1], coefs[2])
    }
    names(start) <- c("coef", "exponent")
    if(trace == TRACE.VARMOD.DETAILS)
        cat(sprintf("iter %d  RSS:   ", iter), names(start), "\n")
    start
}
residmod.power0 <- function(exponent, parent.x, parent.y, abs.resids,
                            weights, trace, parent, prev.residmod, iter, ...)
{
    if(exponent != 1) # TODO allow this?
        stop0("the exponent argument cannot be used with varmod.method=\"power0\"")
    parent.fit <- predict(parent)
    check.that.most.are.positive(
        parent.fit, "predict(parent)", "varmod.method=\"power0\"", "nonpositive")
    parent.fit[parent.fit < 0] <- 0 # force negative values to zero
    data <- data.frame(abs.resids, apply.exponent(parent.fit, exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    start <- estimate.power0.start.values(prev.residmod, abs.resids, data, weights, trace, iter)
    nls.wrapper(abs.resids~coef * FITTED^exponent,
                data, start, weights, abs.resids, trace)
}
residmod.lm <- function(exponent, parent.x, parent.y, abs.resids,
                        weights, trace, parent, ...)
{
    # we use FITTED instead of colnames(parent.y) because we have applied exponent
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    lm(abs.resids~., data=data, weights=weights, y=TRUE)
}
residmod.rlm <- function(exponent, parent.x, parent.y, abs.resids,
                        weights, trace, parent, ...)
{
    # we use FITTED instead of colnames(parent.y) because we have applied exponent
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    library(MASS)
    mod <- MASS::rlm(abs.resids~., data=data, weights=weights, method="MM")
    # make model data available for plotmo and plotmor
    mod$y             <- abs.resids
    mod$fitted.values <- predict(mod)
    mod
}
residmod.earth <- function(exponent, minspan, parent.x, parent.y, abs.resids,
                           weights, trace, parent, ...)
{
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    earth(abs.resids~., data=data, weights=weights,
          keepxy=TRUE, trace=trace, minspan=minspan, ...)
}
gam.package.name <- function() # should we use the gam or the mgcv package?
{
    mgcv.package.loaded <- "package:mgcv" %in% search()

    if(mgcv.package.loaded && "package:gam" %in% search()) {
        # prevent downstream strange error messages
        stop0("cannot use varmod.method=\"gam\" when both the ",
               "\"gam\" and \"mgcv\" package are loaded")
    }
    if(mgcv.package.loaded)
        return("mgcv")

    "gam" # note that we prefer the gam package if neither is loaded
}
residmod.gam <- function(exponent, parent.x, parent.y, abs.resids,
                         weights, trace, parent, ...)
{
    data <- data.frame(abs.resids, apply.exponent(predict(parent), exponent))
    colnames(data) <- c("abs.resids", "FITTED")
    package.name <- gam.package.name()
    if(package.name == "gam") {
        if(trace >= 1) {
            printf("gam::gam(abs.resids~s(FITTED), data=data)\n")
            printf("x=TRUE, y=TRUE\n")
        }
        library(gam)
        residmod <- gam::gam(abs.resids~s(FITTED), data=data, weights=weights,
                           x=TRUE, y=TRUE)
    } else { # package.name "mgcv"
        if(trace >= 1)
            printf("mgcv::gam(abs.resids~s(FITTED), data=data)\n")
        library(mgcv)
        residmod <- mgcv::gam(abs.resids~s(FITTED), data=data, weights=weights)
        residmod$data <- data
    }
    residmod
}
residmod.x.lm <- function(exponent, parent.x, parent.y, abs.resids,
                          weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument cannot be used with varmod.method=\"x.lm\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    lm(abs.resids~., data=data, weights=weights, y=TRUE)
}
residmod.x.rlm <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument cannot be used with varmod.method=\"x.rlm\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    library(MASS)
    mod <- MASS::rlm(abs.resids~., data=data, weights=weights, method="MM", y.ret=TRUE)
    # make model data available for plotmo and plotmor
    mod$y             <- abs.resids
    mod$fitted.values <- predict(mod)
    mod
}
residmod.x.earth <- function(exponent, minspan, parent.x, parent.y, abs.resids,
                             weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument cannot be used with varmod.method=\"x.earth\"")
    data <- data.frame(abs.resids, parent.x)
    colnames(data) <- c("abs.resids", colnames(parent.x))
    earth(abs.resids~., data=data, weights=weights,
          keepxy=TRUE, trace=trace, minspan=minspan, ...)
}
residmod.x.gam <- function(exponent, parent.x, parent.y, abs.resids,
                           weights, trace, ...)
{
    if(exponent != 1)
        stop0("the exponent argument cannot be used with varmod.method=\"x.gam\"")
    if(ncol(parent.x) != 1)
        stop0("varmod.method=\"x.gam\" is not allowed when x has more than one column")

    RHS <- parent.x[,1] # needed so we can use s() in call to gam
    data <- data.frame(abs.resids=abs.resids, RHS=RHS)
    colnames(data) <- c("abs.resids", "RHS")

    package.name <- gam.package.name()
    if(package.name == "gam") {
        library(gam)
        if(trace >= 1) {
            printf("gam::gam(abs.resids~s(RHS), data=data, ")
            printf("x=TRUE, y=TRUE\n")
        }
        residmod <- gam::gam(abs.resids~s(RHS), data=data, weights=weights,
                           x=TRUE, y=TRUE)
    } else { # package.name "mgcv"
        library(mgcv)
        if(trace >= 1)
            printf("mgcv::gam(abs.resids~s(RHS), data=data)\n")
        residmod <- mgcv::gam(abs.resids~s(RHS), data=data, weights=weights)
        residmod$data <- data # for later access by plotmo etc.
    }
    residmod
}
get.quant <- function(level) # e.g for level=.95 return 1.96
{
    check.level.arg(level, zero.ok=FALSE)
    stopifnot(level > 0 && level < 1)
    level <- 1 - (1 - level) / 2 # .95 becomes .975
    qnorm(level)                 # .975 becomes 1.96
}
predict.abs.residual <- function(object, newdata, level)
{
    # unfortunately needed to get model formulas to work for some varmod methods
    hack.colnames <- function(newdata, method)
    {
        if(!is.null(newdata) &&
           method %in% c("power", "power0", "lm", "rlm", "earth", "gam", "x.gam")) {
            if(NCOL(newdata) != 1) {
                stop0("predict.varmod: NCOL(newdata) must be 1 ",
                      "when method=\"", method, "\" ",
                      "(implementation restriction)")
            }
            newdata <- as.data.frame(newdata)
            if(method == "x.gam")
                colnames(newdata) <- "RHS"
            else
                colnames(newdata) <- "FITTED"
        }
        newdata
    }
    if(is.null(newdata))
        abs.resid <- predict(object$residmod)
    else if(method.uses.fitted.response(object$method)) {
        parent.fit <- predict(object$parent, newdata=newdata)
        parent.fit <- apply.exponent(parent.fit, object$exponent)
        parent.fit <- data.frame(parent.fit)
        parent.fit <- hack.colnames(parent.fit, object$method)
        if(object$method %in% c("power", "power0"))
            parent.fit[parent.fit < 0] <- 0 # force negative values to zero
        abs.resid <- predict(object$residmod, newdata=parent.fit)
        stopifnot(length(abs.resid) == NROW(parent.fit))
    } else {
        newdata <- hack.colnames(newdata, object$method)
        abs.resid <- predict(object$residmod, newdata=newdata)
        stopifnot(length(abs.resid) == NROW(newdata))
    }
    abs.resid <- as.vector(abs.resid)

    # clamp at object$min.sd
    min.abs.resid <- (object$min.sd ^ object$lambda) / LAMBDA.FACTOR
    pmax(abs.resid, min.abs.resid)
}
predict.se <- function(object, newdata, level)
{
    to.sd(predict.abs.residual(object, newdata, level), object$lambda)
}
get.parent.fit <- function(object, newdata)
{
    parent.fit <- predict(object$parent, newdata=newdata)
    check.vec(parent.fit, "parent.fit")
    stopifnot(!is.null(dim(parent.fit))) # check parent.fit is a matrix or dataframe
    parent.fit[,1]
}
predict.pint <- function(object, newdata, level) # newdata allowed
{
    se <- predict.se(object, newdata)
    parent.fit <- get.parent.fit(object, newdata)
    stopifnot(length(parent.fit) == length(se))
    quant <- get.quant(level)
    data.frame(fit = parent.fit,
               lwr = parent.fit - quant * se,
               upr = parent.fit + quant * se)
}
predict.training.int <- function(object, se, newdata, level, int.name)
{
    if(!is.null(newdata))
        stop0("predict.varmod: newdata is not allowed with interval=\"",
              int.name, "\"")
    parent.fit <- get.parent.fit(object, newdata)
    stopifnot(length(se) == length(parent.fit))
    quant <- get.quant(level)
    data.frame(fit = parent.fit,
               lwr = object$meanfit - quant * se,
               upr = object$meanfit + quant * se)
}
predict.varmod <- function(
    object  = stop("no 'object' arg"),
    newdata = NULL,
    type    = c("pint", "training.pint", "training.cint", "se", "abs.residual"),
    level   = .95,
    trace   = FALSE, # currently unused
    ...)
{
    check.classname(object, deparse(substitute(object)), "varmod")
    warn.if.dots.used("predict.varmod", ...)
    switch(match.arg1(type),
        pint = {
            predict.pint(object, newdata, level)
        },
        training.pint = {
            se <- predict.se(object, newdata)
            predict.training.int(object, se, newdata, level, "training.pint")
        },
        training.cint = {
            predict.training.int(object,
                se=sqrt(object$model.var), newdata, level, "training.cint")
        },
        se = {
            if(level != .95)
                stop0("level argument is not allowed with this predict type")
            predict.se(object, newdata)
        },
        abs.residual = {
            if(level != .95)
                stop0("level argument is not allowed with this predict type")
            predict.abs.residual(object, newdata)
        })
}
# Example: if digits=3, then "%.*f" becomes "%.3f"
# Needed because R printf doesn't support * in printf formats
# and we need it to make the digits arg work in printfs

dot.star.to.digits <- function(s, digits)
{
    check.integer.scalar(digits, min=1)
    stopifnot(floor(digits) == digits)
    stopifnot(digits > 0 && digits < 20)
    gsub("\\.\\*", sprintf(".%d", digits), s)
}
# Example:
#       fix.coef.names(coef.names=h(y-123), response.name="y", exponent=.5)
#   returns
#       h(sqrt(y)-123) # the func knows that the special case of exponent=.5 is the sqrt

fix.coef.names <- function(coef.names, response.name, exponent)
{
    if(length(coef.names) == 1)
        return(coef.names) # do nothing if intercept only model
    stopifnot(length(response.name) == 1)
    stopifnot(exponent > 0)
    new.response.name <-
        if(exponent > .33 && exponent < .34)
            sprintf("cbrt(%s)", response.name)
        else if(exponent == .5)
            sprintf("sqrt(%s)", response.name)
        else if(exponent == 1)
            response.name
        else if(exponent == 2)
            sprintf("sq(%s)", response.name)
        else
            sprintf("%s^%.3g", response.name, exponent)
    coef.names <- gsub("FITTED", response.name, coef.names, fixed=TRUE)
    if(exponent == 1)
        coef.names
    else {
        # TODO revisit, will fail if response.name is substring of a token in
        #      coef.names or if response.name="h" and coef.names="h(h-12)"
        gsub(response.name, new.response.name, coef.names, fixed=TRUE)
    }
}
# restore original exponent, it doesn't get scaled like the other coefficients

restore.exponent <- function(coef, org.coef, method)
{
    if(method == "power")
        coef[3] <- org.coef[3]  # exponent is in coef[3]
    else if(method == "power0")
        coef[2] <- org.coef[2]  # exponent is in coef[2]
    coef
}
coef.varmod <- function(object, as.sd=TRUE, ...)
{
    warn.if.dots.used("coef.varmod", ...)
    coef <- coef(object$residmod)
    if(is.null(coef))
        stop0("coef.varmod: cannot get coefficients for \"",
               class(object$residmod)[1], "\" residmod")
    as.sd <- check.boolean(as.sd)
    if(as.sd) {
        org.coef <- coef
        negs <- coef < 0
        coef <- to.sd(abs(coef), object$lambda)
        coef[negs] <- -coef[negs]
        coef <- restore.exponent(coef, org.coef, object$method)
    }
    names(coef) <- fix.coef.names(names(coef),
                                  colnames(object$parent.y), object$exponent)
    coef
}
VARMOD.COEF.TAB.STYLES <- c("standard", "unit")

print.varmod.coef.tab <- function(
    object,
    style  = VARMOD.COEF.TAB.STYLES,
    digits = 2)
{
    style <- match.arg1(style)
    coef <- coef.varmod(object, as.sd=TRUE)

    # if style="unit", normalize coef if possible
    unit <- 1
    if(style == "unit") {
        # choose which coef will be the unit
        if(length(coef) == 1) # method == "const"?
            unit <- abs(coef[1])
        else
            unit <- abs(coef[2])
        if(unit < 1e-3) {
            warning0("coef=", unit, " is very small, forcing style=\"standard\"\n")
            style <- "standard"
            unit <- 1
        }
    }
    org.coef <- coef
    coef <- coef / unit
    coef <- restore.exponent(coef, org.coef, object$method)

    # get stderr
    NAs <- repl(NA, length(coef))
    stderr <- NAs
    if(!is.null(object$iter.stderr)) {
        org.stderr <- object$iter.stderr
        stderr <- to.sd(org.stderr, object$lambda) / unit
        stderr <- restore.exponent(stderr, org.stderr, object$method)
    }
    # get stderr.percent
    abs.coef <- abs(coef)
    stderr.percent <- 100 * stderr / abs.coef
    stderr.percent.as.char <- sprintf("%.0f", stderr.percent)
    stderr.percent.as.char[stderr.percent >= 1e3] <- "big"

    # create a data.frame and print that
    coef.names <- names(coef)
    coef.names <- gsub("`", "", coef.names) # remove backquotes added by lm etc.
    if(style == "unit") {
        coef.tab <- data.frame(
            zapsmall(c(coef, unit, object$min.sd / unit), digits+1),
            # sprintf below so print "NA" not "<NA>"
            c(sprintf("%g", zapsmall(stderr, digits)), " ", " "),
            c(stderr.percent.as.char, " ", " "))
        rownames(coef.tab) <- c(coef.names, "unit", "clamp")
    } else {
        coef.tab <- data.frame(
            zapsmall(coef, digits+1),
            c(sprintf("%g", zapsmall(stderr, digits+1))),
            c(stderr.percent.as.char))
        rownames(coef.tab) <- coef.names
    }
    if(object$method %in% c("earth", "x.earth")) {
        order <- reorder.earth(object$residmod, decomp="anova")
        order <- c(order, nrow(coef.tab)) # add "min.sd" row
        coef.tab <- coef.tab[order,]
    }
    colnames(coef.tab) <- c("coefficients", "iter.stderr", "iter.stderr%")
    print(coef.tab, digits=digits)
    coef.tab
}
print.interval.tab <- function(object, level, digits)
{
    level <- check.level.arg(level, zero.ok=FALSE)
    predict <- predict.varmod(object, type="pint", level=level)
    interval <- predict$upr - predict$lwr
    interval <- interval[order(interval)]
    tab <- data.frame(
        " ", mean(interval),
        " ", interval[1],
        " ", interval[length(interval)],
        " ", interval[length(interval)] / interval[1])
    colnames(tab) <- c(
        " ", "mean",
        " ", "smallest",
        " ", "largest",
        " ", "ratio")
    rownames(tab) <- sprintf("%g%% prediction interval", 100*level)
    print(tab, digits=digits)
    # remove spaces in returned table
    tab[, c("mean", "smallest", "largest", "ratio")]
}
get.model.response.from.formula <- function(formula, data, subset) # extract response from data
{
    call. <- match.call(expand.dots=FALSE)
    mf <- call.[c(1, match(c("formula", "data"), names(call.), 0))]
    mf[[1]] <- as.name("model.frame")
    y <- try(model.response(eval.parent(mf), "any"))
    if(is.try.error(y))
        stop0("cannot get the model response from newdata")
    if(is.null(y)) # probably not needed
        stop0("cannot get the model response from newdata")
    as.matrix(y)
}
percent.inconf <- function(object, level, parent.y, newdata)
{
    predict <- predict.varmod(object, newdata, type="pint", level=level)
    inconf <- parent.y >= predict$lwr & parent.y <= predict$upr
    100 * sum(inconf) / length(inconf)
}
print.inconf.tab <- function(object, parent.y, newdata)
{
    if(NCOL(parent.y) != 1) {
        warning0("multiple response model: the table is for the first response")
        parent.y <- parent.y[,1]
    }
    inconf68 <- percent.inconf(object, .68, parent.y, newdata)
    inconf80 <- percent.inconf(object, .80, parent.y, newdata)
    inconf90 <- percent.inconf(object, .90, parent.y, newdata)
    inconf95 <- percent.inconf(object, .95, parent.y, newdata)

    # .5 below adjusts for rounding in printf %.0f
    lt <- function(x, level) if(x < level-.5) "<" else " "

    tab <- data.frame(
        " ", sprintf("%.0f%%%s", inconf68, lt(inconf68, 68)),
        " ", sprintf("%.0f%%%s", inconf80, lt(inconf80, 80)),
        " ", sprintf("%.0f%%%s", inconf90, lt(inconf90, 90)),
        " ", sprintf("%.0f%%%s", inconf95, lt(inconf95, 95)))

    colnames(tab) <- c(
        " ", "68% ",
        " ", "80% ",
        " ", "90% ",
        " ", "95% ")

    if(is.null(newdata))
        rowname <- "response values in prediction interval"
    else
        rowname <- "newdata in prediction interval"
    rownames(tab) <- rowname
    print(tab)

    # return value is the table but not in string form
    tab <- data.frame(inconf68, inconf80, inconf90, inconf95)
    colnames(tab) <- c("68%", "80%", "90%", "95%")
    rownames(tab) <- rowname
    tab
}
print.varmod <- function(
    x       = stop("no 'x' arg"), # x is a varmod obj
    level   = .95,        # use 0 to not print the interval tabs
    style   = "standard", # one of VARMOD.COEF.TAB.STYLES
    digits  = 2,
    newdata = NULL,
    ...)
{
    check.classname(x, deparse(substitute(x)), "varmod")
    object <- x # minimize confusion with x, the regression input matrix
    remove(x)   # not necessary but prevents mistakes later
    warn.if.dots.used("print.varmod", ...)
    if(!is.null(newdata)) { # if newdata, print just the inconf table
        object$inconf.tab <-
            print.inconf.tab(object,
                             get.model.response(object$parent, newdata), newdata)
        return(invisible(object))
    }
    printf("method \"%s\"", object$method)
    space <- if(object$exponent != 1 || object$lambda != 1) "" else "  "
    if(object$exponent != 1)
        printf("%s  exponent %.3f", space, object$exponent)
    if(object$lambda != 1)
        printf("%s  lambda %g", space, object$lambda)
    printf("%s  min.sd %.3g", space, object$min.sd)
    if(!is.null(object$iter.rsq)) {
        printf("%s  iter.rsq %.3f", space, object$iter.rsq)
        # TODO prints too many digits
        # printf(dot.star.to.digits(", iter.rsq %.*f", digits+1), object$iter.rsq)
    }
    printf("\n\nstddev of predictions%s:\n",
           if(style == "unit") " (scaled by unit)" else "")
    object$coef.tab <- print.varmod.coef.tab(object, style, digits)
    level <- check.level.arg(level, zero.ok=TRUE)
    if(is.specified(level)) {
        printf("\n")
        object$interval.tab <- print.interval.tab(object, level, digits)
        printf("\n")
        object$inconf.tab <- print.inconf.tab(object, object$parent.y, newdata=NULL)
    }
    invisible(object)
}
print.summary.varmod <- function(
    x       = stop("no 'x' arg"), # x is a summary.varmod obj
    level   = x$level,
    style   = x$style,
    digits  = x$digits,
    newdata = x$newdata,
    ...)
{
    check.classname(x, deparse(substitute(x)), "varmod")
    warn.if.dots.used("print.summary.varmod", ...)
    if(is.null(level))
        level <- .95
    if(is.null(style))
        style <- "standard"
    if(is.null(digits))
        digits <- 2
    if(!is.null(newdata)) # if newdata, print just the inconf table
        print.varmod(x, level, style, digits, newdata)
    else {
        my.print.call("Parent model: ", x$parent$call)
        printf("\n")
        print.varmod(x, level, style, digits)
        printf("\nRegression submodel (%s):\n", get.resids.name(x))
        if(!(class(x$residmod)[1] %in% c("lm", "x.lm"))) # TODO check this
            printf("\n")
        print(x$residmod, digits=digits)
    }
    invisible(x)
}
summary.varmod <- function(
    object  = stop("no \"object\" arg"),
    level   = .95,
    style   = "standard", # one of VARMOD.COEF.TAB.STYLES
    digits  = 2,
    newdata = NULL,
    ...)
{
    check.classname(object, deparse(substitute(object)), "varmod")
    warn.if.dots.used("summary.varmod", ...)
    object$level   <- level   # pass level on to print.summary.varmod
    object$style   <- style   # ditto
    object$digits  <- digits  # ditto
    object$newdata <- newdata # ditto
    class(object)  <- c("summary.varmod", "varmod")
    object
}
get.resids.name <- function(object)
{
    if(object$lambda == 1)
        sprintf("Abs Residuals")
    else if(object$lambda == 2)
        sprintf("Squared Residuals")
    else
        sprintf("(Squared Residuals) ^ %.2g", object$lambda / 2)
}
get.varmod.ylab <- function(object, as.sd)
{
    sprintf("ParentMod %s", if(as.sd) "StdDev" else get.resids.name(object))
}
min.sd.line <- function(object, col.min.sd, lwd) # draw horizontal line at min.sd
{
    if(is.specified(col.min.sd)) {
        # TODO need to apply lambda exponent here?
        abline(h=object$min.sd / LAMBDA.FACTOR, col=col.min.sd, lty=2, lwd=lwd)
    }
}
sd.axis <- function(object) # draw righthand axis in standard deviation scale
{
    sd <- LAMBDA.FACTOR * object$abs.resids ^ (1 / object$lambda) # for righthand axis
    pretty.sd <- pretty(range(sd))
    axis(side=4, at=pretty.sd / LAMBDA.FACTOR, labels=pretty.sd, srt=90)
    mtext(get.varmod.ylab(object, as.sd=TRUE), side=4, line=1.5, cex=par("cex"))
}
plot.varmod <- function(
    x          = stop("no 'x' arg"), # a varmod object
    which      = 1:4,
    do.par     = length(which) > 1,
    info       = FALSE,
    cex        = NULL,
    caption    = if(do.par) NULL else "", # NULL for auto
    main       = NULL, # NULL means auto
    col.line   = "lightblue",
    col.min.sd = col.line,
    ...)    # unused, for compat with the generic
{
    check.classname(x, deparse(substitute(x)), "varmod")
    object.name <- substr(deparse(substitute(x)), 1, 40)
    object <- x # minimize confusion with x, the regression input matrix
    remove(x)   # needed else get.plotmo.x gets this x instead of the x matrix
    warn.if.dots.used("plot.varmod", ...)
    info <- check.boolean(info)
    plotmo::check.index(which, "which", 1:4)
    stopifnot(length(do.par) == 1)
    stopifnot(do.par == 0 || do.par == 1 || do.par == 2)
    old.par <- par.for.plot(do.par, length(which), 1, "CAPTION1\nCAPTION2",
                            right.axis=TRUE)
    if(do.par == 1)
        on.exit(par(old.par))
    else if(!is.null(cex)) {
        old.cex <- par("cex")
        on.exit(par(cex=old.cex))
        par(cex=cex)
    }
    if(is.null(cex))
        cex <- get.cex.points(-1, length(object$parent.y))
    if(!is.null(main))
        main <- repl(main, 4) # recycle for up to 4 plots
    ylim <- fix.lim(c(min(object$abs.resids, 0), max(object$abs.resids)))
    parent.fitted <- predict(object$parent)[,1]
    order <- order(parent.fitted)
    col.smooth <- if(info) 2 else 0
    # lwd <- ifelse(length(parent.fitted) > 1000, 3, 2) # npoints matches plotmor
    lwd <- 3
    for(iwhich in seq_along(which)) {
        if(which[iwhich] == 1) {            #--- fitted versus parent fitted ---
            plot(parent.fitted[order], object$abs.resids[order],
                 main=if(is.null(main))
                        sprintf("%s vs Fitted", get.resids.name(object))
                      else
                        main[iwhich],
                 ylim=ylim, pch=20, cex=cex, xlab="Fitted",
                 ylab=get.varmod.ylab(object, as.sd=FALSE))
            min.sd.line(object, col.min.sd, lwd) # horizontal line at min.sd
            sd.axis(object)                # right hand axis in stddev scale
            # fitted values of residual model
            fitted <- predict.varmod(object, type="abs.residual")
            lines(parent.fitted[order], fitted[order], col=col.line, lwd=lwd)
            if(info) {
                # lowess smooth
                smooth <- lowess(parent.fitted[order], object$abs.resids[order], f=.5)
                lines(smooth$x, smooth$y, col=col.smooth, lwd=1)
            }

        } else if(which[iwhich] == 2) {     #--- fitted versus parent first predictor ---
            plotmo::plotmo(object, type="abs.residual",
                ylim=ylim, degree1=1, degree2=0, do.par=FALSE, trace=-1,
                col.response=1, col.degree1=col.line, cex.response=cex,
                lwd.degree1=lwd, col.smooth=col.smooth,
                ylab=get.varmod.ylab(object, as.sd=FALSE),
                main=if(is.null(main))
                        sprintf("%s vs First Predictor", get.resids.name(object))
                     else
                        main[iwhich])
            min.sd.line(object, col.min.sd, lwd) # horizontal line at min.sd
            sd.axis(object)                 # right hand axis in stddev scale

        } else if(which[iwhich] == 3) {     #--- residual plot ---
            plotmor(object$residmod, which=3, npoints=-1, cex.points=cex,
                    xlab=get.varmod.ylab(object, as.sd=FALSE),
                    ylab="VarMod Residuals", info=info, center=FALSE)

        } else if(which[iwhich] == 4) {     #--- model selection graph ---
            if(class(object$residmod)[1] == "earth")
                plot.earth(object$residmod, which=1,
                           main=if(is.null(main)) "VarMod Model Selection" else main[iwhich])
        } else
            stop0("plot.varmod: illegal value %g in \"which\" argument", which[iwhich])
    }
    if(is.null(caption) || any(nchar(caption))) { # show caption?
        trim <- FALSE
        if(is.null(caption)) { # auto caption?
            trim <- TRUE
            parent.call <-
                strip.white.space(paste0(deparse(object$parent$call), collapse=""))
            caption <- sprintf("Variance Model  method=\"%s\"\nParentMod: %s",
                               object$method, parent.call)
        }
        show.caption(caption, trim=trim)
    }
    invisible()
}

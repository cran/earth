# plotmor.R: plot model residuals
#
# "plotres" is already used in several R packages, hence the name "plotmor"

plotmor <- function(
    x           = stop("no 'x' arg"),
    which       = 2:4,
    info        = FALSE,
    delever     = FALSE,  # divide by sqrt(1 - h_ii)
    pearson     = FALSE,  # divide by estimated stderr
    level       = 0,
    versus      = NULL,

    nresponse   = 1,
    npoints     = 1000,
    id.n        = 3,
    labels.id   = rownames(residuals(object, warn=FALSE)),
    center      = TRUE,
    loess.f     = .5,

    env         = NULL,
    do.par      = length(which) > 1,
    xlim        = NULL,
    ylim        = NULL,
    main        = NULL,
    cex.main    = 1,
    caption     = if(do.par) NULL else "",
    xlab        = NULL,
    ylab        = NULL,

    pch         = 20,
    col.line    = "lightblue",
    col.loess   = NULL,
    lwd.loess   = 2,
    col.cv      = "red",
    col.qq      = col.line,
    col.grid    = "lightgray",
    col.points  = 1,
    cex.points  = NULL,
    shade.pints = "mistyrose2",
    shade.cints = "mistyrose4",

    cum.grid    = "percentages",

    ...)
{
    object.name <- substr(deparse(substitute(x)), 1, 40)
    object <- x  # minimize confusion with x, the regression input matrix
    remove(x)    # not necessary but prevents mistakes later
    versus.spec <- versus
    stop.if.dots.used("plotmor", ...)
    plotmo::check.index(which, "which", 1:8)
    stopifnot(length(do.par) == 1)
    stopifnot(do.par == 0 || do.par == 1 || do.par == 2)

    rinfo <- get.rinfo(object, nresponse, delever, pearson)

    fitted <- get.plotmor.fitted(object, nresponse, rinfo)

    if(is.null(cex.points))
        cex.points <- get.cex.points(npoints, length(rinfo$resids))
    if(is.null(labels.id))
        labels.id <- paste(1:length(rinfo$resids))

    check.integer.scalar(npoints, min=-1, null.ok=TRUE)
    check.integer.scalar(id.n, min=0, null.ok=TRUE)
    info <- check.boolean(info)

    if(is.null(versus.spec)) {
        if(do.par && length(which) > 1) {
            old.par <- par.for.plot(do.par, length(which), cex.main, caption)
            if(do.par == 1)
                on.exit(par(old.par))
        }
        plotmor1(object, fitted, fitted, rinfo,
            which, info, delever, pearson, level,
            nresponse, npoints, id.n, labels.id, center, loess.f,
            env, xlim, ylim, main, cex.main, caption,
            xlab=if(is.null(xlab)) "Fitted" else xlab, ylab,
            pch, col.line, col.loess, lwd.loess, col.cv, col.qq, col.grid, col.points,
            cex.points, shade.pints, shade.cints,
            cum.grid)
    } else { # versus arg was specified
        # remove all from which except standard resid and abs resid plots
        which <- which[which %in% c(3,5)]
        if(length(which) == 0) {
            warning0("plot.earth: nothing to plot\n",
                     "(the \"which\" argument is empty after removing certain ",
                     "plots because \"versus\" was specified)")
            return(invisible())
        }
        # get the versus matrix (either earth's x or bx matrix)
        if(is.character(versus.spec) && length(versus.spec) &&
           substr(versus.spec, 1, 1) == "*") {
            # versus.spec begins with "*", so use the earth basis matrix
            versus <- object$bx[,-1, drop=FALSE]     # drop the intercept
            versus.spec <- substring(versus.spec, 2) # drop the "*"
            versus.spec <- plotmo::check.index(versus.spec, "versus",
                                1:NCOL(versus), colnames=colnames(versus))
        } else {
            if(is.null(env))
                env <- get.model.env(object, parent.frame)
            versus <- plotmo::get.plotmo.x.wrapper(object, env, trace=0)
            if(is.null(ncol(versus)))
                versus <- matrix(versus, nrow=length(versus))
            if(is.null(colnames(versus)))
                colnames(versus) <- paste("x", 1:ncol(versus), sep="")
            versus.spec <- plotmo::check.index(versus.spec, "versus",
                                1:NCOL(versus), colnames=colnames(versus))
        }
        nfigs <- length(which) * length(versus.spec)
        if(do.par && nfigs > 1) {
            old.par <- par.for.plot(do.par, nfigs, cex.main, caption)
            if(do.par == 1)
                on.exit(par(old.par))
        }
        for(iplot in seq_along(versus.spec))
            plotmor1(object, versus[, versus.spec[iplot]], fitted, rinfo,
                which, info, delever, pearson, level,
                nresponse, npoints, id.n, labels.id, center, loess.f,
                env, xlim, ylim, main, cex.main, caption,
                xlab=colnames(versus)[versus.spec[iplot]], ylab,
                pch, col.line, col.loess, lwd.loess, col.cv, col.qq, col.grid, col.points,
                cex.points, shade.pints, shade.cints,
                cum.grid)
    }
    if(do.par && (is.null(caption) || any(nchar(caption)))) { # show caption?
        trim <- is.null(caption) # trim if auto-generate the caption
        show.caption(get.caption(which, object, nresponse, object.name, caption))
    }
    invisible()
}
plotmor1 <- function(
    object,
    versus,
    fitted,
    rinfo,

    which,
    info,
    delever,
    pearson,
    level,

    nresponse,
    npoints,
    id.n,
    labels.id,
    center,
    loess.f,

    env, # unused
    xlim,
    ylim,
    main,
    cex.main,
    caption,
    xlab,
    ylab,

    pch,
    col.line,
    col.loess,
    lwd.loess,
    col.cv,
    col.qq,
    col.grid,
    col.points,
    cex.points,
    shade.pints,
    shade.cints,

    cum.grid)
{
    stopifnot(NCOL(versus) == 1)
    iresids <- get.iresids(npoints, rinfo$resids, versus, id.n)
    for(iwhich in seq_along(which)) {
        if(which[iwhich] == 2) {            #--- cumulative distribution ---
            if(length(which) !=1 || is.null(xlim))
                xlim <- range(abs(rinfo$scale * rinfo$resids))
            plot.cum(rinfo, main,
                     xlim=xlim, col=1, col.grid, cum.grid, add=FALSE, jitter=0)

        } else if(which[iwhich] == 4)       #--- QQ plot ---
            plot.qq(rinfo, iresids, main, info,
                    col.points, cex.points, col.qq, pch, id.n, labels.id)

        else {                               #--- various residual plots ---
            plot.resids(object, versus, fitted, rinfo, iresids, id.n,
                which[iwhich], info, level, labels.id, center, loess.f,
                if(length(which)!=1) NULL else xlim,
                if(length(which)!=1) NULL else ylim,
                main, cex.main, xlab, ylab,
                pch, col.line, col.loess, lwd.loess,
                col.cv, col.grid, col.points, cex.points,
                shade.pints, shade.cints)
        }
    }
}
# Get the environment in which the model function was originally called.
# If that is not available, use the environment in which plotmo was called.

get.model.env <- function(object, parent.frame, trace=0)
{
    .Environment <- attr(object$terms, ".Environment")
    if(is.null(.Environment)) {
        env <- parent.frame
        if(trace >= 1)
            printf("Using env parent.frame()\n")
    } else {
        env <- .Environment
        if(trace >= 1)
            printf("Using env attr(object$terms, \".Environment\")\n")
    }
    env
}
# return the response.name (for prepending to caption), if appropriate

get.caption.prefix <- function(which, object, nresponse, caption)
{
    resids <- residuals(object, warn=FALSE)
    if(is.null(residuals))
        stop0("get.caption.prefix: residuals(object) returned NULL")
    if(is.null(ncol(resids)))
        resids <- matrix(resids, ncol=1)
    nresponse <- plotmo::check.index(nresponse, "nresponse", resids, is.col.index=TRUE)
    if((all(which == 1) && NCOL(resids) > 1))
        return(NULL)    # Model Selection graph only, no multi-resp prefix
    colnames <- colnames(resids)
    if(!is.null(colnames) && !is.null(colnames[nresponse]) &&
            !is.na(colnames[nresponse]) && colnames[nresponse] != "")
        colnames[nresponse]
    else if(NCOL(resids) > 1)
        paste("Response", nresponse)
    else
        NULL
}
get.caption <- function(which, object, nresponse=NA, object.name, caption)
{
    if(!is.null(caption))
        return(caption)    # don't modify caption explictly set by user
    caption <- NULL
    if(!all(is.na(nresponse)))
        caption.prefix <- get.caption.prefix(which, object, nresponse, caption)
    if(!is.null(object$ifold)) # object is one fold of a cross-validated model?
        caption <- object.name
    else
        caption <- get.caption.from.call(caption, object)
    if(!is.null(caption.prefix))
        caption <- paste0(caption.prefix, ": ", caption)
    caption
}
plot.resids <- function(
    object,

    versus,
    fitted,
    rinfo,
    iresids,
    id.n,

    which,
    info,
    level,
    labels.id,
    center,
    loess.f,

    xlim,
    ylim,
    main,
    cex.main,
    xlab,
    ylab,

    pch,
    col.line,
    col.loess,
    lwd.loess,
    col.cv,
    col.grid,
    col.points,
    cex.points,
    shade.pints,
    shade.cints)
{
    stopifnot(length(which) == 1)
    info <- check.boolean(info)

    ok <- which %in% c(3,5:8)
    if(!all(ok))
        stop0("which=", which[!ok][1], " is not allowed")

    id.indices <- NULL
    if(which %in% c(3, 4:7))
        id.indices <- get.id.indices(rinfo$scale * rinfo$resids, id.n)

    level <- check.level.arg(level, zero.ok=TRUE)
    if(which %in% (5:8))
        level <- 0 # no pints

    varmod <- NULL
    if(which == 3 && !is.null(object$varmod) && is.specified(level) &&
            (is.specified(shade.pints) || is.specified(shade.cints)))
        varmod <- object$varmod

    pints <- NULL
    cints <- NULL
    if(!is.null(varmod) && is.specified(level)) {
        if(is.specified(shade.pints))
            pints <- predict(object, interval="training.pint", level=level)
        if(is.specified(shade.cints))
            cints <- predict(object, interval="training.cint", level=level)
    }
    resids <- rinfo$scale * rinfo$resids

    if((which %in% (6:8)))
        check.that.most.are.positive(
            versus, "fitted", sprintf("which=%d", which), "nonpositive")

    trans.versus <- trans.versus(versus[iresids], which)

    if(which == 8)
        check.that.most.are.positive(
            abs(resids), "abs(residuals)", sprintf("which=%d", which), "negative")

    trans.resids <- trans.resids(resids[iresids], which)

    xlab <- get.resids.xlab(xlab, which, versus)
    ylab <- get.resids.ylab(ylab, which, rinfo$name)
    main <- get.resids.main(main, xlab, ylab, varmod, level) # title of plot

    col.points <- repl(col.points, length(resids)) # recycle
    cex.points <- repl(cex.points, length(resids))
    pch        <- repl(pch, length(resids))

    ylim <- get.resids.ylim(ylim, fitted, trans.resids, which, info,
                            id.indices, center, pints, cints, rinfo$scale)

    xlim <- get.resids.xlim(xlim, which, trans.versus, ylim) # needs ylim

    plot(trans.versus, trans.resids,
         col=col.points[iresids], cex=cex.points[iresids],
         pch=pch[iresids],
         main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)

    if(!is.null(varmod)) {
        if(is.specified(shade.pints))
            plot.pint.resids(pints, versus, fitted, rinfo$scale, shade.pints)
        if(is.specified(shade.cints))
            plot.pint.resids(cints, versus, fitted, rinfo$scale, shade.cints)
    }
    if(which != 8)
        abline(h=0, lty=1, col=col.grid) # axis

    # replot points (because they are obscured by abline and shaded areas)
    points(trans.versus, trans.resids,
           col=col.points[iresids], cex=cex.points[iresids],
           pch=pch[iresids])

    if(is.null(col.loess)) {
        # must distinguish loess line from regression line
        col.loess <- if(info && which %in% c(5,8)) "red" else col.line
    }
    check.numeric.scalar(lwd.loess, null.ok=TRUE)
    if(is.null(lwd.loess))
        lwd.loess <- if(!is.null(varmod) || length(trans.resids) >= 1000) 3 else 2

    if(which != 8)
        plot.loess(trans.versus, trans.resids, loess.f, col.loess, lwd.loess)

    if(!is.null(varmod) && is.specified(col.cv)) {
        # plot resids of meanfit with col.cv (default red)
        meanfit <- varmod$meanfit - fitted
        meanfit <- rinfo$scale * meanfit
        order <- order(versus)
        lines(trans.versus(versus[order], which),
              trans.resids(meanfit[order], which), col=col.cv)
    }
    if(!is.null(id.indices))
        thigmophobe.labels(
            x=trans.versus(versus, which)[id.indices],
            y=trans.resids(resids, which)[id.indices],
            labels=labels.id[id.indices],
            offset=.33, font=2, cex=.8, col="steelblue4", xpd=NA)

    if(info) {
        col.lm <- col.line
        lwd.lm <- lwd.loess+1
        add.resids.info(which, info, versus, resids, col.lm, lwd.lm)
    }
    else if(!is.null(varmod)) {
        if(is.specified(col.loess) || is.specified(col.cv)) {
            # add legend, else blue and red lines may confuse the user
            # no legend if info=TRUE because not enough space
            legend.txt <- NULL
            legend.col <- NULL
            legend.lwd <- NULL
            legend.lty <- NULL
            if(which != 8 && is.specified(col.loess)) { # loess plotted?
                legend.txt <- "loess"
                legend.col <- col.loess
                legend.lwd <- lwd.loess
                legend.lty <- 1
            }
            if(is.specified(col.cv)) {
                legend.txt <- c(legend.txt, "cross validated oof fit")
                legend.col <- c(legend.col, col.cv)
                legend.lwd <- c(legend.lwd, 1)
                legend.lty <- c(legend.lty, 1)
            }
            legend("bottomleft", legend=legend.txt,
                   col=legend.col, lwd=legend.lwd, lty=legend.lty, bg="white", cex=.8)
        }
    }
}
check.level.arg <- function(level, zero.ok)
{
    if(identical(level, NA) || is.null(level)) # treat NA and NULL as 0
        level <- 0
    check.numeric.scalar(level)
    if(!((zero.ok && level == 0) || level >= 0.5 || level < 1)) {
        msg <- if(zero.ok)
                    "level argument must be zero or between 0.5 and 1"
                else
                    "level argument must be between 0.5 and 1"
        stop0(msg, " (you have level=", level, ")")
    }
    level
}
get.resids.xlim <- function(xlim, which, trans.versus, ylim)
{
    if(is.null(xlim) && which==8) {
        # don't show lower 5% of points
        quant <- quantile(trans.versus, prob=c(.05, 1), na.rm=TRUE)
        min <- quant[1]
        max <- quant[2]
        # extra left margin so slope of linear versus not flattened
        if(min > .2 * ylim[1])
            min <- .2 * ylim[1]
        xlim <- c(min, max)
    } else
        xlim <- fix.lim(range(trans.versus, na.rm=TRUE))
    xlim
}
get.resids.ylim <- function(ylim, fitted, resids, which, info,
                            id.indices, center, pints, cints, scale)
{
    if(!is.null(ylim) && !is.na(ylim)) # user specified a ylim?
        ylim
    else {
        if(!is.null(pints)) {
            pints <- get.pint.resids(pints, fitted, scale)
            min <- min(resids, pints$lwr, na.rm=TRUE)
            max <- max(resids, pints$upr, na.rm=TRUE)
        } else if(!is.null(cints)) {
            cints <- get.pint.resids(cints, fitted, scale)
            min <- min(resids, cints$lwr, na.rm=TRUE)
            max <- max(resids, cints$upr, na.rm=TRUE)
        } else {
            min <- min(resids, na.rm=TRUE)
            max <- max(resids, na.rm=TRUE)
        }
        if(!is.null(id.indices)) { # space for point labels
            range <- abs(max - min)
            min <- min - .03 * (max - min)
            max <- max + .03 * (max - min)
        }
        if(which %in% (5:7))
            min <- 0
        else if(which == 8)
            max <- max + .5  # more space on top, looks better
        else if(which == 3 && center) {
            # want symmetric ylim so can more easily see asymmetry
            if(abs(min) > abs(max))
                max <- -min
            else if(abs(max) > abs(min))
                min <- -max
        }
        if(info) { # space for info text (on top) and density plot (in the bottom)
            max <- max + max * if(!is.null(id.indices)) .2 else .1
            min <- min - max * if(!is.null(id.indices)) .2 else .1
        }
        ylim <- c(min, max)
    }
    stopifnot(is.numeric(ylim))
    stopifnot(length(ylim) == 2)
    fix.lim(ylim)
}
get.pint.resids <- function(pints, fitted, scale)
{
    stopifnot(nrow(pints)   == nrow(fitted))
    stopifnot(length(scale) == nrow(fitted))
    scale * (pints - fitted)
}
plot.pint.resids <- function(pints, versus, fitted, scale, shade.pints)
{
    if(!is.null(pints)) {
        pints <- get.pint.resids(pints, fitted, scale)
        # abscissa must be ordered for polygon to work
        order <- order(versus)
        versus <- versus[order]
        pints <- pints[order,]
        polygon(trans.versus(c(versus, rev(versus)), 0),
                trans.resids(c(pints$lwr, rev(pints$upr)), 0),
                col=shade.pints, border=shade.pints)
    }
}
plot.loess <- function(versus, resids, loess.f, col.loess, lwd.loess)
{
    if(is.specified(col.loess)) {
        check.numeric.scalar(loess.f)
        stopifnot(loess.f > .01 && loess.f < 1)
        # We use lowess rather than loess because loess tends to give ugly warnings,
        # but the argument is named "col.loess" for backward compatibility.
        # na.rm is needed if we take logs of nonpos, see check.that.most.are.positive.
        # That's why we calculate delta explicitly instead of using lowess default.
        delta <- .01 * diff(range(versus, na.rm=TRUE))
        smooth <- lowess(versus, resids, f=loess.f, delta=delta)
        lines(smooth$x, smooth$y, col=col.loess, lwd=lwd.loess)
    }
}
get.resids.xlab <- function(xlab, which, versus)
{
    if(is.null(xlab)) {
        stopifnot(!is.null(colnames(versus)))
        xlab <- colnames(versus)[1]
    }
    if(which %in% (6:8))
        xlab <- sprintf("Log %s", xlab)
    xlab
}
get.resids.ylab <- function(ylab, which, rinfo.name)
{
    if(is.null(ylab))
        ylab <- sprintf("%ss", rinfo.name)
    if(which == 5 || which == 6)
        ylab <- sprintf("Abs %s", ylab)
    else if(which == 7)
        ylab <- sprintf("Cube Root Squared %s", ylab)
    else if(which == 8)
        ylab <- sprintf("Log Abs %s", ylab)
    ylab
}
get.resids.main <- function(main, xlab, ylab, varmod, level) # title of plot
{
    # TODO should really use strwidth for newline calculation
    # The "Fitted" is for back compat, helps with limitations of formula below
    newline <- xlab != "Fitted" && nchar(ylab) + nchar(xlab) > 15
    if(is.null(main))
        main <- sprintf("%s vs%s%s", ylab, if(newline) "\n" else " ", xlab)
    if(!is.null(varmod) && !newline) # two newlines is too many
        main <- sprintf("%s\n%g%% level shaded", main, 100*(level))
    main
}
plot.density.along.the.bottom <- function(x)
{
    den <- density(x, adjust=.5, na.rm=TRUE)
    usr <- par("usr") # xlim1, xlim2, ylim1, ylim2
    den$y <- den$y * .08 * (usr[4] - usr[3]) / (max(den$y) - min(den$y))
    den$y <- den$y + usr[3]
    den.col <- "#909090"
    lines(den, col=den.col)
}
add.resids.info <- function(which, info, versus, resids, col.lm, lwd.lm)
{
    plot.density.along.the.bottom(trans.versus(versus, which))

    lm.text <- ""
    usr <- par("usr") # xlim1, xlim2, ylim1, ylim2
    if(which == 5 || which == 8) { # add linear regression line?
        trans.versus <- trans.versus(versus, which)
        trans.resids <- trans.resids(resids, which)

        # ignore lower 10% of points (very small residuals i.e. very negative logs)
        quant <- quantile(trans.versus, prob=.1, na.rm=TRUE)
        ok <- which(trans.versus > quant)
        trimmed.trans.versus  <- trans.versus[ok]
        trimmed.trans.resids <- trans.resids[ok]

        library(MASS) # for rlm

        # linear regression on 10 bootstrap samples so we can see variance of versus
        # for(i in 1:10) {
        #   j <- sample.int(length(trimmed.trans.versus), replace=TRUE)
        #   trimmed.trans.resids1 <- trimmed.trans.resids[j]
        #   trimmed.trans.fit1 <- trimmed.trans.versus[j]
        #   rlm <- rlm(trimmed.trans.resids1~trimmed.trans.fit1,
        #              method="MM", na.action="na.omit")
        #   abline(coef(rlm), col="lightgray", lwd=.6)
        # }

        # robust linear regression of trans.resids on trans.versus
        # na.omit is needed if some versus (or resids) were nonpos so log(versus) is NA

        rlm <- rlm(trimmed.trans.resids~trimmed.trans.versus,
                   method="MM", na.action="na.omit")
        abline(coef(rlm), col=col.lm, lwd=lwd.lm)
        lm.text <- sprintf("     slope %.2g", coef(rlm)[2])
    }
    # exact=FALSE else get warning "Cannot compute exact p-value with ties"
    cor <- cor.test(versus, abs(resids), method="spearman", exact=FALSE)
    text <- sprintf("spearman  %.2f%s", cor$estimate, lm.text)
    text.on.white(usr[1] + .6 * strwidth(text),
                  usr[4] - 1.2 * strheight(text), text,
                  cex=1, font=2, col="steelblue4", xpd=NA)
}
get.se <- function(object)
{
    if(inherits(object, "earth")) {
        if(is.null(object$varmod))
            stop0("\"pearson\" is not allowed because\n",
                  "the model was not built with varmod.method")
        predict(object, type="earth", interval="se")
    } else if(inherits(object, "lm")) {
        predict <- predict(object, se.versus=TRUE)
        predict$residual.scale # TODO correct?
    } else
       stop0("\"pearson\" is not allowed because\n",
             "get.se doesn't know how to get the stderrs of a \"",
             class(object)[1], "\" object")
}
# get information on the residuals

get.rinfo <- function(object, nresponse, delever, pearson)
{
    resids <- residuals(object, warn=FALSE)
    if(is.null(resids))
        stop0("cannot get residuals for \"", class(object), "\" object")
    if(NCOL(resids) == 1)
        resids <- matrix(resids, ncol=1)
    nresponse <- plotmo::check.index(nresponse, "nresponse", resids, is.col.index=TRUE)
    resids <- resids[, nresponse, drop=FALSE]
    scale <- repl(1, length(resids))
    name <- "Residual"

    delever <- check.boolean(delever)
    if(delever) {
        hat <- hatvalues(object)
        stopifnot(all(hat > 0))
        scale <- scale / sqrt(1 - hat)
        name <- sprintf("Delevered %s", name)
    }
    pearson <- check.boolean(pearson)
    if(pearson) {
        se <- get.se(object)
        stopifnot(length(se) == length(resids))
        scale <- scale / se
        name <- sprintf("Pearson %s", name)
    }
    check.vec(scale, "scale", length(resids))
    check(scale, "scale", "non-positive value", function(x) { x <= 0 })
    list(resids = resids, # raw resids (delever and pearson not applied)
         scale  = scale,  # will be 1 unless delever or pearson set
         name   = name)   # "Residual" or "Delevered Residuals" etc.
}
get.plotmor.fitted <- function(object, nresponse, rinfo)
{
    fitted <- fitted(object)
    if(is.null(fitted))
        stop0("plotmor: object$fitted is NULL")
    if(is.null(ncol(fitted)))
        fitted <- matrix(fitted, ncol=1)
    nresponse <- plotmo::check.index(nresponse, "nresponse", fitted, is.col.index=TRUE)
    fitted <- fitted[, nresponse, drop=FALSE] # TODO correct?
    if(is.null(colnames(fitted)))
        colnames(fitted) <- "Fitted"
    stopifnot(nrow(fitted) == length(rinfo$resids))
    fitted
}
my.log10 <- function(x)
{
    i <- which(x < max(x) / 1e6)
    x[i] <- 1
    x <- log10(x)
    x[i] <- NA
    x
}
trans.versus <- function(versus, which)
{
    if(which %in% (6:8))
        my.log10(versus)
    else
        versus
}
trans.resids <- function(resid, which) # transform the residuals
{
    if(which == 5 || which == 6)
        abs(resid)
    else if(which == 7) {
        # do it in two steps so no problem with negative numbers
        resid <- resid^2
        resid^(1/3)
    } else if(which == 8)
        my.log10(abs(resid))
    else
        resid
}
plot.cum <- function(rinfo, main,
                     xlim, col, col.grid, cum.grid, add, jitter)
{
    show.percents <- function(q)
    {
        is.space.available <- function(i)
        {
            q[i] - q[i-1] > 1.2 * strwidth && q[i+1] - q[i] > 1.2 * strwidth
        }
        show.percent <- function(i, label, ...)
        {
            text.on.white(x=q[i], y=1.02, label, cex1, ...)
        }
        #--- show.percents starts here ---
        cex1 <- .6
        strwidth <- strwidth("25%", cex=cex1)
        show.percent(1, "0%", xpd=NA) # xpd=NA to allow text out of plot region
        if(is.space.available(2)) show.percent(2, "25%")
        show.percent(3, "50%")
        if(is.space.available(4)) show.percent(4, "75%")
        show.percent(5, "90%")
        if(is.space.available(6)) show.percent(6, "95%")
        show.percent(7, "100%", xpd=NA)
    }
    #--- plot.cum starts here ---
    trans.resids <- abs(rinfo$scale * rinfo$resids)
    cum <- ecdf(trans.resids)
    if(jitter > 0.05)
       stop0("\"jitter\" ", jitter , " is too big, try something like jitter=0.01")
    if(jitter > 0)
        environment(cum)$"y" <- jitter(environment(cum)$"y", amount=2 * jitter)
    # col.points=0 gives a finer resolution graph (points are quite big regardless of pch)
    xlim <- fix.lim(xlim)
    xlab <- rinfo$name
    xlab <- sprintf("abs(%ss)", xlab)
    plot.stepfun(cum, add=add,
        main=if(is.null(main)) "Cumulative Distribution" else main,
        xlab=xlab, ylab="Proportion",
        xlim=xlim, col.points=0, col.hor=col, col.vert=col)
    if(is.specified(col.grid) && !add) {
        cum.grid <- match.choices(cum.grid[1],
                                  c("none", "grid", "percentages"))
        if(cum.grid %in% c("grid", "percentages")) {
            # add annotated grid lines, unattractive but useful
            for(h in c(0,.25,.5,.75,.90,.95,1)) # horizontal lines
                abline(h=h, lty=1, col=col.grid)
            q <- quantile(trans.resids, probs=c(0, .25, .50, .75, .9, .95, 1))
            for(v in q)    # vertical lines at 0,25,50,75,90,95,100% quantiles
                abline(v=v, lty=1, col=col.grid)
            if(cum.grid == "percentages")
                show.percents(q)
            plot.stepfun(cum, add=TRUE, verticals=TRUE, # replot data over grid
                xlim=xlim, col.points=0, col.hor=col, col.vert=col)
        }
    }
}
plot.qq <- function(rinfo, iresids, main, info,
                    col.points, cex.points, col.qq, pch, id.n, labels.id)
{
    par(pty="s") # square
    # we figure out the shape of the qq line with all resids but
    # plot only npoints points (selecting them with iresids)
    resids <- rinfo$scale * rinfo$resids
    qq <- qqnorm(resids, main=main, plot.it=FALSE)
    xlim <- ylim <- NULL
    id.indices <- get.id.indices(resids, id.n)
    if(!is.null(id.indices)) {
        # figure out xlim and ylim with extra space for point labels
        min <- min(qq$x)
        max <- max(qq$x)
        xlim <- c(min - .1 * (max - min), max + .1 * (max - min))
        min <- min(qq$y)
        max <- max(qq$y)
        ylim <- c(min - .05 * (max - min), max + .05 * (max - min))
        if(info)
            ylim[1] <- ylim[1] - .1 * (max - min) # space for density plot
    }
    col.points <- repl(col.points, length(resids)) # recycle
    cex.points <- repl(cex.points, length(resids))
    pch        <- repl(pch, length(resids))
    ylab <- rinfo$name
    if(is.null(main))
        main <- sprintf("%s QQ", ylab)
    ylab <- sprintf("%s Quantiles", ylab)
    plot(qq$x[iresids], qq$y[iresids],
         col=col.points[iresids], cex=cex.points[iresids],
         pch=pch[iresids],
         main=main, xlab="Theoretical Quantiles", ylab=ylab,
         xlim=xlim, ylim=ylim)
    if(info)
        plot.density.along.the.bottom(qq$x) # plot density of x along the bottom
    qqline(resids, col=col.qq)
    if(!is.null(id.indices))
        thigmophobe.labels(x=qq$x[id.indices], y=qq$y[id.indices],
             labels=labels.id[id.indices],
             offset=.33, font=2, cex=.8, col="steelblue4", xpd=NA)
    par(pty="m") # back to maximal (the default)
}
# get the subset of residuals we're going to display

get.iresids <- function(npoints, resids, versus, id.n)
{
    check.vec(resids, "residuals", length(versus))
    iresids <- 1:nrow(resids)
    if(npoints > 0 && npoints < nrow(resids)) {
        # take a sample, but make sure it includes the largest resids
        npoints <- min(npoints, nrow(resids))
        nlargest <- min(max(8, id.n), npoints)
        isorted <- order(abs(resids), decreasing=TRUE)
        ikeep <- 1:nlargest
        if(npoints > nlargest)
            ikeep <- c(ikeep, seq(from=nlargest + 1, to=nrow(resids),
                        length.out=npoints - nlargest))
        iresids <- isorted[ikeep]
    }
    # index vector into resids and versus, sorted on versus
    iresids[order(versus[iresids])]
}
# get the indices of the id.n largest resids
# requires a sort so is quite slow

get.id.indices <- function(resids, id.n)
{
    check.integer.scalar(id.n, min=0, null.ok=TRUE)
    id.n <- if(is.null(id.n)) 0 else as.integer(id.n)
    ncases <- if(!is.null(nrow(resids))) nrow(resids) else length(resids)
    if(id.n > ncases)
        id.n <- ncases
    if(id.n > 0)
        sort.list(abs(resids), decreasing = TRUE)[1:id.n]
    else
        NULL
}

# plotmor.R: plot model residuals
#
# "plotres" is already used in several R packages, hence the name "plotmor"

plotmor <- function(
    x           = stop("no 'x' arg"),
    which       = 2:4,
    info        = FALSE,
    student     = FALSE,  # divide by stderr * sqrt(1 - h_ii)
    delever     = FALSE,  # divide by sqrt(1 - h_ii)
    level       = 0,
    versus      = 1,

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
    cex.main    = 1.1,
    caption     = if(do.par) NULL else "",
    xlab        = NULL,
    ylab        = NULL,

    pch         = 20,
    col.loess   = "red",
    lwd.loess   = 1,
    col.cv      = "lightblue",
    col.qq      = "gray",
    col.grid    = "lightgray",
    col.points  = 1,
    cex.points  = NULL,
    shade.pints = "mistyrose2",
    shade.cints = "mistyrose4",

    cum.grid    = "percentages",
    legend.pos  = NULL,

    ...)
{
    object.name <- substr(deparse(substitute(x)), 1, 40)
    object <- x  # minimize confusion with x, the regression input matrix
    remove(x)    # not necessary but prevents mistakes later
    stop.if.dots.used("plotmor", ...)
    plotmo::check.index(which, "which", 1:9)
    stopifnot(length(do.par) == 1)
    stopifnot(do.par == 0 || do.par == 1 || do.par == 2)
    check.integer.scalar(npoints, min=-1, null.ok=TRUE)
    check.integer.scalar(id.n, min=0, null.ok=TRUE)
    check.integer.scalar(nresponse, min=1)
    info <- check.boolean(info)
    rinfo <- get.rinfo(object, nresponse, student, delever,
                if(any(which %in% c(3,5:9))) "plotted as a star"
                else "ignored")
    fitted <- get.plotmor.fitted(object, nresponse, rinfo) # returns a matrix
    if(is.null(cex.points))
        cex.points <- get.cex.points(npoints, length(rinfo$resids))
    if(is.null(labels.id))
        labels.id <- paste(1:length(rinfo$resids))
    if(is.null(env))
        env <- get.model.env(object, parent.frame)

    # get the values we will plot against (by default the fitted values)
    temp <- get.versus.mat(which, versus, object, fitted, env, nresponse)
        which              <- temp$which      # plots we don't want will be removed
        versus.mat         <- temp$versus.mat # either fitted, response, x, or bx
        icolumns           <- temp$icolumns   # desired column indices in versus.mat
        versus.is.response <- temp$versus.is.response
        versus.is.index    <- temp$versus.is.index
        versus.is.leverage <- temp$versus.is.leverage

    nfigs <- length(which) * length(icolumns)
    if(nfigs == 0) {
        warning0("plot.earth: nothing to plot")
        return(invisible())
    }
    if(do.par && nfigs > 1) {
        old.par <- par.for.plot(do.par, nfigs, cex.main, caption)
        if(do.par == 1)
            on.exit(par(old.par))
    }
    for(icolumn in icolumns)
        plotmor1(object, versus.mat[, icolumn], fitted, rinfo,
            which, info, student, delever, level,
            nresponse, npoints, id.n, labels.id, center, loess.f,
            env, xlim, ylim, main, cex.main, caption,
            if(is.null(xlab)) colnames(versus.mat)[icolumn] else xlab,
            ylab, pch, col.loess, lwd.loess,
            col.cv, col.qq, col.grid, col.points,
            cex.points, shade.pints, shade.cints,
            cum.grid, legend.pos,
            versus.is.response, versus.is.index, versus.is.leverage)

    if(do.par && (is.null(caption) || any(nchar(caption)))) { # show caption?
        trim <- is.null(caption) # trim if auto-generate the caption
        show.caption(get.caption(which, object, nresponse, object.name, caption))
    }
    invisible()
}
versus.err.msg <- function()
{
    stopf("%s\n%s\n%s\n%s\n%s\n%s\n%s",
        "versus must be an integer or a string:",
        "  1    fitted (default)",
        "  2    obs numbers",
        "  3    response",
        "  4    leverages",
        "  \"*\"  earth terms",
        "  \"\"   predictors")
}
get.versus.mat <- function(which, versus, object, fitted, env, nresponse)
{
    versus.mat         <- fitted
    icolumns           <- 1
    versus.is.response <- FALSE
    versus.is.index    <- FALSE
    versus.is.leverage <- FALSE
    trim.which         <- FALSE
    got.versus         <- FALSE

    if(is.numeric(versus)) {
        got.versus <- TRUE
        if(length(versus) != 1)
            stop0("illegal versus (length of versus must be 1 when versus is numeric)")
        if(floor(versus) != versus)
            versus.err.msg()
        if(versus == 1)         # versus is fitted
            NULL
        else if(versus == 2)    # versus is index
            versus.is.index <- TRUE
        else if(versus == 3) {  # versus is response
            versus.mat <- plotmo::get.plotmo.y.wrapper(object,
                            env, y.column=nresponse,
                            expected.len=NROW(fitted(object)), trace=0)
            versus.mat <- matrix(versus.mat, ncol=1)
            colnames(versus.mat) <- "Response"
            versus.is.response <- TRUE
        } else if(versus == 4) { # versus is leverage
            versus.is.leverage <- TRUE
            versus.mat <- matrix(hatvalues(object), ncol=1)
            colnames(versus.mat) <- "Leverage"
        } else
            versus.err.msg()
    }
    else if(!is.character(versus))
        versus.err.msg()
    else if(length(versus) == 1 && nchar(versus) &&
            substr(versus, 1, 1) == "*") { # use the earth basis matrix
        got.versus <- TRUE
        trim.which <- TRUE
        if(is.null(object$bx))
            stop0("versus=\"*\" is not allowed for this object (object$bx is NULL)")
        if(NCOL(object$bx) == 1) # intercept only model?
            versus.mat <- object$bx
        else {
            versus.mat <- object$bx[,-1, drop=FALSE] # drop the intercept
            versus <- substring(versus, 2)           # drop the "*"
            icolumns <- plotmo::check.index(versus, "versus",
                                            1:NCOL(versus.mat),
                                            colnames=colnames(versus.mat))
        }
    }
    if(!got.versus) { # user specified x variables
        trim.which <- TRUE
        prefix <- substr(versus, 1, 1)
        # following are needed if versus is a vector
        if(any(prefix == "*"))
            stop0("\"*\" is not allowed in this context in the versus argument\n",
                  "      Your versus argument is ", quote.with.c(versus))
        versus.mat <- plotmo::get.plotmo.x.wrapper(object, env, trace=0)
        if(is.null(ncol(versus.mat)))
            versus.mat <- matrix(versus.mat, nrow=length(versus.mat))
        if(is.null(colnames(versus.mat)))
            colnames(versus.mat) <- paste("x", 1:ncol(versus.mat), sep="")
        icolumns <- plotmo::check.index(versus, "versus",
                                        1:NCOL(versus.mat),
                                        colnames=colnames(versus.mat))
    }
    if(trim.which) {
        # remove all entries from which except standard resid and abs resid plots
        which <- which[which %in% c(3,5)]
        if(length(which) == 0)
            warning0("the \"which\" argument is empty after removing certain ",
                     "plots because \"versus\" was specified")
    }
    list(which              = which,
         versus.mat         = versus.mat, # either fitted, response, x, or bx
         icolumns           = icolumns,   # desired column indices in versus.mat
         versus.is.response = versus.is.response,
         versus.is.index    = versus.is.index,
         versus.is.leverage = versus.is.leverage)
}
plotmor1 <- function(
    object,
    versus1, # what we plot along the x axis, a vector
    fitted,
    rinfo,

    which,
    info,
    student,
    delever,
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
    col.loess,
    lwd.loess,
    col.cv,
    col.qq,
    col.grid,
    col.points,
    cex.points,
    shade.pints,
    shade.cints,

    cum.grid,
    legend.pos,
    versus.is.response,
    versus.is.index,
    versus.is.leverage)
{
    stopifnot(NCOL(versus1) == 1)
    iresids <- get.iresids(npoints, rinfo$resids, versus1, id.n, versus.is.leverage)
    if(is.null(col.loess))
        col.loess <- "red"
    for(iwhich in seq_along(which)) {
        if(which[iwhich] == 2) {            #--- cumulative distribution ---
            if(length(which) !=1 || is.null(xlim))
                xlim <- range(abs(rinfo$scale * rinfo$resids), na.rm=TRUE)
            plot.cum(rinfo, main,
                     xlim=xlim, col=1, col.grid, cum.grid, add=FALSE, jitter=0)

        } else if(which[iwhich] == 4)       #--- QQ plot ---
            plot.qq(rinfo, iresids, main, info,
                    col.points, cex.points, col.qq, pch, id.n, labels.id, col.loess)

        else {                               #--- various residual plots ---
            plot.resids(object, versus1, fitted, rinfo, iresids, id.n,
                which[iwhich], info, student, level, labels.id, center, loess.f,
                if(length(which)!=1) NULL else xlim,
                if(length(which)!=1) NULL else ylim,
                main, cex.main, xlab, ylab,
                pch, col.loess, lwd.loess,
                col.cv, col.grid, col.points, cex.points,
                shade.pints, shade.cints, legend.pos,
                versus.is.response, versus.is.index, versus.is.leverage)
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

    versus1,
    fitted,
    rinfo,
    iresids,
    id.n,

    which,
    info,
    student,
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
    col.loess,
    lwd.loess,
    col.cv,
    col.grid,
    col.points,
    cex.points,
    shade.pints,
    shade.cints,

    legend.pos,
    versus.is.response,
    versus.is.index,
    versus.is.leverage)
{
    stopifnot(length(which) == 1)
    info <- check.boolean(info)

    ok <- which %in% c(3,5:9)
    if(!all(ok))
        stop0("which=", which[!ok][1], " is not allowed")

    id.indices <- NULL
    if(which %in% c(3, 4:8))
        id.indices <- get.id.indices(rinfo$scale * rinfo$resids, id.n,
                                     if(versus.is.leverage) hatvalues(object) else NULL)

    level <- check.level.arg(level, zero.ok=TRUE)
    if(which %in% (5:9))
        level <- 0 # no pints

    varmod <- NULL
    if(which == 3 && !is.null(object$varmod) && is.specified(level) &&
            (is.specified(shade.pints) || is.specified(shade.cints)))
        varmod <- object$varmod

    pints <- NULL
    cints <- NULL
    if(!is.null(varmod) && is.specified(level)) {
        if(is.specified(shade.pints))
            pints <- predict(object, interval="pint", level=level)
        if(is.specified(shade.cints))
            cints <- predict(object, interval="cint", level=level)
    }
    resids <- rinfo$scale * rinfo$resids

    if((which %in% 7:9))
        check.that.most.are.positive(
            versus1, "fitted", sprintf("which=%d", which), "nonpositive")

    # TODO following is redundant after above check?
    if(which %in% 7:9) # abs(resids) must be strictly positive to take their log
        check.that.most.are.positive(
            abs(resids), "abs(residuals)", sprintf("which=%d", which), "zero")

    trans.versus <- trans.versus(versus1[iresids], which)
    trans.resids <- trans.resids(resids[iresids], which)

    xlab <- get.resids.xlab(xlab, which, versus1, versus.is.index)
    ylab <- get.resids.ylab(ylab, which, rinfo$name)
    main <- get.resids.main(main, xlab, ylab, varmod, level) # title of plot

    col.points <- repl(col.points, length(resids)) # recycle
    cex.points <- repl(cex.points, length(resids))
    pch        <- repl(pch, length(resids))

    ylim <- get.resids.ylim(ylim, fitted, trans.resids, which, info, student,
                            id.indices, center, pints, cints, rinfo$scale,
                            versus.is.leverage)

    x <- if(versus.is.index) 1:length(trans.versus) else trans.versus

    xlim <- get.resids.xlim(xlim, which, x, trans.versus, ylim, versus.is.leverage)

    plot(x, trans.resids,
         col=col.points[iresids], cex=cex.points[iresids],
         pch=pch[iresids],
         main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)

    if(!is.null(varmod) && !versus.is.leverage) {
        if(is.specified(shade.pints))
            plot.pint.resids(pints, versus1, fitted, rinfo$scale, shade.pints, versus.is.index)
        if(is.specified(shade.cints))
            plot.pint.resids(cints, versus1, fitted, rinfo$scale, shade.cints, versus.is.index)
    }
    if(which != 9)
        abline(h=0, lty=1, col=col.grid) # axis

    # Removed because cooks levels probably not meaningful for non constant variance
    # if(versus.is.leverage && student)
    #   show.cooks.levels(object, resids)

    # replot points (because they are obscured by abline and shaded areas)
    points(x, trans.resids,
           col=col.points[iresids], cex=cex.points[iresids],
           pch=pch[iresids])

    # plot points with leverage one as stars
    plot.bad.leverage.as.star(x, rinfo, iresids, cex.points,
                              col=if(is.specified(col.loess)) col.loess else "red")

    check.numeric.scalar(lwd.loess, null.ok=TRUE)
    if(is.null(lwd.loess))
        lwd.loess <- 1

    if(which != 9)
        plot.loess(x, trans.resids, loess.f, col.loess, lwd.loess, rinfo$scale[iresids])

    if(!is.null(varmod) && !is.null(object$cv.oof.fit.tab) && is.specified(col.cv))
        plot.oof.meanfit(object$cv.oof.fit.tab, fitted, versus1, rinfo,
                         which, col.cv, versus.is.index)

    # TODO implement id.indices for versus.is.index
    if(!is.null(id.indices) && !versus.is.index)
        plotrix::thigmophobe.labels(
            x=trans.versus(versus1, which)[id.indices],
            y=trans.resids(resids, which)[id.indices],
            labels=labels.id[id.indices],
            offset=.33, font=2, cex=.8, xpd=NA,
            col=if(is.specified(col.loess)) col.loess else "red")

    if(info)
        add.resids.info(which, info, versus1, resids,
                        col.lm="lightblue", lwd.lm=lwd.loess+1,
                        versus.is.response, versus.is.index, versus.is.leverage)

    else if(!is.null(varmod) && (is.null(legend.pos) || !all(is.na(legend.pos))) &&
            (is.specified(col.loess) || is.specified(col.cv))) {
        # add legend, else red and blue may confuse the user
        # no legend if info=TRUE because not enough space

        if(is.null(legend.pos))
            legend.pos <- "bottomleft"
        legend.txt <- NULL
        legend.col <- NULL
        legend.lwd <- NULL
        legend.lty <- NULL
        if(which != 9 && is.specified(col.loess)) { # loess plotted?
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
        legend(legend.pos, legend=legend.txt,
               col=legend.col, lwd=legend.lwd, lty=legend.lty, bg="white", cex=.8)
    }
}
check.level.arg <- function(level, zero.ok)
{
    if(identical(level, NA) || is.null(level)) # treat NA and NULL as 0
        level <- 0
    check.numeric.scalar(level)
    if(!((zero.ok && level == 0) || level >= .5 || level < 1)) {
        msg <- if(zero.ok)
                    "level argument must be zero or between 0.5 and 1"
                else
                    "level argument must be between 0.5 and 1"
        stop0(msg, " (you have level=", level, ")")
    }
    level
}
get.resids.xlim <- function(xlim, which, x, trans.versus, ylim, versus.is.leverage)
{
    if(is.null(xlim)) { # auto xlim?
        if(which==9) {
            # don't show lower 5% of points
            quant <- quantile(trans.versus, prob=c(.05, 1), na.rm=TRUE)
            min <- quant[1]
            max <- quant[2]
            # extra left margin so slope of linear fit not flattened
            if(min > .2 * ylim[1])
                min <- .2 * ylim[1]
            xlim <- c(min, max)
        } else if(versus.is.leverage)
            xlim <- c(0, 1.1 * max(x, na.rm=TRUE)) # extra room for labels on hi lev points
        else
            xlim <- range(x, na.rm=TRUE)
    }
    stopifnot(is.numeric(xlim))
    stopifnot(length(xlim) == 2)
    fix.lim(xlim)
}
get.resids.ylim <- function(ylim, fitted, resids, which, info, student,
                            id.indices, center, pints, cints, scale,
                            versus.is.leverage)
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
        range <- abs(max - min)
        if(!is.null(id.indices)) { # space for point labels
            min <- min - .03 * range
            max <- max + .03 * range
        }
        if(which %in% (5:8))
            min <- 0
        else if(which == 3 && center) {
            # want symmetric ylim so can more easily see asymmetry
            if(abs(min) > abs(max))
                max <- -min
            else if(abs(max) > abs(min))
                min <- -max
        } else if(which == 9)
            max <- max + .5  # more space on top, looks better

        # if(versus.is.leverage && student)
        #     max <- max + .1 * range # space for cooks distance legend

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
plot.pint.resids <- function(pints, versus, fitted, scale, shade.pints, versus.is.index)
{
    if(!is.null(pints)) {
        pints <- get.pint.resids(pints, fitted, scale)
        # abscissa must be ordered for polygon to work
        order <- order(versus)
        versus <- versus[order]
        pints <- pints[order,]
        x <- if(versus.is.index) c(1:length(versus), length(versus):1)
             else                trans.versus(c(versus, rev(versus)), 0)
        polygon(x, trans.resids(c(pints$lwr, rev(pints$upr)), 0),
                col=shade.pints, border=shade.pints)
    }
}
# show.cooks.levels <- function(object, resids)
# {
#   # based on code in stats::plot.lm.R
#   leverage <- hatvalues(object)
#   p <- length(coef(object))
#   cook.levels <- c(0.5, 1.0)
#   usr <- par("usr")
#     leverage.range <- range(leverage, na.rm=TRUE) # though should never have NA
#   x <- seq.int(0, 1, length.out=101)
#   col <- "darkgray"
#   for(cook.level in cook.levels) {
#       cl <- sqrt(cook.level * p *(1 - x) / x)
#       lines(x,  cl, col=col)
#       lines(x, -cl, col=col)
#   }
#   # don't use bottomleft like plot.lm because we may plot the density there
#   legend("topleft", legend="Cook's distance", lty=1, col=col, bty="n", bg="white")
#   xmax <- min(0.99, usr[2])
#   ymult <- sqrt(p * (1 - xmax) / xmax)
#   axis(4, at=c(-sqrt(rev(cook.levels)) * ymult, sqrt(cook.levels)*ymult),
#        labels=paste(c(rev(cook.levels), cook.levels)),
#        mgp=c(.25,.25,0), las=2, tck=0,
#        cex.axis=1, col.axis=col)
# }

# Plot points with leverage one as stars.  We plot them on
# the axis, which is arguably not correct but still useful.

plot.bad.leverage.as.star <- function(x, rinfo, iresids, cex.points, col)
{
    which <- which(is.na(rinfo$scale[iresids]))
    if(length(which) > 0) {
        # use .8 on the cex else the stars get displayed as very big
        points(x[which], 0,
               col=col, cex=.8 * cex.points[iresids], pch=8) # pch 8 is a star
        # add label if possible (not possible if not all points are plotted, see npoints)
        if(length(iresids) == length(rinfo$scale)) {
            label <- which(is.na(rinfo$scale))
            text.on.white(x[which], 0, label, cex=.8, adj=-.5, col=col, font=2, xpd=NA)
        }
    }
}
plot.loess <- function(versus, resids, loess.f, col.loess, lwd.loess, scale)
{
    if(is.specified(col.loess)) {
        check.numeric.scalar(loess.f)
        stopifnot(loess.f > .01 && loess.f < 1)

        # na.rm is needed if we take logs of nonpos, see check.that.most.are.positive.
        # That's why we calculate delta explicitly instead of using lowess default.
        delta <- .01 * diff(range(versus, na.rm=TRUE))

        # Replace points with NA scale with 0 (else lowess stops at the NA).
        # Zero is appropriate because the points are 0 resids with leverage 1.
        resids[which(is.na(scale))] <- 0

        # We use lowess rather than loess because loess tends to give ugly warnings,
        # but the argument is named "col.loess" for backward compatibility.
        smooth <- lowess(versus, resids, f=loess.f, delta=delta)
        lines(smooth$x, smooth$y, col=col.loess, lwd=lwd.loess)
    }
}
get.resids.xlab <- function(xlab, which, versus1, versus.is.index)
{
    if(is.null(xlab)) {
        stopifnot(!is.null(colnames(versus1)))
        xlab <- colnames(versus1)[1]
    }
    if(which %in% (7:9))
        xlab <- sprintf("Log %s", xlab)
    if(versus.is.index)
        xlab <- sprintf("%s index", xlab)
    xlab
}
get.resids.ylab <- function(ylab, which, rinfo.name)
{
    if(is.null(ylab))
        ylab <- sprintf("%ss", rinfo.name)
    if(which == 5)
        ylab <- sprintf("Abs %s", ylab)
    else if(which == 6)
        ylab <- sprintf("Sqrt Abs %s", ylab)
    else if(which == 7)
        ylab <- sprintf("Abs %s", ylab)
    else if(which == 8)
        ylab <- sprintf("Cube Root Squared %s", ylab)
    else if(which == 9)
        ylab <- sprintf("Log Abs %s", ylab)
    ylab
}
get.resids.main <- function(main, xlab, ylab, varmod, level) # title of plot
{
    # TODO should really use strwidth for newline calculation
    # The "Fitted" helps with limitations of formula below and retains back compatibility
    newline <- xlab != "Fitted" && xlab != "Fitted index" && xlab != "Response" &&
               nchar(ylab) + nchar(xlab) > 15

    if(length(grep("Studentized", ylab)) > 0 ||
       length(grep("Delevered",   ylab)) > 0)
        newline <- TRUE

    if(is.null(main))
        main <- sprintf("%s%svs %s", ylab, if(newline) "\n" else " ", xlab)
    if(!is.null(varmod) && !newline) # two newlines is too many
        main <- sprintf("%s\n%g%% level shaded", main, 100*(level))

    main
}
# plot resids of oof meanfit with col.cv (default lightblue)

plot.oof.meanfit <- function(cv.oof.fit.tab, fitted, versus1,
                             rinfo, which, col.cv, versus.is.index)
{
    # mean of each row of oof.fit.tab
    meanfit <- apply(cv.oof.fit.tab,  1, mean)
    meanfit <- rinfo$scale * (meanfit - fitted)
    order <- order(versus1)
    trans.versus1 <- trans.versus(versus1[order], which)
    x <- if(versus.is.index) 1:length(trans.versus1) else trans.versus1
    lines(x, trans.resids(meanfit[order], which), col=col.cv)
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
add.resids.info <- function(which, info, versus1, resids, col.lm, lwd.lm,
                            versus.is.response, versus.is.index, versus.is.leverage)
{
    trans.versus <- trans.versus(versus1, which)
    x <- if(versus.is.index) 1:length(trans.versus) else trans.versus
    plot.density.along.the.bottom(x)
    if(!versus.is.leverage) {
        lm.text <- ""
        slope.text <- ""
        if(which == 5 || which == 9) # add linear regression line?
            slope.text <- add.rlm.line(which, info, versus1, resids, col.lm, lwd.lm,
                                   versus.is.index)
        # exact=FALSE else get warning "Cannot compute exact p-value with ties"
        cor.abs <- cor.test(versus1, abs(resids), method="spearman", exact=FALSE)
        if(versus.is.response) {
            cor <- cor.test(versus1, resids, method="spearman", exact=FALSE)
            text <- sprintf("spearman abs  %.2f   resids %.2f\n%s",
                             cor.abs$estimate, cor$estimate, slope.text)
        } else
            text <- sprintf("spearman abs  %.2f%s", cor.abs$estimate, slope.text)
        usr <- par("usr") # xlim1, xlim2, ylim1, ylim2
        text.on.white(usr[1] + strwidth("x", font=2),
                      usr[4] - strheight(text, font=2) - .5 * strheight("X", font=2),
                      text,
                      cex=1, adj=c(0, 0), font=2, col=1, xpd=NA)
    }
}
add.rlm.line <- function(which, info, versus1, resids, col.lm, lwd.lm, versus.is.index)
{
    trans.resids <- trans.resids(resids, which)
    trans.versus <- trans.versus(versus1, which)
    x <- if(versus.is.index) 1:length(trans.versus) else trans.versus
    if(which == 9) {
        # ignore lower 10% of points (very small residuals i.e. very neg logs)
        quant <- quantile(trans.versus, prob=.1, na.rm=TRUE)
        ok <- which(x > quant)
        x <- x[ok]
        trans.resids <- trans.resids[ok]
    }
    # # regression on 10 bootstrap samples so we can see variance of versus1
    # for(i in 1:10) {
    #   j <- sample.int(length(x), replace=TRUE)
    #   trans.resids1 <- trans.resids[j]
    #   trimmed.trans.fit1 <- x[j]
    #   rlm <- MASS::rlm(trans.resids1~trimmed.trans.fit1,
    #                    method="MM", na.action="na.omit")
    #   abline(coef(rlm), col="lightgray", lwd=.6)
    # }

    # robust linear regression of trans.resids on x
    # na.omit is needed if some versus1 (or resids) were nonpos so log(versus1) is NA
    rlm <- MASS::rlm(trans.resids~x, method="MM", na.action="na.omit")
    abline(coef(rlm), col=col.lm, lwd=lwd.lm)
    sprintf("    slope %.2g", coef(rlm)[2])
}
get.plotmor.fitted <- function(object, nresponse, rinfo) # returns a matrix
{
    fitted <- fitted(object)
    if(is.null(fitted))
        stop0("plotmor: object$fitted is NULL")
    if(is.null(ncol(fitted)))
        fitted <- matrix(fitted, ncol=1)
    nresponse <- plotmo::check.index(nresponse, "nresponse", fitted, is.col.index=TRUE)
    fitted <- fitted[, nresponse, drop=FALSE]
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
    if(which %in% (7:9))
        my.log10(versus)
    else
        versus
}
trans.resids <- function(resid, which) # transform the residuals
{
    if(which == 5)
        abs(resid)
    else if(which == 6)
        sqrt(abs(resid))
    else if(which == 7)
        abs(resid)
    else if(which == 8) {
        # do it in two steps so no problem with negative numbers
        resid <- resid^2
        resid^(1/3)
    } else if(which == 9)
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
    # TODO check what happens here if NA in trans.resids (leverage==1)
    cum <- ecdf(trans.resids)
    if(jitter > .05)
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
            q <- quantile(trans.resids, probs=c(0, .25, .50, .75, .9, .95, 1), na.rm=TRUE)
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
                    col.points, cex.points, col.qq, pch, id.n, labels.id, col.loess)
{
    par(pty="s") # square
    # we figure out the shape of the qq line with all resids but
    # plot only npoints points (selecting them with iresids)
    resids <- rinfo$scale * rinfo$resids
    # TODO check what happens here if NA in trans.resids (leverage==1)
    qq <- qqnorm(resids, main=main, plot.it=FALSE)
    xlim <- ylim <- NULL
    id.indices <- get.id.indices(resids, id.n)
    if(!is.null(id.indices)) {
        # figure out xlim and ylim with extra space for point labels
        min <- min(qq$x, na.rm=TRUE)
        max <- max(qq$x, na.rm=TRUE)
        xlim <- c(min - .1 * (max - min), max + .1 * (max - min))
        min <- min(qq$y, na.rm=TRUE)
        max <- max(qq$y, na.rm=TRUE)
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
        plotrix::thigmophobe.labels(
            x=qq$x[id.indices], y=qq$y[id.indices],
            labels=labels.id[id.indices],
            offset=.33, font=2, cex=.8, xpd=NA,
            col=if(is.specified(col.loess)) col.loess else "red")
    par(pty="m") # back to maximal (the default)
}
# get the subset of residuals we're going to display, ordered on versus

get.iresids <- function(npoints, resids, versus, id.n, versus.is.leverage)
{
    check.vec(resids, "residuals", length(versus))
    iresids <- 1:nrow(resids)
    # note that to keep things simple, if leverage we always plot all points
    if(npoints > 0 && npoints < nrow(resids) && !versus.is.leverage) {
        # take a sample, but make sure it includes the biggest resids
        npoints <- min(npoints, nrow(resids))
        nbiggest <- min(max(8, id.n), npoints)
        isorted <- order(abs(resids), decreasing=TRUE)
        ikeep <- 1:nbiggest
        if(npoints > nbiggest)
            ikeep <- c(ikeep, seq(from=nbiggest + 1, to=nrow(resids),
                        length.out=npoints - nbiggest))
        iresids <- isorted[ikeep]
    }
    # index vector into resids and versus, sorted on versus
    iresids[order(versus[iresids])]
}
# get the indices of the id.n biggest resids
# requires a sort so is quite slow

get.id.indices <- function(resids, id.n, hatvalues=NULL)
{
    check.integer.scalar(id.n, min=0, null.ok=TRUE)
    id.n <- if(is.null(id.n)) 0 else as.integer(id.n)
    ncases <- if(!is.null(nrow(resids))) nrow(resids) else length(resids)
    if(id.n > ncases)
        id.n <- ncases
    id.indices <- NULL
    if(id.n > 0) {
        id.indices <- sort.list(abs(resids), decreasing=TRUE, na.last=TRUE)[1:id.n]
        if(!is.null(hatvalues)) {
            # add the worst hatvalues
            ibad <- order(hatvalues, decreasing=TRUE)[1:id.n]
            # ibad <- ibad[hatvalues[ibad] > .5]
            id.indices <- unique(c(id.indices, ibad))
        }
    }
    id.indices
}

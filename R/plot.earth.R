# plot.earth.R: plotting routines for the earth package
#
# Comments containing "$$" mark known issues.
# Stephen Milborrow Mar 2007 Petaluma
#
# $$ allow model comparison using newdata instead of data used to build model
# $$ add lty.grsq and lty.rsq vector arguments to plot.earth.models
# $$ add legend entries argument i.e. vector of strings
# $$ add an x range option to model plot
# $$ add code to not print percentages if too close to each other in plot.cum
# $$ allow col.vline=col.grsq in plot.earth.models
# $$ add ycolumn to plot.earth.models

#--------------------------------------------------------------------------------------------
plot.earth <- function(
    x = stop("no 'x' arg"),     # an earth object (called x for consistency with generic)
    which       = 1:4,          # which plots to plot: 1 model 2 cumul 3 residuals 4 qq

    ycolumn     = 1,            # which column of response to use if multiple col response

    caption = if(do.par) NULL else "",
                                # overall caption
                                #   "string"  string
                                #   ""        no caption
                                #   NULL      generate a caption from x$call

    col.rsq     = "lightblue",  # color of rsq line
    col.loess   = col.rsq,      # color of residuals plot loess line
    col.qq      = col.rsq,      # color of qq line
    col.grid    = "grey",       # color of grid lines in cumulative distribution plot
    col.vline   = "grey",       # color of vertical line at best model in model plot
    lty.vline   = 3,            # color of vertical line at best model in model plot
    col.legend  = 1,            # legend (inside plot area) of model plot, if 0 no legend
    col.npreds  = 1,            # if 0, don't show "nbr of preds" in model plot

    nresiduals  = 1000,         # max number of residuals to plot, -1 for all
                                # default limits to 1000 to minimize overplotting
                                # the largest few residuals are always included

    cum.grid   = "percentages", # "none"        no grid on cumulative distribution plot
                                # "grid"        add grid
                                # "percentages" add grid and percentage labels

    rlim   = c(-1,-1),          # min.max values on rsq and grsq plot
                                # special value min=-1 means use min(rsqVec[-1])
                                # special value max=-1 means use max(rsqVec)

    jitter      = 0,            # non zero val allows overlaid plots to be visible

    id.n        = 3,            # number of residuals to be labelled
    labels.id   = rownames(residuals(x)),  # residual names

    legend.pos  = NULL,         # NULL means auto, else specify c(x,y) in user coords
    do.par      = TRUE,         # call par() as appropriate

    main        = NULL,         # par() settings, listed so caller can override our use
    pch         = 1,
    ...)                        # extra args passed to plotting and method funcs
{
    get.iresiduals <- function(nresiduals, residuals, fitted.values)
    {
        if(any(!is.finite(residuals)))
            warning1("non finite residuals")
        stopifnot(nrow(residuals) == nrow(fitted.values))
        iresiduals = 1:nrow(residuals)
        if(nresiduals > 0 && nresiduals < nrow(residuals)) {
            # take a sample, but make sure it includes the largest residuals

            nresiduals = min(nresiduals, nrow(residuals))
            nlargest = min(max(3, id.n), nresiduals)
            isorted = order(abs(residuals), decreasing=TRUE)
            ikeep = 1:nlargest
            if(nresiduals > nlargest)
                ikeep = c(ikeep, seq(from=nlargest + 1, to=nrow(residuals),
                            length.out=nresiduals - nlargest))
            iresiduals <- isorted[ikeep]
        }
        # index vector into residuals and fitted.values, sorted on fitted.values
        iresiduals[order(object$fitted.values[iresiduals])]
    }
    get.id.indices <- function()  # get the indices of the id.n largest residuals
    {
        id.n <- if(is.null(id.n)) 0 else as.integer(id.n)
        if(id.n < 0)
            stop1("negative 'id.n'")
        ncases <- nrow(object$residuals)
        if(id.n > ncases)
            id.n <- ncases
        if(id.n > 0)
            sort.list(abs(residuals), decreasing = TRUE)[1:id.n]
        else
            NULL
    }
    plot.residuals <- function(fitted.values, residuals)    # residuals versus fitted
    {
        plot.loess <- function()
        {
            do.loess = (col.loess != 0)
            if(length(unique(object$fitted.values[,1])) < 10)  # prevent warnings from loess
                do.loess = FALSE
            if(do.loess) {
                fitted.values <- fitted.values
                loessModel <- loess(residuals ~ fitted.values, degree=2)
                y.loess <- predict(loessModel, newdata=fitted.values)
                lines(fitted.values, y.loess, col=col.loess, ...)
            }
        }
        if(is.null(main))
            main="Residuals vs Fitted"
        plot(fitted.values, residuals, main=main,xlab="Fitted",ylab="Residuals",pch=pch,...)
        abline(h=0, lty=3, col=col.grid)
        plot.loess()
        if(!is.null(id.indices)) {
            mean = mean(fitted.values)
            text(fitted.values[id.indices], residuals[id.indices], labels.id[id.indices],
                 pos=ifelse(fitted.values[id.indices] < mean, 4, 2), cex=.8)
        }
        NULL
    }
    plot.qq <- function(residuals)  # qqplot of residuals
    {
        par(pty="s")                # square
        if(is.null(main))
            main="Normal Q-Q"
        qq <- qqnorm(residuals, main=main,
                xlab="Theoretical Quantiles", ylab="Residual Quantiles", pch=pch, ...)
        qqline(residuals, col=col.qq, ...)
        if(!is.null(id.indices))
            text(qq$x[id.indices], qq$y[id.indices], labels.id[id.indices],
                 pos=ifelse(qq$x[id.indices] < 0, 4, 2), cex=.8)
        par(pty="m")                # back to maximal (the default)
        NULL
    }
    # return the response.name (for prepending to caption), if appropriate
    get.caption.prefix <- function(residuals, ycolumn, caption)
    {
        if (!is.null(caption))
            return(NULL)                # don't modify caption explictly set by user
        colnames <- colnames(residuals)
        if(!is.null(colnames) && !is.null(colnames[1]) && colnames[1] != "")
            colnames[1]
        else if(NCOL(residuals) > 1)
            paste("Response", ycolumn)
        else
            NULL
    }
    # plot.earth starts here
    object <- x
    plot.earth.prolog(object, deparse(substitute(x)))
    check.index.vec("which", which, 1:4)
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    nfigs <- sum(show)
    if(nfigs == 0) {
        warning1("nothing to plot")
        return(invisible())
    }
    must.trim.caption <- is.null(caption)
    check.index.vec("ycolumn", ycolumn, object$residuals,
                    check.empty = TRUE, use.as.col.index=TRUE)
    caption.prefix = get.caption.prefix(object$residuals, ycolumn, caption)
    caption <- get.caption.from.call(caption, object)
    if(!is.null(caption.prefix))
        caption <- paste(caption.prefix, ": ", caption, sep="")
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        on.exit(par(old.par))
        nrows <- ceiling(sqrt(nfigs))
        par(mfrow=c(nrows, nrows))
        par(mar = c(4, 4, 2, 1))    # small margins and text to pack figs in
        par(mgp = c(1.6, 0.6, 0))   # flatten axis elements
        par(cex = 0.7)
        make.space.for.caption(caption)
    }
    nresponses <- NCOL(object$residuals)
    if(nresponses > 1) {
        # already called check.index.vec(ycolumn) above, so no need to call again
        object$fitted.values <- object$fitted.values[, ycolumn, drop=FALSE]
        object$residuals <- object$residuals[, ycolumn, drop=FALSE]
    }
    if(show[1]) {
        rlim <- get.rlim(object, rlim, 1, col.rsq)
        plot.earth.model(object, col.grsq=1, col.rsq=col.rsq, lty.rsq=1,
                col.npreds=col.npreds, col.vline=col.vline, lty.vline=lty.vline, col.vseg=0,
                col.legend=col.legend, rlim=rlim, add=FALSE, do.par=do.par,
                max.nterms=length(object$rss.per.subset),
                max.npreds=max(get.nused.preds.per.subset(object$dirs, object$prune.terms)),
                jitter=jitter, legend.pos=legend.pos, nresponses,
                main=main, ...)
    }
    if(any(show[2:4])) {
        if(NROW(object$fitted.values) != nrow(object$residuals))
            stop1("NROW(object$fitted.values) ",  NROW(object$fitted.values),
                " != nrow(object$residuals) ", nrow(object$residuals))
        iresiduals <- get.iresiduals(nresiduals, object$residuals, object$fitted.values)
        residuals <- object$residuals[iresiduals]
        fitted.values <- object$fitted.values[iresiduals]
    }
    if(show[2])
        plot.cum(residuals, main, xlim=range(abs(object$residuals)),
                 col=1, col.grid, cum.grid, add=FALSE, ...)
    if(any(show[3:4])) {
        id.indices <- get.id.indices()
        if(is.null(labels.id))
            labels.id <- paste(iresiduals)
        else
            labels.id <- labels.id[iresiduals]
    }
    if(show[3])
        plot.residuals(fitted.values, residuals)
    if(show[4])
        plot.qq(residuals)
    show.caption(caption, trim=must.trim.caption)
    invisible()
}

#--------------------------------------------------------------------------------------------
plot.earth.models <- function(  # compare earth models by plotting them
    x           = stop("no 'x' arg"), # list of earth objects, can just be one object
                                      # called x for consistency with the generic
    which       = c(1:2),       # which plots to plot: 1 model 2 cumul
    caption = "",               # overall caption
                                #   "string"  string
                                #   ""        no caption
                                #   NULL      generate a caption from objects[1]$call

    rlim        = c(0,1),       # min.max values on rsq and grsq plot, same as plot.earth
    jitter      = 0,            # non zero val allows overlaid plots to be visible
                                # All the col arguments: if 0 don't plot this graph element
    col.grsq    = discrete.plot.cols(length(x)),
    col.rsq     = 0,            # 0 means don't superimpose rsq plot
                                # NULL means use col.grsq or col.rsq in npreds plot
                                # else specify a color vector for npred plot

    col.npreds  = 0,            # same usage as col.npred but for nbr of preds plot

    col.cum     = NULL,         # same usage as col.npred but for cumul dist plot
    col.vline   = "grey",
    lty.vline   = 3,
    col.legend  = 1,            # 0 for no legend
    legend.pos  = NULL,         # NULL means auto, else specify c(x,y)
    legend.text = NULL,         # vector of strings to use as legend text, NULL for auto
    do.par      = TRUE,         # call par() as appropriate

                                # par() settings, listed so caller can override our use
    main        = "Model Comparison",
    ...)
{
    do.legend <- function(xmax)
    {
        lty <- NULL
        col <- NULL
        if(is.null(legend.text))
            legend <- NULL
        else
            legend <- rep(legend.text, length.out=length(col.rsq) + length(col.grsq))
        maxchars <- 20
        args <- get.arg.strings(objects, maxchars)
        if(col.rsq[1] != 0) {       # RSq plotted?
            col <- c(col, col.rsq)
            lty <- c(lty, rep(2, length.out=length(col)))
            if(is.null(legend.text))
                for(imodel in seq_along(objects))
                    legend <- c(legend, paste(imodel, args[[imodel]]))
        }
        if(col.grsq[1] != 0) {      # GRSq plotted?
            col <- c(col, col.grsq)
            lty <- c(lty, rep(1, length.out=length(col)))
            if(is.null(legend.text))
                for(imodel in seq_along(objects))
                    legend <- c(legend, paste(imodel, args[[imodel]]))
        }
        len.legend <- 2 * maxchars * strwidth("x", "figure")
        cex1 <- 1
        if(is.null(legend.pos)) {
            ypos <- strheight("") * 1.1 * (length(legend) + 2)  # 2 for top and bot border
            xpos <- max(1.5, xmax / 3)
            # If legend is too big relative to figure, then reduce char size and move left.
            # This is just a hack that seems to work most of the time.
            if(len.legend > .9) {
                cex1 <- .7
                xpos <- max(1.5, xmax / 6)
                ypos <- ypos * .8
            } else if(len.legend > .6) {
                cex1 <- .8
                xpos <- max(1.5, xmax / 6)
                ypos <- ypos * .8
            }
        } else {                    # user specified legend position
            if(length(legend.pos) < 2)
                stop1("length(legend.pos) < 2")
            xpos = legend.pos[1]
            ypos = legend.pos[2]
            if(xpos < 0 || xpos > xmax || ypos < 0 || ypos > 1.1)
                warning1("out of range legend.pos")
            if(len.legend > .6)
                cex1 <- .8
        }
        legend(x=xpos, y=ypos, bg="white", legend=legend, col=col, lty=lty, cex=cex1)
    }
    # plot.earth.models starts here
    objects <- x
    if(!is.list(objects))       # note that is.list returns TRUE for a single object
        stop1("'x' is not an \"earth\" object or a list of \"earth\" objects")
    if(typeof(objects[[1]]) != "list") # if user specified just one object, convert to list
        objects <- list(objects)
    for(imodel in seq_along(objects))
        plot.earth.prolog(objects[[imodel]], paste("objects[[", imodel, "]]", sep=""))
    check.index.vec("which", which, 1:2)
    show <- rep(FALSE, 2)
    show[which] <- TRUE
    nfigs <- sum(show)
    if(nfigs == 0) {
        warning1("nothing to plot")
        return(invisible())
    }
    if(is.null(col.rsq))
        col.rsq <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(is.null(col.npreds))
        col.npreds <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(is.null(col.cum))
        col.cum <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(show[1] && col.grsq[1] == 0 && col.rsq[1] == 0)
        stop1("both col.grsq[1] and col.rsq[1] are zero")
    if(show[2] && is.null(col.cum))
        stop1("col.cum is NULL, and unable to use col.grsq or col.rsq instead")
    nmodels <- length(objects)
    col.grsq   <- rep(col.grsq,   length.out=nmodels)
    col.rsq    <- rep(col.rsq,    length.out=nmodels)
    col.npreds <- rep(col.npreds, length.out=nmodels)
    col.cum    <- rep(col.cum,    length.out=nmodels)
    col.vline  <- rep(col.vline,  length.out=nmodels)
    lty.vline  <- rep(lty.vline,  length.out=nmodels)
    must.trim.caption <- is.null(caption)
    caption <- get.caption.from.call(caption, objects[1])
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        on.exit(par(old.par))
        if(nfigs != 1)
            par(mfrow=c(2, 2))
        par(mar = c(4, 4, 2, 3))    # small margins and text
        par(mgp = c(1.6, 0.6, 0))   # flatten axis elements
        par(cex = 0.7)
        make.space.for.caption(caption)
    }
    max.npreds <- 1
    max.nterms <- 1
    for(imodel in seq_along(objects)) {
        object <- objects[[imodel]]
        rlim <- range(rlim, get.rlim(object, rlim, col.grsq[imodel], col.rsq[imodel]))
        max.npreds <- max(max.npreds,
              get.nused.preds.per.subset(object$dirs, object$prune.terms))
        max.nterms <- max(max.nterms, length(object$rss.per.subset))
    }
    if(show[1]) {
        for(imodel in seq_along(objects))
            plot.earth.model(objects[[imodel]],
                col.grsq    = col.grsq[imodel],
                col.rsq     = col.rsq[imodel],
                lty.rsq     = 2,
                col.npreds  = col.npreds[imodel],
                col.vline   = col.vline[imodel],
                lty.vline   = lty.vline[imodel],
                col.vseg    = col.grsq[imodel],
                col.legend  = 0,                    # we plot our own legend
                rlim        = rlim,
                add         = (imodel > 1),
                do.par      = do.par,
                max.nterms  = max.nterms,
                max.npreds  = max.npreds,
                jitter      = jitter,
                legend.pos  = NULL,
                nresponses  = NCOL(object[[imodel]]$residuals),
                main        = if(imodel > 1) "" else main,
                ...)
        if(col.legend != 0 && length(objects) > 1 && !show[2])
            do.legend(max.nterms)
    }
    if(show[2]) {
        multiple.responses <- FALSE
        xlim = c(0,0)
        for(object in objects)
            if(NCOL(object$residuals) > 1) {
                multiple.responses <- TRUE
                xlim = range(xlim, abs(object$residuals[,1]))
            } else
                xlim = range(xlim, abs(object$residuals))
        for(imodel in seq_along(objects))
            plot.cum(
                residuals = if(NCOL(objects[[imodel]]$residuals) > 1)
                                objects[[imodel]]$residuals[,1]
                            else
                                objects[[imodel]]$residuals,
                main      = if(imodel > 1)              ""
                            else if(multiple.responses) "Cumul Distrib (response 1)"
                            else                        "Cumulative Distribution",
                xlim      = xlim,
                col       = if(length(col.cum) > 1)         col.cum[imodel]
                            else if(col.grsq[imodel] != 0)  col.grsq[imodel]
                            else                            col.rsq[imodel],
                col.grid  = 0,
                cum.grid  = "none",
                add       = (imodel > 1),
                pch       = 20,
                ...)
        if(col.legend != 0 && length(objects) > 1)
            do.legend(xmax=xlim[2])
        }
    show.caption(caption, trim=must.trim.caption)
    invisible()
}

#--------------------------------------------------------------------------------------------
plot.earth.model <- function(   # show prune results and cumul distribution
    object,
    col.grsq,
    col.rsq,
    lty.rsq,
    col.npreds,
    col.vline,
    lty.vline,
    col.vseg,
    col.legend,
    rlim,
    add,
    do.par,
    max.nterms,
    max.npreds,
    jitter,         # allows overlaid plots to be visible
    legend.pos,     # NULL means auto, else c(x,y)
    nresponses,     # ncol(y)
    main,           # par() settings
    ...)
{
    scale1 <- function(x, Min, Max)
    {
        return((x-Min)/(Max-Min))
    }
    do.legend <- function()
    {
        lty <- 1
        col <- 1
        legend <- NULL
        if(col.grsq != 0)
            legend <- "GRSq"
        if(col.rsq != 0) {
            # figure out if grsq plot obscures rsq plot so legend can say so

            i <- rsqVec >= rlim[1] & rsqVec <= rlim[2]
            nobscured <- sum(abs(rsqVec[i] - grsqVec[i]) <
                                (rlim[2] - rlim[1]) / 100)
            if(nobscured > .8 * sum(i))
                legend <- c(legend, "RSq (obscured)")
            else
                legend <- c(legend, "RSq")
            col <- c(col, col.rsq)
            lty <- c(lty, 1)
        }
        if(col.npreds != 0) {
            legend <- c(legend, "Nbr preds")
            col <- c(col, col.npreds)
            lty <- c(lty, 3)
        }
        cex1 <- .8
        if(is.null(legend.pos)) {
            xpos <- max(1.5, max.nterms/2.5)
            ypos <- cex1 * strheight("") * 1.2 * (length(legend)+2)  # 2 for top + bot border
        } else {
            if(length(legend.pos) < 2)
                stop1("length(legend.pos) < 2")
            xpos = legend.pos[1]
            ypos = legend.pos[2]
            if(xpos < 0 || xpos > max.nterms || ypos < 0 || ypos > 1.1)
                warning1("out of range legend.pos")
        }
        legend(x=xpos, y=ypos, bg="white", legend=legend, col=col, lty=lty, cex=cex1)
    }
    left.axis <- function()
    {
        Pretty <- pretty(c(rlim[1], rlim[2]))
        if(length(Pretty) >= 2) # test to prevent "no locations are finite" error
            axis(side=2, at=scale1(Pretty, rlim[1], rlim[2]), lab=Pretty, srt=90)
        if(col.rsq != 0)        #$$ mtext needs cex=par("cex"), not sure why
            mtext("GRSq  RSq", side=2, line=1.6, cex=par("cex"))
        else
            mtext("GRSq", side=2, line=1.6, cex=par("cex"))
    }
    right.axis <- function()
    {
        Pretty <- pretty(c(0, max.npreds))
        if(length(Pretty) >= 2) # test to prevent "no locations are finite" error
            axis(side=4, at=scale1(Pretty, 0, max.npreds), lab=Pretty, srt=90)
        mtext("Number of used predictors", side=4, line=1.6, cex=par("cex"))
    }
    plot.nused.preds <- function()  # plot nbr of used predictors
    {                               # nothing actually plotted if col.npreds=0
        nused.preds <- get.nused.preds.per.subset(object$dirs, object$prune.terms)
        nused.preds.vec <- scale1(nused.preds, 0, max.npreds)
        if(jitter > 0)      # 2*jitter seems to work better relative to jitter on GRSq
            nused.preds.vec <- jitter(nused.preds.vec, amount=2*jitter)
        lines(nused.preds.vec, type="l", lty=3, col=col.npreds)
        NULL
    }
    # plot.earth.model starts here
    plot.earth.prolog(object, deparse(substitute(object)))
    warn.if.dots.used("plot.earth.model", ...)
    if(is.null(main)) {
        if(nresponses > 1)
            main <- "Model Selection (all responses)"
        else
            main <- "Model Selection"
    }
    if(is.null(object$prune.terms)) {       # no prune data?
        if(!add)
            plot(c(0,1), col=0, xlab="", ylab="")
        legend(x=1, y=1, bty="n", legend=c("No model selection data", "",
                "Run update.earth() to generate", "model selection data"))
        return(NULL)
    }
    if(jitter > 0.05)
        stop1("'jitter' ", jitter , " is too big, try something like jitter=0.01")
    if(!add) {
        if(do.par && col.npreds != 0)
            make.space.for.right.axis()
        # set up so vertical scale is 0..1, horizontal is 0..max.nterms
        plot(0:max.nterms, (0:max.nterms)/max.nterms,
            type="n", main=main, xlab="Number of terms", ylab="", yaxt="n")
        left.axis()
        if(col.npreds != 0)
            right.axis()
    }
    plot.nused.preds()

    # plot vertical line at best model
    abline(v=length(object$selected.terms), col=col.vline, lty=lty.vline)

    # plot a colored marker at the top of the above line
    points(x=length(object$selected.terms), y=1.02, col=col.vseg, pch=6)

    rsqVec  <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
    if(jitter > 0)
        rsqVec  <- jitter(rsqVec, amount=jitter)
    lines(scale1(rsqVec,  rlim[1], rlim[2]), col=col.rsq, lty=lty.rsq)

    grsqVec <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
    if(jitter > 0)
        grsqVec <- jitter(grsqVec, amount=jitter)
    lines(scale1(grsqVec, rlim[1], rlim[2]), col=col.grsq)

    if(col.legend != 0)
        do.legend()
    NULL
}

#--------------------------------------------------------------------------------------------
plot.cum <- function(           # plot cumulative distribution of absolute residuals
    residuals,
    main,
    xlim,
    col,
    col.grid,
    cum.grid,
    add,
    ...)
{
    if(is.null(main))
        main <- "Cumulative Distribution"
    abs.residuals <- abs(residuals)
    cum <- ecdf(abs.residuals)
    # col.points=0 gives a finer resolution graph (points are quite big regardless of pch)
    plot.stepfun(cum, add=add, main=main, xlab="abs(Residuals)", ylab="Proportion",
            xlim=xlim, col.points=0, col.hor=col, col.vert=col, ...)
    if(col.grid != 0 && !add) {
        ngrid <- match.choices(cum.grid[1], c("none", "grid", "percentages"), "cum.grid")
        if(ngrid >= 2) {
            # add annotated grid lines, unattractive but useful
            for(h in c(0,.05,.10,.25,.5,.75,.90,.95,1)) # horizontal lines
                abline(h=h, lty=1, col=col.grid)
            q <- quantile(abs.residuals)
            for(v in q)    # vertical lines at 0,25,50,75,100% quantiles
                abline(v=v, lty=1, col=col.grid)
            if(ngrid >= 3) {
                cex1 <- .7
                text(x=0,    y=1.02, "0%",   offset=0, cex=cex1)
                text(x=q[2], y=1.02, "25%",  offset=0, cex=cex1)
                text(x=q[3], y=1.02, "50%",  offset=0, cex=cex1)
                text(x=q[4], y=1.02, "75%",  offset=0, cex=cex1)
                text(x=q[5], y=1.02, "100%", offset=0, cex=cex1)
            }
            plot.stepfun(cum, add=TRUE, verticals=TRUE, # replot data over grid
                xlim=xlim, col.points=0, col.hor=col, col.vert=col, ...)
        }
    }
    NULL
}

#--------------------------------------------------------------------------------------------
# check rlim specified by user, and convert special values in rlim to actual vals

get.rlim <- function(object, rlim, col.grsq, col.rsq)
{
    if(length(rlim) != 2)
        stop1("length(rlim) != 2")
    if(rlim[2] <= rlim[1] && rlim[2] >= 0)
        stop1("rlim[2] <= rlim[1]")
    if(rlim[1] < -1 || rlim[1] >  1 || rlim[2] < -1 || rlim[2] >  1)
        stop1(paste(
            "Bad 'rlim' c(", rlim[1], ",", rlim[2], ")\n",
            "Legal settings are from 0 to 1, with special values:\n",
            "  rlim[1]=-1 means use min(RSq) or min(GRSq) excluding intercept)\n",
            "  rlim[2]=-1 means use max(RSq) or max(GRSq)\n", sep=""))

    if(rlim[1] < 0 || rlim[2] < 0) {
        grsq <- NULL
        if(col.grsq != 0)
            grsq <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
        rsq <- NULL
        if(col.rsq != 0)
            rsq <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
        if(rlim[1] < 0)
            rlim[1] <- min(grsq[-1], rsq[-1])
        if(rlim[2] < 0)
            rlim[2] <- max(grsq, rsq)

    }
    rlim
}

#--------------------------------------------------------------------------------------------
plot.earth.prolog <- function(object, object.name)
{
    check.classname(object, object.name, "earth")
    if(is.null(object$selected.terms))
        stop1("cannot plot because $selected.terms component is NULL")
    if(length(object$selected.terms) < 2)
        stop1("cannot plot because there is only an intercept term")
    NULL
}

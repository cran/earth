# plot.earth.R: plotting routines for the earth package
#
# Stephen Milborrow Mar 2007 Petaluma
#
# TODO allow model comparison using newdata instead of data used to build model
# TODO add nresponse to plot.earth.models

plot.earth <- function(
    x           = stop("no 'x' arg"),
    which       = 1:4,
    info        = FALSE,
    delever     = FALSE,
    pearson     = FALSE,
    level       = 0,
    versus      = NULL,

    nresponse   = 1,
    npoints     = 1000,
    id.n        = 3,
    labels.id   = rownames(residuals(object, warn=FALSE)),
    center      = TRUE,
    loess.f     = .5,

    do.par      = length(which) > 1,
    xlim        = NULL,
    ylim        = NULL,

    main        = NULL,
    cex.main    = 1.2,
    caption     = if(do.par) NULL else "", # NULL for auto
    xlab        = NULL,
    ylab        = NULL,

    pch         = 20,
    col.line    = "lightblue",
    col.loess   = NULL, # NULL for auto, col.line unless regression line also plotted
    lwd.loess   = NULL, # NULL for auto
    col.cv      = "red",
    col.qq      = col.line,
    col.grid    = "lightgray",
    col.points  = 1,
    cex.points  = NULL,
    shade.pints = "mistyrose2",
    shade.cints = "mistyrose4",

    cum.grid    = "percentages",

    col.rsq       = NA, # deprecated, use col.line instead
    col.residuals = NA, # deprecated, use col.points instead
    nresiduals    = NA, # deprecated, use npoints instead

    # following passed to plot.model.selection

    legend.pos          = NULL, # NULL means auto, NA means no legend
    cex.legend          = NULL, # NULL means auto
    col.grsq            = 1,
    col.infold.rsq      = 0,
    col.mean.infold.rsq = 0,
    col.mean.oof.rsq    = "palevioletred",
    col.npreds          = if(is.null(object$cv.oof.rsq.tab)) 1 else 0,
    col.oof.labs        = 0,
    col.oof.rsq         = "mistyrose2",
    col.oof.vline       = col.mean.oof.rsq,
    col.pch.cv.rsq      = 0,
    col.pch.max.oof.rsq = 0,
    col.sel.grid        = 0,
    col.vline           = col.grsq,
    col.vseg            = 0,
    lty.grsq            = 1,
    lty.npreds          = 2,
    lty.rsq             = 5,
    lty.vline           = 3,
    col.legend          = NA, # deprecated, use legend.pos=NA for no legend

    ...) # unused
{
    object.name <- substr(deparse(substitute(x)), 1, 40)
    object <- x # minimize confusion with x, the regression input matrix
    remove(x)   # not necessary but prevents mistakes later
    plot.earth.prolog(object, object.name)
    stop.if.dots.used("plot.earth", ...)
    col.line   <- check.deprecated(col.line, missing(col.line), col.rsq)
    col.points <- check.deprecated(col.points, missing(col.points), col.residuals)
    npoints    <- check.deprecated(npoints, missing(npoints), nresiduals)
    plotmo::check.index(which, "which", 1:8)
    if(is.null(caption) || any(nchar(caption))) { # show caption?
        trim.caption <- is.null(caption) # trim if auto-generate the caption
        caption <- get.caption(which, object, nresponse, object.name, caption)
    }
    env <- get.model.env(object, parent.frame)
    stopifnot(length(do.par) == 1)
    stopifnot(do.par == 0 || do.par == 1 || do.par == 2)
    plotmor.do.par <- do.par
    if(is.null(versus)) {
        plotmor.do.par <- FALSE # we do.par here
        old.par <- par.for.plot(do.par, length(which), cex.main, caption)
        if(do.par == 1)
            on.exit(par(old.par))
    } else { # versus was specified
        # remove all from which except standard resid and abs resid plots
        which <- which[which %in% c(3,5)]
        if(length(which) == 0) {
            warning0("plot.earth: nothing to plot\n",
                     "(the \"which\" argument is empty after removing certain ",
                     "plots because \"versus\" was specified)")
            return(invisible())
        }
        do.par <- FALSE # plotmor will do.par
    }
    which1 <- which(which == 1)
    if(any(which1)) {
        stopifnot(inherits(object, "earth"))
        plot.model.selection.wrapper(object,
            do.par, xlim, ylim, main, col.line,
            legend.pos, cex.legend,
            col.grsq, col.infold.rsq, col.mean.infold.rsq, col.mean.oof.rsq,
            col.npreds, col.oof.labs, col.oof.rsq, col.oof.vline,
            col.pch.cv.rsq, col.pch.max.oof.rsq, col.sel.grid, col.vline,
            col.vseg, lty.grsq, lty.npreds, lty.rsq, lty.vline, col.legend)
        which <- which[which != 1]
    }
    if(any(which > 1))
        plotmor(
            object,
            which,
            info,
            delever,
            pearson,
            level,
            versus,

            nresponse,
            npoints,
            id.n,
            labels.id,
            center,
            loess.f,

            env,
            do.par = plotmor.do.par,
            xlim   = if(any(which1)) NULL else xlim,
            ylim   = if(any(which1)) NULL else ylim,
            main,
            cex.main,
            caption = caption,
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

            cum.grid,

            ...)

    if(do.par && any(nchar(caption))) # show caption?
        show.caption(caption, trim=trim.caption)
    invisible()
}
plot.earth.prolog <- function(object, object.name)
{
    check.classname(object, object.name, "earth")
    if(is.null(object$selected.terms))
        stop0(object.name, " has no $selected.terms field")
    if(is.null(object$fitted.values))
        stop0(object.name, " has no $fitted.values field.\n",
              "       Use keepxy=TRUE in the call to earth.")
    if(is.null(object$residuals)) # probably a model from object$cvlist
        stop0(object.name, " has no $residuals field.\n",
              "       Use keepxy=TRUE in the call to earth.")
}
get.earth.legend.cex <- function(legend.text, min.width=.4, min.cex=.4)
{
    longest.text <- legend.text[which.max(strwidth(legend.text))]
    longest.text <- paste0("AAAAAA ", longest.text) # incorporate line on left of legend
    # reduce cex.legend until legend fits, but not more than to min.cex
    cex <- .8
    while((width <- max(strwidth(longest.text, units="figure", cex=cex))) > min.width &&
            cex > min.cex)
        cex <- cex - .1
    cex
}
plot.earth.models <- function(
    x            = stop("no 'x' arg"),
    which        = c(1:2),
    caption      = "",
    jitter       = 0,
    col.grsq     = discrete.plot.cols(length(x)),
    lty.grsq     = 1,
    col.line     = 0,
    lty.rsq      = 5,
    col.vline    = col.grsq,
    lty.vline    = 3,
    col.npreds   = 0,
    lty.npreds   = 2,
    col.sel.grid = 0,
    ylim         = c(0,1),
    col.legend   = 1,
    cex.legend   = NULL,
    legend.pos   = NULL,
    legend.text  = NULL,
    col.cum      = NULL,
    do.par       = TRUE,
    main         = "Model Comparison",
    cex.main     = 1.2,
    ...)
{
    warn.if.dots.used("plot.earth.models", ...)
    objects <- x
    if(!is.list(objects))       # note that is.list returns TRUE for a single object
        stop0("\"x\" is not an \"earth\" object or a list of \"earth\" objects")
    # check for a common error, using plot.earth.models(mod1, mod2) instead
    # of plot.earth.models(list(mod1, mod2)) instead
    if(inherits(which, "earth"))
        stop0("use plot.earth.models(list(model1, model2)), ",
              "not plot.earth.models(model1, model2)")
    if(typeof(objects[[1]]) != "list") # if user specified just one object, convert to list
        objects <- list(objects)
    for(imodel in seq_along(objects))
        plot.earth.prolog(objects[[imodel]], paste0("objects[[", imodel, "]]"))
    plotmo::check.index(which, "which", 1:2)
    show <- to.logical(which, 4)
    if(length(which) == 0) {
        warning0("plot.earth.models: nothing to plot (the \"which\" argument is empty)")
        return(invisible())
    }
    if(is.null(col.line))
        col.line <- if(is.null(col.grsq)) col.line else col.grsq
    if(is.null(col.npreds))
        col.npreds <- if(is.null(col.grsq)) col.line else col.grsq
    if(is.null(col.cum))
        col.cum <- if(is.null(col.grsq)) col.line else col.grsq
    if(show[1] && col.grsq[1] == 0 && col.line[1] == 0)
        stop0("both col.grsq[1] and col.line[1] are zero")
    if(show[2] && is.null(col.cum))
        stop0("col.cum is NULL, and unable to use col.grsq or col.line instead")
    nmodels <- length(objects)
    col.grsq   <- repl(col.grsq, nmodels)
    lty.grsq   <- repl(lty.grsq, nmodels)
    col.line   <- repl(col.line, nmodels)
    lty.rsq    <- repl(lty.rsq, nmodels)
    col.npreds <- repl(col.npreds, nmodels)
    lty.npreds <- repl(lty.npreds, nmodels)
    col.cum    <- repl(col.cum, nmodels)
    col.vline  <- repl(col.vline, nmodels)
    lty.vline  <- repl(lty.vline, nmodels)
    stopifnot(length(do.par) == 1)
    stopifnot(do.par == 0 || do.par == 1 || do.par == 2)
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        if(do.par == 1)
            on.exit(par(old.par))
        if(length(which) > 1)
            par(mfrow=c(2, 2))
        par(mar = c(4, 4, 2, 3))        # small margins to pack figs in
        if(length(which) == 1 && which == 1 && is.specified(col.npreds))
            par(mar = c(4, 4, 2, 4))    # space for right axis
        par(mgp = c(1.6, 0.6, 0))       # flatten axis elements
        par(cex.main = cex.main)
        make.space.for.caption()
    }
    max.npreds <- 1
    max.nterms <- 1
    for(imodel in seq_along(objects)) {
        object <- objects[[imodel]]
        ylim <- range(ylim,
                      get.model.selection.ylim(object, ylim, col.grsq[imodel], col.line[imodel]))
        max.npreds <- max(max.npreds,
                          get.nused.preds.per.subset(object$dirs, object$prune.terms))
        max.nterms <- max(max.nterms, length(object$rss.per.subset))
    }
    if(show[1]) {
        for(imodel in seq_along(objects))
            plot.model.selection(
                object              = objects[[imodel]],

                do.par              = do.par,
                xlim                = NULL,
                ylim                = ylim,
                main                = if(imodel > 1) "" else main,

                col.line            = col.line[imodel],

                legend.pos          = NA, # we plot our own legend
                cex.legend          = cex.legend,

                col.grsq            = col.grsq[imodel],
                col.infold.rsq      = 0,
                col.mean.infold.rsq = 0,
                col.mean.oof.rsq    = 0,
                col.npreds          = col.npreds[imodel],
                col.oof.labs        = 0,
                col.oof.rsq         = 0,
                col.oof.vline       = 0,
                col.pch.cv.rsq      = 0,
                col.pch.max.oof.rsq = 0,
                col.sel.grid        = col.sel.grid,
                col.vline           = col.vline[imodel],
                col.vseg            = col.grsq[imodel],
                lty.grsq            = lty.grsq[imodel],
                lty.npreds          = lty.npreds[imodel],
                lty.rsq             = lty.rsq,
                lty.vline           = lty.vline[imodel],

                add                 = (imodel > 1),
                nresp               = NCOL(object[[imodel]]$residuals),
                max.nterms          = max.nterms,
                max.npreds          = max.npreds,
                jitter              = if(imodel>1) jitter else 0)

        if(is.specified(col.legend) && length(objects) > 1 && !show[2])
            plot.earth.models.legend(objects, min.width=.4,
                legend.text, legend.pos, cex.legend, col.legend,
                col.line, lty.rsq, col.grsq, lty.grsq)
    }
    if(show[2]) {
        multiple.responses <- FALSE
        xlim <- c(0,0)
        for(object in objects) {
            stopifnot(!is.null(object$residuals)) # already checked in plot.earth.prolog
            if(NCOL(object$residuals) > 1) {
                multiple.responses <- TRUE
                xlim <- range(xlim, abs(object$residuals[,1]))
            } else
                xlim <- range(xlim, abs(object$residuals))
        }
        for(imodel in seq_along(objects)) {
            rinfo <- get.rinfo(objects[[imodel]], 1, FALSE, FALSE)
            plot.cum(
                rinfo    = rinfo,
                main     = if(imodel > 1)              ""
                           else if(multiple.responses) "Cumul Distrib (response 1)"
                           else                        "Cumulative Distribution",
                xlim     = xlim,
                col      = if(length(col.cum) > 1)                 col.cum[imodel]
                           else if(is.specified(col.grsq[imodel])) col.grsq[imodel]
                           else                                    col.line[imodel],
                col.grid = 0,
                cum.grid = "none",
                add      = (imodel > 1),
                jitter   = if(imodel == 1) 0 else jitter)
        }
        if(is.specified(col.legend) && length(objects) > 1)
            plot.earth.models.legend(objects, min.width=.5,
                legend.text, legend.pos, cex.legend, col.legend,
                col.line, lty.rsq, col.grsq, lty.grsq)
    }
    show.caption(get.caption.from.call(caption, objects[1]))
    invisible()
}
plot.earth.models.legend <- function(
    objects,
    min.width,
    legend.text,
    legend.pos,
    cex.legend,
    col.legend,
    col.line,
    lty.rsq,
    col.grsq,
    lty.grsq)
{
    lty <- NULL
    col <- NULL
    if(is.null(legend.text)) {
        if(is.null(names(objects))) {
            args <- get.arg.strings(objects, maxchars=20)
            legend.text <- character(length=length(objects))
            for(imodel in seq_along(objects))
                legend.text[imodel] <- paste(imodel, args[[imodel]])
        } else
            legend.text <- names(objects)
    } else
        legend.text <- repl(legend.text, length(objects))
     if(col.line[1] != 0) {       # RSq plotted?
        col <- c(col, col.line)
        lty <- c(lty, repl(lty.rsq, length(col)))
        if(col.grsq[1] != 0)
            legend1 <- paste("RSq", legend.text)
    }
    if(col.grsq[1] != 0) {      # GRSq plotted?
        col <- c(col, col.grsq)
        lty <- c(lty, repl(lty.grsq, length(col)))
        if(col.line[1] != 0)
            legend.text <- c(legend1, paste("GRSq", legend.text))
    }
    if(is.null(cex.legend))
        cex.legend <- get.earth.legend.cex(legend.text, min.width=min.width)
    if(is.null(legend.pos))
        xpos <- "bottomright"
    else { # user specified legend position
        xpos <- legend.pos[1]
        ypos <- if(length(legend.pos) > 1) legend.pos[2] else NULL
    }
    elegend(x=xpos, y=ypos, bg="white",
            legend=legend.text, col=col, lty=lty, cex=cex.legend,
            # y offset allows vertical lines to be visible below legend
            inset=c(.02, .04))
}

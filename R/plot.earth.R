# plot.earth.R: plotting routines for the earth package
#
# Stephen Milborrow Mar 2007 Petaluma
#
# TODO allow model comparison using newdata instead of data used to build model
# TODO add an x range option to plot.model.selection
# TODO add nresponse to plot.earth.models

plot.earth <- function(
    x                 = stop("no 'x' arg"),
    which             = 1:4,
    nresponse         = 1,
    caption           = if(do.par) NULL else "",
    col.grsq          = 1,
    lty.grsq          = 1,
    col.rsq           = "lightblue",
    lty.rsq           = 5,
    col.vline         = col.grsq,
    lty.vline         = 3,
    col.npreds        = if(is.null(x$cv.oof.rsq.tab)) 1 else 0,
    lty.npreds        = 2,
    col.mean.oof.rsq  = "palevioletred",
    col.oof.rsq       = "mistyrose2",
    col.oof.vline     = col.mean.oof.rsq,
    col.oof.labs      = 0,
    col.pch.max.oof.rsq = 0,
    col.pch.cv.rsq      = 0,
    col.mean.infold.rsq = 0,
    col.infold.rsq      = 0,
    col.sel.grid      = 0,
    ylim              = c(-1,-1),
    col.legend        = 1,
    cex.legend        = NULL,
    legend.pos        = NULL,
    col.cum.grid      = "lightgray",
    cum.grid          = "percentages",
    id.n              = 3,
    labels.id         = rownames(residuals(x, warn=FALSE)),
    col.residuals     = 1,
    col.loess         = col.rsq,
    nresiduals        = 1000,
    col.qq            = col.rsq,
    do.par            = TRUE,
    main              = NULL,
    pch               = 1,
    rlim              = NA,
    col.grid          = NA,
    ...)
{
    get.iresiduals <- function(nresiduals, residuals, fitted.values)
    {
        if(any(!is.finite(residuals)))
            warning0("non finite residuals")
        stopifnot(nrow(residuals) == nrow(fitted.values))
        iresiduals <- 1:nrow(residuals)
        if(nresiduals > 0 && nresiduals < nrow(residuals)) {
            # take a sample, but make sure it includes the largest residuals

            nresiduals <- min(nresiduals, nrow(residuals))
            nlargest <- min(max(3, id.n), nresiduals)
            isorted <- order(abs(residuals), decreasing=TRUE)
            ikeep <- 1:nlargest
            if(nresiduals > nlargest)
                ikeep <- c(ikeep, seq(from=nlargest + 1, to=nrow(residuals),
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
            stop0("negative 'id.n'")
        ncases <- nrow(object$residuals)
        if(id.n > ncases)
            id.n <- ncases
        if(id.n > 0)
            sort.list(abs(residuals), decreasing = TRUE)[1:id.n]
        else
            NULL
    }
    plot.residuals <- function(fitted.values, residuals, col.residuals) # residuals versus fitted
    {
        plot.loess <- function()
        {
            if(is.naz(col.loess))
                return()
            # We use lowess rather than loess because loess tends to give ugly warnings,
            # but the argument is named "col.loess" for backward compatibility.
            # For less smoothing (so we can better judge inflection points),
            # we use a value for f lower than the default 2/3. Also iter=0 is
            # best for lowess with binary responses, so says Harrell 2.4.6.
            smooth <- lowess(fitted.values, residuals, f=.5,
                             iter=if(length(unique(residuals)) > 2) 3 else 0)
            lines(smooth$x, smooth$y, col=col.loess, ...)
        }
        if(is.null(main))
            main="Residuals vs Fitted"
        plot(fitted.values, residuals, main=main,
             xlab="Fitted", ylab="Residuals", pch=pch, col=col.residuals, ...)
        abline(h=0, lty=3, col=col.cum.grid)
        plot.loess()
        if(!is.null(id.indices))
            thigmophobe.labels(x=fitted.values[id.indices], y=residuals[id.indices],
                labels=labels.id[id.indices],
                offset=.33, font=2, cex=.8, col="steelblue4", xpd=NA)
    }
    plot.qq <- function(residuals)  # qqplot of residuals
    {
        par(pty="s")                # square
        if(is.null(main))
            main="Normal Q-Q"
        qq <- qqnorm(residuals, main=main,
                xlab="Theoretical Quantiles", ylab="Residual Quantiles",
                pch=pch, col=col.residuals, ...)
        qqline(residuals, col=col.qq, ...)
        if(!is.null(id.indices))
            thigmophobe.labels(x=qq$x[id.indices], y=qq$y[id.indices],
                 labels=labels.id[id.indices],
                 offset=.33, font=2, cex=.8, col="steelblue4", xpd=NA)
        par(pty="m")                # back to maximal (the default)
    }
    # return the response.name (for prepending to caption), if appropriate
    get.caption.prefix <- function(show, residuals, nresponse, caption)
    {
        if(!is.null(caption))
            return(NULL)    # don't modify caption explictly set by user
        if(all(show == c(TRUE, FALSE, FALSE, FALSE)) && NCOL(residuals) > 1)
            return(NULL)    # Model Selection graph only, no multi-resp prefix
        colnames <- colnames(residuals)
        if(!is.null(colnames) && !is.null(colnames[nresponse]) &&
                !is.na(colnames[nresponse]) && colnames[nresponse] != "")
            colnames[nresponse]
        else if(NCOL(residuals) > 1)
            paste("Response", nresponse)
        else
            NULL
    }
    #--- plot.earth starts here ---
    object <- x
    object.name <- substr(deparse(substitute(x)), 1, 40)
    plot.earth.prolog(object, object.name)
    if(length(rlim) != 1 || !is.na(rlim)) {
        warning0("rlim is deprecated.  Please use ylim instead.")
        if(!missing(ylim))
            stop0("rlim and ylim both specified.  Please use just ylim.")
    }
    rlim <- ylim
    if(length(col.grid) != 1 || !is.na(col.grid)) {
        warning0("col.grid is deprecated.  Please use col.cum.grid instead.")
        if(!missing(col.cum.grid))
            stop0("col.grid and col.cum.grid both specified.  Please use just col.cum.grid.")
        col.cum.grid <- col.grid
    }
    check.index.vec("which", which, 1:4)
    show <- to.logical(which, 4)
    nfigs <- sum(show)
    if(nfigs == 0) {
        warning0("plot.earth: nothing to plot")
        return(invisible())
    }
    must.trim.caption <- is.null(caption)
    check.index.vec("nresponse", nresponse, object$residuals,
                    check.empty = TRUE, use.as.col.index=TRUE)
    caption.prefix <- get.caption.prefix(show, object$residuals, nresponse, caption)
    if(!is.null(object$ifold)) # object is one fold of a cross-validated model?
        caption <- object.name
    else
        caption <- get.caption.from.call(caption, object)
    if(!is.null(caption.prefix))
        caption <- paste0(caption.prefix, ": ", caption)
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        on.exit(par(old.par))
        nrows <- ceiling(sqrt(nfigs))
        par(mfrow=c(nrows, nrows))
        par(mar = c(4, 4, 2, 1))        # small margins to pack figs in
        if(nfigs == 1 && which[1] == 1 && !is.naz(col.npreds))
            par(mar = c(4, 4, 2, 4))    # space for right axis
        par(mgp = c(1.6, 0.6, 0))       # flatten axis elements
        make.space.for.caption(caption)
    }
    nresp <- NCOL(object$residuals)
    if(nresp > 1) {
        # already called check.index.vec(nresponse) above, so no need to call again
        object$fitted.values <- object$fitted.values[, nresponse, drop=FALSE]
        object$residuals <- object$residuals[, nresponse, drop=FALSE]
    }
    if(show[1]) {
        rlim <- get.model.selection.ylim(object, rlim, 1, col.rsq,
                    col.mean.oof.rsq, col.oof.rsq, col.mean.infold.rsq, col.infold.rsq)
        plot.model.selection(object, col.grsq, lty.grsq, col.rsq, lty.rsq,
            col.mean.oof.rsq, col.oof.rsq, col.oof.vline, col.oof.labs, col.sel.grid,
            col.pch.max.oof.rsq, col.pch.cv.rsq, col.mean.infold.rsq, col.infold.rsq,
            col.npreds, lty.npreds, col.vline, lty.vline, col.vseg=0,
            col.legend, cex.legend, legend.pos, rlim, add=FALSE, do.par,
            max.nterms=length(object$rss.per.subset),
            max.npreds=max(1, get.nused.preds.per.subset(object$dirs, object$prune.terms)),
            jitter=0, nresp,
            main, ...)
    }
    if(show[2])
        plot.cum(object$residuals, main, xlim=range(abs(object$residuals)),
                 col=1, col.cum.grid, cum.grid, add=FALSE, jitter=0, ...)
    if(show[3] || show[4]) {
        if(NROW(object$fitted.values) != nrow(object$residuals))
            stop0("NROW(object$fitted.values) ",  NROW(object$fitted.values),
                  " != nrow(object$residuals) ", nrow(object$residuals))
        iresiduals <- get.iresiduals(nresiduals, object$residuals, object$fitted.values)
        residuals <- object$residuals[iresiduals]
        fitted.values <- object$fitted.values[iresiduals]
        col.residuals <- rep(col.residuals, length.out=length(object$residuals)) # recycle
        col.residuals <- col.residuals[iresiduals]
        id.indices <- get.id.indices()
        if(is.null(labels.id))
            labels.id <- paste(iresiduals)
        else
            labels.id <- labels.id[iresiduals]
    }
    if(show[3])
        plot.residuals(fitted.values, residuals, col.residuals)
    if(show[4])
        plot.qq(residuals)
    show.caption(caption, trim=must.trim.caption)
    invisible()
}

plot.earth.models <- function(
    x            = stop("no 'x' arg"),
    which        = c(1:2),
    caption      = "",
    jitter       = 0,
    col.grsq     = discrete.plot.cols(length(x)),
    lty.grsq     = 1,
    col.rsq      = 0,
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
    rlim         = NA,
    ...)
{
    plot.legend <- function(min.width)
    {
        lty <- NULL
        col <- NULL
        legend <- rep(legend.text, length.out=length(objects)) # may be NULL
        if(is.null(legend)) {
            if(is.null(names(objects))) {
                args <- get.arg.strings(objects, maxchars=20)
                legend <- character(length=length(objects))
                for(imodel in seq_along(objects))
                    legend[imodel] <- paste(imodel, args[[imodel]])
            } else
                legend <- names(objects)
        }
        if(col.rsq[1] != 0) {       # RSq plotted?
            col <- c(col, col.rsq)
            lty <- c(lty, rep(lty.rsq, length.out=length(col)))
            if(col.grsq[1] != 0)
                legend1 <- paste("RSq", legend)
        }
        if(col.grsq[1] != 0) {      # GRSq plotted?
            col <- c(col, col.grsq)
            lty <- c(lty, rep(1, length.out=length(col)))
            if(col.rsq[1] != 0)
                legend <- c(legend1, paste("GRSq", legend))
        }
        if(is.null(cex.legend))
            cex.legend <- get.cex.legend(legend, min.width=min.width)
        if(is.null(legend.pos))
            xpos <- "bottomright"
        else { # user specified legend position
            xpos <- legend.pos[1]
            ypos <- if(length(legend.pos) > 1) legend.pos[2] else NULL
        }
        elegend(x=xpos, y=ypos, bg="white", legend=legend, col=col, lty=lty, cex=cex.legend,
                inset=c(.02, .04)) # y offset allows vertical lines to be visible below legend
    }
    #--- plot.earth.models starts here ---
    objects <- x
    if(!is.list(objects))       # note that is.list returns TRUE for a single object
        stop0("'x' is not an \"earth\" object or a list of \"earth\" objects")
    # check for a common error, using plot.earth.models(mod1, mod2) instead
    # of plot.earth.models(list(mod1, mod2)) instead
    if(inherits(which, "earth"))
        stop0("use plot.earth.models(list(model1, model2)), ",
              "not plot.earth.models(model1, model2)")
    if(typeof(objects[[1]]) != "list") # if user specified just one object, convert to list
        objects <- list(objects)
    for(imodel in seq_along(objects))
        plot.earth.prolog(objects[[imodel]], paste0("objects[[", imodel, "]]"))
    if(length(rlim) != 1 || !is.na(rlim)) {
        warning0("rlim is deprecated.  Please use ylim instead.")
        if(!missing(ylim))
            stop0("rlim and ylim both specified.  Please use just ylim.")
    }
    rlim <- ylim
    check.index.vec("which", which, 1:2)
    show <- to.logical(which, 4)
    nfigs <- sum(show)
    if(nfigs == 0) {
        warning0("plot.earth.models: nothing to plot")
        return(invisible())
    }
    if(is.null(col.rsq))
        col.rsq <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(is.null(col.npreds))
        col.npreds <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(is.null(col.cum))
        col.cum <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(show[1] && col.grsq[1] == 0 && col.rsq[1] == 0)
        stop0("both col.grsq[1] and col.rsq[1] are zero")
    if(show[2] && is.null(col.cum))
        stop0("col.cum is NULL, and unable to use col.grsq or col.rsq instead")
    nmodels <- length(objects)
    col.grsq   <- rep(col.grsq,   length.out=nmodels)
    lty.grsq   <- rep(lty.grsq,   length.out=nmodels)
    col.rsq    <- rep(col.rsq,    length.out=nmodels)
    lty.rsq    <- rep(lty.rsq,    length.out=nmodels)
    col.npreds <- rep(col.npreds, length.out=nmodels)
    lty.npreds <- rep(lty.npreds, length.out=nmodels)
    col.cum    <- rep(col.cum,    length.out=nmodels)
    col.vline  <- rep(col.vline,  length.out=nmodels)
    lty.vline  <- rep(lty.vline,  length.out=nmodels)
    must.trim.caption <- is.null(caption)
    caption <- get.caption.from.call(caption, objects[1])
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        on.exit(par(old.par))
        if(nfigs > 1)
            par(mfrow=c(2, 2))
        par(mar = c(4, 4, 2, 3))        # small margins to pack figs in
        if(nfigs == 1 && which[1] == 1 && !is.naz(col.npreds))
            par(mar = c(4, 4, 2, 4))    # space for right axis
        par(mgp = c(1.6, 0.6, 0))       # flatten axis elements
        make.space.for.caption(caption)
    }
    max.npreds <- 1
    max.nterms <- 1
    for(imodel in seq_along(objects)) {
        object <- objects[[imodel]]
        rlim <- range(rlim, get.model.selection.ylim(object, rlim, col.grsq[imodel], col.rsq[imodel]))
        max.npreds <- max(max.npreds,
                          get.nused.preds.per.subset(object$dirs, object$prune.terms))
        max.nterms <- max(max.nterms, length(object$rss.per.subset))
    }
    if(show[1]) {
        for(imodel in seq_along(objects))
            plot.model.selection(objects[[imodel]],
                col.grsq          = col.grsq[imodel],
                lty.grsq          = lty.grsq[imodel],
                col.rsq           = col.rsq[imodel],
                lty.rsq           = lty.rsq,
                col.mean.oof.rsq  = 0,
                col.oof.rsq       = 0,
                col.oof.vline     = 0,
                col.oof.labs      = 0,
                col.sel.grid      = col.sel.grid,
                col.pch.max.oof.rsq = 0,
                col.pch.cv.rsq      = 0,
                col.mean.infold.rsq = 0,
                col.infold.rsq      = 0,
                col.npreds        = col.npreds[imodel],
                lty.npreds        = lty.npreds[imodel],
                col.vline         = col.vline[imodel],
                lty.vline         = lty.vline[imodel],
                col.vseg          = col.grsq[imodel],
                col.legend        = 0,   # we plot our own legend
                cex.legend        = cex.legend,
                legend.pos        = NULL,
                rlim              = rlim,
                add               = (imodel > 1),
                do.par            = do.par,
                max.nterms        = max.nterms,
                max.npreds        = max.npreds,
                jitter            = if(imodel>1) jitter else 0,
                nresp             = NCOL(object[[imodel]]$residuals),
                main              = if(imodel > 1) "" else main,
                ...)
        if(!is.naz(col.legend) && length(objects) > 1 && !show[2])
            plot.legend(min.width=.4)
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
        for(imodel in seq_along(objects))
            plot.cum(
                residuals    = if(NCOL(objects[[imodel]]$residuals) > 1)
                                 objects[[imodel]]$residuals[,1]
                               else
                                 objects[[imodel]]$residuals,
                main         = if(imodel > 1)                ""
                                 else if(multiple.responses) "Cumul Distrib (response 1)"
                               else                          "Cumulative Distribution",
                xlim         = xlim,
                col          = if(length(col.cum) > 1)            col.cum[imodel]
                               else if(!is.naz(col.grsq[imodel])) col.grsq[imodel]
                               else                               col.rsq[imodel],
                col.cum.grid = 0,
                cum.grid     = "none",
                add          = (imodel > 1),
                jitter       = if(imodel == 1) 0 else jitter,
                pch          = 20,
                ...)
        if(!is.naz(col.legend) && length(objects) > 1)
            plot.legend(min.width=.5)
    }
    show.caption(caption, trim=must.trim.caption)
    invisible()
}

plot.cum <- function( # plot cumulative distribution of absolute residuals
    residuals,
    main,
    xlim,
    col,
    col.cum.grid,
    cum.grid,
    add,
    jitter,
    ...)
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
    if(is.null(main))
        main <- "Cumulative Distribution"
    abs.residuals <- abs(residuals)
    cum <- ecdf(abs.residuals)
    if(jitter > 0.05)
       stop0("'jitter' ", jitter , " is too big, try something like jitter=0.01")
    if(jitter > 0)
        environment(cum)$"y" <- jitter(environment(cum)$"y", amount=2 * jitter)
    # col.points=0 gives a finer resolution graph (points are quite big regardless of pch)
    plot.stepfun(cum, add=add, main=main, xlab="abs(Residuals)", ylab="Proportion",
                 xlim=xlim, col.points=0, col.hor=col, col.vert=col, ...)
    if(!is.naz(col.cum.grid) && !add) {
        ngrid <- match.choices(cum.grid[1], c("none", "grid", "percentages"), "cum.grid")
        if(ngrid >= 2) {
            # add annotated grid lines, unattractive but useful
            for(h in c(0,.25,.5,.75,.90,.95,1)) # horizontal lines
                abline(h=h, lty=1, col=col.cum.grid)
            q <- quantile(abs.residuals, probs=c(0, .25, .50, .75, .9, .95, 1))
            for(v in q)    # vertical lines at 0,25,50,75,90,95,100% quantiles
                abline(v=v, lty=1, col=col.cum.grid)
            if(ngrid >= 3)
                show.percents(q)
            plot.stepfun(cum, add=TRUE, verticals=TRUE, # replot data over grid
                xlim=xlim, col.points=0, col.hor=col, col.vert=col, ...)
        }
    }
}

plot.earth.prolog <- function(object, object.name)
{
    check.classname(object, object.name, "earth")
    if(is.null(object$selected.terms))
        stop0(object.name, " has no $selected.terms field")
    if(is.null(object$fitted.values))
        stop0(object.name, " has no $fitted.values field.  Use keepxy=TRUE?")
    if(is.null(object$residuals)) # probably a model from object$cvlist
        stop0(object.name, " has no $residuals field.  Use keepxy=TRUE?")
}

get.cex.legend <- function(legend.text, min.width=.4, min.cex=.4)
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

# plot.model.selection.R:
#
# Stephen Milborrow Jun 2011 Berea

plot.model.selection <- function(
    object,

    do.par,
    xlim,
    ylim,
    main,

    col.line            = "lightblue",

    legend.pos          = NULL, # NULL means auto, NA means no legend
    cex.legend          = NULL,

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

    add,
    nresp,   # ncol(y)
    max.nterms,
    max.npreds,
    jitter,         # allows overlaid plots to be visible

    ...)
{
    possibly.issue.cv.warning <- function()
    {
        if((!identical(col.mean.oof.rsq, "palevioletred") && !identical(col.mean.oof.rsq, 0)) ||
           (!identical(col.oof.rsq, "mistyrose2")         && !identical(col.oof.rsq, 0))      ||
           !identical(col.oof.labs, 0)        ||
           !identical(col.pch.max.oof.rsq, 0) ||
           !identical(col.pch.cv.rsq, 0)      ||
           !identical(col.mean.infold.rsq, 0) ||
           !identical(col.infold.rsq, 0)) {
            # user specifed a cross-validation argument, check that data is available
            if(is.null(object$cv.list))
                warning0("no cross-validation data because nfold not used in original call to earth")
            else if(is.null(object$cv.oof.rsq.tab))
                warning0("cannot plot cross-validation data because ",
                         "the earth model was not built with keepxy=TRUE")
        }
    }
    scale1 <- function(x, Min, Max)
    {
        return((x-Min)/(Max-Min))
    }
    left.axis <- function()
    {
        pretty <- pretty(c(ylim[1], ylim[2]))
        axis(side=2, at=scale1(pretty, ylim[1], ylim[2]), labels=pretty, srt=90)
        text <- ""
        if(is.specified(col.grsq))
            text <- "GRSq"
        if(is.specified(col.line) ||
                is.specified(col.oof.rsq) || is.specified(col.mean.oof.rsq) ||
                is.specified(col.infold.rsq) || is.specified(col.mean.infold.rsq))
            text <- paste0(text, "   RSq")
        # TODO mtext needs cex=par("cex"), not sure why
        mtext(text, side=2, line=2, cex=par("cex"))
    }
    right.axis <- function()
    {
        pretty <- pretty(c(0, max.npreds))
        axis(side=4, at=scale1(pretty, 0, max.npreds), labels=pretty, srt=90)
        mtext("Number of used predictors", side=4, line=1.6, cex=par("cex"))
    }
    plot.sel.grid <- function() # plot the grid
    {
        if(!is.specified(col.sel.grid))
            return()
        col <- col.sel.grid[1]
        abline(v=0:nterms.on.horiz.axis, col=col)
        for(v in seq(-1,1,by=if(abs(ylim[2]-ylim[1]) <= .5) .05 else .1))
            abline(h=scale1(v, ylim[1], ylim[2]), col=col)
    }
    plot.infold.rsqs <- function() # plot rsq's measured on the in-fold data
    {
        if(!is.specified(col.infold.rsq))
            return()
        # recycle col.infold.rsq so can use different colors for different folds
        col.infold.rsq <- repl(col.infold.rsq, length(object$cv.list))
        for(ifold in seq_along(object$cv.list)) {
            infold.rsq <- object$cv.infold.rsq.tab[ifold,]
            if(jitter > 0)
                infold.rsq  <- jitter(infold.rsq, amount=jitter)
            scaled.rsq <- scale1(infold.rsq,  ylim[1], ylim[2])
            lines(scaled.rsq, col=col.infold.rsq[ifold], lty=1)
        }
    }
    plot.oof.rsqs <- function() # plot rsq's measured on the out-of-fold data
    {
        if(!is.specified(col.oof.rsq))
            return()
        # recycle col.oof.rsq so can use different colors for different folds
        col.oof.rsq <- repl(col.oof.rsq, length(object$cv.list))
        if(is.specified(col.oof.labs))
            col.oof.labs <- repl(col.oof.labs, length(object$cv.list))
        for(ifold in seq_along(object$cv.list)) {
            oof.rsq <- object$cv.oof.rsq.tab[ifold,]
            if(jitter > 0)
                oof.rsq  <- jitter(oof.rsq, amount=jitter)
            scaled.rsq <- scale1(oof.rsq,  ylim[1], ylim[2])
            lines(scaled.rsq, col=col.oof.rsq[ifold], lty=1)
            if(is.specified(col.oof.labs)) {
                oof.rsq <- oof.rsq[!is.na(oof.rsq)] # truncate NAs
                text(length(oof.rsq)+.2, scaled.rsq[length(oof.rsq)],
                     substr(names(object$cv.list)[ifold], 5, 15),
                     cex=.6, col=col.oof.labs[ifold], xpd=NA)
            }
        }
        if(is.specified(col.pch.max.oof.rsq) || is.specified(col.pch.cv.rsq)) {
            for(ifold in seq_along(object$cv.list)) {
                oof.rsq <- object$cv.oof.rsq.tab[ifold,]
                scaled.rsq <- scale1(oof.rsq,  ylim[1], ylim[2])
                # show the max oof.rsq for this fold
                nterms <- which.max(oof.rsq)
                points(nterms, scale1(oof.rsq,  ylim[1], ylim[2])[nterms],
                       pch=1, col=col.pch.max.oof.rsq)
                # show the position of the cv.rsq's
                nterms <- length(object$cv.list[[ifold]]$selected.terms)
                points(nterms, scale1(oof.rsq,  ylim[1], ylim[2])[nterms],
                       pch=20, col=col.pch.cv.rsq, cex=.5)
            }
        }
    }
    plot.nused.preds <- function()  # plot nbr of used predictors
    {                               # nothing actually plotted if col.npreds=0
        nused.preds <- get.nused.preds.per.subset(object$dirs, object$prune.terms)
        nused.preds.vec <- scale1(nused.preds, 0, max.npreds)
        if(jitter > 0)  # 2*jitter seems to work better relative to jitter on GRSq
            nused.preds.vec <- jitter(nused.preds.vec, amount=2*jitter)
        else {
            # nudge max value to prevent overplot of maximum RSq(s)
            max <- max(nused.preds.vec)
            nused.preds.vec[nused.preds.vec == max] <- max + max / 150
        }
        lines(nused.preds.vec, type="l", col=col.npreds, lty=lty.npreds)
    }
    plot.vertical.line.at.best.mean.oof.rsq <- function()
    {
        if(!is.specified(col.mean.oof.rsq) || !is.specified(col.oof.vline))
            return()
        v <- which.max(mean.oof.rsq.per.subset)
        # possibly nudge right to prevent overplot of grsq.line
        if(v == length(object$selected.terms))
            v <- v + nterms.on.horiz.axis / 150
        # possibly nudge to prevent overplot of grid
        if(is.specified(col.sel.grid))
            v <- v + nterms.on.horiz.axis / 150
        # lty="23" with palevioletred subjectively looks about same
        # as lty=3 with black used for grsq line
        abline(v=v, col=col.oof.vline, lty="23")
    }
    plot.mean.infold.rsq <- function()
    {
        if(!is.specified(col.mean.infold.rsq))
            return()
        lines(scale1(mean.infold.rsq.per.subset, ylim[1], ylim[2]),
              col=col.mean.infold.rsq, lwd=lwd)
    }
    plot.mean.oof.rsq <- function()
    {
        if(!is.specified(col.mean.oof.rsq))
            return()
        lines(scale1(mean.oof.rsq.per.subset, ylim[1], ylim[2]),
              col=col.mean.oof.rsq, lwd=lwd)
    }
    plot.rsq <- function()
    {
        if(jitter > 0)
            rsq.vec  <- jitter(rsq.vec, amount=jitter)
        lines(scale1(rsq.vec,  ylim[1], ylim[2]), col=col.line, lty=lty.rsq)
    }
    plot.grsq <- function()
    {
        if(jitter > 0)
            grsq.vec <- jitter(grsq.vec, amount=jitter)
        lines(scale1(grsq.vec, ylim[1], ylim[2]), col=col.grsq, lwd=lwd)
    }
    plot.vertical.line.at.best.grsq <- function()
    {
        if(!is.specified(col.vline))
            return()
        v <- length(object$selected.terms)
        # possibly nudge to prevent overplot of grid
        if(is.specified(col.sel.grid))
            v <- v + nterms.on.horiz.axis / 150
        abline(v=v, col=col.vline, lty=lty.vline)
        # possibly plot a colored marker at the top of the above line
        if(is.specified(col.vseg))
            points(x=v, y=1.02, col=col.vseg, pch=6)
    }
    plot.legend <- function()
    {
        # return TRUE if "over" lines obscure "under" lines
        is.obscured <- function(under, over)
        {
            len <- min(length(under), length(over))
            under <- under[1:len]
            over  <- over[1:len]
            i <- under >= ylim[1] & under <= ylim[2]
            i[is.na(under) | is.na(over)] <- FALSE # ignore NAs
            nobscured <- sum(abs(under[i] - over[i]) < (ylim[2] - ylim[1]) / 100)
            nobscured > .8 * sum(i)
        }
        update.legend <- function(text, col=1, lty=1, lwd=1, vert=FALSE)
        {
            if(is.null(legend.text)) { # first time?
                if(text == "")         # spacer between entries?
                    return()           # ignore space when first entry
                legend.text <<- text   # note <<- not <-
                legend.col  <<- col
                legend.lty  <<- lty.as.char(lty)
                legend.lwd  <<- lwd
                legend.vert <<- vert
            } else {
                legend.text <<- c(legend.text, text)
                legend.col  <<- c(legend.col, col)
                legend.lty  <<- c(legend.lty, lty.as.char(lty))
                legend.lwd  <<- c(legend.lwd, lwd)
                legend.vert <<- c(legend.vert, vert)
            }
        }
        #--- plot.legend starts here
        # The is.obscured code assumes that plot order is rsq, mean.oof.rsq, grsq.
        # Obscuring of or by infold.rsq is not yet handled.
        if(!is.null(legend.pos) && all(is.na(legend.pos)))
            return()
        legend.text <- legend.col <- legend.lty <- legend.lwd <- legend.vert <- NULL
        full.model <- if(show.cv.data) " (full model)" else ""
        if(is.specified(col.grsq))
            update.legend(paste0("GRSq", full.model), lwd=lwd)
        if(is.specified(col.vline))
            update.legend("selected model", col.vline, lty.vline, vert=TRUE)
        if(is.specified(col.line)) {
            RSq.string <- if(show.cv.data) "RSq (full model)" else "RSq"
            if(is.specified(col.grsq) && is.obscured(rsq.vec, grsq.vec))
                text <- paste0(RSq.string, full.model, " (obscured)")
            else if(is.specified(col.mean.oof.rsq) &&
                    is.obscured(rsq.vec, mean.oof.rsq.per.subset))
                text <- paste0(RSq.string, " (obscured)")
            else
                text <- RSq.string
            update.legend(text, col.line, lty.rsq)
        }
        added.space <- FALSE
        # We draw the infold legend above the oof legend because the infold
        # curves are usually above the oof curves.
        if(is.specified(col.mean.infold.rsq)) {
            text <- "mean in-fold RSq"
            update.legend("", 0) # dummy entry to leave a vertical space
            added.space <- TRUE
            update.legend(text, col.mean.infold.rsq, lwd=lwd)
        }
        if(is.specified(col.infold.rsq)) {
            if(!added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend("in-fold RSq", col.infold.rsq[1])
        }
        if(is.specified(col.mean.oof.rsq)) {
            if(is.specified(col.grsq) && is.obscured(mean.oof.rsq.per.subset, grsq.vec))
                text <- "mean out-of-fold RSq (obscured)"
            else
                text <- "mean out-of-fold RSq"
            update.legend("", 0) # dummy entry to leave a vertical space
            added.space <- TRUE
            update.legend(text, col.mean.oof.rsq, lwd=lwd)
            if(is.specified(col.oof.vline))
                update.legend("max mean out-of-fold RSq", col.oof.vline, "23", vert=TRUE)
        }
        if(is.specified(col.oof.rsq)) {
            if(!added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend("out-of-fold RSq", col.oof.rsq[1])
        }
        if(is.specified(col.npreds)) {
            if(added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend(paste0("nbr preds", full.model), col.npreds, lty.npreds)
        }
        if(is.null(cex.legend))
            cex.legend <- get.earth.legend.cex(legend.text)
        if(is.null(legend.pos)) { # auto?
            xpos <- if(max.nterms <= 2) "topleft" else "bottomright"
            elegend(x=xpos, bg="white", legend=legend.text, col=legend.col,
                    inset=c(.02, .04), # y offset allows vertical lines to be visible below legend
                    lty=legend.lty, cex=cex.legend, lwd=legend.lwd, xpd=NA, vert=legend.vert)
            if(max.nterms == 1)
                text(.5, .4, "intercept-only model")
        } else { # user specified legend position
            xpos <- legend.pos[1]
            ypos <- if(length(legend.pos) > 1) legend.pos[2] else NULL
            elegend(x=xpos, y=ypos, bg="white", legend=legend.text, col=legend.col,
                    inset=c(.02, .03),
                    lty=legend.lty, cex=cex.legend, lwd=legend.lwd, xpd=NA, vert=legend.vert)
        }
    }
    #--- plot.model.selection starts here ---
    plot.earth.prolog(object, substr(deparse(substitute(x)), 1, 40))
    stop.if.dots.used("plot.model.selection", ...)
    if(is.null(object$prune.terms)) {       # no prune data?
        if(!add)
            plot(c(0,1), col=0, xlab="", ylab="", main="Model Selection")
        legend(x=1, y=1, bty="n", legend=c("No model selection data", "",
               "Run update.earth() to generate", "model selection data"))
        return(NULL)
    }
    if(is.null(main)) {
        if(nresp > 1)
            main <- "Model Selection (all responses)"
        else
            main <- "Model Selection"
    }
    if(jitter > 0.05)
        stop0("\"jitter\" ", jitter , " is too big, try something like jitter=0.01")
    if(!is.specified(lty.grsq))
        col.grsq <- 0
    if(!is.specified(lty.rsq))
        col.line <- 0
    if(!is.specified(lty.npreds))
        col.npreds <- 0
    if(!is.specified(lty.vline))
        col.vline <- 0
    possibly.issue.cv.warning()
    if(is.null(object$cv.oof.rsq.tab)) # if no cv data available, force no display of cv data
        col.mean.oof.rsq <- col.oof.rsq <- col.mean.infold.rsq <- col.infold.rsq <- 0
    show.cv.data <- is.specified(col.mean.oof.rsq)  || is.specified(col.oof.rsq) ||
                    is.specified(col.mean.infold.rsq) || is.specified(col.infold.rsq)
    show.non.cv.data <- is.specified(col.grsq) || is.specified(col.line) || is.specified(col.npreds)
    if(is.null(col.npreds)) # by default, show npreds if not show cv data
        col.npreds <- if(show.cv.data) 0 else 1
    rsq.vec  <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
    grsq.vec <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
    mean.oof.rsq.per.subset <- NULL
    if(is.specified(col.mean.oof.rsq))
        mean.oof.rsq.per.subset <- object$cv.oof.rsq.tab[nrow(object$cv.oof.rsq.tab),]
    if(is.specified(col.mean.infold.rsq))
        mean.infold.rsq.per.subset <- object$cv.infold.rsq.tab[nrow(object$cv.infold.rsq.tab),]
    lwd <- if(show.cv.data) 2 else 1 # want fat non-cv lines if plotting cv data
    nterms.on.horiz.axis <- max.nterms
    if(show.cv.data && !show.non.cv.data)
        nterms.on.horiz.axis <-
            min(nterms.on.horiz.axis, get.max.terms.of.fold.models(object))
    if(!add) {
        if(is.specified(col.npreds)) {
            # make space for right axis
            old.mar <- par("mar")
            if(old.mar[4] < 3.5) {
                mar <- old.mar
                mar[4] <- 3.5
                on.exit(par(mar=old.mar))
                par(mar=mar)
            }
        }
        xlim <- get.model.selection.xlim(object, xlim,
                        mean.oof.rsq.per.subset, col.mean.oof.rsq, col.oof.vline)
        # set up so vertical scale is 0..1, horizontal is 0..nterms.on.horiz.axis
        plot(0:nterms.on.horiz.axis, (0:nterms.on.horiz.axis)/nterms.on.horiz.axis,
            type="n", main=main, xlab="Number of terms", ylab="", yaxt="n", xlim=xlim)
        left.axis()
        if(is.specified(col.npreds))
            right.axis()
        plot.sel.grid()
    }
    # note: if you change the plot order here, modify is.obscured code in plot.legend
    plot.infold.rsqs()
    plot.oof.rsqs()
    plot.nused.preds()
    plot.vertical.line.at.best.grsq()
    plot.vertical.line.at.best.mean.oof.rsq()
    plot.rsq()
    plot.mean.infold.rsq()
    plot.mean.oof.rsq()
    plot.grsq()
    plot.legend()
}
plot.model.selection.wrapper <- function(
    object,

    do.par,
    xlim,
    ylim,
    main,

    col.line            = "lightblue",

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

    col.legend          = NA) # deprecated, use legend.pos=NA for no legend
{
    if(!all(is.na(col.legend))) {
        warning0("col.legend is deprecated.  Please use legend.pos=NA for no legend.")
        legend.pos = NA
    }
    ylim <- get.model.selection.ylim(object, ylim, 1, col.line,
                col.mean.oof.rsq, col.oof.rsq, col.mean.infold.rsq, col.infold.rsq)

    plot.model.selection(
        object,

        do.par,
        xlim,
        ylim,
        main,

        col.line,

        legend.pos,
        cex.legend,

        col.grsq,
        col.infold.rsq,
        col.mean.infold.rsq,
        col.mean.oof.rsq,
        col.npreds,
        col.oof.labs,
        col.oof.rsq,
        col.oof.vline,
        col.pch.cv.rsq,
        col.pch.max.oof.rsq,
        col.sel.grid,
        col.vline,
        col.vseg,
        lty.grsq,
        lty.npreds,
        lty.rsq,
        lty.vline,

        add        = FALSE,
        nresp      = NCOL(object$residuals),
        max.nterms = length(object$rss.per.subset),
        max.npreds = max(1,
                     get.nused.preds.per.subset(object$dirs, object$prune.terms)),
        jitter     = 0)
}
get.model.selection.xlim <- function(object, xlim,
    mean.oof.rsq.per.subset, col.mean.oof.rsq, col.oof.vline)
{
    if(is.null(xlim)) { # not specified by the user?
        # length(object$selected.terms) is nbr selected terms
        # nrow(object$prune.terms) is nk
        # which.max(mean.oof.rsq.per.subset) is max mean out-of-fold RSq
        xmax <- 2 * length(object$selected.terms)
        # if cross-validation vert line is plotted, include that too
        # following "if" matches that in plot.vertical.line.at.best.mean.oof.rsq
        if(!is.null(mean.oof.rsq.per.subset) &&
           is.specified(col.mean.oof.rsq) && is.specified(col.oof.vline)) {
            xmax <- max(xmax, which.max(mean.oof.rsq.per.subset))
        }
        xlim <- c(0, min(xmax + 3, nrow(object$prune.terms)))
    }
    xlim
}
# check ylim specified by user, and convert special values in ylim to actual vals

get.model.selection.ylim <- function(object, ylim, col.grsq, col.line,
                                     col.mean.oof.rsq=0, col.oof.rsq=0,
                                     col.mean.infold.rsq=0, col.infold.rsq=0)
{
    get.fold.min.max <- function()
    {
        min <- Inf
        max <- -Inf
        if(!is.null(object$cv.oof.rsq.tab) &&
                (is.specified(col.mean.oof.rsq) || is.specified(col.oof.rsq))) {
            # will be plotting oof.rsq, so must adjust axis limits for that
            min <- min(object$cv.oof.rsq.tab[,-1], na.rm=TRUE) # -1 to ignore intercept-only model
            max <- max(object$cv.oof.rsq.tab[,-1], na.rm=TRUE)
            # prevent outrageous axis scales caused by wayward cross-validation results
            max <- min(max, 2 * max(rsq))   # 2 is arb
            min <- max(min, -3)             # -3 is arb
        }
        if(!is.null(object$cv.infold.rsq.tab) &&
                (is.specified(col.mean.infold.rsq) || is.specified(col.infold.rsq))) {
            min <- min(min, object$cv.infold.rsq.tab[,-1], na.rm=TRUE)
            max <- max(max, object$cv.infold.rsq.tab[,-1], na.rm=TRUE)
            max <- min(max, 2 * max(rsq))
            min <- max(min, -3)
        }
        list(min=min, max=max)
    }
    #--- get.model.selection.ylim starts here ---
    if(is.null(ylim))
        ylim <- c(-1, -1)
    if(length(ylim) != 2)
        stop0("length(ylim) != 2")
    if(ylim[2] <= ylim[1] && ylim[2] != -1)
        stop0("ylim[2] <= ylim[1]")
    if(ylim[1] < -1 || ylim[1] >  1 || ylim[2] < -1 || ylim[2] >  1)
        stop0(paste0(
              "illegal \"ylim\" c(", ylim[1], ",", ylim[2], ")\n",
              "Legal settings are from -1 to 1, with special values:\n",
              "  ylim[1]=-1 means use min(RSq) or min(GRSq) excluding intercept)\n",
              "  ylim[2]=-1 means use max(RSq) or max(GRSq)\n"))
    if(ylim[1] == -1 || ylim[2] == -1) {
        grsq <- NULL
        if(is.specified(col.grsq))
            grsq <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
        rsq <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
        temp <- get.fold.min.max()
        if(!is.specified(col.line))
            rsq <- NULL
        if(ylim[1] == -1) {
            ylim[1] <- min(grsq[-1], rsq[-1], temp$min, na.rm=TRUE)
            # small model, treat specially so user sees context
            if(length(object$rss.per.subset) <= 3)
                ylim[1] <- min(0, ylim[1])
            ylim[1] <- max(-1, ylim[1]) # clamp minimum ylim at -1
        }
        if(ylim[2] == -1)
            ylim[2] <- max(grsq, rsq, temp$max, na.rm=TRUE)
    }
    # following code gives a decent y axis even with an intercept-only model
    if(abs(ylim[1] - ylim[2]) < 1e-6)
        ylim[2] <- ylim[1] + 1
    ylim
}
get.max.terms.of.fold.models <- function(object)
{
    tab <- object$cv.oof.rsq.tab
    stopifnot(!is.null(tab))
    stopifnot(nrow(tab) > 1)
    max.terms <- 0
    for(i in 1:(nrow(tab)-1)) # -1 to skip last (summary) row
        max.terms <- max(max.terms, sum(!is.na(tab[i,])))
    max.terms
}

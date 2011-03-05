# plotmo.R: plot the model response when varying one or two predictors
#
# Comments containing "TODO" mark known issues.
# Stephen Milborrow Sep 2006 Cape Town
#
# TODO would like to add a data argument, but NA handling seems insoluble
# TODO bug: ycolumn sometimes refers to original y, sometime to expanded y
# TODO are n=* correct in calls to try(eval.parent())?
# TODO weights are ignored
# TODO persp cex handling is not right, because of special handling needed for persp
# TODO with ylim specified persp can plot points outside ylim (xpd=FALSE doesn't help) -- why?
# TODO add code to reduce number of repeated warnings if NAs in x
# TODO allow user to specify x range on a per predictor basis
# TODO the factor code could be simplified maybe
# TODO add factor handling to degree2 plots
# TODO show persp axes even when box=FALSE
# TODO get.plotmo.x should allow allow unnamed cols in x if get from object$x or object$call$x
# TODO would like to make plotmo faster
# TODO would like to add loess of reponse option
# TODO use ndegree1 to limit the number points plotted for "func" argument
#
# TODO following causes error in model.frame: invalid type (list) for variable 'trees[, -3]'
#      a <- earth(trees[,3] ~ as.matrix(trees[,-3])); plotmo(a)
#    work around:
#      a <- earth(trees[,-3], trees[,3]); plotmo(a)
#    or:
#      Girth <- trees[,1]; Height <- trees[,2]; Volume <- trees[,3]
#      a <- earth(Volume ~ Girth + Height); plotmo(a, trace=TRUE)
#
# TODO revisit warning: 'newdata' had 10 rows but variable(s) found have 31 rows
#    e.g. newdata <- as.data.frame(trees[1:10,-3]); predict(a, newdata=newdata)
#
# TODO Fix this i.e. if predictor in munged formula is a numeric then delete it:
# dat <- data.frame(x1=1:20, x2 = sample(1:20), y = rnorm(20)+1:20)
# fm <- lm(y ~ x1+x2+sin(x1*x2/10), dat)
# plotmo(fm, se=2)
# Error in terms.formula(formula, data = data) : invalid model formula in ExtractVars

#------------------------------------------------------------------------------
# Notes for plotmo()
#
# If ylim==NULL we do an initial pass across all plots (both degree1 and degree2)
# without actually plotting to calculate the min and max y limits (actually z limits
# on degree2 plots).  If col.response=TRUE then we also incorporate
# the response min,max.
#
# The name "plotmo" was chosen because it is short, pronounceable as a word,
# yet unlikely to conflict with names in other packages or user code.

plotmo <- function(
    object      = stop("no 'object' arg"),
    degree1     = TRUE,     # index vector specifying main effect plots to include
    degree2     = TRUE,     # index vector specifying interaction plots to include
    ycolumn     = 1,        # which column of resp to use
    type        = "response", # type passed on to predict

    caption = if(do.par) NULL else "",
                            # overall caption
                            #   "string"  string
                            #   ""        no caption
                            #   NULL      generate a caption from object$call

    ylim        = NULL,     # NULL all graphs have same vertical limits
                            # NA   each graph has its own ylim
                            # else specify c(ymin,ymax)

    clip        = TRUE,     # plot only values in range of response of original data

    inverse.func = NULL,    # apply to y before plotting, default NULL uses identity

    col.response = 0,       # color of response values, 0 to not plot response
    pch.response = 1,       # plot character for col.response points
    jitter.response = 0,   # passed on to jitter, NA means no jitter

    trace       = FALSE,    # trace operation

    grid.func   = median,   # func applied to x columns to calc plot grid

    grid.levels = NULL,     # list specifying which factor levels to use plot grid, default is first

                            # Following arguments are for degree1 plots
    ndegree1    = 500,      # max nnbr of points in degree1 plot (ndegree=-1 for nrow(x))
    lty.degree1 = 1,        # line type for degree1 plots
    col.degree1 = 1,        # colour for degree1 plots

    se          = 0,        # plot std err lines at se*stderr, use 0 for none
    col.shade   = "lightblue", # colour for se shading
    col.se      = 0,        # colour for se lines, 0 means no lines just shading
    lty.se      = 2,        # line type for se lines, 0 means no lines just shading

    func        = NULL,     # superimpose func(x) if func is not null
    col.func    = "pink",   # color of func values in degree1 plots
    pch.func    = 20,       # plot character for col.func points

    nrug        = 0,        # number of points in degree1 rug (nrug=-1 for all)

                            # Following arguments are for degree2 plots
    type2       = "persp",  # degree2 plot type: for options see draw.plot2()
    ngrid       = 20,       # grid side length for degree2 plots
    col.persp   = "lightblue",  # color of persp surface
    col.image   = grey(0:9/10), # colors for image() plot. 1 unused, is for clipped vals

                            # par() settings, listed so caller can override our use
    do.par      = TRUE,     # call par() global settings as appropriate
    main        = NULL,     # title of each plot, default is predictor names
    theta       = NA,       # rotation arg for persp(), NA means let this func choose
    phi         = 30,       # phi arg for persp()
    shade       = 0.5,      # shade arg for persp()
    ticktype    = "simple", # ticks on persp() plot, one of: simple detailed
    xlab        = "",       # add x axis labels, use NULL for auto, only applies to degree1
    ylab        = "",       # vertical axis label, default none gives more plottable area
    cex         = NULL,
    cex.lab     = 1,        # used as cex.lab and cex.axis
    ...)                    # extra args passed to plotting funcs
{
    get.ylims <- function() # check ylim arg and calculate min,max limits of y axis
    {
        if(!(is.null(ylim) || is.na(ylim[1]) || length(ylim) == 2))
            stop1("'ylim' must be one of:\n",
                "  NULL         all graphs have same vertical axes\n",
                "  NA           each graph has its own vertical axis\n",
                "  c(ymin,ymax) specify y axis min and max")
        if(length(ylim) == 2 && ylim[2] <= ylim[1])
            stop1("ylim[2] ", ylim[2], " is not greater than ylim[1] ", ylim[1])
        ylims <- c(NA, NA)
        if(is.null(ylim)) { # user wants same vertical ylims for all graphs?
            if(nsingles)    # if so, get ylims=c(miny, maxy) by calling with draw.plot=FALSE
                ylims <- plot1(
                        object, degree1, ylim, ycolumn, type, clip, col.response,
                        pch.response, jitter.response, inverse.func, grid.func, grid.levels,
                        ndegree1, lty.degree1, col.degree1, se, lty.se,
                        col.se, col.shade, func, col.func, pch.func, nrug,
                        draw.plot=FALSE, x, y, Singles, ylims, func.arg, pred.names,
                        ntrace, inverse.func.arg, clip.limits, nfigs,
                        main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
            if(npairs)
                ylims <- plot2(object, degree2, ylim, ycolumn, type, clip,
                            col.response, pch.response, jitter.response, inverse.func,
                            grid.func, grid.levels, type2, ngrid, col.persp, col.image,
                            draw.plot=FALSE, x, y, Pairs, ylims, func.arg, pred.names,
                            ntrace, inverse.func.arg, clip.limits, nfigs, nsingles, npairs,
                            do.par, main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
            if(!is.na.or.zero(col.response))
                ylims <- range(ylims, y, finite=TRUE)
            if(trace > 0)
                cat("\ninitialized ylims", ylims, "\n")
            if(any(!is.finite(ylims)))
                stop1("ylim argument is NULL but cannot generate ",
                    "ylim internally from predicted values")
            if(ntrace > 1)
                ntrace <- ntrace - 1 # prevent repeated trace msgs in plot1 and plot2
        }
        ylims   # calculated min,max vertical axis ylims
    }
    # return the response.name (for prepending to caption), if appropriate
    get.caption.prefix <- function(y, ycolumn, caption)
    {
        if(!is.null(caption))
            return(NULL)             # don't modify caption explictly set by user
        colnames <- colnames(y)
        if(!is.null(colnames) && !is.null(colnames[1]) &&
           !is.na(colnames[1]) && colnames[1] != "")
            colnames[1]
        # TODO following is not quite right: should test against nbr of cols
        # from predict so response 1 is labelled uniformly
        else if(ycolumn > 1)
            paste("Response", ycolumn)
        else
            NULL
    }
    # get object's glm family if there is one, returns a string or NULL
    get.glm.family <- function(object, ycolumn)
    {
        family <- NULL
        if(!is.null(object$glm.list[[ycolumn]])) # object class is "earth"
            family <- object$glm.list[[ycolumn]]$family$family
        else if(!is.null(object$family))         # object class is "glm" and similar
            family <- object$family$family
        family
    }
    get.clip.limits <- function(object, ycolumn, y) # return a vector with two elements
    {
        # Get glm family if there is one.  Needed because if there is a glm family
        # we can't use y to establish clip.limits (because not in "response" scale)
        # TODO this only works for the common glm cases?
        #      so not e.g. for fam="pois", link="log"
        family <- get.glm.family(object, ycolumn)
        if(!is.null(family) &&
                pmatch(family, c("binomial","quasibinomial"), nomatch=0))
            clip.limits <- c(0,1) # binomial or quasibinomial model
        else
            clip.limits <- range(y, finite=TRUE)
        if(trace > 0)
            cat("clip.limits", clip.limits, "\n")
        clip.limits
    }
    # plotmo starts here
    plotmo.prolog(object, deparse(substitute(object)))
    ntrace <- if(trace > 0) 2 else 0
    func.arg <- deparse(substitute(func))
    inverse.func.arg <- deparse(substitute(inverse.func))

    # check se arguments
    if(is.logical(se) && se) { # allow user to use se=TRUE for compat with termplot
        warning1("plotmo: converted se=TRUE to se=2")
        se <- 2
    }
    if(!is.numeric(se) || se < 0 || se > 9)
        stop1("'se' ", se, " is out of range, range is 0...9") # 9 is arbitrary
    if(!missing(lty.se) && lty.se != 0 && col.se[1] == 0)
        col.se <- 1  # needed if user sets just lty.se but doesn't set col.se
    if(se && (col.se == 0 || lty.se == 0) && col.shade == 0)
        warning1(
          "plotmo: 'se' ignored because (col.se == 0 || lty.se == 0) && col.shade == 0)")
    if(se == 0 && (!missing(col.se) || !missing(lty.se) || !missing(col.shade)))
        warning1("plotmo: se color and linetype arguments ignored because se=0")

    if(!is.na(pmatch(type, "terms")))
        stop1("type=\"terms\" is not allowed by plotmo")

    x <- get.plotmo.x(object, trace)
    warn.if.not.all.finite(x, "'x' after calling get.plotmo.x()")
    if(trace > 0) {
        cat("x[", nrow(x), ",", ncol(x), "]\n", sep="")
        print(head(x, 3))
        cat("is.factor", sapply(x, is.factor), "\n")
    }
    pred.names <- colnames(x)
    y <- get.plotmo.y(object, ycolumn, nrow(x), trace)
    caption.prefix <- get.caption.prefix(y, ycolumn, caption)
    y <- apply.inverse.func(y, ycolumn, trace, inverse.func, inverse.func.arg)
    if(clip)
        clip.limits <- get.clip.limits(object, ycolumn, y)
    if(ndegree1 == -1)
        ndegree1 <- nrow(x)
    else if(ndegree1 < 1 || ndegree1 > 1e5)      # 1e5 is arbitrary
        stop1("illegal ndegree1 ", ndegree1, ", allowed range is 1 to 1e5 or -1")

    # Singles is a vector of indices of predictors for degree1 plots
    Singles <- get.singles(object, x, degree1, pred.names, trace)
    nsingles <- length(Singles)
    if(trace > 0)
        if(length(Singles) > 0)
            cat("singles:", paste(Singles, pred.names[Singles], collapse=", "), "\n")
        else
            cat("no singles\n")

    # Each row of Pairs is the indices of two predictors for a degree2 plot
    Pairs <- get.pairs(object, x, degree2, pred.names, trace)
    npairs <- nrow(Pairs)
    if(trace > 0)
        if(nrow(Pairs) > 0) {
            cat("pairs:\n")
            print(matrix(paste(Pairs, pred.names[Pairs]), ncol=2))
        } else
            cat("no pairs\n")

    nfigs <- nsingles + npairs
    if(nfigs == 0) {
        warning1("plotmo: nothing to plot")
        return(invisible())
    }
    ylims <- get.ylims()    # check ylim arg and calculate min,max limits of y axis
    must.trim.caption <- is.null(caption)
    caption <- get.caption.from.call(caption, object)
    if(!is.null(caption.prefix))
        caption <- paste(caption.prefix, ": ", caption, sep="")
    if(do.par) {
        old.par <- par(no.readonly=TRUE)
        on.exit(par(old.par))
        nrows <- ceiling(sqrt(nfigs))
        par(mfrow=c(nrows, nrows))
        par(cex = 0.7)
        par(mar = c(4, 3, 1.7, 0.5))    # small margins and text to pack figs in
        par(mgp = c(1.6, 0.6, 0))       # flatten axis elements
        if(nrows >= 4) {
            par(mar = c(2, 2, 1.7, 0.5))
            par(cex = 0.6)
        }
        if(is.null(ylab) || nchar(ylab) > 0)
            make.space.for.left.axis()
        if(is.null(xlab) || nchar(xlab) > 0)
            make.space.for.bottom.axis()
        make.space.for.caption(caption)
    }
    if(ntrace > 1)
        ntrace <- ntrace - 1    # prevent repeated trace msgs in plot1 and plot2
    if(nsingles)
        plot1(object, degree1, ylim, ycolumn, type, clip, col.response,
            pch.response, jitter.response, inverse.func, grid.func, grid.levels,
            ndegree1, lty.degree1, col.degree1, se, lty.se,
            col.se, col.shade, func, col.func, pch.func, nrug,
            draw.plot=TRUE, x, y, Singles, ylims, func.arg, pred.names,
            ntrace, inverse.func.arg, clip.limits, nfigs,
            main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
    if(npairs) {
        if(do.par) {                                 # use degree1 settings, but tweaked
                if(pmatch(type2, "persp", 0) != 1) { # contour or image plot?
                    make.space.for.bottom.axis()
                    make.space.for.left.axis()
                } else if(pmatch(ticktype, "simple", 0) == 1) {
                    par(mar = c(.5, 0.3, 1.7, 0.1))
                    par(mgp = c(2, 0.6, 0))
                } else {                            # ticktype="detailed"
                    par(mar = c(1, 0.3, 1.7, 0.1))
                    par(mgp = c(2, 0.6, 0))         # TODO this doesn't work
                }
            make.space.for.caption(caption)
        }
        plot2(object, degree2, ylim, ycolumn, type, clip,
            col.response, pch.response, jitter.response, inverse.func,
            grid.func, grid.levels, type2, ngrid, col.persp, col.image,
            draw.plot=TRUE, x, y, Pairs, ylims, func.arg, pred.names,
            ntrace, inverse.func.arg, clip.limits, nfigs, nsingles, npairs,
            do.par, main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
    }
    show.caption(caption, trim=must.trim.caption)
    invisible()
}

# Plot degree one plots i.e. main effects
# This is an auxilary function for plotmo so could be local to plotmo,
# and we could avoid a long argument list.  However, it seems a bit
# clearer to make it a global level function.

plot1 <- function(
    # copy of args from plotmo, some have been tweaked slightly
    object, degree1, ylim, ycolumn, type, clip,
    col.response, pch.response, jitter.response,
    inverse.func, grid.func, grid.levels,
    ndegree1, lty.degree1, col.degree1, se, lty.se, col.se, col.shade,
    func, col.func, pch.func, nrug,

    # args generated in plotmo, draw.plot=FALSE means get ylims but don't actually plot
    draw.plot, x, y, Singles, ylims, func.arg, pred.names,
    ntrace, inverse.func.arg, clip.limits, nfigs,

    # copy of par args from plotmo
    main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
{
    draw.func <- function() # draw the func argument, if specified
    {
        if(exists.and.not.null(func.arg, "function", "func")) {
            if(ntrace) {
                cat("\nApplying 'func' arg to\n")
                print(head(xwork, 3))
            }
            y.func <- func(xwork)
            y.func <- check.and.print.y(y.func, paste("func=", func.arg, sep=""),
                                        1, nrow(xwork), ntrace)  # ycolumn always 1
            points(xwork[,ipred], y.func, col=col.func, pch=pch.func)
        }
        NULL
    }
    draw.plot1 <- function()
    {
        draw.factor.se <- function()   # draw std err bands and lines for a factor predictor
        {
            y.se.lower1 <- split(y.se.lower, xwork[,ipred])
            y.se.upper1 <- split(y.se.upper, xwork[,ipred])
            for(ilev in seq_along(levels(xwork[,ipred]))) {
                min <- min(y.se.lower1[[ilev]])
                max <- max(y.se.upper1[[ilev]])
                if(!is.na.or.zero(col.shade))
                    polygon(c(ilev - .4, ilev - .4, ilev + .4, ilev + .4),
                        c(min, max, max, min), col=col.shade, lty=0, border=NA)
                if(lty.se != 0 && !is.na.or.zero(col.se)) {
                    segments(ilev -.4, min, ilev + .4, min, lty=lty.se, col=col.se)
                    segments(ilev -.4, max, ilev + .4, max, lty=lty.se, col=col.se)
                }
            }
            NULL
        }
        draw.numeric.se <- function() # draw std err bars and lines for a numeric predictor
        {
            if(!is.na.or.zero(col.shade))
                polygon(c(xwork[,ipred], rev(xwork[,ipred])),
                        c(y.se.lower, rev(y.se.upper)),
                        col=col.shade, lty=0, border=NA)
            if(lty.se != 0 && !is.na.or.zero(col.se)) {
                lines(xwork[,ipred], y.se.lower, lty=lty.se, col=col.se)
                lines(xwork[,ipred], y.se.upper, lty=lty.se, col=col.se)
            }
            NULL
        }
        get.main1 <- function(main, isingle, nfigs, degree1, pred.names, ipred)
        {
            main <- ""
            # show degree1 plot numbers in headers if plotting all predictors
            if(nfigs > 1 && !is.specified(degree1))
                main <- paste(isingle, " ", sep="")
            paste(main, pred.names[ipred])
        }
        # draw.plot1 starts here
        if(is.null(ylim.org))           # same ylim for each graph?
            ylim <- ylims
        else if(is.na(ylim.org[1])) {   # each graph has its own ylim?
            ylim <- range1(y.predict, finite=TRUE)
            if(!is.null(y.se.lower))
                ylim <- range1(y.predict, y.se.lower, y.se.upper, finite=TRUE)
            if(any(!is.finite(ylim)))
                stop1("ylim argument to plotmo is NA but cannot generate ",
                    "ylim internally from predicted values (predictor \"",
                    pred.names[ipred], "\")")
        }
        if(nrug && (is.null(ylim.org) || is.na(ylim.org[1])))
            ylim[1] <- ylim[1] - .05 * (ylim[2] - ylim[1])
        main <- get.main(main, isingle, get.main1,
                         nfigs, degree1, pred.names, ipred)
        if(is.null(xlab))
            xlab <- pred.names[ipred]
        levnames <- levels(xwork[,ipred])
        plot.levnames <- is.fac && length(levnames) <= 20
        plot(xwork[,ipred], y.predict, type="n",    # calls boxplot if is.fac
                main=main, xlab=xlab, ylab=ylab, ylim=ylim,
                xaxt=if(plot.levnames) "n" else "s",
                cex.lab=cex.lab, cex.axis=cex.lab)
        if(is.fac) {
            if(!is.null(y.se.lower))
                draw.factor.se()
            if(!is.na.or.zero(col.response))
                points(jitter(as.numeric(x[,ipred]), factor=.5), y,
                       pch=pch.response, col=col.response, cex=cex)
            plot(xwork[,ipred], y.predict, add=TRUE,
                 xaxt=if(plot.levnames) "n" else "s",
                 cex.lab=cex.lab, cex.axis=cex.lab)
            if(plot.levnames) # plot level names vertically along the x axis
                mtext(abbreviate(levnames, minlength=6, strict=TRUE),
                      side=1, at=1:length(levnames),
                      cex=par("cex") * cex.lab, las=2, line=0.5)
        } else {    # not is.fac
            if(!is.null(y.se.lower))
                draw.numeric.se()
            draw.func()
            if(!is.na.or.zero(col.response))
                points(jitter(x[,ipred], factor=jitter.response),
                       jitter(y,         factor=jitter.response),
                       pch=pch.response, col=col.response, cex=cex)
            lines(x=xwork[,ipred], y=y.predict, lty=lty.degree1, col=col.degree1)
        }
        if(nrug)
            rug(jitter(as.numeric(x[irug, ipred])))
        NULL
    }
    # plot1 starts here
    if(ntrace)
        cat("\n--plot1(draw.plot=", draw.plot, ") ntrace ", ntrace, "\n", sep="")
    check.grid.levels(grid.levels, pred.names)
    ylim.org <- ylim
    xgrid <- data.frame(matrix(0, ndegree1, ncol(x), byrow=TRUE))
    for(ipred in 1:ncol(x))
        if(is.factor(x[,ipred])) {
            lev <- get.level(x[,ipred], pred.names[ipred], grid.levels)
            xgrid[,ipred] <- lev
            if(draw.plot && !is.null(grid.levels)) {
                cat("Note: setting", pred.names[ipred], "in the degree1 grid to ");
                print(lev, max.levels=0)
            }
        } else
            xgrid[,ipred] <- grid.func(x[,ipred], na.rm=TRUE)
    warn.if.not.all.finite(xgrid, "'xgrid' for degree1 plots")
    stopif(is.null(pred.names))
    colnames(xgrid) <- pred.names
    if(!draw.plot && ntrace >= 0) {
        # show the grid values, must do some finangling for a nice display
        cat("\n")
        row <- xgrid[1,,drop=FALSE]
        names(row) <- c(paste("grid:   ", names(row)[1]), names(row)[-1])
        rownames(row) <- ""
        print(row)
        cat("\n")
    }
    if(draw.plot && nrug) {
        if(nrug < 0)
            nrug <- nrow(x)
        if(nrug > nrow(x))
            nrug <- nrow(x)
        irug <- as.integer(seq(from=1, to=nrow(x), length.out=nrug))
    }
    for(isingle in seq_along(Singles)) {
        ipred <- Singles[isingle] # ipred is the predictor index i.e. col in model mat
        # following happens with lm if you do e.g. ozone1$doy <- NULL after using ozone1
        # TODO I am not sure if this is enough to always catch such errors
        if(ipred > NCOL(x))
            stop1("bad index (missing column in x?)")
        xwork <- xgrid
        is.fac <- is.factor(x[,ipred])
        if(is.fac) {
            # TODO this is inefficient because
            # (i)  it calls sort and unique
            # (ii) there are lots of repeated rows in xwork
            xwork[,ipred] <- rep(sort(unique(x[,ipred])), length=ndegree1)
        } else {
            xrange <- range1(x[,ipred], finite=TRUE)
            xwork[,ipred] <- seq(from=xrange[1], to=xrange[2], length=ndegree1)
        }
        y.predict <- plotmo.predict(object, xwork, type=type, se.fit=FALSE,
                                    pred.names, ipred, 0, ntrace>1)
        y.predict <- check.and.print.y(y.predict, "plotmo.predict",
                                       ycolumn, nrow(xwork), ntrace>1)
        y.se.lower <- NULL
        y.se.upper <- NULL
        if(!is.na.or.zero(se)) {
            if(ntrace > 1)
                cat("begin se handling, ")
            rval <- plotmo.predict(object, xwork, type=type, se.fit=TRUE,
                                   pred.names, ipred, 0, ntrace>1)
            if(typeof(rval) == "list" && !is.null(rval$se.fit)) {
                rval$se.fit <- check.and.print.y(rval$se.fit, "predict with se=TRUE",
                                                 ycolumn, nrow(xwork), ntrace>1)
                y.se.lower <- y.predict - se * rval$se.fit
                y.se.lower <- apply.inverse.func(y.se.lower, ycolumn, ntrace>1,
                                                 inverse.func, inverse.func.arg)
                y.se.upper <- y.predict + se * rval$se.fit
                y.se.upper <- apply.inverse.func(y.se.upper, ycolumn, ntrace>1,
                                                 inverse.func, inverse.func.arg)
            } else if(ntrace > 1)
                cat("no standard errs because is.null(rval$se.fit)\n")
            if(ntrace > 1)
                cat("end se handling\n")
        }
        y.predict <- apply.inverse.func(y.predict, ycolumn,
                                        ntrace>1, inverse.func, inverse.func.arg)
        if(clip) {
            y.lt <- y.predict < clip.limits[1]
            y.gt <- y.predict > clip.limits[2]
            if(all(y.lt) || all(y.gt)) # TODO revisit for better error handling
                warning1("plotmo: Can't clip y to ",
                    clip.limits[1], " to ", clip.limits[2])
            else
                y.predict[y.lt | y.gt] <- NA
        }
        ylims <- range1(ylims, y.predict, y.se.lower, y.se.upper, finite=TRUE)
        if(draw.plot)
            draw.plot1()
    }
    ylims
}

# Plot degree two plots i.e. first order interactions
# Auxilary function for plotmo

plot2 <- function(
    # copy of args from plotmo, some have been tweaked slightly
    object, degree2, ylim, ycolumn, type, clip,
    col.response, pch.response, jitter.response,
    inverse.func, grid.func, grid.levels, type2, ngrid, col.persp, col.image,

    # args generated in plotmo, draw.plot=FALSE means get ylims but don't actually plot
    draw.plot, x, y, Pairs, ylims, func.arg, pred.names,
    ntrace, inverse.func.arg, clip.limits, nfigs, nsingles, npairs,

    # copy of par args from plotmo
    do.par, main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
{
    draw.plot2 <- function(type2 = c("persp", "contour", "image"))
    {
        plot.persp <- function()
        {
            get.theta <- function()
            {
                get.diag.val <- function(diag1, diag2) # return first non NA along diag
                {
                    vals <- y.predict[diag1, diag2]
                    (vals[!is.na(vals)])[1] # return first non NA in vals
                }
                theta1 <- theta
                if(is.na(theta)) {      # no user specified theta?
                    # rotate graph so highest point is farthest
                    theta1 <- -35
                    imax <- which.max(c(
                            get.diag.val(ngrid:1, ngrid:1),
                            get.diag.val(1:ngrid, ngrid:1),
                            get.diag.val(1:ngrid, 1:ngrid),
                            get.diag.val(ngrid:1, 1:ngrid)))
                    if(length(imax))   # length>0 unless entire diag is NA
                        theta1 <- theta1 + switch(imax, 0, 90, 180, 270)
                }
                theta1                  # theta arg for persp()
            }
            # plot.persp starts here
            theta1 <- get.theta()
            if(is.null(cex))            # do this only if user didn't specify a cex
                cex1 <- par("cex")      # persp needs an explicit cex arg
            else
                cex1 <- cex
            if(ntrace > 1)
                cat(pred.names[i1], ":", pred.names[i2], " theta ", theta1,
                        " ylim ", ylim, " cex ", cex1, " phi ", phi, "\n", sep="")
            persp(x1, x2, y.predict, main=main, theta=theta1, phi=phi, col=col.persp,
                ticktype=ticktype, shade=shade, zlim=ylim, cex=cex1,
                zlab=ylab, xlab=pred.names[i1], ylab=pred.names[i2],
                cex.lab=1.1 * cex.lab, cex.axis=cex.lab)
            NULL
        }
        plot.response.sites <- function()
        {
            if(!is.na.or.zero(col.response))
                points(jitter(x[,i1], factor=jitter.response),
                       jitter(x[,i2], factor=jitter.response),
                       pch=pch.response, col=col.response, cex=cex)
            NULL
        }
        plot.contour <- function()
        {
            contour(x1, x2, y.predict, main=main,
                xlab=pred.names[i1], ylab=pred.names[i2],
                cex.lab=1.1 * cex.lab, cex.axis=cex.lab,
                labcex=.8 * cex.lab * par("cex"), ...)
            plot.response.sites()
        }
        plot.image <- function()
        {
            image(x1, x2, y.predict, main=main, col=col.image,
                xlab=pred.names[i1], ylab=pred.names[i2],
                cex.lab=1.1 * cex.lab, cex.axis=cex.lab, ...)
            plot.response.sites()
        }
        get.main2 <- function(main, imain, nfigs, degree2, ipair, pred.names, i1, i2)
        {
            main <- ""
            # show degree2 plot numbers in headers if plotting all predictors
            if(nfigs > 1 && !is.specified(degree2))
                main <- paste(ipair, " ", sep="")
            paste(main, pred.names[i1], ":", pred.names[i2], sep="")
        }
        # draw.plot2 starts here
        main <- get.main(main, nsingles+ipair, get.main2,
                         nfigs, degree2, ipair, pred.names, i1, i2)
        if(is.null(ylim))           # same ylim for each graph?
            ylim <- ylims
        else if(is.na(ylim[1])) {   # each graph has its own ylim?
            ylim <- range1(y.predict, finite=TRUE)
            if(any(!is.finite(ylim)))
                stop1("ylim argument to plotmo is NA but cannot generate ",
                    "ylim internally from predicted values (predictors \"",
                    pred.names[i1], "\" \"", pred.names[i2], "\")")
        }
        switch(match.arg1(type2),
            plot.persp(),
            plot.contour(),
            plot.image())
        NULL
    }
    # plot2 starts here
    if(ntrace)
        cat("\n--plot2(draw.plot=", draw.plot, ")\n", sep="")
    check.grid.levels(grid.levels, pred.names)
    ngrid <- min(ngrid, nrow(x), na.rm=TRUE)
    # each row of xgrid is identical, each row is the col
    # median for that row of x (or selected level for factors)
    grid.func1 <- function(x)
        if(is.factor(x)) get.level(x, NULL, NULL) else grid.func(x, na.rm=TRUE)
    xgrid.row <- lapply(x, grid.func1)  # one row of xgrid
    for(ipred in 1:ncol(x))             # fix factor levels TODO combine with above?
        if(is.factor(x[,ipred])) {
            lev <- get.level(x[,ipred], pred.names[ipred], grid.levels)
            xgrid.row[[ipred]] <- lev
            if(draw.plot && !is.null(grid.levels)) {
                cat("Note: setting", pred.names[ipred], "in the degree2 grid to ");
                print(lev, max.levels=0)
            }
        }
    xgrid <- data.frame(xgrid.row)
    xgrid[1:ngrid^2, ] <- xgrid.row
    warn.if.not.all.finite(xgrid, "'xgrid' for degree2 plots")
    colnames(xgrid) <- pred.names
    xranges <- sapply(x, range1)
    stopifnot(npairs > 0)
    ignored.factors.msg <- NULL # reminder message to user
    for(ipair in 1:npairs) {
        i1 <- Pairs[ipair,1]    # index of first predictor
        i2 <- Pairs[ipair,2]    # index of second predictor
        #TODO for now, factors can't be plotted
        if(is.factor(xgrid[,i1]) || is.factor(xgrid[,i2])) {
            if(is.factor(xgrid[,i1]))
                ignored.factors.msg <- c(ignored.factors.msg, pred.names[i1])
            if(is.factor(xgrid[,i2]))
                ignored.factors.msg <- c(ignored.factors.msg, pred.names[i2])
            next
        }
        x1 <- seq(from=xranges[1,i1], to=xranges[2,i1], length=ngrid)
        x2 <- seq(from=xranges[1,i2], to=xranges[2,i2], length=ngrid)
        # xwork is a grid of x vals, all x vals are medians (or selected
        # factor level) except at i1 and i2
        xwork <- xgrid
        xwork[, i1] <- rep(x1, ngrid)
        xwork[, i2] <- rep(x2, rep(ngrid, ngrid))
        y.predict <- plotmo.predict(object, xwork, type=type, se.fit=FALSE,
                                    pred.names, i1, i2, ntrace>1)
        y.predict <- apply.inverse.func(y.predict, ycolumn,
                                        ntrace>1, inverse.func, inverse.func.arg)
        y.predict <- check.and.print.y(y.predict, "inverse.func", ycolumn, nrow(xgrid), ntrace>1)
        y.predict <- matrix(y.predict, ncol=ngrid, nrow=ngrid)
        if(clip)
            y.predict[y.predict < clip.limits[1] | y.predict > clip.limits[2]] <- NA
        ylims <- range1(ylims, y.predict, finite=TRUE)
        if(draw.plot)
            draw.plot2(type2)
    }
    # test on draw.plot below prevents same mesage being printed twice
    if(draw.plot && !is.null(ignored.factors.msg))
        cat("Note: plotmo cannot use factors as axes for degree2 plots,\n",
            "      so skipping degree2 plots for: ",
            paste.with.comma(unique(ignored.factors.msg)),
            "\n", sep="")
    ylims
}

get.main <- function(main, imain, main.func, ...)
{
    if(is.null(main))
        main <- main.func(main, imain, ...)
    else if(length(main) > 1) { # user supplied a vector main?
        if(imain > length(main)) {
            if(imain == length(main)+1) # issue warning only once
                warning1("not enough elements in \"main\" ",
                         "(there are more plots than strings in \"main\")")
            main <- main.func(main, imain, ...) # revert to defaults
        } else
            main <- main[imain]
    }
    main
}

range1 <- function(x, ...)
{
    if(is.factor(x))
        c(1, nlevels(x))
    else
        range(x, ...)
}

# sanity check of grid.levels arg
# actual levels will be checked in get.level

check.grid.levels <- function(grid.levels, pred.names)
{
    if(!is.null(grid.levels)) {
        if(!is.list(grid.levels))
            stop1("grid.levels must be a list. ",
                 "Example: grid.levels=list(sex=\"male\")")
        names. <- names(unlist(grid.levels)) # get list element names
        for(name in names.)
            if(length(grep(paste("^", name, "$", sep=""), pred.names)) == 0)
                stop1("illegal predictor name \"", name, "\" in grid.levels")
    }
}

# This returns a factor level by looking for the level associated
# with pred.name in grid.levels.
# Here, grid.levels is a list e.g. list(pclass="2nd", sex="female")
#
# If it finds the pred.name, it returns the specified level in x
# else it returns the first level in x

get.level <- function(x, pred.name, grid.levels)
{
    ilev <- 1                                       # by default use the first level
    lev.names <- levels(x)
    if(!is.null(pred.name) && !is.null(grid.levels)) {
        lev.name <- grid.levels[[pred.name]]
        if(!is.null(lev.name)) {                   # lev.name is in grid.levels?
            ilev <- lev.names == lev.name
            ilev <- which(ilev)
            if(length(ilev) != 1)
                stop1("illegal level \"", lev.name,
                   "\" specified for \"", pred.name,
                   "\" in grid.levels (allowed levels are ",
                   paste.quoted.names(lev.names), ")")
        }
    }
    factor(lev.names, levels=lev.names)[ilev]
}

# check that y is good
# also, if y has multiple columns this returns just the column specified by ycolumn

check.and.print.y <- function(y, msg, ycolumn, expected.len, trace, subset.=NULL)
{
    if(is.null(y))
        stop1(msg,  " returned NULL")
    if(length(y) == 0)
        stop1(msg, " returned a zero length value (length(y) == 0)")
    check.index.vec("ycolumn", ycolumn, y, check.empty = TRUE,
                    use.as.col.index=TRUE, allow.negative.indices=FALSE)
    if(NCOL(y) > 1)
        y <- y[, ycolumn, drop=FALSE]
    # convert dataframe to matrix, needed for test.earth.glm.R test a21
    # TODO revisit
    if(is.data.frame(y))
        y <- as.matrix(y)
    if(NCOL(y) > 1)
        stop1("'ycolumn' specifies more than one column")
    if(NROW(y) > 1 && NCOL(y) > 1) {
        if(trace > 0) {
            cat(msg, ":\n", sep="")
            print(head(y, 3))
        }
        warning1(msg, " returned an unexpected response, it has dimensions ",
            NROW(y), " x ", NCOL(y), " with colnames ",
            if(NCOL(y) > 1) paste.quoted.names(colnames(y)) else "")
    }
    if(any(!is.double(y))) # convert logical or whatever to double
        y <- as.vector(y, mode="numeric")
    if(NROW(y) == 1 && NCOL(y) > 1)
        y <- as.vector(y[,1])
    warn.if.not.all.finite(y, "y")
    if(trace > 0) {
        cat(msg, " returned length ", length(y), sep="")
        if(!is.null(subset.))
            cat(" (before taking subset)")
        try(cat(" min", min(y), "max", max(y)), silent=TRUE)
        cat(" values ")
        for(i in 1:min(10, length(y)))
            cat(y[i], "")
        cat("...\n")
    }
    if(is.null(subset.) && length(y) != expected.len)
        warning1(msg, " returned a response of the wrong length ", length(y),
                 ", expected ", expected.len)
    y
}

apply.inverse.func <- function(y, ycolumn, trace, inverse.func, inverse.func.arg)
{
    if(exists.and.not.null(inverse.func.arg, "function", "inverse.func")) {
        y <- inverse.func(y)
        y <- check.and.print.y(y, paste("inverse.func=", inverse.func.arg, sep=""),
                               ycolumn, length(y), trace)
    }
    y
}

# TRUE if a plot was selected by the user (excluding the default setting)

is.specified <- function(degree)
{
    !is.logical(degree) || (length(degree) > 1)
}

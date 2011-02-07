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
                        pch.response, inverse.func, grid.func, grid.levels,
                        ndegree1, lty.degree1, col.degree1, se, lty.se,
                        col.se, col.shade, func, col.func, pch.func, nrug,
                        draw.plot=FALSE, x, y, Singles, ylims, func.arg, pred.names,
                        ntrace, inverse.func.arg, clip.limits, nfigs,
                        main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
            if(npairs)
                ylims <- plot2(object, degree2, ylim, ycolumn, type, clip,
                            col.response, pch.response, inverse.func,
                            grid.func, grid.levels, type2, ngrid, col.persp, col.image,
                            draw.plot=FALSE, x, y, Pairs, ylims, func.arg, pred.names,
                            ntrace, inverse.func.arg, clip.limits, nfigs, nsingles, npairs,
                            do.par, main, theta, phi, shade, ticktype, xlab, ylab, cex, cex.lab, ...)
            if(col.response != 0)
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
    if(!missing(lty.se) && lty.se != 0 && col.se == 0)
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
            pch.response, inverse.func, grid.func, grid.levels,
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
            col.response, pch.response, inverse.func,
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
    object, degree1, ylim, ycolumn, type, clip, col.response, pch.response,
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
                if(col.shade != 0)
                    polygon(c(ilev - .4, ilev - .4, ilev + .4, ilev + .4),
                        c(min, max, max, min), col=col.shade, lty=0, border=NA)
                if(lty.se != 0 && col.se != 0) {
                    segments(ilev -.4, min, ilev + .4, min, lty=lty.se, col=col.se)
                    segments(ilev -.4, max, ilev + .4, max, lty=lty.se, col=col.se)
                }
            }
            NULL
        }
        draw.numeric.se <- function() # draw std err bars and lines for a numeric predictor
        {
            if(col.shade != 0)
                polygon(c(xwork[,ipred], rev(xwork[,ipred])),
                        c(y.se.lower, rev(y.se.upper)),
                        col=col.shade, lty=0, border=NA)
            if(lty.se != 0 && col.se != 0) {
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
            if(col.response != 0)
                points(jitter(as.numeric(x[,ipred]), .5),
                    y, pch=pch.response, col=col.response, cex=cex)
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
            if(col.response != 0)
                points(x[,ipred], y, pch=pch.response, col=col.response, cex=cex)
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
        if(se != 0) {
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
    object, degree2, ylim, ycolumn, type, clip, col.response, pch.response,
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
        plot.contour <- function()
        {
            contour(x1, x2, y.predict, main=main,
                xlab=pred.names[i1], ylab=pred.names[i2],
                cex.lab=1.1 * cex.lab, cex.axis=cex.lab,
                labcex=.8 * cex.lab * par("cex"), ...)
            if(col.response != 0)
                points(x[,i1], x[,i2], pch=pch.response, col=col.response, cex=cex)
            NULL
        }
        plot.image <- function()
        {
            image(x1, x2, y.predict, main=main, col=col.image,
                xlab=pred.names[i1], ylab=pred.names[i2],
                cex.lab=1.1 * cex.lab, cex.axis=cex.lab, ...)
            if(col.response != 0)
                points(x[,i1], x[,i2], pch=pch.response, col=col.response, cex=cex)
            NULL
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

# Given a formula (as string), return a string with the "naked" predictors.
#
# Example: y ~ x9 + ns(x2,4) + s(x3,x4,df=4) + x5:sqrt(x6)
# becomes: y ~ x9 + x2 + x3 + x4 + x5 + x6
# which will later result in a model.matrix with columns x9 x2 x3 x4 x5 x6.
#
# This routine is not infallible but works for the commonly used formulas.

strip.formula.string <- function(form)
{
    gsubi <- function(pat, rep, x) gsub(pat, rep, x, ignore.case=TRUE)

    igrep <- grep("[a-zA-Z._][a-zA-Z._0-9]\\$", form)   # check for "ident$"
    if(length(igrep) > 0) {
        # TODO formula has vars with $, this confuses predict() later, why?
        # they cause "Warning: after calling plotmo.predict, y has the wrong length"
        stop1("plotmo: names with \"$\" are not yet supported\n",
            "The unacceptable formula is ", form)
    }
    form <- strip.white.space(form)
    args <- gsubi(".*~", "", form)                  # extract everything after ~

    # We don't want to mess with anything in [square brackets]
    # Doing that with regular expressions is tricky, so we adopt this approach:
    # change "+" "-" "," in square brackets to #PLUS# #MINUS# #COMMA# to protect them

    args <- gsubi("(\\[.*)\\+(.*\\])", "\\1#PLUS#\\2", args)
    args <- gsubi("(\\[.*)\\-(.*\\])", "\\1#MINUS#\\2", args)
    args <- gsubi("(\\[.*)\\,(.*\\])", "\\1#COMMA#\\2", args)

    args <- gsubi("[-*/:]", "+", args)              # replace - / * : with +

    # next two gsubs allow us to retain "x=x1" but drop "df=99" from "bs(x=x1, df=99)"

    args <- gsubi("\\([a-z._0-9$]+=", "(", args)    # replace "(ident=" with "("
    args <- gsubi("[a-z._0-9$]+=[^,)]+", "", args)  # delete "ident=any"

    # replace ",ident" with ")+f(ident", thus "s(x0,x1)" becomes "s(x0)f(x1)"

    args <- gsubi(",([a-z._])", ")+s(\\1", args)

    args <- gsubi("[a-z._0-9$]*[(]", "", args)      # remove "ident("
    args <- gsubi("[,)][^+-]*", "", args)           # remove remaining ",arg1,arg2)"

    # change #PLUS# etc. back to what they where
    args <- gsubi("#MINUS#", "-", args)
    args <- gsubi("#PLUS#", "+", args)
    args <- gsubi("#COMMA#", ",", args)

    # workaround for "error: invalid type (list) for variable 'trees[, -3]'"
    # for a <- earth(trees[,3] ~ as.matrix(trees[,-3])); plotmo(a)
    #
    # TODO removed because although it fixes that problem we still get
    # Warning: 'newdata' had 10 rows but variable(s) found have 31 rows
    #
    # if(is.list(eval.parent(parse(text=args, n=4))))
    #   args<-paste("as.matrix(", args, ")")

    response <- ""
    if(length(grep("~", form)))                     # if there is a response
        response <- gsubi("~.*", "~", form)         # then extract all before ~

    # FIXED 7 Dec 2007 reported by Joe Retzer
    # collapse possible multiple element response and args into a single string

    strip.white.space(paste(response[1], paste(args, collapse=" "), collapse=" "))
}

# Given the term.labels, return an npairs x 2 matrix specifying
# which predictors are pairs. The elements in the returned matrix are col indices of x.
#
# It works like this: extract substrings from each term.label that look
# like predictor pairs and qualify them as a valid pair if both predictors
# in the pair are in pred.names.
#
# The following combos of x1 and x2 are considered pairs: "x1*x2" "x1:x2" "s(x1,x2)"
#
# This routine is not infallible but works for the commonly used formulas.

get.pairs.from.term.labels <- function(term.labels, pred.names, trace=TRUE)
{
    if(trace > 0)
        cat("term.labels:", term.labels, "\n")
    Pairs <- matrix(0, nrow=0, ncol=2)          # no pairs
    for(i in 1:length(term.labels)) {
        s <- strip.white.space(term.labels[i])
        s <- gsub("[+*/,]", ":", s)             # replace + * / , with :
        s <- gsub("=[^,)]+", "", s)             # delete "=any"

        # get the indices of expressions of the form "ident1:ident2"
        igrep <- gregexpr(
            "[a-zA-Z._][a-zA-Z._0-9$]*:[a-zA-Z._][a-zA-Z._0-9$]*", s)[[1]]

        if(trace > 0)
            cat("considering", s)

        if(igrep[1] > 0) for(i in seq_along(igrep)) {
            # extract the i'th "ident1:ident2" into Pair
            start <- igrep[i]
            stop <- start + attr(igrep, "match.length")[i] - 1
            Pair <- substr(s, start=start, stop=stop)
            Pair <- strsplit(Pair, ":")[[1]]    # Pair is now c("ident1","ident2")
            ipred1 <- which(pred.names == Pair[1])
            ipred2 <- which(pred.names == Pair[2])
            if(trace > 0)
                cat("->", Pair, "at", if(length(ipred1)) ipred1 else NA,
                    if(length(ipred2)) ipred2 else NA)
            if(length(ipred1) == 1 && length(ipred2) == 1 && Pair[1] != Pair[2])
                Pairs <- c(Pairs, ipred1, ipred2)
        }
        if(trace > 0)
            cat("\n")
    }
    unique(matrix(Pairs, ncol=2, byrow=TRUE))
}

get.subset <- function(object, trace) # called by get.plotmo.x.default and get.plotmo.y.default
{
    subset. <- object$subset
    if(is.null(subset.)) {
        # the n=3 takes us to the caller of plotmo
        subset. <- try(eval.parent(object$call$subset, n=3), silent=TRUE)
        if(is.try.error(subset.))
            subset. <- NULL
        #TODO revisit the following, it converts function (x, ...) UseMethod("subset")
        else if(typeof(subset.) == "closure")
            subset. <- NULL
    }
    if(!is.null(subset.) && trace > 0) {
        cat("subset length " , length(subset.), sep="")
        try(cat(" min", min(subset.), "max", max(subset.)), silent=TRUE)
        cat(" values ")
        for(i in 1:min(10, length(subset.)))
            cat(subset.[i], "")
        cat("...\n")
    }
    subset.
}

#------------------------------------------------------------------------------
# plotmo.prolog gets called at the start of plotmo

plotmo.prolog <- function(object, object.name) UseMethod("plotmo.prolog")

plotmo.prolog.default <- function(object, object.name)
{
    # Here we just establish with some sort of plausibility that object is a model obj.
    # The general idea is to let the user know what is going on if plotmo fails later.

    if(!is.list(object))
        stop1("'", object.name, "' is not a model object")
    else if(length(coef(object)) == 1)
        warning1("'", object.name, "' appears to be an intercept only model")

    NULL
}

#------------------------------------------------------------------------------
plotmo.predict <- function(object, newdata, type, se.fit,
                            pred.names, ipred1, ipred2, trace)
{
    if(trace > 0) {
        if(ipred2 == 0)
            cat("\nplotmo.predict for predictor \"", pred.names[ipred1], "\" ", sep="")
        else
            cat("\nplotmo.predict for predictors \"",
                pred.names[ipred1], "\" and \"", pred.names[ipred2], "\" ", sep="")
        if(se.fit)
            cat("se.fit=TRUE ")
        cat("with newdata[", NROW(newdata), ",", NCOL(newdata), "]:\n", sep="")
        print(head(newdata, 3))
    }
    UseMethod("plotmo.predict")
}

plotmo.predict.default <- function(object, newdata, type, se.fit,
                                   pred.names, ipred1, ipred2, trace)
{
    # We do our own error handling to give some context to the error messages.
    # TODO but this means that traceback() is not helpful
    try1 <- try(
        if(se.fit)
            predict(object, newdata=newdata, trace=trace, type=type, se.fit=TRUE)
        else
            predict(object, newdata=newdata, trace=trace, type=type))

    if(is.try.error(try1)) {
        cat("\n")
        stop1("Call to predict.", class(object)[1], " failed")
    }
    try1
}

#------------------------------------------------------------------------------
# Return the data matrix for the given object with the response deleted.
#
# If the model has a call$formula, the columns of the returned matrix are in the same
# order as the predictors in the formula.
#
# The default function tries hard to get x regardless of the model.
# Note that the alternative approach of simply calling the standard
# model.matrix wouldn't get us what we want here because it can return
# columns with headings like "ns(x3,4)" whereas we want the "naked" predictor x3.
#
# The n=2 and n=3 in the calls to eval.parent() take us to the caller of plotmo.
#
# TODO get.plotmo.x.default uses a nasty "formula string manipulation"
#      hack, must be a better way?

get.plotmo.x <- function(object, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.plotmo.x\n\n")
    UseMethod("get.plotmo.x")
}

get.plotmo.x.default <- function(
    object = stop("no 'object' arg"),
    trace  = FALSE)
{
    # get x by calling model.frame() with a stripped formula

    get.x.from.formula <- function(object, trace)
    {
        Call <- object$call
        if(is.null(Call))
            return(NULL)    # error will be reported later
        m <- match(c("formula", "data"), names(Call), 0)
        if(all(m == 0))
            return(NULL)
        Call <- Call[c(1, m)]
        Call[[1]] <- as.name("model.frame")
        # TODO it would be nice to use whatever na handling the original model function
        # used, but there seems to be no general way of knowing what that is.
        # In the meantime the following hack suffices for my purposes.
        Call$na.action <- if(inherits(object, "rpart")) na.pass else na.fail
        form <- Call$formula
        if(is.null(form))
            return(NULL)
        # following "if" is needed for: form <- Volume ~ .; earth(form, data=trees)
        # fixes bug reported by Martin Maechler and Michael Amrein
        if(typeof(form) != "language")
            form <- eval.parent(form, n=3)
        formula.as.string <- format(form)
        stripped.formula <- strip.formula.string(formula.as.string)
        if(trace > 0)
            cat("formula ", formula.as.string, "\n",
                "stripped formula ", stripped.formula, "\n", sep="")
        Call$formula <- parse(text=stripped.formula)[[1]]
        Call$data <- get.formula.data(object, Call$data, FALSE, trace)
        if(trace > 0)
            my.print.call("about to call ", Call)
        if(length(grep(".+~", stripped.formula))) { # has response?
            x <- try(eval.parent(Call, n=3)[,-1])   # then eval without the response
        } else {
            warning1("formula has no response variable, formula is ", stripped.formula)
            x <- try(eval.parent(Call, n=3))
        }
        if(is.try.error(x)) {
            if(length(grep("missing", x)))
                stop1("could not evaluate formula")
            else
                stop1("could not evaluate formula (variables were deleted?)")
        }
        if(NCOL(x) == 1) {
            # if one predictor, model.matrix returns a vec with no colname, so fix it
            x <- data.frame(x)
            colnames(x) <- strip.formula.string(attr(object$terms, "term.labels")[1])
        }
        x
    }
    badx <- function(x, check.colnames)
    {
        is.null(x) || is.try.error(x) || NROW(x) == 0 ||
            (check.colnames && is.null(colnames(x)))
    }
    # get.plotmo.x.default starts here

    if(!is.list(object))
        stop1("get.plotmo.x.default cannot get x matrix --- object is not a model object")

    # look for an x with column names

    try.error.message <- NULL
    # use of brackets rather than $ below prevents incorrect partial matching
    x <- object[["x"]]
    if(!badx(x, TRUE) && trace > 0)
        cat("got x with colnames from object$x\n")
    if(badx(x, TRUE)) {
        x <- get.x.from.formula(object, trace > 0)
        if(!badx(x, TRUE) && trace > 0)
            cat("got x with colnames from object$call$formula\n")
    }
    if(badx(x, TRUE)) {
        x <- try(eval.parent(object$call[["x"]], n=2), silent=TRUE)
        if(!badx(x, TRUE) && trace > 0)
            cat("got x with colnames from object$call$x\n")
        if(is.try.error(x))
            try.error.message <- x
    }
    # if don't have an x with colnames look for one without colnames
    # the call to as.data.frame below will add V1, V2, ... colnames if necessary

    if(badx(x, TRUE)) {
        x <- object[["x"]]
        if(!badx(x, FALSE) && trace > 0)
            cat("got x without colnames from object$x\n")
        if(badx(x, FALSE)) {
            x <- get.x.from.formula(object, trace)
            if(!badx(x, FALSE) && trace > 0)
                cat("got x without colnames from object$call$formula\n")
        }
        if(badx(x, FALSE)) {
            x <- try(eval.parent(object$call[["x"]], n=2), silent=TRUE)
            if(!badx(x, FALSE) && trace > 0)
                cat("got x without colnames from object$call$x\n")
            if(is.try.error(x))
                try.error.message <- x
        }
    }
    if(badx(x, FALSE)) {
        if(trace > 0) {
            cat("Looked unsuccessfully for an x in the following places:\n")
            cat("\nobject$x:\n")
            print(head(object$x, 3))
            cat("\nobject$call$formula:\n")
            dput(object$call$formula)
            cat("\nobject$call$x:\n")
            if(is.null(try.error.message))
                print(head(eval.parent(object$call$x, n=2), 3))
            else
                cat(gsub("Error in ", "", try.error.message[1]))
            cat("\n")
        }
        stop1("get.plotmo.x.default cannot get x matrix --- ",
              "tried object$x, object$call$formula, and object$call$x")
    }
    x <- as.data.frame(x)
    weights <- weights(object)
    if(!is.null(weights) && any(abs(weights - weights[1]) > 1e-8))
        warning1("'weights' are not supported by 'plotmo', ignoring them")

    subset. <- get.subset(object, trace > 0)
    if(!is.null(subset.)) {
        check.index.vec("subset", subset., x, check.empty=TRUE, allow.duplicates=TRUE)
        x <- x[subset., , drop=FALSE]
    }
    x
}

#------------------------------------------------------------------------------
# get.plotmo.y is similar to model.response but can deal with models
# created without a formula
# The n=2 and n=3 in the calls to eval.parent() take us to the caller of plotmo.

get.plotmo.y <- function(
    object = stop("no 'object' arg"),
    ycolumn,            # which column of response to use if response has multiple cols
    expected.len,
    trace)
{
    if(trace > 0)
        cat("\n--get.plotmo.y\n\n")
    UseMethod("get.plotmo.y")
}

get.plotmo.y.default <- function(
    object = stop("no 'object' arg"),
    ycolumn,            # which column of response to use if response has multiple cols
    expected.len,
    trace)
{
    get.y.from.formula <- function(object)
    {
        Call <- object$call
        if(is.null(Call))
            return(NULL)    # error will be reported later
        m <- match(c("formula", "data"), names(Call), 0)
        if(all(m == 0))
            return(NULL)
        Call <- Call[c(1, m)]

        # replace the formula with the stripped formula
        form <- Call$formula
        if(is.null(form))
            return(NULL)
        if(typeof(form) != "language")
            form <- eval.parent(form, n=3)
        formula.as.string <- paste(format(form), collapse=" ")
        stripped.formula <- strip.formula.string(formula.as.string)
        Call$formula <- parse(text=stripped.formula)[[1]]
        if(trace > 0)
            cat("formula ", formula.as.string, "\n",
                "stripped formula ", stripped.formula, "\n", sep="")

        Call$data <- get.formula.data(object, Call$data, TRUE, trace)
        Call[[1]] <- as.name("model.frame")
        Call$na.action <- if(inherits(object, "rpart")) na.pass else na.fail # TODO hack
        stripped.formula <- strip.formula.string(formula.as.string)
        if(trace > 0)
            my.print.call("about to call ", Call)
        Call <- try(eval.parent(Call, n=3))
        if(!is.try.error(Call))
            model.response(Call, type="any")
    }
    bady <- function(y)
    {
        is.null(y) || is.try.error(y)
    }
    # get.plotmo.y.default starts here
    try.error.message <- NULL
    y <- object[["y"]]
    if(!bady(y) && trace > 0)
        cat("got y from object$y\n")
    if(bady(y)) {
        y <- get.y.from.formula(object)
        if(!bady(y) && trace > 0)
            cat("got y from object$call$formula\n")
    }
    if(bady(y)) {
        y <- try(eval.parent(object$call[["y"]], n=2), silent=TRUE)
        if(!bady(y) && trace > 0)
            cat("got y from object$call$y\n")
        if(is.try.error(y))
            try.error.message <- y
    }
    if(bady(y)) {
        if(trace > 0) {
            cat("Looked unsuccessfully for y in the following places:\n")
            cat("\nobject$y:\n")
            print(head(object$y, 3))
            cat("\nobject$call$formula:\n")
            dput(object$call$formula)
            cat("\nobject$call$y:\n")
            if(is.null(try.error.message))
                print(head(eval.parent(object$call$y, n=2), 3))
            else
                cat(gsub("Error in ", "", try.error.message[1]))
            cat("\n")
        }
        stop1("get.plotmo.y.default cannot get y --- ",
              "tried object$call$formula, object$call$y, and object$y")
    }
    subset. <- get.subset(object, trace)
    y <- check.and.print.y(y, "get.plotmo.y", ycolumn, expected.len, trace, subset.)
    if(!is.null(subset.)) {
        check.index.vec("subset", subset., y, check.empty=TRUE, allow.duplicates=TRUE)
        y <- y[subset.]
    }
    y
}

#------------------------------------------------------------------------------
# Return a vector of indices of predictors for degree1 plots
# The indices are col numbers in the x matrix
# The default method simply returns the indices of all predictors
# See also get.singles.earth

get.singles <- function(object, x, degree1, pred.names, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.singles\n\n")
    UseMethod("get.singles")
}

get.singles.default <- function(object, x, degree1, pred.names, trace)
{
    ifirst <- if(pred.names[1]=="(Intercept)") 2 else 1 # delete intercept, if any
    if(is.character(degree1) && !is.na(pmatch(degree1, "all")))
        ifirst:length(pred.names) # degree1 = "all"
    else {
        check.index.vec("degree1", degree1, pred.names)
        (ifirst:length(pred.names))[degree1]
    }
}

get.formula.data <- function(object, data.name, get.y, trace)
{
    data.is.good <- function(...)
    {
        # the length test is necessary for lm which saves x as an
        # empty list if its x arg is FALSE, don't know why
        good <- !is.null(x) && length(x)
        if(good && trace > 0)
            cat("get.formula.data: got",
                if(get.y) "y" else "x", "from", paste(..., sep=""), "\n")
        good
    }
    xname <- if(get.y) "y" else "x"
    # use of brackets rather than $ below prevents incorrect partial matching
    x <- object[["data"]]
    if(!data.is.good("object$data")) {
        x <- object[[xname]]
        if(!data.is.good("object$", xname)) {
            if(!is.null(data.name)) {
                x <- eval.parent(data.name, n=4)
                if(!data.is.good("data passed to original call to ", class(object)))
                    stop1("the original data \"", data.name,
                          "\" is no longer available",
                          if(inherits(object, "earth"))
                                " (use keepxy=TRUE)"
                          else if(inherits(object, "lm"))
                                paste(" (use ", xname, "=TRUE)", sep="")
                          else  "")
            }
        }
    }
    x
}

#------------------------------------------------------------------------------
# Each row of the returned Pairs matrix is the indices of two predictors
# for a degree2 plot
# The indices are col numbers in the x matrix
# See also get.pairs.earth

get.pairs <- function(object, x, degree2, pred.names, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.pairs\n\n")
    UseMethod("get.pairs")
}

# Predictors x1 and x2 are considered paired if they appear in the formula
# in forms such as x1:x2 or I(x1*x2) or s(x1,x2)

get.pairs.default <- function(object, x, degree2, pred.names, trace=FALSE)
{
    if(is.character(degree2) && !is.na(pmatch(degree2, "all"))) {
        used.vars <- get.singles.default(object, x, "all", pred.names, trace)
        if(length(used.vars) == 0)
            return(matrix(0, nrow=0, ncol=2)) # no pairs
        col1 <- rep(used.vars, times=length(used.vars))
        col2 <- rep(used.vars, each=length(used.vars))
        Pairs <- cbind(col1, col2)
        Pairs <- Pairs[col1 != col2, , drop=FALSE]
        return(unique(t(apply(Pairs, 1, sort)))) # remove duplicate pairs
    }
    Pairs <- matrix(0, nrow=0, ncol=2)   # no pairs
    term.labels <- NULL
    if(!is.null(object$call$formula)) {
        form <- object$call$formula
        if(typeof(form) != "language")
            form <- eval.parent(form, n=2)
        form <- as.formula(form)
        data <- get.formula.data(object, object$call$data, FALSE, trace)
        terms <- terms(form, data=data)
    }
    if(!is.null(terms))
        term.labels <- attr(terms, "term.labels")
    if(!is.null(term.labels))
        Pairs <- get.pairs.from.term.labels(term.labels, pred.names, trace)
    else {
        if(trace > 0)
            cat("no degree2 plots because no $call$formula$term.labels\n")
        if(is.specified(degree2))
            warning1("'degree2' specified but no degree2 plots ",
                     "(because no $call$formula$term.labels)")
    }
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("'degree2' specified but no degree2 plots")
    if(nrow(Pairs) > 0) {
        check.index.vec("degree2", degree2, Pairs)
        Pairs <- Pairs[degree2, , drop=FALSE]
    }
    Pairs
}

#------------------------------------------------------------------------------
# define method function because don't want bogus warnings from plotmo.prolog.default
plotmo.prolog.bagEarth <- function(object) NULL

get.pairs.bagEarth <- function(object, x, degree2, pred.names, trace)
{
    pairs <- matrix(0, nrow=0, ncol=2)
    for(i in 1:length(object$fit))     # TODO could probably be vectorized
        pairs <- rbind(pairs, get.pairs(object$fit[[i]], x, degree2, pred.names, trace))
    pairs <- unique(pairs)
    pairs[order(pairs[,1], pairs[,2]),]
}
#------------------------------------------------------------------------------

get.singles.rpart <- function(object, x, degree1, pred.names, trace)
{
    # get all the variables that are actually used in the tree
    irow <- as.integer(row.names(object$frame))
    var.names <- character(length=max(irow))
    var.names[irow] <- as.character(object$frame$var)
    ivar <- charmatch(var.names, pred.names)
    is.split <- !is.na(ivar) & ivar > 0 # same as var.names != "<leaf>" & var.names !=""
    if(sum(is.split) == 0)
        stop1("the rpart tree has no splits")
    ivar <- unique(ivar[is.split])
    if(is.character(degree1) && !is.na(pmatch(degree1, "all")))
        degree1 <- 1:length(ivar)
    check.index.vec("degree1", degree1, ivar)
    ivar[degree1]
}
# We consider rpart variables paired if one is the direct parent of the
# other in the tree.

get.pairs.rpart <- function(object, x, degree2, pred.names, trace)
{
    if(is.character(degree2) && !is.na(pmatch(degree2, "all"))) {
        # TODO this may be redundant for rpart models
        used.vars <- get.singles.rpart(object, x, "all", pred.names, trace)
        if(length(used.vars) == 0)
            return(matrix(0, nrow=0, ncol=2)) # no pairs
        col1 <- rep(used.vars, times=length(used.vars))
        col2 <- rep(used.vars, each=length(used.vars))
        Pairs <- cbind(col1, col2)
        Pairs <- Pairs[col1 != col2, , drop=FALSE]
        return(unique(t(apply(Pairs, 1, sort)))) # remove duplicate pairs
    }
    irow <- as.integer(row.names(object$frame))
    var.names <- character(length=max(irow))
    var.names[irow] <- as.character(object$frame$var)
    ivar <- charmatch(var.names, pred.names)
    # following is the same as var.names != "<leaf>" & var.names !=""
    is.split <- !is.na(ivar) & ivar > 0
    if(sum(is.split) == 0)
        stop1("the rpart tree has no splits")
    Pairs <- NULL
    for(i in 1:length(ivar)) {
        if(is.split[i]) {
            left <- 2 * i
            if(is.split[left] && ivar[i] != ivar[left])
                Pairs <- c(Pairs, ivar[i], ivar[left])
            right <- left + 1
            if(is.split[right] && ivar[i] != ivar[right])
                Pairs <- c(Pairs, ivar[i], ivar[right])
        }
    }
    if(is.null(Pairs))
        Pairs <- matrix(0, nrow=0, ncol=2)
    else
        Pairs <- matrix(Pairs, ncol=2, byrow=TRUE)
    Pairs <- unique(t(apply(Pairs, 1, sort))) # remove duplicate pairs
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("'degree2' specified but no degree2 plots")
    if(nrow(Pairs) > 0) {
        check.index.vec("degree2", degree2, Pairs)
        Pairs <- Pairs[degree2, , drop=FALSE]
    }
    Pairs
}

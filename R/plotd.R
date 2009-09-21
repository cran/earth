# plotd.R: plot densities of class conditional predicted values
#
# TODO: allow newdata so can plot not only with the training data
# TODO: allow freq arg for histograms

plotd <- function(obj,      # obj is a model object
    hist   = FALSE,         # FALSE to use density(), TRUE to use hist()
    type   = "response",    # passed on to predict
    nresponse = NULL,       # which response, for multiple response models, NULL for all
    dichot = FALSE,
    trace  = FALSE,
    xlim   = NULL,          # NULL means auto
    ylim   = NULL,          # NULL means auto
    jitter = FALSE,
    main   = NULL,          # graph caption
                            #   "string"  string
                            #   ""        no caption
                            #   NULL      generate a caption from x$call
    xlab = "Predicted Value",
    ylab = if(hist) "Count" else "Density",
    lty  = 1,               # linetypes for the plotted lines
    col  = c("grey70", 1, "lightblue", "brown", "pink", 2, 3, 4), # cols for plotted lines
    fill = if(hist) col[1] else 0, # fill color for first hist/density plot
    breaks = "Sturges",     # following passed on to hist, only used if hist=TRUE
    labels = FALSE,
    kernel = "gaussian",   # following passed on to density, only used if hist=FALSE
    adjust = 1,
    zero.line = FALSE,
    legend = TRUE,          # TRUE to draw a legend
    legend.names = NULL,    # NULL means auto, else specify a vector of strings,
    legend.pos = NULL,      # NULL means auto, else specify c(x,y) in user coords
    legend.cex = .8,        # cex for legend
    legend.bg = "white",    # bg for legend
    legend.extra = FALSE,   # print number in each class in legend
    vline.col = 0,          # color of vertical line, use NULL or 0 for no line
    vline.thresh = .5,      # horizontal position of vertical line
    vline.lty = 1,          # lty of vertical line
    vline.lwd = 1,          # lwd of vertical line
    err.thresh=vline.thresh, # thresh for "error areas"
    err.col=0,              # col shading of "error areas"
    err.border=0,
    err.lwd=1,
    xaxt = "s",
    yaxt = "s",
    xaxis.cex = 1,
    sd.thresh = 0.01,
    ...)                    # passed on to predict
{
    # add histogram labels --- lifted from plot.histogram and
    # tweaked to (i) use cex and (ii) not draw zero counts
    add.labels <- function(x, labels, cex)
    {
        if((is.logical <- is.logical(labels) && labels) || is.character(labels)) {
            stopif(is.null(x$counts)) # plotd hist supports only counts, not densities
            if(is.logical) {
                labels <- format(x$counts)
                labels[x$counts == 0] <- ""
            }
            text(x$mids, x$counts, labels=labels, adj=c(0.5, -0.5), cex=cex)
        }
    }
    # plotd starts here

    if(typeof(obj) != "list")
        stop1("'", deparse(substitute(obj)), "' is not a model obj")

    yhat.per.class <- get.yhat.per.class(obj, type, nresponse, dichot, trace, ...)
    nclasses <- length(yhat.per.class)

    # get densities

    densities <- NULL
    for(iclass in 1:nclasses)
        if(!hist)
            densities[[iclass]] <- density(yhat.per.class[[iclass]],
                                           kernel=kernel, adjust=adjust)
        else {
            densities[[iclass]] <- hist(yhat.per.class[[iclass]],
                                        breaks=breaks, plot=FALSE)
            # need x and y components so hist can be treated uniformly with density
            densities[[iclass]]$x <- densities[[iclass]]$breaks
            densities[[iclass]]$y <- densities[[iclass]]$counts
        }

    # get x limits of plot

    if(is.null(xlim)) {
        min1 <- Inf
        max1 <- -Inf
        for(iclass in 1:nclasses) {
             min1 <- min(min1, densities[[iclass]]$x)
             max1 <- max(max1, densities[[iclass]]$x)
        }
        xlim <- c(min1, max1)
    }
    if(length(xlim) != 2)
        stop1("length(xlim) != 2")
    xspan <- xlim[2] - xlim[1]

    # sanity check the ranges of each class, issue warnings if need be

    degenerate <- logical(nclasses)
    for(iclass in 1:nclasses) {
        range <- range(yhat.per.class[[iclass]])
        if(sd(yhat.per.class[[iclass]]) < sd.thresh) {
            warning1("standard deviation of \"", names(yhat.per.class)[iclass],
                    "\" density is ", sd(yhat.per.class[[iclass]]),
                    ",  density is degenerate?")
            degenerate[iclass] <- TRUE
        }
    }
    # add jitter and get ymax for plot

    ymax <- 1e-6
    if(is.logical(jitter) && jitter)
        jitter <- xspan / 100
    for(iclass in 1:nclasses) {
        if(jitter) {
            if(hist)
                densities[[iclass]]$breaks <-
                    densities[[iclass]]$breaks + iclass * jitter
            else
                densities[[iclass]]$x <- densities[[iclass]]$x + iclass * jitter
        }
        ymax <- max(ymax, densities[[iclass]]$y)
    }
    if((is.logical(labels) && labels) || is.character(labels))
        ymax <- 1.1 * ymax # hack to make space for labels
    if(!is.null(ylim)) {
        if(length(ylim) != 2)
            stop1("length(ylim) != 2")
        ymax <- ylim[2]
        if(ymax <= 0)
            stop1("ylim[2] <= 0")
        if(ylim[1] != 0)
            warning1("ignoring ylim[1], treating it as 0")
    }
    # expand lty and other arguments if necessary

    if(length(lty) < nclasses)
        lty <- rep(lty, nclasses)
    if(length(col) < nclasses)
        col <- rep(col, nclasses)

    # main title

    main.org <- main
    if(is.null(main))
        main <- paste(deparse(substitute(obj)), " ", paste(type, collapse=" "),
                      if(missing(nresponse)) ""
                      else paste(" nresp=", paste(nresponse, collapse=","), sep=""),
                      if(missing(dichot)) ""
                      else paste(" dichot=", dichot, sep=""),
                      sep="")

    main <- show.caption(main, if(is.null(main.org)) trim=.5 else trim=0, show=FALSE)

    # we draw our own x axis for type="class"
    # xlims are wrong for histograms if a density is degenerate hence the test
    # TODO weird behaviour in hist? hist gives a 0 lower x val if density degenerate
    draw.own.axis <- hist && xaxt != "n" && !is.na(pmatch(type, "class")) &&
                     !any(degenerate)
    if(draw.own.axis)
        xaxt <- "n"

    # plot the first graph

    ifirst <- 1 # index of first non-degenerate class, 1 if all degenerate
    for(iclass in 1:nclasses)
        if(!degenerate[iclass]) {
            ifirst <- iclass
            break
        }
    if(hist) { # plot.histogram
        plot(densities[[ifirst]], xlim=xlim, ylim=c(0, ymax),
             main=main, xlab=xlab, ylab=ylab, lty=lty[ifirst],
             border=col[ifirst], xaxt=xaxt, yaxt=yaxt,
             col=if(ifirst==1) fill else 0) # fill color
        add.labels(densities[[ifirst]], labels, legend.cex)
    } else {   # plot.density
        plot(densities[[ifirst]], xlim=xlim, ylim=c(0, ymax), col=col[ifirst],
             main=main, xlab=xlab, ylab=ylab, lty=lty[ifirst],
             zero.line=zero.line, xaxt=xaxt, yaxt=yaxt)
        if(ifirst == 1 && fill != 0)
            polygon(densities[[1]], col=fill, border=col[1])
    }
    if(draw.own.axis) {
        at <- xlim[1]:xlim[2]
        # hack to adjust leftmost label TODO not quite right, but ok
        at[1] <- at[1] + strwidth(names(yhat.per.class)[1], "user")
        mtext(names(yhat.per.class), side=1, at=at, font=2, adj=1, cex=xaxis.cex)
    }
    # optional error region shading

    if(!hist && any(err.col != 0))
        add.err.col(densities, err.thresh, err.col, err.border, err.lwd)

    # overlay the graphs

    for(iclass in 1:nclasses)
        if(!degenerate[iclass]) {
            if(!hist)  # lines.density
                lines(densities[[iclass]], col=col[iclass], lty=lty[iclass])
            else {     # lines.histogram
                lines(densities[[iclass]], col=NULL, 
                      border=col[iclass], lty=lty[iclass])
                add.labels(densities[[iclass]], labels, legend.cex)
            }
        }
    # optional vertical line at vline.thresh

    if(!is.null(vline.col) && vline.col != 0)
        abline(v=vline.thresh, col=vline.col, lty=vline.lty, lwd=vline.lwd)

    # Redo optional error region shading if it has borders, because the
    # borders go on top of the other plotted lines.

    if(any(err.border != 0))
        add.err.col(densities, err.thresh, err.col, err.border, err.lwd)

    # optional legend

    if(legend)
        draw.legend(densities, degenerate, yhat.per.class, ymax,
                    hist, xlim, col, fill, lty, legend.names,
                    legend.pos, legend.cex, legend.bg, legend.extra)
    invisible(yhat.per.class)
}

is.numlog <- function(x)
    is.numeric(x) || is.logical(x)

# get the original observed response (it's needed to determine correct classes)

get.observed.response <- function(obj)
{
    if(!is.null(obj$call$formula)) {
        # get y from formula and data used in original call
        data <- get.update.arg(NULL, "data", obj, FALSE)
        Call <- obj$call
        m <- match(c("formula", "data", "na.action", "subset"), names(Call), 0)
        mf <- Call[c(1, m)]
        mf[[1]] <- as.name("model.frame")
        mf <- eval.parent(mf, n=3)      # n=3 takes us to the caller of plotd
        y <- model.response(mf, "any")  # "any" means factors are allowed
        if(NCOL(y) == 1 && is.numlog(y)) {
            # turn into a matrix so we have the column name
            names(y) <- NULL # we don't need row names
            y <- as.matrix(y)
            colnames(y) <- colnames(mf)[attr(obj$terms,"response")]
        }
    } else if(inherits(obj, "lda")) { # hack for "lda", get grouping arg
        y <- eval.parent(obj$call[[3]], n=3) # n=3 takes us to the caller of plotd
        # sanity check
        if(NCOL(y) != 1 || length(y) < 3 || (!is.numeric(y) && !is.factor(y)))
            stop1("cannot get \"grouping\" argument from obj$call")
    } else
        y <- get.update.arg(NULL, "y", obj, trace=FALSE, reeval=FALSE)

    if(inherits(obj, "lda"))
        y <- as.factor(y)  # to make plotd handle response appropriately
    y
}
# special handling for MASS lda predicted response

get.lda.yhat <- function(yhat, type, trace)
{
    if(trace)
        cat("\nSpecial handling of \"type\" argument for lda object\n")

    switch(match.choices(type,
                         c("response", "ld", "class", "posterior"), "type"),
           yhat$x,              # response (default)
           yhat$x,              # ld
           yhat$class,          # class
           yhat$posterior)      # posterior
}
# return the predictions for each class, in a list with each class named

get.yhat.per.class <- function(obj, type, nresponse, dichot, trace, ...)
{
    check.min <- function(x, ...)
    {
        len <- length(x)
        if(len == 0)
            warning1("no occurrences of ", paste(..., sep=""),
                     " in the observed response")
        else if(len < 3) # 3 is arbitrary
            warning1("only ", len, " occurrences of ", paste(..., sep=""),
                     " in the observed response")
    }
    get.binary.class.names <- function(yhat, fitted.values, last.resort)
    {
        if(!is.null(colnames(yhat)))
            c(paste("not", colnames(yhat)), colnames(yhat))
        else if(!is.null(colnames(fitted.values)))
            c(paste("not", colnames(fitted.values)), colnames(fitted.values))
        else
            last.resort
    }
    get.class.names <- function(y, yhat, fitted.values)
    {
        ynames <- paste("response", 1:ncol(yhat), sep="")
        if(length(colnames(yhat)) == ncol(yhat))
             class.names <- colnames(yhat)
        else if(length(colnames(y)) == ncol(y))
            class.names <- colnames(y)
        else if(length(fitted.values(y)) == fitted.values(y))
            class.names <- colnames(y)
        else
            class.names <- ynames
        # fill in missing names, if necessary
        which. <- which(class.names == "")
        if(length(which.))
            class.names[which.] <- ynames[which.]
        class.names
    }
    # return names1 but with yhat column names prefixed if necessary
    get.prefixed.names <- function(names1, yhat, nresponse)
    {
        stopifnot(NCOL(yhat) == 1)
        if(!is.null(colnames(yhat)) && !is.null(nresponse))
            names1 <- paste(colnames(yhat)[1], names1, sep=": ")
        names1
    }
    # determine the threshold to split classes, a bit of a hack
    get.thresh <- function(y, yname)
    {
        thresh <- 0
        ymin <- min(y)
        if(ymin == 1)
            thresh <- 1
        if(!is.null(colnames(y)))
            yname <- colnames(y)[1]
        if(ymin == thresh) {
            text.le <- sprintf("%s == %g", yname, thresh)
            text.gt <- sprintf("%s != %g",  yname, thresh)
        } else {
            text.le <- sprintf("%s <= %g", yname, thresh)
            text.gt <- sprintf("%s > %g",  yname, thresh)
        }
        list(thresh=thresh, text.le=text.le, text.gt=text.gt)
    }
    trace.yhat.per.class <- function(yhat.per.class, iclass, nchar)
    {
        description <- sprintf("predicted.response.per.class[%-*s]",
                               nchar, names(yhat.per.class)[iclass])
        if(trace > 1)
            print.matrix.info(description, yhat.per.class[[iclass]],
                              details=TRUE, all.rows=TRUE)
        else {
            cat(description, ": ", sep="")
            yhat1 <- yhat.per.class[[iclass]]
            names(yhat1) <- NULL
            cat(yhat1[1:min(6,length(yhat))])
            if(6 < length(yhat))
                cat(" ...")
        }
        cat("\n")
    }
    cannot.plot <- function(...)
        stop1("cannot plot this kind of response with type=\"",
              type, "\"\n", ...,
              "Additional information:\n  class(observed)=", class(y),
              if(nlevs > 0) paste(" nlevels(observed)=", nlevs, sep="") else "",
              " ncol(observed)=", NCOL(y),
              if(!is.null(colnames(y))) 
                  sprintf(" colnames(observed) %s", paste(colnames(y), collapse=" ")) 
              else 
                  "",
              "\n  class(predicted)=", class(yhat),
              " ncol(predicted)=", NCOL(yhat),
              if(!is.null(colnames(yhat))) 
                  sprintf(" colnames(response) %s", paste(colnames(yhat), collapse=" ")) 
              else 
                  "")

    trace.response.type <- function(observed, predicted, ...) {
        if(trace)
            cat("\nObserved response: ", observed, "\n",
                "Predicted response type=\"", type, "\": ", predicted, "\n",
                "Grouping criterion: ", ..., "\n\n", sep="")
    }
    # get.yhat.per.class starts here
    # nomeclature: y is the observed response, yhat is the predicted response
    y <- get.observed.response(obj)
    if(trace)
        print.matrix.info("y", y, "\nObserved response",
                          details=TRUE, all.rows=trace>1)
    if(!is.character(type))
        stop1("type of \"type\" is not character")
    if(!is.na(pmatch(type, "terms")))
        stop1("type=\"terms\" is not allowed by plotd")
    yhat <- predict(obj, type=type, ...)
    if(inherits(obj, "lda"))
        yhat <- get.lda.yhat(yhat, type, trace)
    if(trace)
        print.matrix.info("yhat", yhat, "\nPredicted response",
                          details=TRUE, all.rows=trace>1)
    if(!is.null(nresponse)) {
        if(!is.numeric(nresponse) || nresponse < 1 || nresponse > NCOL(yhat))
            stop1("illegal nresponse argument, ",
                  "allowed range for this model and type is 1 to ", NCOL(yhat))
        if(NCOL(yhat) > 1) {
            row.names(yhat) = NULL  # needed for fda type="hier" because dup rownames
            yhat <- yhat[, nresponse, drop=FALSE] # drop=FALSE to retain column name
            if(is.data.frame(yhat)) # TODO needed for fda type="hier", why?
                yhat <- as.matrix(yhat)
            if(trace)
                print.matrix.info("yhat", yhat,
                    paste("Predicted response after selecting nresponse", nresponse),
                    details=TRUE, all.rows=trace>1)
        }
    }
    yhat.per.class <- list() # will put per-class predicted vals in here
    ylevs <- levels(y) # will be NULL if y is not a factor
    nlevs <- 0
    if(NCOL(yhat) == 1) {
        #---single column yhat--------------------------------------------------
        if(is.factor(y) && (is.numlog(yhat) || is.factor(yhat))) {
            nlevs <- nlevels(y)
            if(!is.numlog(yhat) && !is.factor(yhat))
                cannot.plot()
            stopifnot(length(ylevs) > 1)
            if(!is.na(pmatch(type, "class")))
                dichot <- FALSE # no dichot for type="class"
            if(length(ylevs) == 2 || dichot) {
                ylev1 = ylevs[1]
                if(length(ylevs) == 2) {
                    other.level.name <- paste(ylevs[2])
                    observed.string <- "two-level factor"
                } else {
                    other.level.name <- paste("not", ylev1)
                    observed.string <- "multi-level factor, but dichot"
                }
                trace.response.type(observed.string, "numeric or logical vector",
                    "CLASS1 predicted[observed == ", ylev1,
                    "], CLASS2 predicted[observed == ", other.level.name, "]")
                yhat.per.class[[1]] <- yhat[y == ylev1]
                check.min(yhat.per.class[[1]], ylev1)
                yhat.per.class[[2]] <- yhat[y != ylev1]
                check.min(yhat.per.class[[2]], other.level.name)
                names(yhat.per.class) <- get.prefixed.names(c(ylev1, other.level.name), 
                                                            yhat, nresponse)
            } else {
                trace.response.type("multi-level factor", "numeric or logical vector",
                                "predicted[observed == level] for ",
                                length(ylevs), " levels")
                for(iclass in 1:length(ylevs)) {
                    lev <- ylevs[iclass]
                    yhat.per.class[[iclass]] <- yhat[y == lev]
                    check.min(yhat.per.class[[iclass]], lev)
                }
                names(yhat.per.class) <- get.prefixed.names(ylevs, yhat, nresponse)
            }
        } else if(NCOL(y) == 2 && is.numlog(y[,1]) && is.numlog(y[,2])) {
            trace.response.type("numeric or logical vector", "two-column numeric",
                            "CLASS1 observed[,1] <= observed[,2], ",
                            "CLASS2 observed[,1] > observed[,2]")
            # split into two classes based on relative sizes of columns of y
            yhat.per.class[[1]] <- yhat[y[,1] <= y[,2]]
            check.min(yhat.per.class[[1]], "observed[,1] <= observed[,2]")
            yhat.per.class[[2]] <- yhat[y[,1] > y[,2]]
            check.min(yhat.per.class[[2]], "observed[,1] > observed[,2]")
            names(yhat.per.class) <-
                get.binary.class.names(yhat, obj$fitted.values, c("FALSE", "TRUE"))
        } else if(NCOL(y) == 1 && is.numlog(y)) {
            th <- get.thresh(y, "response")
            trace.response.type("numeric or logical vector",
                                "numeric or logical vector",
                                "CLASS1 ", th$text.le, ", CLASS2 ", th$text.gt)
            yhat.per.class[[1]] <- yhat[y <= th$thresh]
            check.min(yhat.per.class[[1]], th$text.le)
            yhat.per.class[[2]] <- yhat[y > th$thresh]
            check.min(yhat.per.class[[2]], th$text.gt)
            names(yhat.per.class) <- NULL
            if(th$thresh == 0)
                names(yhat.per.class) <-
                    get.binary.class.names(yhat, obj$fitted.values,
                                           c(th$text.le, th$text.gt))
            else
                names(yhat.per.class) <- c(th$text.le, th$text.gt)
        } else
            cannot.plot()
  } else {
        #---multiple column yhat------------------------------------------------
        if(!is.numeric(yhat))
            cannot.plot()
        if(NCOL(y) == 1 && is.null(ylevs))
            ylevs <- as.numeric(names(table(y))) # use numeric levels like a factor
        nlevs <- length(ylevs)
        if(NCOL(y) == 1 && nlevs == NCOL(yhat)) {
            if(is.factor(y))
                trace.response.type("factor",
                    "multicolumn numeric, ncol(predicted) == nlevels(observed)",
                    "observed==level for each level in observed response")
            else
                trace.response.type("factor",
                    "multicolumn numeric, ncol(predicted) == nbr.of.unique.vals.in.observed",
                    "observed==val for each unique val in observed response")
            for(iclass in 1:ncol(yhat)) {
                lev <- ylevs[iclass]
                yhat.per.class[[iclass]] <- yhat[y == lev, iclass]
                check.min(yhat.per.class[[iclass]], lev)
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop1("no occurrences of ", lev,
                          " in the observed response")
            }
            if(!is.null(colnames(yhat)))
                names(yhat.per.class) <- colnames(yhat)
            else
                names(yhat.per.class) <- ylevs
#         } else if(NCOL(y) == 1) { # nlevs != NCOL(yhat))
#           trace.response.type("factor",
#                   "multicolumn numeric; ncol(predicted) != nlevels(observed)",
#                   "each column of predicted response is a group")
#             for(iclass in 1:ncol(yhat))
#                 yhat.per.class[[iclass]] <- yhat[, iclass]
#             if(!is.null(colnames(yhat)))
#                 names(yhat.per.class) <- colnames(yhat)
#             else
#                 names(yhat.per.class) <-  paste(type, "[,", 1:ncol(yhat), "]", sep="")
        } else if(is.numeric(y) && NCOL(y) == NCOL(yhat)) {
            th <- get.thresh(y, "response")
            trace.response.type("multicolum numeric",
                "multicolumn numeric with same number of columns as observed response",
                th$text.gt, "for each column of observed response")
            for(iclass in 1:ncol(yhat)) {
                yhat.per.class[[iclass]] <- yhat[y[,iclass] > th$thresh, iclass]
                check.min(yhat.per.class[[iclass]], th$text.gt)
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop1("no occurrences of ", th$text.le,
                          " in the observed response")
            }
            names(yhat.per.class) <- get.class.names(y, yhat, obj$fitted.values)
        } else
            cannot.plot("Remedy: use the \"nresponse\" arg to select ",
                "just one column of the predicted response\n")
    }
    nchar <- max(nchar(names(yhat.per.class)))

    for(iclass in 1:length(yhat.per.class)) {
        # density needs numeric
        yhat.per.class[[iclass]] <- as.numeric(yhat.per.class[[iclass]])
        if(trace)
            trace.yhat.per.class(yhat.per.class, iclass, nchar)
    }
    if(trace)
        cat("\n")
    yhat.per.class
}
draw.legend <- function(densities, degenerate, yhat.per.class, ymax,
                    hist, xlim, col, fill, lty, legend.names,
                    legend.pos, legend.cex, legend.bg, legend.extra)
{
    get.legend.pos <- function()
    {
        # take a stab at positioning the legend correctly --
        # on left or right, away from the highest peak
        pos <- c(0,ymax)
        pos[1] <- xlim[1] # place on left side of graph
        max.left <- 0
        max.right <- 0
        xmid <- xlim[1] + (xlim[2] - xlim[1])/2
        for(iclass in 1:nclasses) {
            if(!degenerate[iclass]) {
                den <- densities[[iclass]]
                if(hist)
                    den$x <- den$x[-1]
                x.left <- (den$x >= xlim[1]) & (den$x <= xmid)
                if(sum(x.left))
                    max.left <- max(max.left, den$y[x.left])
                x.right <- (den$x > xmid) & (den$x <= xlim[2])
                if(sum(x.right))
                    max.right <- max(max.right, den$y[x.right])
            }
        }
        if(max.right < max.left)
            pos[1] <- xlim[1] + (xlim[2] - xlim[1]) / 2.1 # slightly to left of center
        pos
    }
    # draw.legend starts here

    nclasses <- length(yhat.per.class)
    if(is.null(legend.pos))
        legend.pos <- get.legend.pos()
    if(length(legend.pos) != 2)
        stop1("length(legend.pos) != 2")
    if(is.null(legend.names))
        legend.names <- names(yhat.per.class)
    if(length(legend.names) < nclasses) {
        warning1("length ", length(legend.names), " of legend.names ",
             "is less than the number ", nclasses, " of classes")
        legend.names <- rep(legend.names, nclasses)
    }
    else for(iclass in 1:nclasses)
        if(degenerate[iclass])
            legend.names[iclass] <- paste(legend.names[iclass], "(not plotted)")
    lwd <- rep(1, nclasses)
    # if the first histogram is filled in, then make its legend lwd bigger
    if(fill==col[1] && fill != "white" && fill != 0)
        lwd[1] <- 4

    if(legend.extra)
        legend.names <- paste(legend.names, " (", sapply(yhat.per.class, length),
                              " cases)", sep="")

    legend(x=legend.pos[1], y=legend.pos[2], legend=legend.names,
           cex=legend.cex, bg=legend.bg, lty=lty, lwd=lwd, col=col)
}
# shade the "error areas" of the density plots

add.err.col <- function(densities, thresh, col, border, lwd) 
{
    den1 <- densities[[1]]
    den2 <- densities[[2]]
    # is reducible error area to the left or to the right?
    # set iden=1 if to the left, iden=2 if to the right
    iden <- den1$y[den1$x >= thresh][1] > den2$y[den2$x >= thresh][1]
    if(is.na(iden)) { # no overlap between classes?
        warning1("no overlap between (first two) classes, ignoring \"err.col\" argument")
        return(NULL)
    }
    iden <- if(iden) 2 else 1
    if(length(col) < 2)
        col[2] <- col[1]
    if(length(col) < 3)
        col[3] <- col[iden]
    if(length(border) < 2)
        border[2] <- border[1]
    if(length(border) < 3)
        border[3] <- border[iden]
    if(length(lwd) < 2)
        lwd[2] <- lwd[1]
    if(length(lwd) < 3)
        lwd[3] <- lwd[iden]
    if(col[1] != 0 || border[1] != 0) {
        # left side of threshold
        matches <- den2$x < thresh
        if(sum(matches)) {
            x <- c(den2$x[matches])
            y <- c(den2$y[matches])
            len <- length(x)
            x[len] <- thresh # close possible tiny gap
            x[len+1] <- thresh
            y[len+1] <- 0
            polygon(x, y, col=col[1], border=border[1], lwd=lwd[1])
        }
    }
    if(col[2] != 0 || border[2] != 0) {
        # right side of threshold
        matches <- den1$x > thresh
        if(sum(matches)) {
            x <- den1$x[matches]
            y <- den1$y[matches]
            x[1] <- thresh # close possible tiny gap
            len <- length(x)
            x[len+1] <- thresh
            y[len+1] <- 0
            polygon(x, y, col=col[2], border=border[2], lwd=lwd[2])
        }
    }
    if(col[3] != 0 || border[3] != 0) {
        if(iden == 1) {
            # reducible error, left side of threshold
            # get indices i1 of den1 and i2 of den2 where den1 crosses den2
            i2 = length(den2$x)
            for (i1 in length(den1$x):1) {
                while(i2 > 1 && den2$x[i2] > den1$x[i1])
                    i2 <- i2 - 1
                if(den1$x[i1] <= thresh && den1$y[i1] >= den2$y[i2])
                    break
            }
            i1 <- den1$x <= thresh & (1:length(den1$x)) >= i1
            i2 <- den2$x <= thresh & (1:length(den2$x)) >= i2
            if(sum(i1) && sum(i2)) {
                # reverse i1 so polygon ends where it starts
                i1 <- rev(((1:length(den1$x))[i1]))
                # close tiny x gap to left of threshhold line
                den1$x[i1][1] = thresh 
                den2$x[i2][sum(i2)] = thresh
                polygon(x = c(den1$x[i1], den2$x[i2]),
                        y = c(den1$y[i1], den2$y[i2]),
                        col=col[3], border=border[3], lwd=lwd[3])
           }
        } else {
            # reducible error, right side of threshold
            # get indices i1 of den1 and i2 of den2 where den1 crosses den2
            i2 = 1
            for (i1 in 1:length(den1$x)) {
                while(i2 < length(den2$x) && den2$x[i2] < den1$x[i1])
                    i2 <- i2 + 1
                if(den1$x[i1] >= thresh && den2$y[i2] >= den1$y[i1])
                    break
            }
            i1 <- den1$x >= thresh & (1:length(den1$x)) < i1
            i2 <- den2$x >= thresh & (1:length(den2$x)) < i2
            if(sum(i1) && sum(i2)) {
                # reverse i2 so polygon ends where it starts
                i2 <- rev(((1:length(den2$x))[i2])) 
                polygon(x = c(den1$x[i1], den2$x[i2]),
                        y = c(den1$y[i1], den2$y[i2]),
                        col=col[3], border=border[3], lwd=lwd[3])
            }
        }
    }
}
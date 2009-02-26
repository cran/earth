# plotd.R: plot densities of class conditional predicted values
#
# TODO: allow newdata so can plot not only with the training data

plotd <- function(object,   # object is a model object
    hist = FALSE,           # FALSE to use density(), TRUE to use hist()
    type="response",        # passed on to predict
    thresh = .5,            # passed to predict (only used if type="class")
    xlim = NULL,            # NULL means auto
    main = NULL,            # graph caption
                            #   "string"  string
                            #   ""        no caption
                            #   NULL      generate a caption from x$call
    xlab="Predicted Value",
    ylab=if(hist) "Count" else "Density",
    lty = 1,                # linetypes for the plotted lines
    col = c("grey60", 1, "lightblue", "brown", "pink", 2, 3, 4), # cols for plotted lines
    borders = NULL,         # border colors for hist plot, NULL means use col
    labels = FALSE,         # labels parameter for hist plot, only TRUE or FALSE
    zero.line = FALSE,      # parameter for density plot, only used if hist=FALSE
    legend = TRUE,          # TRUE to draw a legend
    legend.names = NULL,    # NULL means auto, else specify a vector of strings,
    legend.pos = NULL,      # NULL means auto, else specify c(x,y) in user coords
    legend.cex = .8,        # cex for legend
    legend.bg = "white",    # bg for legend
    legend.extra = FALSE,   # print number in each class in legend
    vline.col = 0,          # color of vertical line, use NULL or 0 for no line
    vline.thresh = thresh,  # position of vertical line
    vline.lty = 1,          # line type of vertical line
    ...)                    # passed on to density or hist
{
    if(typeof(object) != "list")
        stop1("'", deparse(substitute(object)), "' is not a model object")
    if(is.null(coef(object)))
        warning1("'", deparse(substitute(object)),
            "' doesn't look like a model object, because the coef component is NULL")
    if(is.numeric(labels))
        labels <- as.logical(labels)
    if(!is.logical(labels) || length(labels) != 1)
        stop1("labels parameter must be a logical")

    # issue an error message if a dots argument is not in allowed.args
    # allowed.args is the list of arguments accepted by density() or hist()
    # but excluding the ones we have subsumed for our own use and those
    # that hist uses only for plotting

    if(hist) # we call hist() below
        allowed.args <- c("breaks", "freq", "probability", "include.lowest",
                          "right", "density", "labels", "nclass")
    else # we call density() below
        allowed.args <- c("bw", "adjust", "kernel", "weights", "window",
                          "width", "give.Rkern", "n", "from", "to", "cut")

    stop.if.dot.arg.not.allowed("plotd", allowed.args, ...)

    predictions.per.class <- get.predictions.per.class(object, type, thresh)
    nclasses <- length(predictions.per.class)

    # as.numeric is needed for logical or factor responses
    # because density and hist want a numeric arg
    density1 <- function(x, ...)
        if(hist)
            hist(as.numeric(x), plot=FALSE, ...)
        else
            density(as.numeric(x), ...)

    densities <- lapply(predictions.per.class, density1, ...)

    # get x limits of density plot

    if(is.null(xlim)) {
        xlim <- range(predictions.per.class)
        xspan <- xlim[2] - xlim[1]
        # extend the limits slightly to take into
        # account the "edges" of the density
        xlim[1] <- xlim[1] - xspan / 10
        xlim[2] <- xlim[2] + xspan / 10
    }
    if(length(xlim) != 2)
        stop("length(xlim) != 2")
    xspan <- xlim[2] - xlim[1]
    if(xspan < 0.01)
        warning1("x axis has a small range ", xspan,
                 " (not enough variation in model response?)")

    # get y limits of density plot

    ymax <- -1
    for(iclass in 1:nclasses) {
        if(hist) # need y component so hist can be treated like density below
            densities[[iclass]]$y <- densities[[iclass]]$counts
        ymax <- max(ymax, densities[[iclass]]$y)
    }
    # expand lty and other arguments if necessary

    if(length(lty) < nclasses)
        lty <- rep(lty, nclasses)
    if(length(col) < nclasses)
        col <- rep(col, nclasses)
    if(is.null(borders))
        borders <- col
    if(length(borders) < nclasses)
        borders <- rep(borders, nclasses)
    main <- show.caption(get.caption.from.call(main, object), 
                         if(is.null(main)) trim=.5 else trim=0, show=FALSE)

    # plot the first graph

    if(hist) # plot.hist
        plot(densities[[1]], xlim=xlim, ylim=c(0, ymax), col=col[1],
             main=main, xlab=xlab, ylab=ylab, lty=lty[1], 
             border=borders[1], labels=labels)
    else # plot.density
        plot(densities[[1]], xlim=xlim, ylim=c(0, ymax), col=col[1],
             main=main, xlab=xlab, ylab=ylab, lty=lty[1], zero.line=zero.line)

    # overlay the remaining graphs

    stopifnot(nclasses >= 2)
    for(iclass in 2:nclasses)
        if(hist) # lines.histogram
            lines(densities[[iclass]], col=NULL, lty=lty[iclass], 
                  border=borders[iclass], labels=labels)
        else # lines.density
            lines(densities[[iclass]], col=col[iclass], lty=lty[iclass])

    # optional vertical line at vline.thresh

    if(!is.null(vline.col))
        abline(v=vline.thresh, col=vline.col, lty=vline.lty)

    # optional legend

    if(legend) {
        if(is.null(legend.pos)) {
            # take a stab at positioning the legend correctly --
            # on left or right, away from the highest peak
            legend.pos <- c(0,0)
            if(max(densities[[1]]$y) < max(densities[[nclasses]]$y))
                legend.pos[1] <- xlim[1] # place on left side of graph
            else # slightly to the left of center
                legend.pos[1] <- xlim[1] + xspan / 2.1
            legend.pos[2] <- ymax
        }
        if(length(legend.pos) != 2)
            stop1("length(legend.pos) != 2")
        if(is.null(legend.names))
            legend.names <- names(predictions.per.class)
        if(length(legend.names) < nclasses) {
            warning1("length ", length(legend.names), " of legend.names ",
                 "is less than the number ", nclasses, " of classes")
            legend.names <- rep(legend.names, nclasses)
        }
        if(legend.extra)
            legend.names <- paste(legend.names, " (",
                                  sapply(predictions.per.class, length),
                                  " cases)", sep="")

        legend(x=legend.pos[1], y=legend.pos[2], legend=legend.names, bg=legend.bg,
               lty=lty, col=col, cex=legend.cex)
    }
    invisible(predictions.per.class)
}

# get the original observed response (it's needed to determine correct classes)

get.observed.response <- function(object)
{
    if(is.null(object$call$formula))
        y <- get.update.arg(NULL, "y", object, trace=FALSE, reeval=FALSE)
    else {
        # get y from formula and data used in original call to earth
        data <- get.update.arg(NULL, "data", object, FALSE)
        Call <- object$call
        m <- match(c("formula", "data"), names(Call), 0)
        mf <- Call[c(1, m)]
        mf[[1]] <- as.name("model.frame")
        mf$na.action <- na.fail
        mf <- eval.parent(mf) # TODO eval.parent correct?
        y <- model.response(mf, "any")  # "any" means factors are allowed
    }
    y
}

# return the predictions for each class, in a list

MIN.NBR.PER.CLASS <- 3  # 3 is arbitrary

get.predictions.per.class <- function(object, type, thresh)
{
    check.min <- function(x, ...)
    {
        len <- length(x)
        if(len < MIN.NBR.PER.CLASS) {
            if(len == 0) # prevent problem with density later
                stop1("no occurences of ", paste(..., sep=""), 
                      " in the observed response y")
            warning1("only ", len, " occurences of ", paste(..., sep=""),
                    " in the observed response y")
        }
    }
    get.binary.class.names <- function(yhat, fitted.values)
    {
        stopifnot(NCOL(yhat) == 1)
        stopifnot(NCOL(fitted.values) == 1)
        if(!is.null(colnames(yhat)))
            c(paste("not", colnames(yhat)), colnames(yhat))
        else if(!is.null(colnames(fitted.values)))
            c(paste("not", colnames(fitted.values)), colnames(fitted.values))
        else if(is.logical(yhat))
            c("FALSE", "TRUE")
        else
            c("zero", "non-zero")
    }
    cannot.use <- function(...) 
    {
        stop1("plotd does not support this kind of model\n       Reason: ", ...)
    }
    # get.predictions.per.class starts here
    # nomeclature: y is the observed response, yhat is the predicted response
    y <- get.observed.response(object)
    if(!is.character(type))
        stop1("type of \"type\" is not character")
    if(!is.na(pmatch(type, "terms")))
        stop1("type=\"terms\" is not allowed by plotd")
    yhat <- predict(object, type=type, thresh=thresh) # predicted response
    yhat.per.class <- list() # will put per-class predicted vals in here
    if(NCOL(yhat) == 1) {
        #---single column yhat--------------------------------------------------
        if(is.factor(y)) {
            # single column yhat, factor y
            levels <- levels(y)
            if(length(levels) < 2)
                stop1("length(levels) < 2")
            if(is.factor(yhat)) {
                # factor yhat, split into nlevels classes
                for(iclass in 1:length(levels)) {
                    level <- levels[iclass]
                    yhat.per.class[[iclass]] <- yhat[y == level]
                    len <- length(yhat.per.class[[iclass]])
                    check.min(yhat.per.class[[iclass]], level)
                    }
                names(yhat.per.class) <- levels
            } else {
                # numeric or logical yhat for a factor y
                # split into two classes i.e. the first level versus the rest
                # get here for logistic glm models (but not for earth-glm models)
                if(!is.numeric(yhat) && !is.numeric(yhat))
                    cannot.use("ncol(yhat)==1, is.factor(y), yhat unrecognized")
                level1 <- levels[1]
                yhat.per.class[[1]] <- yhat[y == level1]
                check.min(yhat.per.class[[1]], level1)
                yhat.per.class[[2]] <- yhat[y != level1]
                check.min(yhat.per.class[[2]], "\"not ", level1, "\"")
                names(yhat.per.class) <- c(level1, paste("not", level1))
            }
        } else if(NCOL(y) == 2) {
            # single column yhat, two column y, must be a binomial pair model
            # split into two classes based on relative sizes of y
            if(!is.numeric(y[,1]) || !is.numeric(y[,2]))
                cannot.use("ncol(yhat)==1, ncol(y)==1, ",
                           "!is.numeric(y[,1]) || !is.numeric(y[,2]")
            yhat.per.class[[1]] <- yhat[y[,1] <= y[,2]]
            check.min(yhat.per.class[[1]], "y[,1] <= y[,2]")
            yhat.per.class[[2]] <- yhat[y[,1] > y[,2]]
            check.min(yhat.per.class[[2]], "y[,1] > y[,2]")
            names(yhat.per.class) <- get.binary.class.names(yhat, object$fitted.values)
        } else if(NCOL(y) == 1) {
            # single column yhat, single column y, numeric or logical
            # split into two classes based on zero and non-zero y
            if(!is.numeric(y) && !is.logical(y))
                cannot.use("ncol(yhat)==1, ncol(y)==1, !is.numeric(y) && !is.logical(y)")
            yhat.per.class[[1]] <- yhat[y <= 0]
            check.min(yhat.per.class[[1]], "y <= 0")
            yhat.per.class[[2]] <- yhat[y > 0]
            check.min(yhat.per.class[[2]], "y > 0")
            names(yhat.per.class) <- get.binary.class.names(yhat, object$fitted.values)
        } else
            stop1("plotd does not support this kind of model")
  } else {
        #---multiple column yhat------------------------------------------------
        # split into ncol(yhat) classes, assuming y is a factor or 
        # multiple column numeric response

        if(!is.numeric(yhat))
            cannot.use("ncol(yhat)>1, yhat not numeric")
        if(is.factor(y) && NCOL(y) == 1) {
            levels <- levels(y)
            for(iclass in 1:ncol(yhat)) {
                level <- levels[iclass]
                yhat.per.class[[iclass]] <- yhat[y == level, iclass]
                check.min(yhat.per.class[[iclass]], level)
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop1("no occurences of ", level, 
                          " in the observed response")
            }
            names(yhat.per.class) <- levels(y)
        } else if(is.numeric(y) && NCOL(y) == NCOL(yhat)) {
            for(iclass in 1:ncol(yhat)) {
                yhat.per.class[[iclass]] <- yhat[y[,iclass] > thresh, iclass]
                check.min(yhat.per.class[[iclass]], "y > thresh")
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop1("no occurences of y <= thresh ", thresh,
                          " in the observed response y")
            }
            # jump through a few hoops to set the names of yhat.per.class
            ynames <- paste("y", 1:ncol(yhat), sep="")
            names(yhat.per.class) <- 
                if(length(colnames(yhat)) == ncol(yhat))
                     colnames(yhat)
                else if(length(colnames(y)) == ncol(y))
                    colnames(y)
                else 
                    ynames
             which. <- which(colnames(y) == "")
             if(length(which.))
                 names(yhat.per.class)[which.] <- ynames[which.]
        } else
            cannot.use("invalid y for ncol(yhat)>1")
    }
    yhat.per.class
}

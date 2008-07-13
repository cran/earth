# mars.to.earth.R: convert an mda:mars object to an earth object
#
# Stephen Milborrow Mar 2007 Forden, Wales
#
# The differences between mda:mars and earth objects are:
#
#   1. mars returns the MARS basis matrix in $x;
#      earth returns it in $bx.
#      There is no $x component of earth.
#
#   2. after the forward pass, earth discards lin dep terms in
#      bx, dirs and cuts
#
#   3. mars returns $all.terms; earth doesn't
#      Unneeded because of 2 above.
#
#   4. mars returns $lenb; earth doesn't.
#      Unneeded because of 2: lenb == nrow(cuts) == nrow(dirs)
#
#   4. mars$factor == earth$dirs (i.e. factor renamed to dirs).
#      In general, model$factor (sometimes called factors) is not
#      treated uniformly in the R code, so there seems to be no
#      compelling need to names dirs factor.
#      Note that this is not the same as model$terms$factors,
#      which is treated uniformly (but means something different).
#      Also earth$dirs can have a value of 2, for lin dep terms.
#
#   5. the formal arguments to mars and earth differ, thus $call differs
#
#   6. earth objects can be created through the formula interface and
#      if so will have a $terms field (doesn't apply to the conversion below)
#
#   7. earth objects have some extra components
#
#   8. mars normalizes the wp arg to len 1; earth normalizes the wp len
#      equal to the number of cols in y (so an all 1s wp argument is
#      equivalent to no wp argument).

mars.to.earth <- function(object=stop("no 'object' arg"))
{
    check.classname(object, deparse(substitute(object)), "mars")
    oldcall <- object$call
    newcall <- object$call
    newcall[[1]] <- as.name("earth")
    if(!is.null(object$call$prune) && !eval.parent(object$call$prune))
        newcall$pmethod <- "none"   # prune=FALSE was specified in the original call
    newcall$prune <- NULL
    if(!is.null(object$call$trace.mars) && eval.parent(object$call$trace.mars))
        newcall$trace <- 4          # trace.mars=TRUE was specified in original call
    newcall$trace.mars <- NULL
    y <- eval.parent(object$call$y)
    # convert vector y to ncases x 1 matrices so can access uniformly below.
    if(is.null(dim(y)))
        dim(y) <- c(length(y), 1)
    nresp  <- ncol(y)               # number of responses
    ncases <- nrow(y)               # number of cases
    weights.used <- FALSE
    if(!is.null(object$call$wp)) {
        newcall$wp <- object$call$wp
        object$call$wp <- NULL      # prevent partial match to "w" below
        weights.used <- TRUE
    }
    if(!is.null(object$call[["w"]]) && !is.null(eval.parent(object$call[["w"]]))) {

        warning("the \"w\" argument was used in the original call to mda::mars\n",
                "although mda::mars actually ignores the \"w\" argument")
        newcall$weights <- object$call$w
        weights.used <- TRUE
    }
    newcall$w <- NULL
    newcall$forward.step <- NULL
    newcall$prev.fit <- NULL
    if(!is.null(dim(residuals)))
        dim(residuals) <- c(ncol(y), nrow(y)) # convert vector to ncases x 1 matrix

    nselected <- length(object$selected.terms)
    residuals <- object$residuals
    penalty <- object$penalty

    # Renumber selected.terms.  Needed because earth drops terms from cuts and
    # dirs that are not in all.terms (whereas mars does not).

    selected <- rep(NA, nrow(object$factor))
    selected[object$all.terms] <- FALSE
    selected[object$selected.terms] <- TRUE
    selected <- selected[!is.na(selected)]
    selected.terms <- (1:length(selected))[selected]

    # Fill in the [1] and [nselected] elements of rss.per.subset and gcv.per.subset.
    # This is enough for print.earth() and summary.earth() etc. to work.
    # You can fill in all the elements by calling update.earth() later.
    # We don't call update.earth() now because minor differences between pruning
    # pass implementations could conceivably change selected.terms.

    ntermsVec <- rep(NA, length(object$all.terms))
    ntermsVec[1] <- 1                                # intercept
    ntermsVec[nselected] <- nselected                # nterms of selected model

    rss.per.subset <- rep(NA, length(object$all.terms))
    rss.per.subset[1] <- sum(colSums((y - colMeans(y)) ^ 2)) # null RSS
    rss.per.subset[nselected] <- sum(residuals^2)            # RSS of selected model
    rss <- rss.per.subset[nselected]                         # RSS of selected model

    gcv.per.subset <- get.gcv(rss.per.subset, ntermsVec, penalty, ncases)
    gcv <- gcv.per.subset[nselected]                         # GCV of selected model

    rss.per.response  <- vector(mode="numeric", length=nresp)
    rsq.per.response  <- vector(mode="numeric", length=nresp)
    gcv.per.response  <- vector(mode="numeric", length=nresp)
    grsq.per.response <- vector(mode="numeric", length=nresp)
    for(iresp in 1:nresp) {
        rss.per.response[iresp]  <- sum(residuals[,iresp]^2)
        rss.null                 <- sum((y[,iresp] - mean(y[,iresp]))^2)
        rsq.per.response[iresp]  <- get.rsq(rss.per.response[iresp], rss.null)
        gcv.null                 <- get.gcv(rss.null, 1, penalty, ncases)
        gcv.per.response[iresp]  <- get.gcv(rss.per.response[iresp],
                                            nselected, penalty, ncases)
        grsq.per.response[iresp] <- get.rsq(gcv.per.response[iresp], gcv.null)
    }
    pred.names <- generate.colnames(object$factor)
    term.names <- get.earth.term.name(1:nrow(object$factor),
                                      object$factor, object$cuts, pred.names, NULL)
    dimnames(object$factor) <- list(term.names, pred.names)
    dimnames(object$cuts)   <- list(term.names, pred.names)
    colnames(object$x) <- term.names[selected.terms]
    rownames(object$coefficients)  <- term.names[selected.terms]
    response.names <- generate.colnames(object$coefficients, is.y.arg=TRUE, xname=NULL)
    colnames(object$fitted.values) <- response.names
    colnames(object$residuals)     <- response.names
    colnames(object$coefficients)  <- response.names
    dirs <- object$factor[object$all.terms, , drop=FALSE]

    rval <- structure(list(
        bx                = object$x,
        dirs              = dirs,
        cuts              = object$cuts[object$all.terms, , drop=FALSE],
        selected.terms    = selected.terms,
        prune.terms       = NULL, # init later if you want by calling update.earth()
        rss               = rss,
        rsq               = get.rsq(rss, rss.per.subset[1]),
        gcv               = gcv,
        grsq              = get.rsq(gcv, gcv.per.subset[1]),
        rss.per.response  = rss.per.response,
        rsq.per.response  = rsq.per.response,
        gcv.per.response  = gcv.per.response,
        grsq.per.response = grsq.per.response,
        rss.per.subset    = rss.per.subset,
        gcv.per.subset    = gcv.per.subset,
        fitted.values     = object$fitted.values,
        residuals         = residuals,
        coefficients      = object$coefficients,
        penalty           = object$penalty,
        namesx            = colnames(dirs),
        call              = newcall),
    class = "earth")

    if(weights.used) { # wp or w args used in original call?
        # mars and earth normalize wp differently, see header comments
        # TODO there is probably a better way of handling this

        warning1("w or wp were used in the original call to mars.\n",
                 "         Running update.earth to conform mars ",
                 "use of weights to earth.\n")

        rval <- update(object=rval)
    }
    else if(!isTRUE(all.equal(object$gcv, rval$gcv)))
        warning1("the original mars GCV is              ", object$gcv,
                 "\n         ",
                 "but the GCV recalculated for earth is ", rval$gcv, "\n")

    my.print.call("Converted ", oldcall)
    my.print.call("to        ", newcall)

    rval
}

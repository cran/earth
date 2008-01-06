# mars.to.earth.R: convert an mda:mars object to an earth object
#
# Stephen Milborrow Mar 2007 Forden, Wales
#
# The differences between mda:mars and earth objects are:
#
#   1. mars returns bx in $x; earth returns bx in $bx with only
#      the selected.terms.
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
#   7. earth objects have a few extra components
#
# $$ Note that the w parameter is actually ignored by the mda:mars routines - a bug?
#
# Earth is much faster than mars on large models.
#
#
# $$ need to add support for multiple responses (easy but haven't had time to do it)

mars.to.earth <- function(object=stop("no 'object' arg"))
{
    check.classname(object, deparse(substitute(object)), "mars")
    call <- object$call
    call[[1]] <- as.name("earth")
    if(!is.null(object$call$prune) && !eval(object$call$prune, sys.parent()))
        call$pmethod <- "none"  # prune=FALSE was specified in the original call
    call$prune <- NULL
    if(!is.null(object$call$trace.mars) && eval(object$call$trace.mars, sys.parent()))
        call$trace <- 4         # trace.mars=TRUE was specified in the original call
    call$trace.mars <- NULL
    if(!is.null(object$call$wp) && !is.null(eval(object$call$wp, sys.parent())))
        warning("the mars 'wp' argument is not supported by earth()")
    call$wp <- NULL
    if(!is.null(object$call$w))
        call$weights <- object$call$w
    call$w <- NULL
    call$forward.step <- NULL
    call$prev.fit <- NULL

    y <- eval.parent(object$call$y)
    if(NCOL(y) != 1)
        stop1("'y' has more than one column (multiple responses are not yet supported)")
    nselected <- length(object$selected.terms)

    # Fill in the [1] and [nselected] elements of rss.per.subset and gcv.per.subset.
    # This is enough for print.earth() and summary.earth() etc. to work.
    # You can fill in all the elements by calling update.earth() later.
    # We don't call update.earth() now because minor differences between pruning
    # pass implementations could conceivably change selected.terms.

    ntermsVec <- rep(NA, length(object$all.terms))
    ntermsVec[1] = 1                                # intercept
    ntermsVec[nselected] = nselected                # nterms of selected model

    rss.per.subset <- rep(NA, length(object$all.terms))
    rss.per.subset[1] <- sum((y - mean(y))^2)               # null RSS
    rss.per.subset[nselected] <- sum(object$residuals^2)    # RSS of selected model

    gcv.per.subset <- get.gcv(rss.per.subset, ntermsVec, object$penalty, length(y))

    rss <- rss.per.subset[nselected]                # RSS of selected model
    gcv <- gcv.per.subset[nselected]                # GCV of selected model

    if(!isTRUE(all.equal(object$gcv, gcv)))         # should never happen
        warning("the original mars GCV ", object$gcv,
                " is not equal to the re-calculated GCV ", gcv)

    # Renumber selected.terms, needed because we drop terms from cuts and
    # dirs that are not in all.terms (whereas mars does not)

    selected <- rep(NA, nrow(object$factor))
    selected[object$all.terms] <- FALSE
    selected[object$selected.terms] <- TRUE
    selected <- selected[!is.na(selected)]
    selected.terms <- (1:length(selected))[selected]

    # Add names (actually, the only names we need for plotting etc. are names for dirs)

    pred.names <- colnames(object$factor)
    term.names <- get.earth.term.name(1:nrow(object$factor),
                                      object$factor, object$cuts, pred.names)
    dimnames(object$factor) <- list(term.names, pred.names)
    dimnames(object$cuts)   <- list(term.names, pred.names)

    # Show the new call if trace.mars was enabled in the original mars object.

    if(!is.null(call$trace) && eval.parent(call$trace))
        cat("Converted to", strip.white.space(format(call)), "\n")

    structure(list(
        bx                = object$x,
        dirs              = object$factor[object$all.terms, , drop=FALSE],
        cuts              = object$cuts[object$all.terms, , drop=FALSE],
        selected.terms    = selected.terms,
        prune.terms       = NULL, # init later if you want by calling update.earth()
        rss               = rss,
        rsq               = get.rsq(rss, rss.per.subset[1]),
        gcv               = gcv,
        grsq              = get.rsq(gcv, gcv.per.subset[1]),
        rss.per.response  = rss,
        rsq.per.response  = get.rsq(rss, rss.per.subset[1]),
        gcv.per.response  = gcv,
        grsq.per.response = get.rsq(gcv, gcv.per.subset[1]),
        rss.per.subset    = rss.per.subset,
        gcv.per.subset    = gcv.per.subset,
        fitted.values     = object$fitted.values,
        residuals         = object$residuals,
        coefficients      = object$coefficients,
        penalty           = object$penalty,
        call              = call),
    class = "earth")
}

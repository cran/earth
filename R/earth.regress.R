# earth.regress.R: earth.regress is used only for testing Regress in earth.c

earth.regress <- function(
    x       = stop("no \"x\" arg"), # NAs are not allowed in x or y
    y       = stop("no \"y\" arg"),
    weights = NULL,               # case weights
    used.cols = NULL)
{
    # following copied from header of earth.fit
    # expand factors, convert to double matrix with col names
    env <- parent.frame()
    x <- expand.arg(x, env)
    y <- expand.arg(y, env, is.y.arg=TRUE, deparse(substitute(y)))
    if(nrow(x) == 0)
        stop("no \"x\" values")
    if(ncol(x) == 0)    # this happens for example for earth(Volume~Volume,data=trees)
        stop("no \"x\"")
    if(nrow(x) != nrow(y))
        stop("nrow(x) ", nrow(x), " != nrow(y) ", nrow(y))
    if(!all(is.double(x)))
        stop("non double entries in \"x\" argument")
    if(!all(is.double(y)))
        stop("non double entries in \"y\" argument")
    if(!is.null(weights) && length(weights) != nrow(x))
        stop0("length(weights) ", length(weights), " != nrow(x) ", nrow(y))

    # add intercept to x
    colnames. <- colnames(x)
    x <- cbind(rep(1,nrow(x)), x)
    colnames(x) <- c("(Intercept)", colnames.)

    nresp <- ncol(y)
    ncols <- ncol(x)
    ncases <- nrow(x)

    if(is.null(weights))
        weights <- rep(1, ncases)
    else {
        meanw <- check.weights(weights, nrow(x), "weights")
        weights <- standardize.weights(weights, meanw)
    }
    if(is.null(used.cols)) {
        used.cols <- rep(TRUE, ncols)
        coefficients <- matrix(1.0, nrow=ncol(x), ncol=nresp)
    } else {
        if(!is.logical(used.cols))
            stop("used.cols is not logical")
        if(length(used.cols) != ncol(x)-1)     # -1 for intercept added above
            stop("length(used.cols) != ncol(x)")
        check.index.vec("used.cols", used.cols, x, check.empty=TRUE, use.as.col.index=FALSE)
        used.cols <- c(TRUE, used.cols)         # add intercept
        coefficients <- matrix(1.0, nrow=ncol(x) - sum(!used.cols), ncol=nresp)
    }
    rownames(coefficients) <- colnames(x)[used.cols]
    colnames(coefficients) <- colnames(y)

    rval <- .C("RegressR",
        coefficients = coefficients, # double  Betas[]     out: nUsedCols * nResp
        residuals = matrix(1.0, nrow=ncases, ncol=nresp),
                                # double       Residuals[] out: nCases * nResp
        rss = double(1),        # double       *pRss       out: RSS, summed over all nResp
        diags = double(ncols),  # double       Diags[]     out:
        rank = integer(1),      # int          *pnRank     out: nbr of indep cols in x
        pivots = integer(ncols),# int          iPivots[]   out: nCols, can be NULL
        as.double(x),           # const double x[]         in: nCases x nCols
        as.double(y),           # const double y[]         in: nCases x nResp
        as.double(weights),     # const double weights[]   in: nCases x 1, normalized weights
        as.integer(ncases),     # const int    *pnCases    in: number of rows in x and in y
        as.integer(nresp),      # const int    *pnResp     in: number of cols in y
        as.integer(ncols),      # int          *pnCols     in: number of columns in x
        as.integer(used.cols),  # const bool   UsedCols[]) in: specifies used columns in x
        PACKAGE="earth")

    rval$fitted.values <- y - rval$residuals
    rval$call <- match.call(expand.dots=TRUE)
    rval
}

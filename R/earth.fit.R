# earth.fit.R:
#
# Functions are in alphabetical order.


# returns the number of arguments to the user's "allowed" function

check.allowed.arg <- function(allowed) # check earth's "allowed" argument
{
    len <- 0
    if(!is.null(allowed)) {
        allowed.func.needs <- paste0(
            "  The 'allowed' function needs the following arguments ",
            "(but namesx and first are optional):\n      ",
            paste.collapse(c("degree", "pred", "parents", "namesx", "first")))

        if(!identical(typeof(allowed), "closure"))
            stop0("your 'allowed' argument is not a function")
        names. <- names(formals(allowed))
        len <- length(names.)
        if(len < 3 || len > 5)
            stop0("your 'allowed' function does not have the correct number of arguments\n",
                  allowed.func.needs)

        if(names.[1] != "degree" || names.[2] != "pred" || names.[3] != "parents" ||
           (len >= 4 && names.[4] != "namesx") || (len >= 5 && names.[5] != "first")) {
              stop0(allowed.func.needs,
                "\n  You have:\n      ", paste.collapse(names.))
        }
    }
    len
}
check.weights <- function(w, wname, expected.len, tweak.zero.weights) # invoked for both wp and weights
{
    check.vec(w, wname, expected.len)
    check(w, wname, "negative value", function(x) { x < 0 })
    if(tweak.zero.weights)
        w <- tweak.zero.weights(w, wname)
    w
}
convert.linpreds.to.logical <- function(linpreds, npreds, x)
{
    linpreds <- check.index(linpreds, "linpreds", x,
                            is.col.index=TRUE, allow.empty=TRUE)
    to.logical(linpreds, npreds)
}
# This is called from earth.default or earth.formula, not directly
# because the x and y args must be expanded for factors first.

earth.fit <- function(
    x       = stop("no 'x' argument"), # x and y already processed by model.matrix
    y       = stop("no 'y' argument"), # NAs are not allowed in x or y, an error msg if so
    weights = NULL,         # case weights (row weights)
    wp      = NULL,         # response weights (column weights)
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail
    offset  = NULL,         # offset term in formula
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
    keepxy  = FALSE,
    trace = 0,              # 0 none 1 overview 2 forward 3 pruning
                            # 4 model mats, memory use, more pruning, etc. 5  ...
    glm     = NULL,         # glm parameter from earth.formula or earth.default
    degree  = 1,            # max degree of interaction (1=additive model) (Friedman's mi)

    penalty = if(degree > 1) 3 else 2,
                            # GCV penalty per knot:
                            #   0 penalizes only terms (not knots)
                            #   special case -1 means no penalty (so GRSq==RSq)

                            # Following affect forward pass only, not pruning pass

    nk             = min(200, max(20, 2 * ncol(x))) + 1,
                            # max number of model terms including intercept

    thresh         = 0.001, # used as one of the conditions to stop adding terms in forw pass
                            # stop if RSqDelta<thresh or 1-RSq<thresh

    minspan        = 0,     # consider knots that are minspan apart
                            # special value 0 means use internally calculated min span
                            # negative values specify the number of knots

    endspan        = 0,     # special value 0 means use internally calculated end span

    newvar.penalty = 0,     # penalty for adding a new variable in forward pass

    fast.k         = 20,    # Fast MARS K: 0 means use all terms i.e. no Fast MARS
    fast.beta      = 1,     # Fast MARS ageing coefficient

                            # Following affect pruning only, not forward pass
                            # If you change these, update prune.only.args too!

    linpreds       = FALSE, # index vector specifying which preds must enter linearly

    allowed        = NULL,  # constraint function specifying allowed interactions

    nprune          = NULL, # max nbr of terms (including intercept) in pruned subset

                            # Following will not usually be used by the user

    Object  = NULL,         # if null, recreate earth model from scratch with forward pass
                            # if not null: no forward pass, just pruning pass

    Scale.y            = NULL,  # TRUE to scale y in the forward pass
    Adjust.endspan     = 2,     # for adjusting endspan for interaction terms
    Auto.linpreds      = TRUE,  # assume predictor linear if knot is min predictor value
    Force.weights      = FALSE, # TRUE to force use of weight code in earth.c
    Use.beta.cache     = TRUE,  # TRUE to use beta cache, for speed
    Force.xtx.prune    = FALSE, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    Get.leverages      = NROW(x) < 1e5, # TRUE to get leverages
    Exhaustive.tol     = 1e-10,
    ...)                        # unused
{
    init.global.data()
    check.no.family.arg.to.earth(..., is.null.glm.arg=is.null(glm))
    stop.if.dots(...)
    if(is.logical(trace) && trace) {
        warning0("earth: converted trace=TRUE to trace=4")
        trace <- 4
    }
    on.exit(init.global.data()) # release memory on exit
    Force.xtx.prune <- check.boolean(Force.xtx.prune)
    Use.beta.cache  <- check.boolean(Use.beta.cache)
    Get.leverages   <- check.boolean(Get.leverages)
    env <- parent.frame()
    y.org <- y
    if(!is.null(offset))
        y <- y - offset # apply offset before weights, like stats::lm.wfit
    pmethod <- match.arg1(pmethod, "pmethod")
    print_earth_fit_args(trace, x, y)
    # we do basic parameter checking here but more in ForwardPass in earth.c
    check.integer.scalar(nk, min=1, max=1000)
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop0("illegal 'na.action', only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop0("illegal 'na.action', only na.action=na.fail is allowed")
    na.action <- na.fail
    n.allowed.args <- check.allowed.arg(allowed)
    check.integer.scalar(degree, min=0, max=10)
    check.numeric.scalar(penalty)
    check.numeric.scalar(thresh)
    check.integer.scalar(minspan)
    check.integer.scalar(endspan, min=0)
    check.numeric.scalar(newvar.penalty)
    check.numeric.scalar(fast.k)
    check.numeric.scalar(fast.beta)
    check.integer.scalar(nprune, null.ok=TRUE)
    check.numeric.scalar(Adjust.endspan)
    check.boolean(Auto.linpreds)
    check.boolean(Force.weights)
    check.boolean(Use.beta.cache)
    check.boolean(Force.xtx.prune)
    check.boolean(Get.leverages)
    check.numeric.scalar(Exhaustive.tol)
    stopifnot(!is.null(dim(x))) # should have been converted to mat in earth.default or earth.formula
    stopifnot(!is.null(dim(y))) # ditto
    if(nrow(x) == 0)
        stop0("'x' has no rows")
    if(ncol(x) == 0)    # this happens for example for earth(Volume~Volume,data=trees)
        stop0("'x' has no columns")
    if(nrow(x) != nrow(y))
        stop0("nrow(x) ", nrow(x), " != nrow(y) ", nrow(y))
    # x and y must be double for calls to C functions
    if(!all(is.double(x)))
        stop0("non double entries in 'x' argument")
    if(!all(is.double(y)))
        stop0("non double entries in 'y' argument")
    check.no.na.in.mat(x)
    check.no.na.in.mat(y)
    # Oct 2021: Check for duplicate column names (possibly due to factor expansion)
    if(length(unique(colnames(x))) != length(colnames(x)))
        stop0("Duplicate colname in x (colnames are ",
              paste.with.quotes(colnames(x), maxlen=100), ")")

    glm.arg <- process.glm.arg(glm)
    # determine is.bpairs before applying weights (which can change y values)
    is.bpairs <- !is.null(glm.arg) && is.bpairs(y, glm.arg$family, trace, is.earth=TRUE)
    ret <- get.weights(weights, y, is.bpairs,
                       !is.null(subset), !is.null(glm.arg), Force.weights, trace)
        weights               <- ret$weights
        weights.before.bpairs <- ret$weights.before.bpairs
        use.weights           <- ret$use.weights

    # subset
    if(!is.null(subset)) {
        # duplicates are allowed in subset so user can specify a bootstrap sample
        subset <- check.index(subset, "subset", x, allow.dups=TRUE, allow.zeros=TRUE)
        x <- x[subset, , drop=FALSE]
        y <- y[subset, , drop=FALSE]
        weights <- weights[subset]
        if(!is.null(weights.before.bpairs) && !is.null(glm.arg))
            weights.before.bpairs <- weights.before.bpairs[subset]
        trace1(trace, "%d cases after taking subset\n", nrow(x))
        trace2(trace, "\n")
    }
    # case weights (applied to y)
    yw <- NULL # weighted version of y (also an indicator that weights should be used)
    if(use.weights)
        yw <- sqrt(weights) * y

    # column weights
    if(!is.null(wp)) {
        wp <- check.weights(wp, "wp", ncol(y), tweak.zero.weights=TRUE)
        wp <- sqrt(wp / mean(wp)) # normalize
        # multiply each column of y by its normalized weight
        y <- y  * outer(repl(1, nrow(y)), wp)
        if(!is.null(yw))
            yw <- yw * outer(repl(1, nrow(y)), wp)
    }
    # glm bpairs
    yfrac <- NULL # used only if response is glm bpairs
    if(is.bpairs) {
        trace1(trace,
           "Response columns %s and %s are a binomial pair (%g obs in total)\n",
           colnames(y)[1], colnames(y)[2], sum(y))
        yfrac <- bpairs.yfrac(y, trace) # fraction true
        y <- yfrac
        if(use.weights)
            yw <- sqrt(weights) * y
    }
    Scale.y <- get.Scale.y(Scale.y, y, use.weights, wp, Force.weights)
    maxmem <- maxmem(x, nk, trace)
    possible.gc(maxmem, trace, "before forward.pass")
    if(is.null(Object)) {
        rv <- forward.pass(x, y, yw, weights,
                    trace, degree, penalty, nk, thresh,
                    minspan, endspan, newvar.penalty, fast.k, fast.beta,
                    linpreds, allowed,
                    Scale.y, Adjust.endspan, Auto.linpreds, Use.beta.cache,
                    n.allowed.args, env, maxmem)
        termcond <- rv$termcond
        bx       <- rv$bx
        dirs     <- rv$dirs
        cuts     <- rv$cuts
    } else {
        # no forward pass: get here if update.earth() called me with no forward pass params
        trace1(trace, "Skipped forward pass\n")
        check.classname(Object, substitute(Object), "earth")
        dirs     <- Object$dirs
        cuts     <- Object$cuts
        termcond <- Object$termcond
        bx       <- get.bx(x, seq_len(nrow(dirs)), dirs, cuts) * sqrt(weights) # weight bx
    }
    possible.gc(maxmem, trace, "after forward.pass")
    penalty1 <- penalty

    # add column names to bx (necessary for possible err msg in pruning.pass)
    pred.names <- colnames(x)
    possible.gc(maxmem, trace, "after pruning.pass") # apply(range) in get.earth.term.name uses much mem
    term.names <- get.earth.term.name(seq_len(nrow(dirs)), dirs, cuts, pred.names, x)
    possible.gc(maxmem, trace, "after get.earth.term.name") # release memory used by get.earth.term.name
    colnames(bx) <- term.names
    rv <- pruning.pass(if(use.weights) sqrt(weights) * x else x,
                       if(use.weights) yw else y, bx,
                       pmethod, penalty, nprune,
                       trace, dirs, Force.xtx.prune, Exhaustive.tol)
    rss.per.subset <- rv$rss.per.subset # vector of RSSs for each model (index on subset size)
    gcv.per.subset <- rv$gcv.per.subset # vector of GCVs for each model
    prune.terms    <- rv$prune.terms    # triang mat: each row is a vector of term indices
    selected.terms <- rv$selected.terms # vec of term indices of selected model

    remove(rv) # free memory
    bx <- bx[, selected.terms, drop=FALSE]
    bx <- bx / sqrt(weights) # unweight bx
    colnames(bx) <- term.names[selected.terms]
    dimnames(dirs) <- list(term.names, pred.names)
    dimnames(cuts) <- list(term.names, pred.names)

    nresp <- ncol(y)  # number of responses
    nselected <- length(selected.terms)

    # regress y on bx to get coefficients, fitted.values etc.
    if(use.weights)
        lm.fit <- lm.wfit(bx, y, w=weights, singular.ok=FALSE)
    else
        lm.fit <- lm.fit(bx, y, singular.ok=FALSE)

    # the as.matrix calls are needed if y is a vector
    # so the fitted.values etc. are always matrices (not vectors)
    fitted.values <- as.matrix(lm.fit$fitted.values)
    residuals     <- as.matrix(lm.fit$residuals)
    coefficients  <- as.matrix(lm.fit$coefficients)
    resp.names <- colnames(y)
    colnames(fitted.values) <- resp.names
    colnames(residuals)     <- resp.names
    colnames(coefficients)  <- resp.names
    if(!is.null(offset)) {
        stopifnot(NROW(offset) == NROW(fitted.values))
        fitted.values <- fitted.values + offset
    }
    leverages <- NULL
    if(Get.leverages) {
        # when n >> p, memory use peaks here
        qr <- lm.fit$qr
        remove(lm.fit) # free memory
        leverages <- hatvalues.qr.wrapper(qr, maxmem, trace)
        possible.gc(maxmem, trace, "after hatvalues.qr.wrapper")
    }
    if(!is.null(wp)) {
        tt <- outer(repl(1, nrow(y)), wp)
        fitted.values <- fitted.values / tt # divide each column by its wp
        residuals <- residuals / tt
        y <- y / tt
        coefficients <- coefficients / outer(repl(1, nselected), wp)
    }
    # build glm model(s) if glm argument is not NULL
    glm.list <- NULL    # glm.list is a list of glm models, NULL if none
    glm.coefs <- NULL   # glm.coefs is a nselected x nresponses matrix
    if(!is.null(glm.arg)) {
        y.glm <- y.org
        if(!is.null(subset))
            y.glm <- y.glm[subset, , drop=FALSE]
        glm.list <- earth.glm(bx, y.glm, weights.before.bpairs, na.action, offset,
                              glm.arg, trace, is.bpairs, env)
        glm.coefs <- get.glm.coefs(glm.list=glm.list,
                                   nresp=if(is.bpairs) 1 else ncol(coefficients),
                                   selected.terms, term.names, resp.names)
    }
    # following is for consistency when running test suite on different machines
    rss.per.subset[rss.per.subset < 1e-10] <- 0
    gcv.per.subset[gcv.per.subset < 1e-10] <- 0

    # prepare returned summary statistics

    rss.per.response  <- vector(mode="numeric", length=nresp)
    rsq.per.response  <- vector(mode="numeric", length=nresp)
    gcv.per.response  <- vector(mode="numeric", length=nresp)
    grsq.per.response <- vector(mode="numeric", length=nresp)
    tss.per.response  <- vector(mode="numeric", length=nresp)
    offset.fitted.values <- if(is.null(offset)) fitted.values else fitted.values - offset
    if(nresp == 1 && !use.weights) { # special case to save memory
        rss.per.response[1] <- sos(residuals)
        rsq.per.response[1] <- get.weighted.rsq(y, offset.fitted.values, weights)
        gcv.per.response[1] <- get.gcv(rss.per.response[1], nselected, penalty, nrow(bx))
        tss.per.response[1] <- sos(y - mean(y))
    } else for(iresp in seq_len(nresp)) { # multiple response or weighted model
        y1 <- y[,iresp]
        offset.fitted.values1 <- offset.fitted.values[,iresp]
        rss.per.response[iresp] <- sos(y1 - offset.fitted.values1, weights)
        rsq.per.response[iresp] <- get.weighted.rsq(y1, offset.fitted.values1, weights)
        gcv.per.response[iresp] <- get.gcv(rss.per.response[iresp], nselected, penalty, nrow(bx))
        tss.per.response[iresp] <- sos(y1 - weighted.mean(y1, weights), weights)
    }
    for(iresp in seq_len(nresp)) {
        gcv.null <- get.gcv(tss.per.response[iresp], 1, penalty, nrow(bx))
        grsq.per.response[iresp] <-
            if(!use.weights)
                get.rsq(gcv.per.response[iresp], gcv.null)
            else if(nresp == 1)
                get.rsq(gcv.per.subset[nselected], gcv.per.subset[1])
            else # TODO grsq is wrong when weights and multiple responses
                get.rsq(gcv.per.response[iresp], gcv.null)
    }
    rss  <- rss.per.subset[nselected]
    rsq  <- get.rsq(rss, rss.per.subset[1])
    gcv  <- gcv.per.subset[nselected]
    grsq <- get.rsq(gcv, gcv.per.subset[1])
    rv <- structure(list(   # term 1 is the intercept in all returned data
        rss            = rss,   # RSS, across all responses if y has multiple cols
        rsq            = rsq,   # R-Squared, across all responses
        gcv            = gcv,   # GCV, across all responses
        grsq           = grsq,  # GRSq across all responses

        bx             = bx,    # selected terms only
        dirs           = dirs,  # all terms including unselected: nterms x npreds
        cuts           = cuts,  # all terms including unselected: nterms x npreds
        selected.terms = selected.terms,# row indices into dirs and cuts
        prune.terms    = prune.terms,   # nprune x nprune, each row is vec of term indices

        fitted.values  = fitted.values, # ncases (after subset) x nresp
        residuals      = residuals,     # ncases (after subset) x nresp
        coefficients   = coefficients,  # selected terms only: nselected x nresp

        rss.per.response  = rss.per.response,   # nresp x 1, RSS for each response
        rsq.per.response  = rsq.per.response,   # nresp x 1, RSq for each response
        gcv.per.response  = gcv.per.response,   # nresp x 1, GCV for each response
        grsq.per.response = grsq.per.response,  # nresp x 1, GRSq for each response

        rss.per.subset = rss.per.subset,# nprune x 1, RSS of each model, across all resp
        gcv.per.subset = gcv.per.subset,# nprune x 1, GCV of each model, across all resp

        leverages      = leverages,
        pmethod        = pmethod,
        nprune         = nprune,
        penalty        = penalty,          # copy of penalty argument
        nk             = nk,               # copy of nk argument
        thresh         = thresh,           # copy of thresh argument
        termcond       = termcond,         # reason we terminated the forward pass
        weights        = if(use.weights) weights.before.bpairs else NULL,
        Scale.y        = Scale.y),         # return Scale.y so can pass to earth.cv

    class = "earth")

    if(!is.null(offset))
        rv$offset <- offset

    if(!is.null(glm.list)) {
        rv$glm.list         <- glm.list   # list of glm models, NULL if none
        rv$glm.coefficients <- glm.coefs  # matrix of glm coefs, nselected x nresp
        rv$glm.stats        <- get.glm.stats(glm.list, colnames(rv$fitted.values))
        rv$glm.bpairs       <- if(is.bpairs) c(TRUE, FALSE) else NULL # backwards compat
        rv$glm.yfrac        <- yfrac # fraction true
    }
    rv
}
effective.nbr.of.params <- function(ntermsVec, nknotsVec, penalty)  # for GCV calculation
{
    if(penalty < 0) # special case: term and knots are free so GCV == RSS/ncases
        repl(0, length(ntermsVec))
    else
        ntermsVec + (penalty * nknotsVec)
}
forward.pass <- function(x, y, yw, weights, # must be double, but yw can be NULL
                         trace, degree, penalty, nk, thresh,
                         minspan, endspan, newvar.penalty, fast.k, fast.beta,
                         linpreds, allowed,
                         Scale.y, Adjust.endspan, Auto.linpreds, Use.beta.cache,
                         n.allowed.args, env, maxmem)
{
    if(nrow(x) < 2)
        stop0("the x matrix must have at least two rows")
    stopifnot(nrow(x) == nrow(y))
    if(!is.null(yw)) {
        stopifnot(nrow(y) == nrow(yw))
        stopifnot(ncol(y) == ncol(yw))
    }
    npreds <- ncol(x)
    fullset <- repl(0, nk) # element will be set TRUE if corresponding term used
    linpreds <- convert.linpreds.to.logical(linpreds, npreds, x)
    if(trace >= 2)
        print_linpreds(linpreds, x)
    if(Scale.y) {
        y <- get.scaled.y(y, "y")
        if(!is.null(yw))
            yw <- scale(yw, center=FALSE, scale=attr(y, "scaled:scale"))
    }
    print_scaled_y(trace, Scale.y, y) # prints only if trace >= 5

    # we are careful to initialize the "out" variables for ForwardPassR here
    # in a way that does not require them to be duplicated in ForwardPassR
    fullset  <- as.integer(fullset)
    bx       <- matrix(0, nrow=nrow(x), ncol=nk)
    dirs     <- matrix(0, nrow=nk, ncol=npreds)
    cuts     <- matrix(0, nrow=nk, ncol=npreds)
    termcond <- integer(length=1) # reason we terminated the forward pass

    stopifnot(!is.null(colnames(x))) # ensure we have predictor names
    stopifnot(is.double(x)) # no typecast in .Call below
    stopifnot(is.double(y))
    stopifnot(is.null(weights) || is.double(weights))

    on.exit(.C("FreeEarth", PACKAGE="earth")) # if error or user interrupt, free mem

    .Call("ForwardPassR",
        fullset,                           # out: int FullSet[]
        bx,                                # out: double bx[]
        dirs,                              # out: double Dirs[]
        cuts,                              # out: double Cuts[]
        termcond,                          # out: int*
        x,                                 # in: double x[]
        y,                                 # in: double y[]
        yw,                                # in: double yw[] or NULL
        weights,                           # in: double WeightsArg[] (never NULL)
        as.integer(nrow(x)),               # in: int* nCases
        as.integer(ncol(y)),               # in: int* nResp
        as.integer(npreds),                # in: int* nPreds
        as.integer(degree),                # in: int* nMaxDegree
        as.double(penalty),                # in: double* Penalty
        as.integer(nk),                    # in: int* nMaxTerms
        as.double(thresh),                 # in: double* Thresh
        as.integer(minspan),               # in: int* nMinSpan
        as.integer(endspan),               # in: int* nEndSpan
        as.integer(fast.k),                # in: int* nFastK
        as.double(fast.beta),              # in: double* FastBeta
        as.double(newvar.penalty),         # in: double* NewVarPenalty
        as.integer(linpreds),              # in: int LinPreds[]
        allowed,                           # in: SEXP Allowed
        as.integer(n.allowed.args),        # in: int* nAllowedArgs
        env,                               # in: SEXP Env
        as.double(Adjust.endspan),         # in: double AdjustEndSpan
        as.integer(Auto.linpreds),         # in: int* nAutoLinPred
        as.integer(Use.beta.cache),        # in: int* nUseBetaCache
        as.double(max(trace, 0)),          # in: double* Trace
        colnames(x),                       # in: char* sPredNames[]
        NAOK = TRUE, # we check for NAs etc. internally in C ForwardPass
        PACKAGE="earth")

    fullset  <- as.logical(fullset)
    list(termcond = termcond,
         bx       = bx[, fullset, drop=FALSE],
         dirs     = dirs[fullset, , drop=FALSE],
         cuts     = cuts[fullset, , drop=FALSE])
}
# Used when building the model to name the columns of bx and rows of dirs etc.
# Also called by mars.to.earth.
# Return string like "h(55-x1)*h(x2-58)".
# h represents the hockey stick func.
# If ntermsVec is a vector, this returns a vector of strings.
# x can be NULL (currently only when called from mars.to.earth), it is used
# only for simplifying terms with factor predictors.

get.earth.term.name <- function(ntermsVec, dirs, cuts, pred.names, x, warn.if.dup=TRUE)
{
    get.term.name1 <- function(nterm, dirs, cuts, pred.names, xrange, form1, form2)
    {
        get.name <- function(ipred) # return "name" if possible, else "x[,i]"
        {
            pred.name <- pred.names[ipred]
            if(is.null(pred.name) || anyNA(pred.name))
                paste0("x[,", ipred, "]")
            else
                pred.name
        }
        if(nterm == 1)
            return("(Intercept)")
        s <- ""
        first.fac <- TRUE
        stopifnot(ncol(dirs) > 0)
        for(ipred in seq_len(ncol(dirs)))
            if(dirs[nterm,ipred]) {
                if(!first.fac)
                    s <- paste0(s, "*")
                first.fac <- FALSE
                if(dirs[nterm,ipred] == 2)  # linear predictor?
                    s <- pastef(s, "%s", get.name(ipred))
                else if(dirs[nterm,ipred] == -1)
                    s <- pastef(s, form1, cuts[nterm,ipred], get.name(ipred))
                else if(dirs[nterm,ipred] == 1) {
                    if(cuts[nterm,ipred] == 0 && !is.null(xrange) &&
                            xrange[1, ipred] == 0 && xrange[2, ipred] < 100 &&
                            all(x[,ipred] == floor(x[,ipred]))) # all integer?
                        # simplify to no hinge function, it's a factor
                        s <- pastef(s, "%s", get.name(ipred))
                    else
                        s <- pastef(s, form2, get.name(ipred), cuts[nterm,ipred])
                } else if(dirs[nterm,ipred] != 0)
                    stop0("illegal direction ", dirs[nterm,ipred], " in dirs")
            }
        s
    }
    #--- get.earth.term.name starts here ---
    stopifnot(ncol(dirs) == ncol(x))
    xrange <- NULL      # 1st row is min, 2nd row is max, a column for each pred
    if(!is.null(x))
        xrange <- apply(x, 2, range) # for simplifying "h(ldose-0)" to "ldose"

    # get format strings for sprint later
    ndigits <- getOption("digits")
    if(ndigits <= 7) {      # for back compat with previous versions of earth
        form1 <- "h(%g-%s)" # let %g figure out the nbr of digits
        form2 <- "h(%s-%g)"
    } else {
        form1 <- sprint("h(%%.%dg-%%s)", ndigits) # e.g. "h(%.9g-%s)"
        form2 <- sprint("h(%%s-%%.%dg)", ndigits) # e.g. "h(%s-%.9g)"
    }
    term.names <- sapply(seq_along(ntermsVec), get.term.name1, dirs, cuts,
                         pred.names, xrange, form1, form2)
    # check for duplicated term names
    duplicated <- duplicated(term.names)
    if(warn.if.dup && any(duplicated))
        warning0("duplicate term name \"", term.names[which(duplicated)[1]], "\"\n",
                 "This is usually caused by cuts that are very close to each other\n",
                 "Remedy: use options(digits=NDIGITS), ",
                 "typically NDIGITS has to be at least 7 ",
                 "(currently NDIGITS=", ndigits, ")")
    term.names
}
# get.gcv returns GCVs as defined in Friedman's MARS paper, with an
# extension for penalty < 0

get.gcv <- function(
    rss.per.subset,
    ntermsVec,          # number of MARS regression terms including intercept
    penalty,            # penalty per knot, argument from earth.fit()
    ncases)             # number of cases
{
    stopifnot(length(rss.per.subset) == length(ntermsVec))
    nknotsVec <- get.nknots(ntermsVec)
    nparams <- effective.nbr.of.params(ntermsVec, nknotsVec, penalty)
    ifelse(nparams >= ncases,
           Inf, # ensure that GCVs are non-decreasing as number of terms increases
           rss.per.subset / (ncases * (1 - nparams/ncases)^2))
}
# Return the estimated number of knots
#
# TODO This is not quite correct?  It assumes that each term pair adds one
# knot.  Thus each term adds "half a knot".  But if we have deleted a term
# in a pair then the remaining term should add a knot, not half a knot.

get.nknots <- function(nterms)
{
    (nterms - 1 ) / 2
}
get.Scale.y <- function(Scale.y, y, use.weights, wp, Force.weights)
{
    if(is.null(Scale.y)) # auto Scale.y
        Scale.y <- NCOL(y) == 1 && is.null(wp)
    else { # user specified Scale.y
        Scale.y <- check.boolean(Scale.y)
        if(Scale.y && !is.null(wp))
            stop0("Scale.y=TRUE is not allowed with wp (implementation restriction)")
    }
    Scale.y
}
get.scaled.y <- function(y, yname)
{
    y.scaled <- scale(y)
    # check that scaling was ok
    i <- which(attr(y.scaled, "scaled:scale") == 0)
    if(length(i)) {
        if(ncol(y) > 1)
            warning0("Cannot scale column ", i[1], " of ", yname,
                     " (values are all equal to ", y[1,i], ")\n",
                     "         Use Scale.y=FALSE to silence this warning")
        else
            warning0("Cannot scale ", yname,
                     " (values are all equal to ", y[1,1], ")\n",
                     "         Use Scale.y=FALSE to silence this warning")
        y.scaled <- y # fall back
    }
    y.scaled
}
get.weights <- function(weights, y, is.bpairs, is.subset.arg, is.glm.arg, Force.weights, trace)
{
    n <- nrow(y)
    weights.specified <- !is.null(weights)
    if(is.null(weights))
        weights <- repl(1, n)
    else
        weights <- check.weights(weights, "weights", n, tweak.zero.weights=FALSE)
    weights.before.bpairs <- if(weights.specified || Force.weights) weights else NULL
    if(is.bpairs) {
        # note: all zero rows will have weight zero, which will later
        # be changed to almost-zero in tweak.zero.weights below
        weights <- weights * rowSums(y)
    }
    # use weights only if necessary, because earth is much faster without weights
    use.weights <- check.boolean(Force.weights) ||
                   any(abs(weights - weights[1]) > 1e-8)
    if(!use.weights)
        weights <- repl(1, n)
    if(!is.double(weights)) # paranoia: weights must be double for calls to C functions
        weights <- as.double(weights)

    trace.weights(trace, weights, weights.before.bpairs, is.bpairs,
                  is.glm.arg, weights.specified, use.weights)

    list(weights               = tweak.zero.weights(weights, "weights"),
         weights.before.bpairs = weights.before.bpairs,
         use.weights           = use.weights)
}
# hatvalues() doesn't work on lm.fit objects, so get leverages ourselves
# TODO this hasn't been tested for multiple response models
hatvalues.qr <- function(qr, maxmem, trace)
{
    possible.gc(maxmem, trace, "hatvalues.qr1")
    y <- diag(1, nrow=nrow(qr$qr), ncol=qr$rank)
    possible.gc(maxmem, trace, "hatvalues.qr2")
    qr.qy <- qr.qy(qr, y)^2 # calculate (Q %*% y)^2
    remove(y) # free memory for rowSums
    possible.gc(maxmem, trace, "hatvalues.qr3")
    leverages <- rowSums(qr.qy)
    leverages[leverages < 1e-8]     <- 0 # allow for numerical error
    leverages[leverages > 1 - 1e-8] <- 1 # allow for numerical error
    leverages
}
# this wrapper is because when n >> p we run out of memory in qr.qy
hatvalues.qr.wrapper <- function(qr, maxmem, trace)
{
    # treat memory allocation warning in qr.qy as an error, not a warning
    old.warn <- getOption("warn")
    on.exit(options(warn=old.warn))
    options(warn=2) # treat warnings as errors
    if(trace == 1.5)
        printf("Getting leverages\n")
    leverages <- try(hatvalues.qr(qr, maxmem, trace), silent=TRUE)
    if(is.try.err(leverages)) {
        options(warn=1)
        warning0("Not enough memory to get leverages (but otherwise the model is fine)")
        leverages <- NULL
    }
    leverages
}
init.global.data <- function()
{
    assignInMyNamespace("lamba.global",                        -999)
    assignInMyNamespace("lamba.factor.global",                 -999)
    assignInMyNamespace("prev.coef.global",                    NULL)
    assignInMyNamespace("trace.ncoef.global",                  0)
    assignInMyNamespace("issued.singularities.warning.global", FALSE)
}
# the following calculation of max mem needed by earth is just an approximation
maxmem <- function(x, nk, trace)
{
    sint <- 4       # sizeof int
    sdouble <- 8    # sizeof double
    n <- nrow(x)
    p <- ncol(x)
    x <- n * p * sdouble
    y <- n * sdouble
    bx <- n * nk * sdouble
    # matrices used in the C code
    xorder <- n * p * sint
    bxorth <- bx
    bxorthcenteredt <- bx
    xbx <- n * sdouble
    # TODO what about xUsed and work (used "After forward pass GRSq ...")
    maxmem <- x + y + bx + xorder + bxorth + bxorthcenteredt + xbx
    maxmem <- maxmem / 1024^3 # Bytes to GBytes
    if(trace >= 5 || trace == 1.5 || trace == 1.6)
        printf("maxmem %.1f GB\n", maxmem)
    maxmem # estimated maximum memory used by earth, in GBytes
}
possible.gc <- function(maxmem, trace, msg)
{
    if(maxmem > 1) { # more than 1 GByte of memory required by earth?
        if(trace == 1.5 || trace == 1.6)
            old.memsize <- memory.size()
        gc()
        if(trace == 1.5 || trace == 1.6) {
            printf("memsize %.1f to %.1f max %.1f GB %s\n",
                old.memsize / 1024, memory.size() / 1024,
                memory.size(max=TRUE) / 1024, msg)
        }
    }
}
print_earth_fit_args <- function(trace, x, y)
{
    if(trace >= 4)
        cat("\n")
    if(trace >= 1 && trace < 7) { # don't print matrices when doing very detailed earth.c tracing
        tracex <- if(trace >= 5) 4 else 2 # adjust trace for print_summary
        details <- if(trace >= 4) 2 else if(trace >= 1) -1 else 0
        print_summary(x, "x", tracex, details=details)
        if(details > 1) printf("\n")
        print_summary(y, "y", tracex, details=details)
        if(details > 1) printf("\n")
    }
}
print_linpreds <- function(linpreds, x)
{
    if(any(linpreds != 0)) {
        cat("Linear predictors ")
        colnames. <- colnames(x)
        index <- (1:length(linpreds))[linpreds]
        if(!is.null(colnames.))
            cat(paste(index, "=", colnames.[linpreds], sep="", collapse=" "))
        else
            cat(paste.collapse((1:length(linpreds))[linpreds]))
        cat("\n")
    }
}
print_scaled_y <- function(trace, Scale.y, y)
{
    if(trace >= 5 && trace < 7) { # don't print matrices when doing very detailed earth.c tracing
        if(!Scale.y)
            printf("\nScale.y = FALSE\n\n")
        else {
            printf("\nScale.y = TRUE: yscale %s ycenter %s\n\n",
                   paste.trunc(format(attr(y, "scaled:scale"), digits=5)),
                   paste.trunc(format(attr(y, "scaled:center"), digits=5)))
            tracex <- if(trace >= 5) 4 else 2 # adjust trace for print_summary
            details <- if(trace >= 4) 2 else if(trace >= 1) -1 else 0
            print_summary(y, "yscaled", tracex, details=details)
            if(details > 1) printf("\n")
        }
    }
}
trace.weights <- function(trace, weights, weights.before.bpairs, is.bpairs,
                          is.glm.arg, weights.specified, use.weights)
{
    if(trace < 1)
        return(NULL)
    if(is.bpairs) {
        if(use.weights) {
            printf("weights used by earth internally: %s\n",
                   paste.trunc(sprint("%.4g", weights), collapse=", ", maxlen=50))
            trace2(trace, "weights passed to glm (which will adjust by rowsums): %s\n",
                   if(is.null(weights.before.bpairs)) "NULL"
                   else paste.trunc(sprint("%.4g", weights.before.bpairs),
                                    collapse=", ", maxlen=30))
        } else
            printf("earth and glm: unweighted%s\n",
                   if(weights.specified) " (because all weights equal)" else "")
    } else if(is.glm.arg) {
        if(use.weights)
            printf("earth and glm weights[%d]: %s\n", NROW(weights),
                   paste.trunc(sprint("%.4g", weights), collapse=", ", maxlen=40))
        else if(weights.specified)
            printf("earth and glm: unweighted (because all weights equal)\n")
    } else {
        if(use.weights)
            printf("weights[%d]: %s\n", NROW(weights),
                   paste.trunc(sprint("%.4g", weights), collapse=", ", maxlen=60))
        else if(weights.specified)
            printf("weights: no weights (because all weights equal)\n")
    }
}
tweak.zero.weights <- function(w, wname)
{
    meanw <- mean(w)
    if(meanw == 0)
        stop0("mean of '", wname, "' is zero")
    else if(meanw < 1e-8)
        stop0("mean of '", wname, "' is (almost) zero")
    # TODO fix zero weights (but should maybe do what lm.wfit does and delete cols)
    almost.zero <- meanw / 1e8  # note that 1e8 becomes 1e4 after sqrt later
    w[w < almost.zero] <- almost.zero
    w
}

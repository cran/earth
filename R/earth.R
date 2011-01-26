# earth.R: an implementation of Friedman's Multivariate Adaptive
#          Regression Splines, commonly known as MARS.
#
# This code is derived from code in mda.R by Hastie and Tibshirani.
# Functions are in alphabetical order except for some method
# functions which are at the end
# Comments containing "TODO" mark known issues.
# Stephen Milborrow Mar 2007 Petaluma
#
#-----------------------------------------------------------------------------
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses
#
#-----------------------------------------------------------------------------
# Notes for earth() that didn't make it into the man pages.
#
# --- subset argument (for selecting cases)
#
# All subset handling is done in earth.fit not in earth.formula or
# update.earth.  This is because we want to allow the user to specify
# a subset even when he or she isn't using the formula based approach
# i.e.  using earth.default() and not earth.formula().
#
# --- Eval.model.subsets (and Print.pruning.pass)
#
# This function evaluates subsets of the MARS model.  It returns (in
# prune.terms) the best terms to use for each number of terms i.e. for
# each model size.  See the default function eval.model.subsets() for an
# example Eval.model.subsets.
#
#-----------------------------------------------------------------------------

# This is a list of those formal arguments of earth.fit that can be changed without
# requiring a new forward pass.
# NOTE: if you change the pruning formal arguments in earth.fit(), update this too!

prune.only.args <- c("glm", "trace", "nprune", "pmethod",
                     "Get.crit", "Eval.model.subsets", "Print.pruning.pass",
                     "Force.xtx.prune")

# returns the number of arguments to the user's "allowed" function

check.allowed.arg <- function(allowed) # check earth's "allowed" argument
{
    len <- 0
    if(!is.null(allowed)) {
        allowed.func.needs <- paste(
            "  The \"allowed\" function needs the following arguments ",
            "(but namesx and first are optional):\n      ",
            paste.with.space(c("degree", "pred", "parents", "namesx", "first")),
            sep="")

        if(!identical(typeof(allowed), "closure"))
            stop("your \"allowed\" argument is not a function");
        names. <- names(formals(allowed))
        len <- length(names.)
        if(len < 3 || len > 5)
            stop1("your \"allowed\" function does not have the correct number of arguments\n",
                  allowed.func.needs)

        if(names.[1] != "degree" || names.[2] != "pred" || names.[3] != "parents" ||
           (len >= 4 && names.[4] != "namesx") || (len >= 5 && names.[5] != "first")) {
              stop1(allowed.func.needs,
                "\n  You have:\n      ", paste.with.space(names.))
        }
    }
    len
}

check.and.standarize.wp <- function(wp, expected.len)
{
    meanw <- check.weights(wp, expected.len, "wp")
    sqrt(standardize.weights(wp, meanw))
}

# following is used for "weights" and "wp"

check.weights <- function(w, expected.len, wname)
{
    if(!is.vector(w))
        stop1("\"", wname, "\" is not a vector")
    if(is.logical(w))
        w <- as.numeric(w)
    if(!is.numeric(w))
        stop1("\"", wname, "\" is not numeric")
    if(length(w) != expected.len)
        stop1("\"", wname, "\" has a bad length ",
              length(w), ", expected ", expected.len)
    if(any(is.na(w)))
        stop1("NA in \"", wname, "\"")
    if(any(!is.finite(w)))
        stop1("non finite value in \"", wname, "\"")
    if(any(w < 0))
        stop1("negative value in \"", wname, "\"")
    meanw <- mean(w)
    if(meanw < 1e-8)
        stop1("mean of \"", wname, "\" is 0")
    meanw  # return mean
}

check.which.terms <- function(dirs, which.terms) # ensure which.terms is valid
{
    if(is.null(which.terms))
        stop1("\"which.terms\" is NULL")
    if(length(which.terms) == 0)
        stop1("length(which.terms) == 0")
    if(which.terms[1] != 1)
        stop1("first element of \"which.terms\" must be 1, the intercept term")
    check.index.vec("which.terms", which.terms[1], dirs)
}

convert.linpreds.to.logical <- function(linpreds, npreds, x)
{
    check.index.vec("linpreds", linpreds, x, use.as.col.index=TRUE)
    to.logical(linpreds, npreds)
}

earth <- function(...) UseMethod("earth")

earth.default <- function(
    x       = stop("no 'x' arg"), # NAs are not allowed in x or y, an error msg if so
    y       = stop("no 'y' arg"),
    weights = NULL,         # case weights
    wp      = NULL,         # response column weights
    scale.y = (NCOL(y)==1), # TRUE to scale y in the forward pass
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail
    glm     = NULL,
    trace   = 0,
    keepxy  = FALSE,        # true to retain x and y in returned value
    nfold   = 0,            # cross validation folds, 0 means no cross validation
    stratify = TRUE,        # stratify levels in cross validation folds
    ...)                    # passed on to earth.fit
{
    trace <- check.trace.arg(trace)
    Call <- make.call.generic(match.call(expand.dots=TRUE), "earth")
    if(trace >= 4)
        my.print.call("Call: ", Call)
    if(!is.null(Call$data))
        stop1("\"data\" argument not allowed in earth.default")
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop1("illegal \"na.action\", only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop1("illegal \"na.action\", only na.action=na.fail is allowed")
    if(keepxy) {
        x.org <- x
        y.org <- y
    }
    xname <- deparse(substitute(x))
    if(is.vector(x))
        namesx.org <-  xname    # ensure x always has a name
    else
        namesx.org <- colnames(x)
    namesx <- generate.colnames(x, is.y.arg=FALSE, xname=xname)

    # expand factors, convert to double matrix with col names

    env <- parent.frame()
    x <- expand.arg(x, env, FALSE, xname=xname)
    rownames(x) <- possibly.delete.rownames(x)
    ylevels <- levels(y)          # save levels if any before expanding y
    if(is.logical(y))
        ylevels <- c(FALSE, TRUE)  # needed for predict.earth(type="class")
    if(is.logical(y))
        ylevels <- c(FALSE, TRUE)
    y <- expand.arg(y, env, is.y.arg=TRUE, deparse(substitute(y)))
    rownames(y) <- possibly.delete.rownames(y)

    rval <- earth.fit(x=x, y=y, weights=weights, wp=wp, scale.y=scale.y,
                      subset=subset, na.action=na.action, glm=glm,
                      trace=trace, ...)
    rval$call <- Call
    rval$namesx.org <- namesx.org # name chosen not to alias with rval$x
    rval$namesx <- namesx         # ditto
    rval$levels <- ylevels
    rval$wp <- wp
    if(keepxy) {
        rval$x <- x.org
        # following ruse is needed else vector y's lose their name
        if(is.vector(y.org)) {
            dim(y.org) <- c(nrow=length(y.org), ncol=1)
            colnames(y.org) <- colnames(y)[1]
        }
        rval$y <- y.org
        rval$subset <- subset
        rval$weights <- weights
    }
    # earth.cv will return null unless nfold>1
    # subset parameter was already checked in earth.fit so it's safe to use it here
    cv <- earth.cv(if(is.null(subset)) x else x[subset,,drop=FALSE],
                   if(is.null(subset)) y else y[subset,,drop=FALSE],
                   weights, wp, scale.y, subset, na.action,
                   glm, rval$glm.bpairs, trace, keepxy, nfold, stratify, ...)
    if(!is.null(cv)) {
         rval$cv.list <- cv$list.
         rval$cv.nterms <- cv$nterms
         rval$cv.nvars <- cv$nvars
         rval$cv.rsq.tab <- cv$rsq.tab
         rval$cv.maxerr.tab <- cv$maxerr.tab
         rval$cv.deviance.tab <- cv$deviance.tab
         rval$cv.calib.int.tab <- cv$calib.int.tab
         rval$cv.calib.slope.tab <- cv$calib.slope.tab
         rval$cv.auc.tab <- cv$auc.tab
         rval$cv.cor.tab <- cv$cor.tab
         rval$cv.groups <- cv$groups
    }
    rval
}

earth.formula <- function(
    formula = stop("no 'formula' arg"), # intercept will be ignored
    data    = NULL,
    weights = NULL,         # case weights
    wp      = NULL,         # response column weights
    scale.y = (NCOL(y)==1), # TRUE to scale y in the forward pass
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail
    glm     = NULL,
    trace   = 0,
    keepxy  = FALSE,        # true to retain x and y in returned value
    nfold   = 0,            # cross validation folds, 0 means no cross validation
    stratify = TRUE,        # stratify levels in cross validation folds
    ...)                    # passed on to earth.fit
{
    # get x column names from model frame

    get.namesx <- function(mf)
    {
        t <- attr(mf,"terms")
        iresp <- attr(t,"response")
        namesx <- character(0)
        for(i in 1:ncol(mf))
            if(i != iresp)
                if(!is.null(colnames(mf[[i]])))
                    namesx <- c(namesx, colnames(mf[[i]]))
                else
                    namesx <- c(namesx, colnames(mf)[i])

        namesx
    }
    trace <- check.trace.arg(trace)
    Call <- make.call.generic(match.call(expand.dots=TRUE), "earth")
    if(trace >= 4)
        my.print.call("Call: ", Call)
    if(!is.null(Call$x))
        stop1("\"x\" argument not allowed in earth.formula")
    if(!is.null(Call$y))
        stop1("\"y\" argument not allowed in earth.formula")
    Call2 <- match.call(expand.dots=FALSE)

    # subset and weights handled in earth.fit, so match only on formula, data, and na.action

    m <- match(c("formula", "data", "na.action"), names(Call2), 0)
    mf <- Call2[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    if(!is.null(mf$na.action))
        stop1("\"na.action\" argument is not allowed (it is set internally to na.fail)")
    mf$na.action <- na.fail
    mf <- eval.parent(mf)

    # expand factors in x, convert to double matrix, add colnames

    x <- model.matrix(attr(mf, "terms"), mf)
    rownames(x) <- possibly.delete.rownames(x)
    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]     # silently discard intercept
    else
        warning1("ignored -1 in formula (earth objects always have an intercept)")
    # strip white space for better reading earth formula for e.g. earth(y~I(X-3))
    # because model.matrix inserts spaces around the minus
    if(!is.null(colnames(x)))
        colnames(x) <- strip.white.space(colnames(x))

    # expand factors in y, convert to double matrix, add colnames

    y <- model.response(mf, "any")  # "any" means factors are allowed
    ylevels <- levels(y)            # save levels if any before expanding y
    if(is.logical(y))
        ylevels <- c(FALSE, TRUE)    # needed for predict.earth(type="class")
    terms. <- attr(mf, "terms")
    yname <- NULL                   # yname=NULL means let expand.arg choose name
    if(is.vector(y)) {
        # following ruse is needed else vector y's lose their name
        names. <- names(attr(terms., "dataClasses"))
        yname <- names.[[attr(terms., "response")]]
    }
    env <- parent.frame()
    y <- expand.arg(y, env, is.y.arg=TRUE, xname=yname)
    rownames(y) <- possibly.delete.rownames(y)

    rval <- earth.fit(x=x, y=y, weights=weights, wp=wp, scale.y=scale.y,
                      subset=subset, na.action=na.action, glm=glm,
                      trace=trace, ...)

    rval$namesx.org <- get.namesx(mf)
    rval$namesx <- make.unique(rval$namesx.org)
    rval$levels <- ylevels
    rval$terms <- terms.
    rval$call <- Call
    rval$wp <- wp

    if(keepxy) {
        if(!is.null(data))
            rval$data <- data
        else # OLD: if(trace >= 1)
            warning1("No 'data' argument to earth so ignoring keepxy for 'data'\n")
        rval$subset <- subset
        rval$weights <- weights
    }
    # earth.cv will return null unless nfold>1
    # subset parameter was already checked in earth.fit so it's safe to use it here
    cv <- earth.cv(if(is.null(subset)) x else x[subset,,drop=FALSE],
                   if(is.null(subset)) y else y[subset,,drop=FALSE],
                   weights, wp, scale.y, subset, na.action,
                   glm, rval$glm.bpairs, trace, keepxy,  nfold, stratify, ...)
    if(!is.null(cv)) {
         rval$cv.rsq.tab <- cv$rsq.tab
         rval$cv.maxerr.tab <- cv$maxerr.tab
         rval$cv.deviance.tab <- cv$deviance.tab
         rval$cv.calib.int.tab <- cv$calib.int.tab
         rval$cv.calib.slope.tab <- cv$calib.slope.tab
         rval$cv.auc.tab <- cv$auc.tab
         rval$cv.cor.tab <- cv$cor.tab
         rval$cv.groups <- cv$groups
         rval$cv.nterms <- cv$nterms
         rval$cv.nvars <- cv$nvars
         rval$cv.list <- cv$list.
    }
    rval
}

# This is called from earth.default or earth.formula, not directly
# because the x and y args should be expanded for factors first.

earth.fit <- function(
    x       = stop("no 'x' arg"), # x and y already processed by model.matrix
    y       = stop("no 'y' arg"), # NAs are not allowed in x or y, an error msg if so
    weights = NULL,         # case weights
    wp      = NULL,         # response column weights
    scale.y = (NCOL(y)==1), # TRUE to scale y in the forward pass
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail
    glm = NULL,             # glm parameter from earth.formula or earth.default
    trace = 0,              # 0 none 1 overview 2 forward 3 pruning
                            # 4 model mats, memory use, more pruning, etc. 5  ...

    nk      = max(21, 2 * ncol(x) + 1),
                            # max number of model terms including intercept

    degree         = 1,     # max degree of interaction (1=additive model) (Friedman's mi)

    penalty = if(degree > 1) 3 else 2,
                            # GCV penalty per knot:
                            #   0 penalizes only terms (not knots)
                            #   special case -1 means no penalty (so GRSq==RSq)

                            # Following affect forward pass only, not pruning pass

    thresh         = 0.001, # used as one of the conditions to stop adding terms in forw pass
                            # stop if RSqDelta<thresh or 1-RSq<thresh

    minspan        = 0,     # consider knots that are minspan apart
                            # special value 0 means use internally calculated min span
                            # special value -1 (for back compatibility) means calc minspan
                            #   using the method used by earth versions prior to 2.4

    newvar.penalty = 0,     # penalty for adding a new variable in forward pass

    fast.k         = 20,    # Fast MARS K: 0 means use all terms i.e. no Fast MARS
    fast.beta      = 1,     # Fast MARS ageing coefficient

                            # Following affect pruning only, not forward pass
                            # If you change these, update prune.only.args too!

    linpreds       = FALSE, # index vector specifying which preds must enter linearly

    allowed        = NULL,  # constraint function specifying allowed interactions

    pmethod = "backward",   # for legal values see eval.model.subsets.*
    nprune  = NULL,         # max nbr of terms (including intercept) in prune subset

                            # Following will not usually be used by the user

    Object  = NULL,         # if null, recreate earth model from scratch with forward pass
                            # if not null: no forward pass, just pruning pass

    Get.crit           = get.gcv,            # criterion func for model select during pruning
    Eval.model.subsets = eval.model.subsets, # function used to evaluate model subsets
    Print.pruning.pass = print.pruning.pass, # function used to print pruning pass results
    Force.xtx.prune    = FALSE, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    Use.beta.cache     = TRUE,  # TRUE to use beta cache, for speed
    ...)                        # unused
{
    check.no.family.arg.to.earth(...)
    warn.if.dots.used("earth.fit", ...)
    if(is.logical(trace) && trace) {
        warning1("earth: converted trace=TRUE to trace=4")
        trace <- 4
    }
    y.org <- y

    # For binomial glms we have to drop paired y cols before passing y to
    # to the C earth routines.  The logical vector glm.bpairs keeps track
    # of which cols are used.
    glm.bpairs <- NULL                       # NULL means all columns used
    if(!is.null(glm)) {
        glm <- get.glm.arg(glm)
        glm.bpairs <- get.glm.bpairs(y, glm) # returns NULL if all cols used
    }
    if(trace >= 1) {
        print.matrix.info("x", x, NULL, NULL,       details=(trace >= 4), all.names=(trace >= 2))
        print.matrix.info("y", y, NULL, glm.bpairs, details=(trace >= 4), all.names=(trace >= 2))
    }
    # we do basic parameter checking here but much more in ForwardPass in earth.c
    if(nk < 1)
        stop1("\"nk\" ", nk, " is less than 1")
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop1("illegal \"na.action\", only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop1("illegal \"na.action\", only na.action=na.fail is allowed")
    na.action <- na.fail
    # TODO implementation of case weights is not complete so don't allow case weights
    if(!is.null(weights) &&
            !isTRUE(all.equal(weights, rep(weights[1], length(weights))))) {
        warning1("\"weights\" are not supported by \"earth\", ignoring them")
        weights <- NULL
    }
    n.allowedfunc.args <- check.allowed.arg(allowed)
    stopif(is.vector(x)) # should have been converted to matrix in earth.default or earth.formula
    stopif(is.vector(y)) # ditto
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
    if(!is.null(weights)) # must check weights now so can use if subset below
        meanw <- check.weights(weights, nrow(x), "weights")
    if(!is.null(subset)) {
        # duplicates are allowed so subset can specify a bootstrap sample
        check.index.vec("subset", subset, x,
                        check.empty=TRUE, allow.duplicates=TRUE, allow.zeroes=TRUE)
        x <- x[subset, , drop=FALSE]
        y <- y[subset, , drop=FALSE]
        if(!is.null(weights)) {
            weights <- weights[subset, , drop=FALSE]
            meanw <- mean(weights)
        }
    }
    standardized.weights <- NULL
    if(!is.null(weights))
        standardized.weights <- standardize.weights(weights, meanw)
    if(!is.null(wp)) {
        wp <- check.and.standarize.wp(wp, ncol(y))
        # multiply each column of y by its normalized weight
        y <- y * outer(rep(1, nrow(y)), wp)
    }
    if(!is.null(glm.bpairs))
        y <- y[, glm.bpairs, drop=FALSE]
    if(is.null(Object)) {
        env <- parent.frame()
        rval <- forward.pass(x, y, standardized.weights, scale.y, trace, penalty,
                             nk, degree, linpreds,
                             allowed, thresh, minspan, newvar.penalty,
                             fast.k, fast.beta, Use.beta.cache,
                             n.allowedfunc.args, env)
          bx   <- rval[[1]]
          dirs <- rval[[2]]
          cuts <- rval[[3]]
    } else {
        # no forward pass: get here if update() called me with no forward pass params
        if(trace >= 1)
            cat("Skipped forward pass\n")
        check.classname(Object, deparse(substitute(Object)), "earth")
        bx   <- NULL
        dirs <- Object$dirs
        cuts <- Object$cuts
    }
    rval <- pruning.pass(x, y, standardized.weights, subset,
                         trace, penalty, pmethod, nprune, bx, dirs, cuts, Get.crit,
                         Eval.model.subsets, Print.pruning.pass, Force.xtx.prune)
      bx             <- rval[[1]]
      rss.per.subset <- rval[[2]]   # vector of RSSs for each model (index on subset size)
      gcv.per.subset <- rval[[3]]   # vector of GCVs for each model
      prune.terms    <- rval[[4]]   # triangular mat: each row is a vector of term indices
      selected.terms <- rval[[5]]   # vector of term indices of selected model
    bx <- bx[, selected.terms, drop=FALSE]

    # add names for returned values

    pred.names <- colnames(x)
    term.names <- get.earth.term.name(1:nrow(dirs), dirs, cuts, pred.names, x)
    colnames(bx) <- term.names[selected.terms]
    dimnames(dirs) <- list(term.names, pred.names)
    dimnames(cuts) <- list(term.names, pred.names)

    # Regress y on bx to get fitted.values etc.
    # The as.matrix calls after the call to lm are needed if y is
    # a vector so the fitted.values etc. are always arrays.

    if(is.null(weights))
        lfit <- lm.fit(bx, y, singular.ok=FALSE)
    else
        lfit <- lm.wfit(bx, y, w=weights, singular.ok=FALSE)
    fitted.values <- as.matrix(lfit$fitted.values)
    residuals     <- as.matrix(lfit$residuals)
    coefficients  <- as.matrix(lfit$coefficients)
    response.names <- colnames(y)
    colnames(fitted.values) <- response.names
    colnames(residuals)     <- response.names
    colnames(coefficients)  <- response.names

    nresp <- ncol(y)  # number of responses
    nselected <- length(selected.terms)

    if(!is.null(wp)) {
        tt <- outer(rep(1, nrow(y)), wp)
        fitted.values <- fitted.values / tt
        residuals <- residuals / tt
        y <- y / tt
        coefficients <- coefficients / outer(rep(1, nselected), wp)
    }
    # build glm model(s) if glm argument is not NULL

    glm.list <- NULL    # glm.list is a list of glm models, NULL if none
    glm.coefs <- NULL   # glm.coefs is a nselected x nresponses matrix
    if(!is.null(glm)) {
        y.glm <- y.org
        if(!is.null(subset))
            y.glm <- y.glm[subset, , drop=FALSE]
        glm.list <- earth.glm(bx, y.glm, weights, na.action, glm,
                              trace, glm.bpairs, response.names[1])
        glm.coefs <- get.glm.coefs(glm.list, ncol(coefficients),
                                   selected.terms, term.names, response.names)
    }
    # prepare returned summary statistics

    rss  <- rss.per.subset[nselected]
    rsq  <- get.rsq(rss, rss.per.subset[1])
    gcv  <- gcv.per.subset[nselected]
    grsq <- get.rsq(gcv, gcv.per.subset[1])
    rss.per.response  <- vector(mode="numeric", length=nresp)
    rsq.per.response  <- vector(mode="numeric", length=nresp)
    gcv.per.response  <- vector(mode="numeric", length=nresp)
    grsq.per.response <- vector(mode="numeric", length=nresp)
    for(iresp in 1:nresp) {
        rss.per.response[iresp]  <- sum(residuals[,iresp]^2)
        rss.null                 <- sum((y[,iresp] - mean(y[,iresp]))^2)
        rsq.per.response[iresp]  <- get.rsq(rss.per.response[iresp], rss.null)
        gcv.null                 <- Get.crit(rss.null, 1, penalty, nrow(x))
        gcv.per.response[iresp]  <- Get.crit(rss.per.response[iresp], nselected,
                                             penalty, nrow(x))
        grsq.per.response[iresp] <- get.rsq(gcv.per.response[iresp], gcv.null)
    }
    rval <- structure(list(     # term 1 is the intercept in all returned data
        bx             = bx,    # selected terms only
        dirs           = dirs,  # all terms including unselected: nterms x npreds
        cuts           = cuts,  # all terms including unselected: nterms x npreds
        selected.terms = selected.terms,# row indices into dirs and cuts
        prune.terms    = prune.terms,   # nprune x nprune, each row is vec of term indices

        rss            = rss,   # RSS, across all responses if y has multiple cols
        rsq            = rsq,   # R-Squared, across all responses
        gcv            = gcv,   # GCV, across all responses
        grsq           = grsq,  # GRSq across all responses

        rss.per.response  = rss.per.response,   # nresp x 1, RSS for each response
        rsq.per.response  = rsq.per.response,   # nresp x 1, RSq for each response
        gcv.per.response  = gcv.per.response,   # nresp x 1, GCV for each response
        grsq.per.response = grsq.per.response,  # nresp x 1, GRSq for each response

        rss.per.subset = rss.per.subset,# nprune x 1, RSS of each model, across all resp
        gcv.per.subset = gcv.per.subset,# nprune x 1, GCV of each model, across all resp

        fitted.values  = fitted.values, # ncases (after subset) x nresp
        residuals      = residuals,     # ncases (after subset) x nresp
        coefficients   = coefficients,  # selected terms only: nselected x nresp

        penalty        = penalty),      # copy of penalty argument
    class = "earth")

    if(!is.null(glm.list)) {
        rval$glm.list         <- glm.list   # list of glm models, NULL if none
        rval$glm.coefficients <- glm.coefs  # matrix of glm coefs, nselected x nresp
        rval$glm.bpairs       <- glm.bpairs # will be null unless paired binomial cols
    }
    rval
}

effective.nbr.of.params <- function(ntermsVec, nknotsVec, penalty)  # for GCV calculation
{
    if(penalty < 0) # special case: term and knots are free so GCV == RSS/ncases
        0
    else
        ntermsVec + (penalty * nknotsVec)
}

# this returns the RSS and selected terms for each subset of size 1:nprune

eval.model.subsets <- function(     # this is the default Eval.model.subsets
    bx,         # basis model matrix
    y,          # model response
    weights,
    pmethod,
    nprune,     # max nbr of terms (including intercept) in prune subset, in range 1..nterms
    Force.xtx.prune) # TRUE to always call EvalSubsetsUsingXtx rather than leaps
{
    stopifnot(nprune >= 1 && nprune <= nrow(bx))
    if(Force.xtx.prune ||   # user explicitly asked for xtx subset evaluation
            ncol(y) > 1 ||  # leaps cannot deal with multiple responses
            ncol(bx) <= 2)  # leaps code gives an error for small number of cols
        eval.model.subsets.using.xtx(bx, y, weights, pmethod, nprune, Force.xtx.prune)
    else
        eval.model.subsets.with.leaps(bx, y, weights, pmethod, nprune)
}

eval.model.subsets.with.leaps <- function(
    bx,
    y,
    weights,
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep"),
    nprune)
{
    convert.lopt <- function(lopt)  # convert lopt format to prune.terms format
    {
        # Assignment fills matrix column wise. We want row wise, so
        # take upper triangle and then transpose.

        prune.terms <- matrix(0, nrow=nprune, ncol=nprune)
        prune.terms[upper.tri(prune.terms, diag=TRUE)] <- lopt
        t(prune.terms)
    }
    # eval.model.subsets.with.leaps starts here

    pacify <- pmatch(pmethod[1], "exhaustive", 0) == 1 && nrow(bx) > 40 && nprune > 8
    if(pacify) {                # pruning could take a while so print a reminder
        cat("Pruning...")
        flush.console()
    }
    if(is.null(weights))
        rprune <- leaps:::leaps.setup(x=bx, y=y,
                    force.in=1,        # make sure intercept is in model
                    force.out=NULL,
                    intercept=FALSE,   # we have an intercept so leaps.setup must not add one
                    nvmax=nprune, nbest=1, warn.dep=TRUE)
    else
        rprune <- leaps:::leaps.setup(x=bx, y=y, wt=weights,
                    force.in=1,
                    force.out=NULL,
                    intercept=FALSE,
                    nvmax=nprune, nbest=1, warn.dep=TRUE)

    rprune <- switch(match.arg1(pmethod),
                    leaps:::leaps.backward(rprune),     # "backward"
                    leaps:::leaps.backward(rprune),     # "none"
                    leaps:::leaps.exhaustive(rprune, really.big=TRUE),
                    leaps:::leaps.forward(rprune),
                    leaps:::leaps.seqrep(rprune))
    if(pacify)
        cat("\n")
    rss.per.subset <- as.vector(rprune$ress) # convert from n x 1 mat to vector
    prune.terms <- convert.lopt(rprune$lopt)

    list(rss.per.subset,    # vector of RSSs for each model (index on subset size)
         prune.terms)       # triangular mat: each row is a vector of term indices
}

# This calls the earth.c routine EvalSubsetsUsingXtxR.
# Unlike the leaps code, it can deal with multiple responses (i.e. multiple y columns)

eval.model.subsets.using.xtx <- function(
    bx,
    y,
    weights,
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep"),
    nprune,
    Force.xtx.prune)
{
    bad.pmethod <- function()
    {
        # following reasons must match conditions in eval.model.subsets
        reason <- "unknown"
        if(Force.xtx.prune)
            reason <- "Force.xtx.prune==TRUE"
        else if(ncol(y) > 1)
            reason <- "ncol(y) > 1"
        else if(ncol(bx) <= 2)
            reason <- "ncol(bx) <= 2"
        stop1("You have pmethod=\"", pmethod, "\" ",
            "but only pmethod=\"backward\" or \"none\"\n",
            "is supported by eval.model.subsets.using.xtx\n",
            "(eval.model.subsets.using.xtx was invoked because ", reason, ")\n")
        NULL
    }
    backward <- function(bx, y, weights.)
    {
        ncases <- nrow(bx)
        nterms <- ncol(bx)
        nresp <- ncol(y)
        stopifnot(is.double(bx))
        stopifnot(is.double(y))
        stopifnot(is.null(weights.) || is.double(weights.))
        rval <- .C("EvalSubsetsUsingXtxR",
            prune.terms = matrix(0, nrow=nterms, ncol=nterms),  # double PruneTerms[]
            rss.per.subset = vector(mode="numeric", length=nterms),
            as.integer(ncases),                       # const int *pnCases
            as.integer(nresp),                        # const int *pnResp
            as.integer(nterms),                       # const int *pnMaxTerms
            bx,                                       # const double bx[]
            y,                                        # const double y[]
            if(is.null(weights.)) NULL else weights., # const double weights[]
            DUP = FALSE,
            PACKAGE="earth")

        # above always evaluates all subsets, so trim back to nprune

        list(rval$rss.per.subset[1:nprune], rval$prune.terms[1:nprune, 1:nprune, drop=FALSE])
    }
    # eval.model.subsets.using.xtx starts here

    rprune <- switch(match.arg1(pmethod),
                backward(bx, y, weights),    # "backward"
                backward(bx, y, weights),    # "none"
                bad.pmethod(),
                bad.pmethod(),
                bad.pmethod())
}

forward.pass <- function(x, y, weights,  # must be all double
                         scale.y, trace, penalty, nk, degree, linpreds,
                         allowed, thresh, minspan, newvar.penalty,
                         fast.k, fast.beta, Use.beta.cache,
                         n.allowedfunc.args, env)
{
    npreds <- ncol(x)

    # CHANGED Mar 12, 2008: check.vec() moved to ForwardPass in earth.c
    # because it uses lots of memory even if we do a gc afterwards, not
    # sure why. Added NAOK=TRUE to .C call below to allow this.
    #
    # check.vec("x", x, xvec)
    # check.vec("y", y, yvec)

    stopifnot(nrow(x) == nrow(y))
    fullset <- rep(0, nk)   # element will be set TRUE if corresponding term used

    linpreds <- convert.linpreds.to.logical(linpreds, npreds, x)
    if(trace >= 2)
        print.linpreds(linpreds, x)

    # CHANGED Mar 12, 2008: We want DUP=FALSE in the .C call below to
    # reduce memory useage.  But this isn't permitted if you pass a
    # string array to .C, so we set pred.names=NULL unless tracing
    # or unless we need pred.names for the "allowed" function.

    if(trace >= 2 || n.allowedfunc.args >= 4)
        pred.names <- colnames(x)
    else
        pred.names <- NULL

    nbytes <- 8 * (nk^2 * ncol(x) + (nrow(x) * (3 + 2*nk + ncol(x)/2)))
    if(nbytes > 1e7)  # need more than 10 MBytes to build model?
        gc()

    # CHANGED Jan 22, 2009: scale y for better stability in the forward pass
    if(scale.y) {
        y.scaled <- scale(y)
        # make sure that scaling was ok
        i <- which(attr(y.scaled,"scaled:scale") == 0)
        if(length(i)) {
            if(ncol(y) > 1)
                stop1("cannot scale column ", i[1], " of y (values are all equal to ", y[1,i], ")")
            else
                stop1("cannot scale y (values are all equal to ", y[1,1], ")")
        }
    } else
        y.scaled <- y

    on.exit(.C("FreeR",PACKAGE="earth")) # frees memory if user interupts

    rval <- .C("ForwardPassR",
        fullset = as.integer(fullset),          # out: int FullSet[]
        bx = matrix(0, nrow=nrow(x), ncol=nk),  # out: double bx[]
        dirs = matrix(0, nrow=nk, ncol=npreds), # out: double Dirs[]
        cuts = matrix(0, nrow=nk, ncol=npreds), # out: double Cuts[]
        x,                              # in: const double x[]
        y.scaled,                       # in: const double y[]
        weights,                        # in: const double weightsArg[]
        as.integer(nrow(x)),            # in: const int *pnCases
        as.integer(ncol(y)),            # in: const int *pnResp
        as.integer(npreds),             # in: const int *pnPreds
        as.integer(degree),             # in: const int *pnMaxDegree
        as.integer(nk),                 # in: const int *pnMaxTerms
        as.double(penalty),             # in: const double *pPenalty
        as.double(thresh),              # in: double *pThresh
        as.integer(minspan),            # in: const int *pnMinSpan
        as.integer(fast.k),             # in: const int *pnFastK
        as.double(fast.beta),           # in: const double *pFastBeta
        as.double(newvar.penalty),      # in: const double *pNewVarPenalty
        as.integer(linpreds),           # in: const int LinPreds[]
        allowed,                        # in: const SEXP Allowed
        as.integer(n.allowedfunc.args), # in: const int *pnAllowedFuncArgs
        env,                            # in: const SEXP Env for Allowed
        as.integer(Use.beta.cache),     # in: const int *pnUseBetaCache
        as.integer(trace),              # in: const int *pnTrace
        pred.names,                     # in: const char *sPredNames[]
        NAOK = TRUE, # we check for NAs etc. internally in C ForwardPass
        DUP = !is.null(pred.names),     # see above comment
        PACKAGE="earth")

    fullset <- as.logical(rval$fullset)

    list(rval$bx[, fullset, drop=FALSE],
        rval$dirs[fullset, , drop=FALSE],
        rval$cuts[fullset, , drop=FALSE])
}

# Return a vec which specifies the degree of each term in dirs.
# Each row of dirs specifies one term so we work row-wise in dirs.

get.degrees.per.term <- function(dirs)
{
    if(nrow(dirs) == 1)            # intercept only model?
        return(0)
    degrees <- double(nrow(dirs))
    for(i in seq_along(degrees))
        degrees[i] <- sum(dirs[i,] != 0)
    degrees
}

# Used when building the model to name the columns of bx and rows of dirs etc.
# Also called by mars.to.earth.
# Return string like "h(55-x1)*h(x2-58)".
# h represents the hockey stick func.
# If ntermsVec is a vector, this returns a vector of strings.
# x can be NULL (currently only when called from mars.to.earth), it is used
# only for simplifying terms with factor predictors.

get.earth.term.name <- function(ntermsVec, dirs, cuts, pred.names, x)
{
    get.term.name1 <- function(nterm, dirs, cuts, pred.names, xrange, form1, form2)
    {
        get.name <- function(ipred) # return "name" if possible, else "x[,i]"
        {
            pred.name <- pred.names[ipred]
            if(is.null(pred.name) || is.na(pred.name))
                paste("x[,", ipred, "]", sep="")
            else
                pred.name
        }
        if(nterm == 1)
            return("(Intercept)")
        s <- ""
        first.fac <- TRUE
        stopifnot(ncol(dirs) > 0)
        for(ipred in 1:ncol(dirs))
            if(dirs[nterm,ipred]) {
                if(!first.fac)
                    s <- paste(s, "*", sep="")
                first.fac <- FALSE
                if(dirs[nterm,ipred] == 2)  # linear predictor?
                    s <- pastef(s, "%s", get.name(ipred))
                else if(dirs[nterm,ipred] == -1) {
                    s <- pastef(s, form1, cuts[nterm,ipred], get.name(ipred))
                } else if(dirs[nterm,ipred] == 1) {
                    if(cuts[nterm,ipred] == 0 && !is.null(xrange) &&
                            xrange[1, ipred] == 0 && xrange[2, ipred] < 100 &&
                            x[,ipred] == floor(x[,ipred])) # all integer?
                        # simplify to no hinge function, it's a factor
                        s <- pastef(s, "%s", get.name(ipred))
                    else
                        s <- pastef(s, form2, get.name(ipred), cuts[nterm,ipred])
                } else if(dirs[nterm,ipred] != 0)
                    stop("illegal direction ", dirs[nterm,ipred], " in dirs")
            }
        s
    }
    # --- get.earth.term.name starts here ---
    stopifnot(ncol(dirs) == ncol(x))
    xrange <- NULL      # 1st row is min, 2nd row is max, a column for each pred
    if(!is.null(x))
        xrange <- apply(x, 2, range) # for simplifying "h(ldose-0)" to "ldose"
    # Get format strings for sprintf later
    ndigits <- getOption("digits")
    if(ndigits <= 7) {      # for back compat with previous versions of earth
        form1 <- "h(%g-%s)" # let %g figure out the nbr of digits
        form2 <- "h(%s-%g)"
    } else {
        form1 <- sprintf("h(%%.%dg-%%s)", ndigits) # e.g. "h(%.9g-%s)"
        form2 <- sprintf("h(%%s-%%.%dg)", ndigits) # e.g. "h(%s-%.9g)"
    }
    term.names <- sapply(seq_along(ntermsVec), get.term.name1, dirs, cuts,
                         pred.names, xrange, form1, form2)
    # check for duplicate term names
    duplicated <- duplicated(term.names)
    if(any(duplicated))
        warning1("duplicate term name \"", term.names[which(duplicated)[1]], "\"\n",
                 "This is usually caused by cuts that are very close to each other\n",
                 "Remedy: use options(digits=NDIGITS), ",
                 "typically NDIGITS has to be at least 7 ",
                 "(currently NDIGITS=", ndigits, ")")
    term.names
}

# get.gcv returns GCVs as defined in Friedman's MARS paper, with an
# extension for penalty < 0

get.gcv <- function(    # default Get.crit function
    rss.per.subset,
    ntermsVec,          # number of MARS regression terms including intercept
    penalty,            # penalty per knot, argument from earth.fit()
    ncases)             # number of cases
{
    stopifnot(length(rss.per.subset) == length(ntermsVec))
    nknotsVec <- get.nknots(ntermsVec)
    nparams <- effective.nbr.of.params(ntermsVec, nknotsVec, penalty)

    if(max(nparams, na.rm=TRUE) >= ncases)
        warning1("effective number ", max(nparams, na.rm=TRUE),
                 " of GCV parameters >= number ", ncases, " of cases")

    rss.per.subset / (ncases * (1 - nparams/ncases)^2)
}

# Return the estimated number of knots
#
# TODO This is not quite correct.  It assumes that each term pair adds one
# knot.  Thus each term adds "half a knot".  But if we have deleted a term
# in a pair then the remaining term should add a knot, not half a knot.

get.nknots <- function(nterms)
{
    (nterms - 1 ) / 2
}

get.nterms.per.degree <- function(object, which.terms = object$selected.terms)
{
    check.classname(object, deparse(substitute(object)), "earth")
    check.which.terms(object$dirs, which.terms)
    table(get.degrees.per.term(object$dirs[which.terms, , drop=FALSE]))
}

get.nused.preds.per.subset <- function(dirs, which.terms)
{
    # object was converted from mars? if so, ugly hack to allow plot routines to work
    if(is.null(which.terms))
        which.terms <- matrix(1:ncol(dirs), ncol(dirs), ncol(dirs))

    # allow which.terms to be a vector or matrix
    if(NROW(which.terms) == 1 || NCOL(which.terms) == 1)
        which.terms <- matrix(which.terms, nrow=1,
                            ncol=NROW(which.terms) * NCOL(which.terms == 1))

    nmodels <- NROW(which.terms)
    stopifnot(nmodels > 0)
    nused <- vector(mode="numeric", nmodels)
    for(i in 1:nmodels) {
        check.which.terms(dirs, which.terms)
        nused[i] <- sum(0 != colSums(abs(
                             dirs[which.terms[i,,drop=FALSE], , drop=FALSE])))
    }
    nused
}

get.rsq <- function(rss, rss.null)
{
    1 - rss / rss.null
}

# remove useless(?) "1" "2" "3" ... rownames for x (added by
# model.matrix) so earth.formula x is same as earth.default x,
# and to save memory (although not as important now that R
# hashes strings internally).

possibly.delete.rownames <- function(x)
{
    if(!is.null(rownames(x)) && # decide by looking at first few names
            (is.null(rownames(x)[1]) || rownames(x)[1] == "1") &&
            (is.null(rownames(x)[2]) || rownames(x)[2] == "2") &&
            (is.null(rownames(x)[3]) || rownames(x)[3] == "3"))
        NULL
    else
        rownames(x)
}

print.linpreds <- function(linpreds, x)
{
    if(any(linpreds != 0)) {
        cat("linear predictors ")
        colnames. <- colnames(x)
        index <- (1:length(linpreds))[linpreds]
        if(!is.null(colnames.))
            cat(paste(index, "=", colnames.[linpreds], sep="", collapse=" "))
        else
            cat(paste.with.space((1:length(linpreds))[linpreds]))
        cat("\n")
    }
}

print.pruning.pass <- function(     # this is the default Print.pruning.pass
    trace,
    pmethod,
    penalty,
    nprune,
    selected.terms,
    prune.terms,
    rss.per.subset,
    gcv.per.subset,
    dirs)
{
    nselected <- length(selected.terms)
    prev.grsq <- 0
    if(trace >= 3) {
        cat("Subset size        GRSq     RSq  DeltaGRSq nPreds")
        if(trace >= 4)
            cat("  Terms (col nbr in bx)")
        cat("\n")
        for(iterm in seq_along(rss.per.subset)) {
            grsq <- get.rsq(gcv.per.subset[iterm], gcv.per.subset[1])
            delta.grsq <- grsq - prev.grsq
            prev.grsq <- grsq
            selected <- prune.terms[iterm,]
            selected <- selected[selected != 0]
            cat(if(iterm==nselected) "chosen " else "       ",
                format(iterm, width=4),
                sprintf("%12.4f ", grsq),
                sprintf("%7.4f",   get.rsq(rss.per.subset[iterm], rss.per.subset[1])),
                sprintf("%11.4f ", delta.grsq),
                sprintf("%6d",     get.nused.preds.per.subset(dirs, selected)),
                "  ", sep="")
            if(trace >= 4)
                cat(selected)
        cat("\n")
        }
    cat("\n")
    }
    if(trace >= 1) {
        cat("Prune method \"", pmethod[1], "\" penalty ", penalty,
                " nprune ", nprune, ": selected ", nselected, " of ", sep="")
        selected <- prune.terms[nselected,]
        selected <- selected[selected != 0]
        cat(nrow(dirs), "terms, and",
            get.nused.preds.per.subset(dirs, selected),
            "of", ncol(dirs), "predictors\n")
        cat("After backward pass GRSq",
            format(get.rsq(gcv.per.subset[nselected], gcv.per.subset[1]), digits=4),
            "RSq",
            format(get.rsq(rss.per.subset[nselected], rss.per.subset[1]), digits=4),
            "\n")
    }
    NULL
}

pruning.pass <- function(x, y, weights, subset,
                         trace, penalty, pmethod, nprune, bx, dirs, cuts,
                         Get.crit, Eval.model.subsets, Print.pruning.pass, Force.xtx.prune)
{
    get.nprune <- function()  # convert user's nprune argument to a valid value
    {
        nterms <- nrow(dirs)
        if(is.null(nprune))
            nprune <- nterms
        if(nprune > nterms) {
            warning1("specified \"nprune\" ", nprune,
                    " is greater than the number ", nterms, " of available model terms ",
                    "\nForcing \"nprune\" to ", nterms)
            nprune <- nterms
        }
        if(nprune < 1)
            stop1("\"nprune\" is less than 1")
        nprune
    }
    # pruning.pass starts here

    if(!is.null(bx) && is.null(subset))
        bx.prune <- bx   # use bx created in forward pass
    else                 # else recreate bx with all terms (update invoked earth)
        bx.prune <- get.bx(x, 1:nrow(dirs), dirs, cuts) # all terms
    stopifnot(nrow(bx) == nrow(y))
    nprune <- get.nprune()
    rval <- Eval.model.subsets(bx.prune, y, weights,
                               pmethod, nprune, Force.xtx.prune)
      rss.per.subset <- rval[[1]]   # RSS for each subset (across all responses)
      prune.terms    <- rval[[2]]   # terms in each subset

    # sanity checks here because Eval.model.subsets may be user written
    stopifnot(length(rss.per.subset) == nprune)
    stopifnot(nrow(prune.terms) == nprune)
    stopifnot(ncol(prune.terms) == nprune)

    stopif(any(prune.terms[,1] != 1))   # check intercept column
    gcv.per.subset <- Get.crit(rss.per.subset, 1:nprune, penalty, nrow(bx.prune))
    if(!all(is.finite(rss.per.subset)))
        warning1("earth: non finite RSS in model subsets ",
                 "(see the rss.per.subset returned by earth)")
    else if(!all(is.finite(gcv.per.subset)))
        warning1("earth: non finite GCV in model subsets ",
                 "(see the gcv.per.subset returned by earth)")
    do.prune <- pmatch(pmethod[1], "none", 0) != 1
    selected.terms <- 1:nprune  # all terms
    if(do.prune) {
        # choose the subset which has the lowest GCV in the vector of GCVS
        selected.terms <- prune.terms[which.min(gcv.per.subset),]
        selected.terms <- selected.terms[selected.terms != 0]
    }
    Print.pruning.pass(trace, pmethod, penalty, nprune,
        selected.terms, prune.terms, rss.per.subset, gcv.per.subset, dirs)
    list(bx.prune,
        rss.per.subset,     # vector of RSSs for each model (index on subset size)
        gcv.per.subset,     # vector of GCVs for each model (index on subset size)
        prune.terms,        # triangular mat: each row is a vector of term indices
        selected.terms)     # vector of model terms in best model
}

# return a vector of term numbers, ordered as per the "anova" decomposition

reorder.terms.anova <- function(dirs, cuts)
{
    nterms <- nrow(dirs)
    key.degrees <- get.degrees.per.term(dirs)   # sort first on degree
    first.fac.order <- double(nterms)           # order of first factor
    key.x <- double(nterms)                     # order of preds in factors
    if(nterms > 1)
        for(i in 2:nterms) {                    # start at 2 to skip intercept
            used <- which(dirs[i,] != 0)
            first.fac.order[i] <- used[1]
            key.x[i] <- 1e6 * used[1]           # 1st factor
            if(!is.na(used[2])) {               # 2nd factor if any
                key.x[i] <- key.x[i] + 1e3 * used[2]
                if(!is.na(used[3]))             # 3rd factor if any
                    key.x[i] <- key.x[i] + used[3]
            }
    }
    key.linpreds <- double(nterms)              # put lin pred factors first
    key.cuts <- double(nterms)                  # cut values
    if(nterms > 1)
        for(i in 2:nterms) {
            key.linpreds[i] <- -sum(dirs[i, ] == 2)
            key.cuts[i] <- cuts[i, first.fac.order[i]]
    }
    order(key.degrees, key.linpreds, key.x, key.cuts)
}

# Return an index vector suitable for indexing into object$coefficients
# and ordered using the specified "decomp":
#
# "none"  Order the terms as created during the earth forward pass
#
# "anova" Order the terms using the "anova decomposition"
#         i.e. in increasing order of interaction
#
# The first arg is actually an object but called x for consistency with generic

reorder.earth <- function(
    x           = stop("no 'x' arg"),
    which.terms = x$selected.terms,
    decomp      = c("anova", "none"),
    degree      = 99,       # max degree, 0 returns just the intercept
    min.degree  = 0,
    ...)                    # unused
{
    warn.if.dots.used("reorder.earth", ...)
    if(degree < 0)
        stop1("degree ", degree, " < 0")
    if(min.degree < 0)
        stop1("min.degree ", min.degree, " < 0")
    if(degree < min.degree)
        stop1("degree ", degree, " < min.degree ", min.degree)
    check.which.terms(x$dirs, which.terms)
    dirs <- x$dirs[which.terms, , drop=FALSE]
    new.order <- switch(match.arg1(decomp),
                   reorder.terms.anova(                         # anova
                       dirs, x$cuts[which.terms,,drop=FALSE]),
                   1:length(which.terms))                       # none
    degrees <- get.degrees.per.term(dirs[new.order, , drop=FALSE])
    new.order[degrees >= min.degree & degrees <= degree]
}

# We normalize wp to a length equal to the number of responses
# so wp=const gives same the same results as no wp.
# Note that mda::mars normalizes wp to length 1.

standardize.weights <- function(w, meanw)
{
    # TODO finagle zero weights (but should really do what lm.wfit does and delete columns)
    almost.zero <- meanw / 1e8
    w[w < almost.zero] <- almost.zero
    w / meanw
}

# update.earth is based on update.default but:
#
# a) If a forward pass is needed (i.e. regenerate the earth model
#    from scratch) then it removes any "Object" argument from the call.
#
#    Conversely, if the forward pass is unneeded (i.e. we just need to
#    re-prune the earth model) then it adds an "Object" argument to the call.
#
#    The global character vector prune.only.args says which
#    args are needed only for the pruning pass.
#
#    This default decision to do a forward pass or not can be overriden
#    with the ponly argument.
#
# b) This function also deals appropriately with objects that were or were
#    not created using a formula i.e. were created by earth.formula() or
#    by earth.default().
#
# c) This function retrieves x and y from object$x and object$y if need be
#    and also data, weights, wp, and subset.

update.earth <- function(
    object   = stop("no 'object' arg"),
    formula. = NULL,    # formula. is optional
    ponly    = FALSE,   # force prune only, no forward pass
    ...,                # dots passed on to earth()
    evaluate = TRUE)    # for compatibility with generic update
{
    # update.earth starts here

    check.classname(object, deparse(substitute(object)), "earth")
    Call <- object$call
    stopif(is.null(Call))
    do.forward.pass <- FALSE
    if(!is.null(formula.)) {
        if(is.null(Call$formula))
            stop1("\"formula.\" argument cannot be used on ",
                  "objects created without a formula")
        Call$formula <- update.formula(formula(object), formula.)
        do.forward.pass <- TRUE
    }
    # figure out what trace should be

    this.call <- match.call(expand.dots=TRUE)
    trace <- get.update.arg(this.call$trace, "trace", object,
                            trace1=NULL, "update.earth", print.trace=FALSE)
    if(is.null(trace))
        trace <- 0

    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots) > 0) {
        if(any(is.na(pmatch(names(dots), prune.only.args))))
            do.forward.pass <- TRUE
        # conservative approach: always do forward pass if bpairs
        # argument used even if it hasn't changed
        else if(!is.null((dots$glm)$b) || !is.null((Call$glm)$b)) {
            if(trace >= 1)
                cat("update.earth: forcing forward pass because bpairs argument used\n")
            do.forward.pass <- TRUE
        } else if(!is.null(dots$nfold) || !is.null(Call$nfold)) {
            if(trace >= 1)
                cat("update.earth: forcing forward pass because nfold argument used\n")
            do.forward.pass <- TRUE
        }
        existing <- !is.na(match(names(dots), names(Call)))
        for(i in names(dots)[existing])     # replace existing args
            Call[[i]] <- dots[[i]]
        if(any(!existing)) {                # append new args
            Call <- c(as.list(Call), dots[!existing])
            Call <- as.call(Call)
        }
    }
    if(is.null(Call$formula)) {
        Call$x <- get.update.arg(this.call$x, "x", object, trace)
        Call$y <- get.update.arg(this.call$y, "y", object, trace)
    } else
        Call$data <- get.update.arg(this.call$data, "data", object, trace)

    Call$subset  <- get.update.arg(this.call$subset,  "subset",  object, trace)
    Call$weights <- get.update.arg(this.call$weights, "weights", object, trace)
    Call$wp      <- get.update.arg(this.call$wp,      "wp",      object, trace)
    if(ponly)
        do.forward.pass <- FALSE
    Call$Object <- if(do.forward.pass) NULL else substitute(object)
    if(evaluate)
        eval.parent(Call)
    else
        Call
}

#-----------------------------------------------------------------------------
# Method functions for generics in models.R
# All standard method funcs are supported, I think, or give a warning.
# TODO Some of the functions that give just a warning could possibly give
# a meaningful value.

deviance.earth <- function(object, warn=TRUE, ...)
{
    if(warn && !is.null(object$glm.list))
        warning1("deviance.earth: returning earth (not GLM) deviance")
    object$rss
}

effects.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning1("effects.earth: returning NULL")
    NULL
}

family.earth <- function(object, ...)
{
    family(object$glm.list[[1]])
}

anova.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning1("anova.earth: returning NULL")
    NULL
}

# variable.names.earth returns "name" if possible, else return "x[,i]"

variable.names.earth <- function(object, ..., use.names=TRUE)
{
    warn.if.dots.used("variable.names.earth", ...)
    ipred <- 1:ncol(object$dirs)
    var.name <- NULL
    if(use.names)
        var.name <- colnames(object$dirs)[ipred]
    if(!use.names || is.null(var.name) || is.na(var.name))
        paste("x[,", ipred, "]", sep="")
    else
        var.name
}

case.names.earth <- function(object, ..., use.names=TRUE)
{
    if(is.null(row.names(object$residuals)))
        paste(1:nrow(object$residuals))
    else
        row.names(object$residuals)
}

# type is one of "earth", "deviance", "pearson", "working", "response", "partial"

residuals.earth <- function(object = stop("no 'object' arg"), type = NULL, warn=TRUE, ...)
{
    glm.resids <- function(object, type)
    {
        g <- object$glm.list
        if(is.null(g))
            stop1("residuals.earth: type \"", type, "\" can only be used on earth-glm models")
        colnames <- ""
        for(imodel in seq_along(g)) {
            rval1 <- residuals(g[[imodel]], type)
            if(imodel == 1)
                rval <- rval1
            if(NROW(rval1) != NROW(rval)) # should never happen
                stop1("residuals.earth: glm.list[[", imodel, "]] does ",
                      "not conform to glm.list[[", 1, "]] (residuals have a different length)")
            if(imodel > 1) {
                colnames <- c(colnames)
                rval <- cbind(rval, rval1)
            }
        }
        rval
    }
    # residuals.earth starts here

    warn.if.dots.used("residuals.earth", ...)
    if(warn && is.null(type) && !is.null(object$glm.list))
        warning1("residuals.earth: returning earth (not glm) residuals")

    if(is.null(type))
        type <- "earth"
    itype <- match.choices(type,
                  choices=c("earth", "deviance", "pearson", "working", "response", "partial"),
                  arg.name="type")
    rval <- switch(itype,
        object$residuals,                               # earth
        if(is.null(object$glm.list)) object$residuals   # deviance
            else glm.resids(object, type),
        glm.resids(object, type),                       # pearson
        glm.resids(object, type),                       # working
        glm.resids(object, type),                       # response
        glm.resids(object, type))                       # partial

    if(!is.matrix(rval))
        rval <- matrix(rval, ncol = 1)
    if(type != "partial")
        colnames(rval) <- colnames(object$residuals)
    rownames(rval) <- case.names(object)
    rval
}

resid.earth <- function(object = stop("no 'object' arg"), type = NULL, warn=TRUE, ...)
{
    residuals.earth(object, type, warn, ...)
}

# Fake the AIC by returning the GCV. This is enough for step() to work.

extractAIC.earth <- function(fit, scale = 0, k = 2, warn=TRUE, ...)
{
    if(warn)
        warning1("extractAIC.earth: returning GCV instead of AIC")
    if(scale != 0)
        warning1("extractAIC.earth: ignored scale parameter ", scale)
    if(k != 2)
        warning1("extractAIC.earth: ignored k parameter ", k)
    warn.if.dots.used("extractAIC.earth", ...)
    nterms <- length(fit$selected.terms)
    c(effective.nbr.of.params(nterms, get.nknots(nterms), fit$penalty), fit$gcv)
}

#-----------------------------------------------------------------------------
# Method functions for plotmo.R

# get.singles.earth and get.pairs.earth exist because we need to look at
# the built earth model to determine singles and pairs, because:
#
# (i)  get.singles.default and get.pairs.default return ALL singles and
#      pairs, even if unused
#
# (ii) earth CREATES pairs,  whereas lm will only have pairs if forms such
#      as x1:x2 of x1*x2 appear in the formula

get.singles.earth <- function(object, x, degree1, pred.names, trace)
{
    dataClasses <- attr(object$terms, "dataClasses")
    if(is.character(degree1) && !is.na(pmatch(degree1, "all"))) {
        # user wants all used predictors
        used.vars <- NULL
        selected <- object$selected.terms[reorder.earth(object, decomp="anova")]
        if(length(selected) > 0)
            used.vars <- unique(
                        which(object$dirs[selected, , drop=FALSE] != 0, arr.ind=TRUE)[,2])
        return(used.vars)
    }
    check.classname(object, deparse(substitute(object)), "earth")
    if(any((dataClasses == "factor") | (dataClasses == "ordered"))) {
        # NOV 2008: new code, only use if factors in x
        # TODO this can give extra predictors if variable names alias
        #      e.g. "x" and "x1" are both variable names
        used.colnames <- apply(object$dirs, 2, any1)
        colnames <- colnames(object$dirs)[used.colnames]
        used.preds <- NULL
        for(ipred in seq_along(object$namesx.org)) {
            if(is.factor(x[,ipred])) {
                # This knows how to deal with expanded factor names because
                # it e.g. looks for "^pclass" in "pclass3rd"
                if(length(grep(paste("^", object$namesx.org[ipred], sep=""), colnames)) > 0)
                    used.preds <- c(used.preds, ipred)
            } else {
                # exact match
                if(length(grep(paste("^", object$namesx.org[ipred], "$", sep=""), colnames)) > 0)
                    used.preds <- c(used.preds, ipred)
                used.preds <- c(used.preds, ipred)
            }
        }
        Singles <- unique(used.preds)
    } else {
        # original code, use if no factors in x
        Singles <- NULL
        selected <- object$selected.terms[  # selected is all degree 1 terms
                        reorder.earth(object, decomp="anova", degree=1, min.degree=1)]
        if(length(selected) > 0)
            Singles <- unique(
                        which(object$dirs[selected, , drop=FALSE] != 0, arr.ind=TRUE)[,2])
    }
    if(length(Singles) == 0 && is.specified(degree1))
        warning1("\"degree1\" specified but no degree1 plots")
    if(length(Singles) > 0) {
        check.index.vec("degree1", degree1, Singles)
        Singles <- Singles[degree1]
    }
    Singles # Singles is a vector of indices of predictors for degree1 plots

}

get.pairs.earth <- function(object, x, degree2, pred.names, trace=FALSE)
{
    if(is.character(degree2) && !is.na(pmatch(degree2, "all"))) {
        # user wants all combos of all used predictors
        used.vars <- NULL
        selected <- object$selected.terms[reorder.earth(object, decomp="anova")]
        if(length(selected) > 0)
            used.vars <- unique(
                        which(object$dirs[selected, , drop=FALSE] != 0, arr.ind=TRUE)[,2])
        if(length(used.vars) == 0)
            return(matrix(0, nrow=0, ncol=2)) # no pairs
        col1 <- rep(used.vars, times=length(used.vars))
        col2 <- rep(used.vars, each=length(used.vars))
        Pairs <- cbind(col1, col2)
        Pairs <- Pairs[col1 != col2, , drop=FALSE]
        return(unique(t(apply(Pairs, 1, sort)))) # remove duplicate pairs
    }
    Pairs <- matrix(0, nrow=0, ncol=2)      # no pairs
    selected <- object$selected.terms[      # selected is all degree 2 terms
                    reorder.earth(object, decomp="anova", degree=2, min.degree=2)]
    Pairs <- vector(mode="numeric")
    for(i in selected)                      # append indices of the two preds in term i
        Pairs <- c(Pairs, which(object$dirs[i,] != 0))
    Pairs <- unique(matrix(Pairs, ncol=2, byrow=TRUE))
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("\"degree2\" specified but no degree2 plots")
    if(nrow(Pairs) > 0) { # any pairs?
        check.index.vec("degree2", degree2, Pairs)
        Pairs <- Pairs[degree2, , drop=FALSE]
        if(nrow(Pairs) && any(sapply(x, is.factor))) { # any columns in x are factors?
            # Pairs works off expanded factor names, so replace each name with
            # with index of original variable name
            # TODO this can give wrong results if variable names alias
            #      e.g. if "x" and "x1" are both variable names
            #      this takes the LAST of the matching names so correct with "x" "x1" but not "x1" "x"
            dir.colnames <- colnames(object$dirs)
            prednames <- object$namesx.org
            prednames.hat <- paste("^", prednames, sep="")
            for(i in 1:nrow(Pairs))
                for(j in 1:2) {
                    ipred1 <- 0
                    for(ipred in seq_along(prednames.hat))
                        if(length(grep(prednames.hat[ipred], dir.colnames[Pairs[i, j]])) > 0)
                            ipred1 <- ipred
                    if(ipred1 == 0)
                        stop1("internal error: illegal ipred1 in get.pairs.earth")
                    Pairs[i, j] <- ipred1
                }
            Pairs <- unique(Pairs)  # unique is needed if multiple factors were converted to single predictor
        }
    }
    Pairs
}

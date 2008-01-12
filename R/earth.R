# earth.R: an implementation of Friedman's Multivariate Adaptive
#          Regression Splines, commonly known as MARS.
#
# This code is derived from code in mda.R by Hastie and Tibshirani.
# Comments containing "$$" mark known issues.
# Stephen Milborrow Mar 2007 Petaluma
#
#--------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------
# Notes for earth() that didn't make it into the man pages.
#
# --- subset argument (for selecting cases)
#
# All subset handling is done in earth.default not in earth.formula or
# update.earth.  This is because we want to allow the user to specify
# a subset even when he or she isn't using the formula based approach
# i.e.  using earth.default() directly and not earth.formula().
#
# --- Eval.model.subsets (and Print.pruning.pass)
#
# This function evaluates subsets of the MARS model.  It returns (in
# prune.terms) the best terms to use for each number of terms i.e. for
# each model size.  See the default function eval.model.subsets() for an
# example Eval.model.subsets.
#
# $$ feature: add get.var.importance function
# $$ there are places where efficiency could be improved I think
#
#--------------------------------------------------------------------------------------------

# This is a list of those formal arguments of earth.default that can be changed without
# requiring a new forward pass.
# NOTE: if you change the pruning formal arguments in earth.default(), update this too!

prune.args.global <- c("trace", "pmethod", "ppenalty", "nprune",
                        "Get.crit", "Eval.model.subsets", "Print.pruning.pass",
                        "Force.xtx.prune", "Use.beta.cache")

#--------------------------------------------------------------------------------------------
earth <- function(...) UseMethod("earth")

earth.default <- function(
    x       = stop("no 'x' arg"), # NAs are not allowed in x or y, an error msg if so
    y       = stop("no 'y' arg"),
    weights = NULL,         # not yet supported
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail

    penalty = if(degree > 1) 3 else 2,
                            # GCV penalty per knot:
                            #   0 penalizes only terms (not knots)
                            #   special case -1 means no penalty (so GRSq==RSq)

    trace = 0,              # 0 none 1 overview 2 forward 3 pruning 4 more pruning 5 ...

    keepxy  = FALSE,        # true to retain x and y (before taking subset) in returned value

                            #------------------------------------------------------------
                            # Following affect forward pass only, not pruning pass

    nk      = max(21, 2 * NCOL(x) + 1),
                            # max number of model terms including intercept

    degree         = 1,     # max degree of interaction (1=additive model) (Friedman's mi)

    linpreds       = FALSE, # index vector specifying which preds that must enter linearly

    allowed        = NULL,  # constraint function specifying allowed interactions

    thresh         = 0.001, # used as one of the conditions to stop adding terms in forw pass
                            # stop if RSqDelta<thresh or 1-RSq<thresh

    minspan        = 1,     # consider knots that are minspan apart
                            # special value 0 means use internally calculated min span

    newvar.penalty = 0,     # penalty for adding a new variable in forward pass

    fast.k         = 20,    # Fast MARS K: 0 means use all terms i.e. no Fast MARS
    fast.beta      = 1,     # Fast MARS ageing coefficient

                            #------------------------------------------------------------
                            # Following affect pruning only, not forward pass
                            # If you change these, update prune.args.global too!

    pmethod = "backward",   # for legal values see eval.model.subsets.*
    ppenalty = penalty,     # like penalty, but for the pruning pass
    nprune  = NULL,         # max nbr of terms (including intercept) in prune subset

    Object  = NULL,         # if null, recreate earth model from scratch with forward pass
                            # if not null: no forward pass, just pruning pass

    Get.crit           = get.gcv,            # criterion func for model select during pruning
    Eval.model.subsets = eval.model.subsets, # function used to evaluate model subsets
    Print.pruning.pass = print.pruning.pass, # function used to print pruning pass results
    Force.xtx.prune    = FALSE, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    Use.beta.cache     = TRUE,  # TRUE to use beta cache, for speed
    ...)                    # unused
{
    forward.pass <- function()
    {
        check.vec <- function(xname, x, xvec)
        {
            if(any(is.na(x)))
                stop1("NA in '", xname, "'")
            else if(any(is.na(xvec)))
                stop1("NA in '", xname, "' after conversion to numeric")
            else if(any(!is.finite(xvec)))
                stop1("non finite value in '", xname, "' after conversion to numeric")
        }
        # forward.pass starts here
        npreds <- ncol(x)
        xvec <- as.vector(x, mode="numeric") # convert ints etc. to double for C code
        yvec <- as.vector(y, mode="numeric")
        check.vec("x", x, xvec)
        check.vec("y", y, yvec)
        stopifnot(nrow(x) == nrow(y))
        fullset <- rep(0, nk)   # element will be set TRUE if corresponding term used

        linpreds1 = rep(0, npreds)
        check.index.vec("linpreds", linpreds, x, use.as.col.index=TRUE)
        linpreds1[linpreds] = TRUE

        if(is.null(colnames(x))) {
           # Generate predictor names (empty strings).  Needed because there seems to be
            # no easy way to pass a C NULL to a C function (passing R NULL doesn't do it).
            pred.names <- rep("", npreds)
        } else
            pred.names <- colnames(x)
        on.exit(.C("FreeR",PACKAGE="earth")) # frees memory if user interupts

        rval <- .C("ForwardPassR",
            fullset = as.integer(fullset),          # out: int FullSet[]
            bx = matrix(0, nrow=nrow(x), ncol=nk),  # out: double bx[]
            dirs = matrix(0, nrow=nk, ncol=npreds), # out: double Dirs[]
            cuts = matrix(0, nrow=nk, ncol=npreds), # out: double Cuts[]
            xvec,                           # in: const double x[]
            yvec,                           # in: const double y[]
            as.integer(nrow(x)),            # in: const int *pnCases
            as.integer(ncol(y)),            # in: const int *pnResponses
            as.integer(npreds),             # in: const int *pnPreds
            as.integer(degree),             # in: const int *pnMaxDegree
            as.integer(nk),                 # in: const int *pnMaxTerms
            as.double(penalty),             # in: const double *pPenalty
            as.double(thresh),              # in: double *pThresh
            as.integer(minspan),            # in: const int *pnMinSpan
            as.integer(fast.k),             # in: const int *pnFastK
            as.double(fast.beta),           # in: const double *pFastBeta
            as.double(newvar.penalty),      # in: const double *pNewVarPenalty
            as.integer(linpreds1),          # in: const int LinPreds[]
            allowed,                        # in: const SEXP Allowed
            parent.frame(),                 # in: const SEXP Env
            as.integer(Use.beta.cache),     # in: const int *pnUseBetaCache
            as.integer(trace),              # in: const int *pnTrace
            as.character(pred.names),       # in: const char *sPredNames[]
            PACKAGE="earth")

        fullset <- as.logical(rval$fullset)

        list(rval$bx[, fullset, drop=FALSE],
            rval$dirs[fullset, , drop=FALSE],
            rval$cuts[fullset, , drop=FALSE])
    }
    pruning.pass <- function()
    {
        get.nprune <- function()  # convert user's nprune argument to something we can use
        {
            nterms <- nrow(dirs)
            if(is.null(nprune))
                nprune <- nterms
            if(nprune > nterms) {
                warning1("specified 'nprune' ", nprune,
                        " is greater than the number of available model terms ", nterms,
                        ", forcing 'nprune' to ", nterms)
                nprune <- nterms
            }
            if(nprune < 1)
                stop1("'nprune' is less than 1")
            nprune
        }
        # pruning.pass starts here

        # get bx.prune and y.prune, with subset handling
        if(!is.null(bx) && is.null(subset))
            bx.prune <- bx   # use bx created in forward pass
        else                 # else recreate bx with all terms (update invoked me)
            bx.prune <- get.bx(x.org, 1:nrow(dirs), dirs, cuts) # all terms

        y.prune <- y.org
        stopifnot(nrow(bx) == nrow(y))
        if(!is.null(subset)) {
            check.index.vec("subset", subset, y.prune, check.empty=TRUE)
            y.prune <- y.prune[subset, , drop=FALSE]
            bx.prune <- bx.prune[subset, , drop=FALSE]
        }
        nprune <- get.nprune()
        rval <- Eval.model.subsets(bx.prune, y.prune, pmethod, nprune, Force.xtx.prune)
          rss.per.subset <- rval[[1]]   # RSS for each subset
          prune.terms    <- rval[[2]]   # terms in each subset
        stopifnot(length(rss.per.subset) == nprune)
        stopifnot(nrow(prune.terms) == nprune)
        stopifnot(ncol(prune.terms) == nprune)
        stopif(any(prune.terms[,1] != 1))   # check intercept column
        gcv.per.subset <- Get.crit(rss.per.subset, 1:nprune, ppenalty, nrow(bx.prune))
        if(!all(is.finite(rss.per.subset)))
            warning1("non finite RSS in model subsets (see rss.per.subset)")
        else if(!all(is.finite(gcv.per.subset)))
            warning1("non finite GCV in model subsets (see gcv.per.subset)")
        do.prune <- pmatch(pmethod[1], "none", 0) != 1
        selected.terms <- 1:nprune  # all terms
        if(do.prune) {
            # choose the subset which has the lowest GCV in the vector of GCVS
            selected.terms <- prune.terms[which.min(gcv.per.subset),]
            selected.terms <- selected.terms[selected.terms != 0]
        }
        Print.pruning.pass(trace, pmethod, ppenalty, nprune,
            selected.terms, prune.terms, rss.per.subset, gcv.per.subset, dirs)
        list(bx.prune,
            y.prune,
            rss.per.subset,     # vector of RSSs for each model (index on subset size)
            gcv.per.subset,     # vector of GCVs for each model (index on subset size)
            prune.terms,        # triangular mat: each row is a vector of term indices
            selected.terms)     # vector of model terms in best model
    }
    # earth.default starts here

    warn.if.dots.used("earth.default", ...)

    if(trace >= 4) {           # show the call?
        cat("Call: ")
        # $$ there must be a better way of doing this
        if (length(paste(substitute(x))) > 100 || length(paste(substitute(y))) > 100)
            cat("too long to display")
        else
            cat(strip.white.space(format(match.call(expand.dots=TRUE), "earth")))

        cat("\n\n")
    }
    if(nk < 1)
        stop1("'nk' ", nk, " is less than 1")
    if(!identical(na.action, na.fail))
        stop1("illegal 'na.action', only na.action=na.fail is currently allowed")
    if(!is.null(weights) && !isTRUE(all.equal(weights, rep(weights[1], length(weights)))))
        warning1("'weights' are not yet supported by 'earth', ignoring them")
    if(is.logical(trace) && trace) {
        warning1("converted trace=TRUE to trace=4")
        trace <- 4
    }
    # FIXED 30 Oct 2007: changed as.matrix to data.matrix because as.matrix(x)
    # converts x to a character matrix if x has a factor column which
    # creates problems later if update.earth is invoked.
    x <- data.matrix(x)   # if x is a vec this converts it to an ncases x 1 matrix
    y <- data.matrix(y)   # ditto for y
    x.org <- x
    y.org <- y
    if(nrow(x) == 0)
        stop1("no 'x' values")
    if(ncol(x) == 0)    # this happens for example for earth(Volume~Volume,data=trees)
        stop1("no 'x'")
    if(nrow(x) != nrow(y))
        stop1("nrow(x) ", nrow(x), " != nrow(y) ", nrow(y))
    if(!is.null(subset)) {
        check.index.vec("subset", subset, x, check.empty=TRUE)
        x <- x[subset, , drop=FALSE]
        y <- y[subset, , drop=FALSE]
    }
    if(!is.null(allowed)) {
        if(typeof(allowed) != "closure")
            stop("the given \"allowed\" argument is not a function");
        names. <- names(formals(allowed))
        if (length(names.) != 3)
            stop("the given \"allowed\" function does not have 3 arguments")
        if (!identical(names., c("degree", "pred", "parents")))
            stop("the \"allowed\" function needs the following arguments ",
                 paste.quoted.names(c("degree", "pred", "parents")), "\n",
                 "You have ", paste.quoted.names(names.))
    }
    if(is.null(Object)) {
        rval <- forward.pass()
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
    rval <- pruning.pass()
      bx             <- rval[[1]]
      y              <- rval[[2]]
      rss.per.subset <- rval[[3]]   # vector of RSSs for each model (index on subset size)
      gcv.per.subset <- rval[[4]]   # vector of GCVs for each model
      prune.terms    <- rval[[5]]   # triangular mat: each row is a vector of term indices
      selected.terms <- rval[[6]]   # vector of term indices of selected model
    bx <- bx[, selected.terms, drop=FALSE]

    # add names for returned values

    pred.names <- colnames(x)
    term.names <- get.earth.term.name(1:nrow(dirs), dirs, cuts, pred.names)
    colnames(bx)   <- term.names[selected.terms]
    rownames(dirs) <- term.names
    colnames(dirs) <- pred.names
    rownames(cuts) <- term.names
    colnames(cuts) <- pred.names

    # Regress y on bx to get fitted.values etc.
    # The as.matrix calls after the call to lm are needed if y is
    # a vector so the fitted.values etc. are always arrays.

    lfit <- lm.fit(bx, y, singular.ok=FALSE)
       fitted.values <- as.matrix(lfit$fitted.values)
       residuals     <- as.matrix(lfit$residuals)
       coefficients  <- as.matrix(lfit$coefficients)

    # prepare returned summary statistics

    nselected <- length(selected.terms)
    rss  <- rss.per.subset[nselected]
    rsq  <- get.rsq(rss, rss.per.subset[1])
    gcv  <- gcv.per.subset[nselected]
    grsq <- get.rsq(gcv, gcv.per.subset[1])
    nresponses <- ncol(y)
    rss.per.response  <- vector(mode="numeric", length=nresponses)
    rsq.per.response  <- vector(mode="numeric", length=nresponses)
    gcv.per.response  <- vector(mode="numeric", length=nresponses)
    grsq.per.response <- vector(mode="numeric", length=nresponses)
    for(iresponse in 1:nresponses) {
        rss.per.response[iresponse]  <- sum(residuals[,iresponse]^2)
        rss.null                     <- sum((y - mean(y[,iresponse]))^2)
        rsq.per.response[iresponse]  <- get.rsq(rss.per.response[iresponse], rss.null)
        gcv.null                     <- Get.crit(rss.null, 1, ppenalty, nrow(x))
        gcv.per.response[iresponse]  <- Get.crit(rss.per.response[iresponse], nselected,
                                                 ppenalty, nrow(x))
        grsq.per.response[iresponse] <- get.rsq(gcv.per.response[iresponse], gcv.null)
    }
    rval <- structure(list(             # term 1 is the intercept in all returned data
        bx             = bx,            # selected terms only
        dirs           = dirs,          # all terms including unselected: nterms x npreds
        cuts           = cuts,          # all terms including unselected: nterms x npreds
        selected.terms = selected.terms,# row indices into dirs and cuts
        prune.terms    = prune.terms,   # nprune x nprune, each row is vec of term indices

        rss            = rss,           # RSS, across all responses if y has multiple cols
        rsq            = rsq,           # R-Squared, across all responses
        gcv            = gcv,           # GCV, across all responses
        grsq           = grsq,          # GRSq across all responses

        rss.per.response  = rss.per.response,   # nresponses x 1, RSS for each response
        rsq.per.response  = rsq.per.response,   # nresponses x 1, RSq for each response
        gcv.per.response  = gcv.per.response,   # nresponses x 1, GCV for each response
        grsq.per.response = grsq.per.response,  # nresponses x 1, GRSq for each response

        rss.per.subset = rss.per.subset,# nprune x 1, RSS of each model, across all responses
        gcv.per.subset = gcv.per.subset,# nprune x 1, GCV of each model, across all responses

        fitted.values  = fitted.values, # ncases (after subset) x nresponses
        residuals      = residuals,     # ncases (after subset) x nresponses
        coefficients   = coefficients,  # selected terms only: nselected x nresponses

        ppenalty       = ppenalty,      # copy of ppenalty argument
        call           = make.call.generic(
                            strip.dots.from.call(match.call(expand.dots=FALSE)), "earth")),
    class = "earth")
    if(keepxy) {
        rval$x <- x.org
        rval$y <- y.org
        rval$subset <- subset
    }
    rval
}

# NAs are not allowed.
# If there are NAs in the data the call to model.matrix below issues
# the message "Error in na.fail.default: missing values in object"

earth.formula <- function(
    formula = stop("no 'formula' arg"), # intercept will be ignored
    data,
    ...)                                # args same as earth.default, but no x and y
{
    call <- match.call(expand.dots=FALSE)
    if(!is.null(call$na.action) && !identical(na.fail, eval(call$na.action, sys.parent())))
        stop1("illegal 'na.action', only na.action=na.fail is allowed")
    # subset, weights, na.action handled elsewhere, so match only on formula and data
    m <- match(c("formula", "data"), names(call), 0)
    model.frame <- call[c(1, m)]
    model.frame[[1]] <- as.name("model.frame")
    model.frame$na.action <- na.fail    # confirms the default
    model.frame <- eval.parent(model.frame)
    x <- model.matrix(attr(model.frame, "terms"), model.frame)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint)
        x <- x[, -xint, drop=FALSE]     # silently discard intercept
    else
        warning1("ignored -1 in formula (earth objects always have an intercept)")
    y <- model.response(model.frame, "numeric")
    rval <- earth.default(x=x, y=y, ...)
    rval$terms <- attr(model.frame, "terms")
    # $$ would like to strip unused args from call here, use strip.dots.from.call?
    rval$call <- make.call.generic(match.call(), "earth")
    rval
}

# update.earth is based on update.default but:
#
# a) If a forward pass is needed (i.e. regenerate the earth model
#    from scratch) then it removes any "Object" argument from the call.
#
#    Conversely, if the forward pass is unneeded (i.e. we just need to
#    re-prune the earth model) then it adds an "Object" argument to the call.
#
#    The global character vector prune.args.global says which
#    args are needed only for the pruning pass.
#
# b) This function also deals appropriately with objects that were or were
#    not created using a formula i.e. were created by earth.formula() or
#    by calling earth.default() directly.
#
# c) This function retrieves x and y from object$x and object$y if need be

update.earth <- function(
    object   = stop("no 'object' arg"),
    formula. = NULL,                        # formula. is optional
    ...,                                    # args the same as earth()
    evaluate = TRUE)
{
    check.classname(object, deparse(substitute(object)), "earth")
    call <- object$call
    stopif(is.null(call))
    do.forward.pass <- FALSE
    if(!is.null(formula.)) {
        if(is.null(call$formula))
            stop1("'formula.' argument cannot be used on objects created without a formula")
        call$formula <- update.formula(formula(object), formula.)
        do.forward.pass <- TRUE
    }
    this.call <- match.call(expand.dots=FALSE)
    dots <- this.call$...
    if(length(dots) > 0) {
        if(any(is.na(pmatch(names(dots), prune.args.global))))
            do.forward.pass <- TRUE
        existing <- !is.na(match(names(dots), names(call)))
        for(i in names(dots)[existing])     # replace existing args
            call[[i]] <- dots[[i]]
        if(any(!existing)) {                # append new args
            call <- c(as.list(call), dots[!existing])
            call <- as.call(call)
        }
    }
    # Which x should we use? The precedence is [1] the x parameter, if any,
    # in this call to update [2] the $x in the earth object (which exists
    # if keepxy=TRUE was used the original call to earth) [3] the x found
    # by eval.parent().

    this.call <- match.call(expand.dots=TRUE)
    if(is.null(this.call$x)) {  # no x argument in this call to update?
        x <- try(eval.parent(object$x), silent=TRUE)
        if(!is.null(x) && class(x) != "try-error")
            call$x <- x         # use object$x
    }
        if(is.null(this.call$y)) {   # same as above, but for y
                y <- try(eval.parent(object$y), silent=TRUE)
                if(!is.null(y) && class(y) != "try-error")
                        call$y <- y
    }
    call$Object <- if(do.forward.pass) NULL else substitute(object)
    if(evaluate)
        eval.parent(call)
    else
        call
}

eval.model.subsets.with.leaps <- function(
    bx,
    y,
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

    pacify <- pmatch(pmethod[1], "exhaustive", 0) == 1 && nrow(bx) > 30 && nprune > 8
    if(pacify) {                # pruning could take a while so print a reminder
        cat("Pruning...")
        flush.console()
    }
    rprune <- leaps.setup(x=bx, y=y,
             force.in=1,        # make sure intercept is in model
             force.out=NULL,
             intercept=FALSE,   # we have an intercept so leaps.setup must not add one
             nvmax=nprune, nbest=1, warn.dep=TRUE)

    rprune <- switch(match.arg1(pmethod),
                leaps.backward(rprune),     # "backward"
                leaps.backward(rprune),     # "none"
                leaps.exhaustive(rprune, really.big=TRUE),
                leaps.forward(rprune),
                leaps.seqrep(rprune))
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
    backward <- function(bx, y)
    {
        ncases <- nrow(bx)
        nterms <- ncol(bx)
        nresponses <- ncol(y)
        rval <- .C("EvalSubsetsUsingXtxR",
            prune.terms = matrix(0, nrow=nterms, ncol=nterms),  # double PruneTerms[]
            rss.per.subset = vector(mode="numeric", length=nterms),
            as.integer(ncases),                         # const int *pnCases
            as.integer(nresponses),                     # const int *pnresponses
            as.integer(nterms),                         # const int *pnMaxTerms
            as.double(as.vector(bx, mode="numeric")),   # const double bx[]
            as.double(as.vector(y, mode="numeric")),    # const double y[]
            PACKAGE="earth")

        # above always evaluates all subsets, so trim back to nprune

        list(rval$rss.per.subset[1:nprune], rval$prune.terms[1:nprune, 1:nprune, drop=FALSE])
    }
    # eval.model.subsets.using.xtx starts here

    rprune <- switch(match.arg1(pmethod),
                backward(bx, y),    # "backward"
                backward(bx, y),    # "none"
                bad.pmethod(),
                bad.pmethod(),
                bad.pmethod())
}

# This returns the RSS and selected terms for each subset of size 1:nprune

eval.model.subsets <- function(     # this is the default Eval.model.subsets
    bx,         # basis model matrix
    y,          # model response
    pmethod,
    nprune,     # max nbr of terms (including intercept) in prune subset, in range 1..nterms
    Force.xtx.prune) # TRUE to always call EvalSubsetsUsingXtx rather than leaps
{
    stopifnot(nprune >= 1 && nprune <= nrow(bx))
    if(Force.xtx.prune ||   # user explicityly asked for xtx subset evaluation
            ncol(y) > 1 ||  # leaps cannot deal with multiple responses
            ncol(bx) <= 2)  # leaps code gives an error for small number of cols, so avoid
        eval.model.subsets.using.xtx(bx, y, pmethod, nprune, Force.xtx.prune)
    else
        eval.model.subsets.with.leaps(bx, y, pmethod, nprune)
}

print.pruning.pass <- function(     # this is the default Print.pruning.pass
    trace,
    pmethod,
    ppenalty,
    nprune,
    selected.terms,
    prune.terms,
    rss.per.subset,
    gcv.per.subset,
    dirs)
{
    nselected <- length(selected.terms)
    if(trace >= 3) {
        cat("Subset size         GCV         GRSq      RSq  nPreds")
        if(trace >= 4)
            cat("  Terms (index into bx)")
        cat("\n")
        for(iterm in seq_along(rss.per.subset)) {
            selected <- prune.terms[iterm,]
            selected <- selected[selected != 0]
            cat(if(iterm==nselected) "chosen " else "       ",
                format(iterm, width=4),
                sprintf("%12g ",   gcv.per.subset[iterm]),
                sprintf("%12.4f ", get.rsq(gcv.per.subset[iterm], gcv.per.subset[1])),
                sprintf("%8.4f",   get.rsq(rss.per.subset[iterm], rss.per.subset[1])),
                sprintf("%8d",     get.nused.preds.per.subset(dirs, selected)),
                "  ", sep="")
            if(trace >= 4)
                cat(selected)
        cat("\n")
        }
    cat("\n")
    }
    if(trace >= 1) {
        cat("Prune method \"", pmethod[1], "\" ppenalty ", ppenalty,
                " nprune ", nprune, ": selected ", nselected, " of ", sep="")
        selected <- prune.terms[nselected,]
        selected <- selected[selected != 0]
        cat(nrow(dirs), "terms, and",
            get.nused.preds.per.subset(dirs, selected),
            "of", ncol(dirs), "predictors\n")
        cat("GRSq:",
            format(get.rsq(gcv.per.subset[nselected], gcv.per.subset[1]), digits=4),
            "RSq:",
            format(get.rsq(rss.per.subset[nselected], rss.per.subset[1]), digits=4),
            "\n")
    }
    NULL
}

# The first arg is actually an object but called x for consistency with generic

print.earth <- function(x, digits=getOption("digits"), ...)
{
    form <- function(x, pad) sprintf("%-*s", digits+pad, format(x, digits=digits))

    warn.if.dots.used("print.earth", ...)

    cat("Selected", length(x$selected.terms),
        "of", nrow(x$dirs), "terms, and",
        get.nused.preds.per.subset(x$dirs, x$selected.terms),
        "of", ncol(x$dirs), "predictors")

    nlinpreds <- sum(x$dirs[x$selected.terms,] == 2)
    if (nlinpreds == 1)
        cat(" (", nlinpreds, " linear predictor)", sep="")
    else if (nlinpreds > 1)
        cat(" (", nlinpreds, " linear predictors)", sep="")

    nterms.per.degree <- get.nterms.per.degree(x, x$selected.terms)

    cat("\nNumber of terms at each degree of interaction:", nterms.per.degree)

    cat(switch(length(nterms.per.degree),
            " (intercept only model)",
            " (additive model)"),
        "\n", sep="")

    nresponses <- NCOL(x$coefficients)

    if(nresponses > 1) {
        cat("\n")
        for(iresponse in 1:nresponses)
            cat("Response ", iresponse, "  ",
                "GCV: ",  form(x$gcv.per.response[iresponse], 10),
                "RSS: ",  form(x$rss.per.response[iresponse], 10),
                "GRSq: ", form(x$grsq.per.response[iresponse], 8),
                "RSq: ",  form(x$rsq.per.response[iresponse], 0),
                "\n", sep="")
        cat("All         ")
    }
    cat("GCV: ",  form(x$gcv, 10),
        "RSS: ",  form(x$rss, 10),
        "GRSq: ", form(x$grsq, 8),
        "RSq: ",  form(x$rsq, 0),
        "\n", sep="")

    invisible(x)
}

summary.earth <- function(  # returns a superset, not a summary in the strict sense
    object = stop("no 'object' arg"),
    digits = getOption("digits"),
    ...)                    # extra args passed on to format.earth e.g. decomp="none"
{
    rval <- object
    rval$string <- format.earth(object, digits, ...)
    rval$digits <- digits   # allows us to pass digits arg on to print.summary.earth
    class(rval) <- "summary.earth"
    rval
}

# The first arg is actually an object but called x for consistency with generic

print.summary.earth <- function(
    x = stop("no 'x' arg"),     # "summary.earth" object
    digits,
    ...)
{
    warn.if.dots.used("print.summary.earth", ...)
    cat("Call:\n")
    call. = x$call
    if(is.null(x$terms)) {
        # earth was called with the x,y interface i.e. not the formula interface
        # don't print x or y if they are too long
        # $$ there must be a better way of doing this
        x. <- x$call$x
        if (length(paste(substitute(x.))) > 100)
            call.$x = paste("[", NROW(call.$x), ",", NCOL(call.$x),
                            "] too long to display", sep="")
        y. <- x$call$y
        if (length(paste(substitute(y.))) > 100)
            call.$y = paste("[", NROW(call.$y), ",", NCOL(call.$y),
                            "] too long to display", sep="")
    }
    dput(call.)
    cat("\n")
    nresponses <- NCOL(x$coefficients)
    for(iresponse in 1:nresponses)
        cat(if(nresponses == 1)
                "Expression:\n"
            else
                paste("Response ", iresponse, " expression:\n", sep=""),
            x$string[iresponse], sep="")
    cat("\nNumber of cases: ", length(x$residuals), "\n", sep="")
    print.earth(x, if(missing(digits)) x$digits else digits)
    invisible(x)
}

get.rsq <- function(rss, rss.null)
{
    1 - rss / rss.null
}

# Return the estimated number of knots
#
# $$ This is not quite correct.  It assumes that each term pair adds one
# knot.  Thus each term adds "half a knot".  But if we have deleted a term
# in a pair then the remaining term should add a knot, not half a knot.

get.nknots <- function(nterms)
{
    (nterms - 1 ) / 2
}

effective.nbr.of.params <- function(ntermsVec, nknotsVec, penalty)  # for GCV calculation
{
    if(penalty < 0) # special case: term and knots are free so GCV == RSS/ncases
        0
    else
        ntermsVec + (penalty * nknotsVec)
}

# get.gcv returns GCVs as defined in Friedman's MARS paper, with an
# extension for penalty < 0

get.gcv <- function(    # default Get.crit function
    rss.per.subset,
    ntermsVec,          # number of MARS regression terms including intercept
    penalty,            # penalty per knot, argument from earth.default()
    ncases)             # number of cases
{
    stopifnot(length(rss.per.subset) == length(ntermsVec))
    nknotsVec <- get.nknots(ntermsVec)
    nparams <- effective.nbr.of.params(ntermsVec, nknotsVec, penalty)

    if(max(nparams, na.rm=TRUE) >= ncases)
        warning1("GCV effective number of parameters ", max(nparams, na.rm=TRUE),
            " >= number of cases ", ncases)

    rss.per.subset / (ncases * (1 - nparams/ncases)^2)
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

# Return a vec which specifies the degree of each term in dirs.
# Each row of dirs specifies one term so we work row-wise in dirs.

get.degrees.per.term <- function(dirs)
{
    if (nrow(dirs) == 1)            # intercept only model?
        return(0)
    degrees <- double(nrow(dirs))
    for (i in seq_along(degrees))
        degrees[i] = sum(dirs[i,] != 0)
    degrees
}

get.nterms.per.degree <- function(object, which.terms = object$selected.terms)
{
    check.classname(object, deparse(substitute(object)), c("earth", "summary.earth"))
    check.which.terms(object$dirs, which.terms)
    table(get.degrees.per.term(object$dirs[which.terms, , drop=FALSE]))
}

# Return string like "h(55-x1)*h(x2-58)".
# h represents the hockey stick func.
# If ntermsVec is a vector, this returns a vector of strings.

get.earth.term.name <- function(ntermsVec, dirs, cuts, pred.names)
{
    get.term.name1 <- function(nterm, dirs, cuts, pred.names)
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
                if(dirs[nterm,ipred] == 2)  # predictor enters linearly?
                    s <- pastef(s, "%s", get.name(ipred))
                else if(dirs[nterm,ipred] == -1)
                    s <- pastef(s, "h(%g-%s)", cuts[nterm,ipred], get.name(ipred))
                else
                    s <- pastef(s, "h(%s-%g)", get.name(ipred), cuts[nterm,ipred])
            }
        s
    }
    sapply(seq_along(ntermsVec), get.term.name1, dirs, cuts, pred.names)
}

model.matrix.earth <- function(     # returns bx
    object      = stop("no 'object' arg"),
    x           = NULL,
    subset      = NULL,
    which.terms = NULL,
    ...)
{
    warn.if.dots.used("model.matrix.earth", ...)
    check.classname(object, deparse(substitute(object)), "earth")
    if(is.null(x) && is.null(subset) && is.null(which.terms))
        return(object$bx)
    x <- get.earth.x(object, x)
    if(is.null(which.terms))
        which.terms <- object$selected.terms
    if(!is.null(subset)) {
        check.index.vec("subset", subset, x, check.empty=TRUE)
        x <- x[subset, , drop=FALSE]
    }
    get.bx(x, which.terms, object$dirs, object$cuts)
}

get.bx <- function(x, which.terms, dirs, cuts)
{
    stopifnot(all(dirs[1,] == 0))   # intercept term dirs must all be 0
    check.index.vec("which.terms", which.terms, dirs)
    stopifnot(which.terms[1] == 1)
    stopifnot(ncol(x) > 0)
    bx <- matrix(0, nrow=nrow(x), ncol=length(which.terms))
    ibx <- 1
    for(iterm in which.terms) {
        temp1 <- 1
        for(ipred in 1:ncol(x))
            if(dirs[iterm, ipred] == 2)  # predictor enters linearly?
                temp1 <- temp1 * x[, ipred]
            else if(dirs[iterm, ipred] != 0) {
                temp2 <- dirs[iterm, ipred] * (x[, ipred] - cuts[iterm, ipred])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            }
        bx[, ibx] <- temp1
        ibx <- ibx + 1
    }
    colnames(bx) <- rownames(dirs[which.terms,])
    bx
}

get.earth.x <- function(    # returns x matrix
    object      = stop("no 'object' arg"),
    data        = NULL)     # can be a dataframe, matrix, or vector
{
    # allow data to be a vector, if its length matchs ncol(original data)
    convert.data <- function(data)
    {
        if(is.vector(data))
            if(length(data) == ncol(object$dirs))
                data <- matrix(data, nrow=1, ncol=ncol(object$dirs))
            else
                stop1("length(data) ", length(data),
                    " is not equal to the number of predictors ", ncol(object$dirs),
                    " in 'object'")
        data
    }
    # get.earth.x starts here
    if(is.null(object$terms)) {
        # object was created without calling earth.formula

        if(is.null(data))
            x <- eval.parent(object$call$x)
        else
            x <- convert.data(data)
        # FIXED 30 Oct 2007: changed as.matrix to data.matrix
        x <- data.matrix(x)   # to cope with dataframes
        if(is.null(colnames(x)))
            colnames(x) <- colnames(object$dirs)
    } else {
        # object was created by calling earth.formula

        Terms <- delete.response(object$terms)
        if(is.null(data)) {
            object$call$subset <- NULL
            data <- model.frame(object)
        }
        else {
            data <- convert.data(data)
            if(is.null(colnames(data)))
                colnames(data) <- colnames(object$dirs)
            data <- as.data.frame(data)
            data <- model.frame(Terms, data)
            classes <- attr(Terms, "dataClasses")
            if(!is.null(classes))
                .checkMFClasses(classes, data)
        }
        x <- model.matrix(Terms, data)
        xint <- match("(Intercept)", colnames(x), nomatch=0)
        if(xint > 0)
            x <- x[, -xint, drop=FALSE] # silently discard intercept
    }
    # check that the column names of data match predictor names used to build the object

    if(!is.data.frame(data)) {
        data.names <- colnames(data)
        dirs.names <- colnames(object$dirs)
        for(iname in seq_along(dirs.names)) {
            if(!is.null(data.names[iname]) && data.names[iname] != "" &&
                    !is.null(dirs.names[iname]) && dirs.names[iname] != "" &&
                    data.names[iname] != dirs.names[iname]) {
                warning1("the variable names in 'data' do not match those in 'object'",
                        "\n  data:        ", paste.with.space(data.names),
                        "\n  object$dirs: ", paste.with.space(dirs.names))
                break
                }
        }
    }
    if(NCOL(x) != NCOL(object$dirs))
        warning1("ncol(x) ", NCOL(x), " does not match the number of cols ",
            NCOL(object$dirs), " of the 'object' input matrix")
    if(NROW(x) == 0)
        stop1("empty model matrix")
    x
}

predict.earth <- function(
    object  = stop("no 'object' arg"),
    newdata = NULL,
    ...,
    type = c("response", "terms"))  # "terms" returns just the additive terms!
                                    # and just the first response if more than one
{
    get.response <- function()
    {
        if(is.null(newdata))
            object$fitted.values
        else
            model.matrix.earth(object, get.earth.x(object, newdata)) %*% object$coefficients
    }
    get.terms <- function()         # returns just enough for termplot to work
    {
        if(is.null(newdata))
            bx <- object$bx
        else
            bx <- model.matrix(object, get.earth.x(object, newdata))
        dirs <- object$dirs[object$selected.terms, ]
        # retain only additive terms
        additive.terms <- get.degrees.per.term(dirs) == 1
        bx <- bx[, additive.terms]
        dirs <- dirs[additive.terms, ]
        coefs <- object$coefficients[additive.terms, 1, drop=FALSE]
        additive.preds <- colSums(abs(dirs)) != 0
        dirs <- dirs[, additive.preds]
        var.names <- variable.names(object, use.names=TRUE)[additive.preds]
        termMat <- matrix(0, nrow=nrow(bx), ncol=ncol(dirs))
        colnames(termMat) <- var.names
        if(ncol(bx) > 1)
            for(ipred in 1:ncol(dirs))
                for(iterm in 1:ncol(bx))
                    if(dirs[iterm, ipred] != 0)
                        termMat[, ipred] = termMat[, ipred] + coefs[iterm] * bx[, iterm]
        termMat
    }
    # predict.earth starts here
    warn.if.dots.used("predict.earth", ...)
    switch(match.arg1(type),
        get.response(),
        get.terms())
}

check.which.terms <- function(dirs, which.terms) # ensure which.terms is valid
{
    if(is.null(which.terms))
        stop1("'which.terms\' is NULL")
    if(length(which.terms) == 0)
        stop1("length(which.terms) == 0")
    if(which.terms[1] != 1)
        stop1("first element of 'which.terms' must be 1, the intercept term")
    check.index.vec("which.terms", which.terms[1], dirs)
}

factors.present <- function(object)
{
    dataClasses <- attr(object$terms, "dataClasses")
    return(any(dataClasses == "factor" | dataClasses == "ordered"))
}

# return a list of term numbers, ordered as per the "anova" decomposition

anova.decomp <- function(dirs, cuts)
{
    nterms = nrow(dirs)
    key.degrees = get.degrees.per.term(dirs)    # sort first on degree
    first.fac.order <- double(nterms)           # order of first factor
    key.x <- double(nterms)                     # order of preds in factors
    if (nterms > 1)
        for (i in 2:nterms)      {              # start at 2 to skip intercept
            used = which(dirs[i,] != 0)
            first.fac.order[i] <- used[1]
            key.x[i] = 1e6 * used[1]            # 1st factor
            if (!is.na(used[2])) {              # 2nd factor if any
                key.x[i] = key.x[i] + 1e3 * used[2]
                if (!is.na(used[3]))            # 3rd factor if any
                    key.x[i] = key.x[i] + used[3]
            }
    }
    key.linpreds <- double(nterms)              # put lin pred factors first
    key.cuts <- double(nterms)                  # cut values
    if (nterms > 1)
        for (i in 2:nterms) {
            key.linpreds[i] = -sum(dirs[i, ] == 2)
            key.cuts[i] = cuts[i, first.fac.order[i]]
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
    degree      = 99,       # 0 returns just the intercept
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
                   anova.decomp(dirs, x$cuts[which.terms,,drop=FALSE]), # anova
                   1:length(which.terms))                               # none
    degrees <- get.degrees.per.term(dirs[new.order, , drop=FALSE])
    new.order[degrees >= min.degree & degrees <= degree]
}

get.coef.width <- function(coefs, digits)   # get print width for earth coeffs
{
    if(length(coefs) > 0)
        max(nchar(format(abs(coefs), digits=digits)))
    else
        10 # arbitrary width if no coefs
}

format.one.response <- function(
    iresponse,      # response index i.e. column in y matrix
    object,         # "earth" object
    digits,
    use.names,      # use predictor names, else "x[,1]" etc
    add.labels,     # add comments labelling each term
    decomp)         # see reorder.earth for legal decomp values
{
    # get.width returns the width for printing elements of the earth expression.
    # This is used to keep things lined up without too much white space.
    # This returns the widest of all possible printed elements.

    get.width <- function()
    {
        if(length(which.terms) == 1)
            return(10)  # return arbitrary width for intercept only model
        used.dirs <- dirs[which.terms, , drop=FALSE]
        # used.preds is a logical index vector which selects used x predictors
        used.preds <- apply(used.dirs, 2, any)
        # as.list is needed so format treats each cut independently
        max(nchar(var.names[used.preds]),
            nchar(format(as.list(cuts[which.terms, used.preds]), digits=digits)))
    }
    dirs <- object$dirs
    cuts <- object$cuts
    new.order <- reorder.earth(object, decomp=decomp)
    coefs <- object$coefficients[new.order, iresponse]
    which.terms <- object$selected.terms[new.order]
    check.which.terms(object$dirs, which.terms)
    coef.width <- get.coef.width(coefs[-1], digits)
    var.names <- variable.names(object, use.names=use.names)
    width <- get.width()
    s <- ""     # result goes into this string
    s <- pastef(s, "  %.*g %s\n", digits=digits, coefs[1], if(add.labels) " # 1" else"")
    iterm <- 2
    while(iterm <= length(which.terms)) {
        isel.term <- which.terms[iterm]
        dir <- dirs[isel.term, , drop=FALSE]
        cut <- cuts[isel.term, , drop=FALSE]
        coef <- coefs[iterm]
        if(coef < 0)
            s <- pastef(s, "  - %s ",
                  format(-coef, justify="left",w=coef.width,digits=digits,format="%g"))
        else
            s <- pastef(s, "  + %s ",
                  format(coef, justify="left",w=coef.width,digits=digits,format="%g"))
        npreds <- ncol(cuts)
        for(ipred in 1:npreds)
            if(dir[ipred] == -1)
                s <- pastef(s, "* pmax(0, %s - %*s) ",
                        format(cut[ipred], width=width, digits=digits),
                        width, var.names[ipred])
            else if(dir[ipred] == 1)
                s <- pastef(s, "* pmax(0, %*s - %s) ",
                        width=width, var.names[ipred],
                        format(cut[ipred], width=width, digits=digits))
            else if(dir[ipred] == 2) # linear predictor
                s <- pastef(s, "* %-*s %*s            ",
                        width=width, var.names[ipred],
                        width=width, "")
            else if(dir[ipred] != 0)
                stop1("illegal dir in 'dirs'")
        if(add.labels)
            s <- pastef(s, " # %d", iterm)
        s <- pastef(s, "\n")
        iterm <- iterm + 1
    }
    s
}

# format.earth returns a string representing the earth model.
# There is one term per line. Each term (except the intercept)
# is made up of a coefficent which multiplies one or more hockey stick funcs.
# The result looks like this:
#
#         23.208244
#         +  5.7459616 * pmax(0,  Girth -   12.9)
#         -  2.8664516 * pmax(0,   12.9 -  Girth)
#         + 0.71833643 * pmax(0, Height -     76)
#
# If there are multiple responses, this function returns a vector of strings
#
# decomp argument: see reorder.earth()
#
# The first arg is actually an object but called x for consistency with generic
#
# $$ would be nice to add an option to first list bx terms and then combine them
# $$ would be nice to add an option to print in term importance order
# $$ would be simpler using print with print.gap rather than my hand made alignment?

format.earth <- function(
    x           = stop("no 'x' arg"),   # "earth" object, called x for consist with generic
    digits      = getOption("digits"),
    use.names   = TRUE,                 # use predictor names, else "x[,1]" etc
    add.labels  = FALSE,                # add comments labelling each term
    decomp      = "anova",              # see reorder.earth for legal decomp values
    ...)                                # unused
{
    check.classname(x, deparse(substitute(a)), "earth")
    warn.if.dots.used("format.earth", ...)
    nresponses <- NCOL(x$coefficients)
    s <- vector(mode = "character", length=nresponses)
    for(iresponse in 1:nresponses)
        s[iresponse] <- format.one.response(iresponse,x,digits,use.names,add.labels,decomp)
    s
}

#--------------------------------------------------------------------------------------------
# This func is commented out because it doesn't really belong in the earth package.
# It seems to work but I haven't done complete testing.
# Example: a <- lm(Volume ~ ., data = trees); cat(format(a))
#
# # Return a string representing the linear model.
# # The result looks like this:
# #
# #   -58
# #   +  4.71 * Girth
# #   + 0.339 * Height
# #
# # The first arg is actually an object but called x for consistency with generic
#
# format.lm <- function(
#   x           = stop("no 'x' arg"),    # "lm" object
#   digits      = getOption("digits"),
#   use.names   = TRUE,
#   add.labels  = FALSE)                 # add comments labelling each term
# {
#   get.name <- function(ipred) # return "name" if possible, else "x[,i]"
#   {
#       pred.name <- variable.names(x)
#       if(!is.null(pred.name))
#           pred.name <- pred.name[ipred+1] # +1 to skip intercept
#       if(is.null(pred.name) || is.na(pred.name) || !use.names)
#           paste("x[,", ipred, "]", sep="")
#       else
#           pred.name
#   }
#   check.classname(x, deparse(substitute(x)), "lm")
#   coefs = coef(x)
#   s <- sprintf("  %.*g%s\n", digits=digits, coefs[1], if(add.labels) " # 1" else "")
#   coefs <- coefs[-1]                  # drop intercept $$ should only do if intercept
#   coef.width <- get.coef.width(coefs, digits)
#   for(ipred in seq_along(coefs)) {
#       coef <- coefs[ipred]
#       if(coef < 0)
#           s <- pastef(s, "  - %s ",
#                   format(-coef, justify="left", w=coef.width, digits=digits, format="%g"))
#       else
#           s <- pastef(s, "  + %s ",
#                   format(coef, justify="left", w=coef.width, digits=digits, format="%g"))
#       s <- pastef(s, "* %s", get.name(ipred))
#       if(add.labels)
#           s <- pastef(s, " # %d", ipred)
#       s <- pastef(s, "\n")
#   }
#   s
# }

#--------------------------------------------------------------------------------------------
# Method functions for generics in models.R
# All standard method funcs are supported, I think, except anova, effects, offset, family.

deviance.earth <- function(object, ...) object$rss
case.names.earth <- function(object, ...) rownames(object$residuals)

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

# Fake the AIC by returning the GCV. This is enough for step() to work.

extractAIC.earth <- function(fit, scale = 0, k = 2, ...)
{
    warning1("extractAIC.earth: returning GCV instead of AIC")
    if(scale != 0)
        warning1("extractAIC.earth: ignored scale parameter ", scale)
    if(k != 2)
        warning1("extractAIC.earth: ignored k parameter ", k)
    warn.if.dots.used("extractAIC.earth", ...)
    nterms <- length(fit$selected.terms)
    c(effective.nbr.of.params(nterms, get.nknots(nterms), fit$ppenalty), fit$gcv)
}

#--------------------------------------------------------------------------------------------
# Method functions for plotmo.R
# We look at the built model to determine singles and pairs -- this is a
# better approach for earth than that used in get.singles.default and
# get.pairs.default because it knows which predictors are unused and it
# knows which predictors were paired up during model building.

get.singles.earth <- function(object, degree1, pred.names, trace)
{
    Singles <- NULL
    dataClasses <- attr(object$terms, "dataClasses")
    if(any((dataClasses == "factor") | (dataClasses == "ordered")))
        stop1("a predictor has class 'factor'") # model.matrix cols are messed up by factors
    else {
        # selected is all degree 1 terms
        selected <- object$selected.terms[
                        reorder.earth(object, decomp="anova", degree=1, min.degree=1)]
        if(length(selected) > 0)
            Singles <- unique(
                        which(object$dirs[selected, , drop=FALSE] != 0, arr.ind=TRUE)[,2])
    }
    if(length(Singles) == 0 && is.specified(degree1))
        warning1("'degree1' specified but no degree1 plots")
    if(length(Singles) > 0) {
        check.index.vec("degree1", degree1, Singles)
        Singles <- Singles[degree1]
    }
    Singles # Singles is a vector of indices of predictors for degree1 plots
}

get.pairs.earth <- function(object, x, degree2, pred.names, trace=FALSE)
{
    Pairs <- matrix(0, nrow=0, ncol=2)      # no pairs
    selected <- object$selected.terms[      # selected is all degree 2 terms
                    reorder.earth(object, decomp="anova", degree=2, min.degree=2)]
    Pairs <- vector(mode="numeric")
    for(i in selected)                      # append indices of the two preds in term i
        Pairs <- c(Pairs, which(object$dirs[i,] != 0))
    Pairs <- unique(matrix(Pairs, ncol=2, byrow=TRUE))
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("'degree2' specified but no degree2 plots")
    if(nrow(Pairs) > 0) {
        check.index.vec("degree2", degree2, Pairs)
        Pairs <- Pairs[degree2, , drop=FALSE]
    }
    Pairs
}

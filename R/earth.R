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
# --- Eval.model.subsets
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
                     "Get.crit", "Eval.model.subsets",
                     "Force.xtx.prune")

# returns the number of arguments to the user's "allowed" function

check.allowed.arg <- function(allowed) # check earth's "allowed" argument
{
    len <- 0
    if(!is.null(allowed)) {
        allowed.func.needs <- paste0(
            "  The \"allowed\" function needs the following arguments ",
            "(but namesx and first are optional):\n      ",
            paste.with.space(c("degree", "pred", "parents", "namesx", "first")))

        if(!identical(typeof(allowed), "closure"))
            stop("your \"allowed\" argument is not a function")
        names. <- names(formals(allowed))
        len <- length(names.)
        if(len < 3 || len > 5)
            stop0("your \"allowed\" function does not have the correct number of arguments\n",
                  allowed.func.needs)

        if(names.[1] != "degree" || names.[2] != "pred" || names.[3] != "parents" ||
           (len >= 4 && names.[4] != "namesx") || (len >= 5 && names.[5] != "first")) {
              stop0(allowed.func.needs,
                "\n  You have:\n      ", paste.with.space(names.))
        }
    }
    len
}
check.weights <- function(w, wname, expected.len) # invoked for both "wp" and "weights"
{
    check.vec(w, wname, expected.len)
    check(w, wname, "negative value", function(x) { x < 0})
    meanw <- mean(w)
    if(meanw < 1e-8)
        stop0("mean of \"", wname, "\" is 0")
    # TODO finagle zero weights (but should really do what lm.wfit does and delete cols)
    almost.zero <- meanw / 1e8  # note that 1e8 becomes 1e4 after sqrt later
    w[w < almost.zero] <- almost.zero
    w
}
check.which.terms <- function(dirs, which.terms) # ensure which.terms is valid
{
    if(is.null(which.terms))
        stop0("\"which.terms\" is NULL")
    if(length(which.terms) == 0)
        stop0("length(which.terms) == 0")
    if(which.terms[1] != 1)
        stop0("first element of \"which.terms\" must be 1, the intercept term")
    if(NCOL(which.terms) > 1) {
        for(i in 1:NCOL(which.terms))
            plotmo::check.index(which.terms[,i], "which.terms", dirs,
                                allow.zeroes=TRUE, allow.dups=TRUE)
    } else
            plotmo::check.index(which.terms, "which.terms", dirs,
                                allow.zeroes=TRUE, allow.dups=TRUE)
}
convert.linpreds.to.logical <- function(linpreds, npreds, x)
{
    linpreds <- plotmo::check.index(linpreds, "linpreds", x,
                                    is.col.index=TRUE, allow.empty=TRUE)
    to.logical(linpreds, npreds)
}
earth <- function(...) UseMethod("earth")

earth.default <- function(
    x               = stop("no 'x' arg"), # NAs are not allowed in x or y, an error msg if so
    y               = stop("no 'y' arg"),
    weights         = NULL,         # case weights
    wp              = NULL,         # response column weights
    subset          = NULL,         # which rows in x to use
    na.action       = na.fail,      # only legal value is na.fail
    keepxy          = FALSE,        # true to retain x, y, etc in returned value
    trace           = 0,
    glm             = NULL,
    ncross          = 1,            # number of cross-validations, ignored unless nfold>0
    nfold           = 0,            # number of folds per cross-validation
    stratify        = TRUE,         # stratify levels in cross-validation folds
    varmod.method   = "none",       # estimate cross-validation pred intervals
    varmod.exponent = 1,            # power transform applied to fitted response
    varmod.conv     = 1,            # max mean percent coef change for IRLS iterations
    varmod.clamp    = .1,           # min.sd predicted by varmod is varmod.clamp * mean(sd)
    varmod.minspan  = -5,           # minspan for varmod call to earth
    Scale.y         = (NCOL(y)==1), # TRUE to scale y in the forward pass
    ...)                            # passed on to earth.fit
{
    trace <- check.trace.arg(trace)
    env <- parent.frame() # the environment from which earth was called
    Call <- make.call.generic(match.call(expand.dots=TRUE), "earth")
    if(trace >= 4)
        my.print.call("Call: ", Call)
    if(!is.null(Call$data))
        stop0("\"data\" argument not allowed in earth.default")
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop0("illegal \"na.action\", only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop0("illegal \"na.action\", only na.action=na.fail is allowed")
    keepxy <- check.boolean(keepxy)
    if(keepxy) {
        x.org <- x
        y.org <- y
    }
    xname <- deparse(substitute(x))
    namesx <- generate.colnames(x, is.y.arg=FALSE, xname=xname)
    namesx.org <- namesx

    # the "if" saves memory when x is already in canonical form
    if(!is.matrix(x) || !is.double(x[,1]) || !good.colnames(x)) {
        # expand factors, convert to double matrix with column names
        x <- expand.arg(x, env, FALSE, xname=xname)
        rownames(x) <- possibly.delete.rownames(x)
    }
    ylevels <- get.ylevels(y, glm)
    y <- expand.arg(y, env, is.y.arg=TRUE, deparse(substitute(y)))
    rownames(y) <- possibly.delete.rownames(y)

    rval <- earth.fit(x=x, y=y, weights=weights, wp=wp, subset=subset,
                      na.action=na.action, keepxy=keepxy, trace=trace, glm=glm,
                      Scale.y=Scale.y, ...)

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
    }
    # earth.cv will return null unless nfold > 1
    # subset parameter was already checked in earth.fit so it's safe to use it here
    cv <- earth.cv(rval,
                   if(is.null(subset)) x else x[subset,,drop=FALSE],
                   if(is.null(subset)) y else y[subset,,drop=FALSE],
                   subset, weights, na.action, keepxy, trace, glm,
                   ncross, nfold, stratify,
                   varmod.method, varmod.exponent, varmod.conv, varmod.clamp, varmod.minspan,
                   Scale.y, env, ...)
    if(!is.null(cv)) {
         rval$cv.list            <- cv$list.
         rval$cv.nterms          <- cv$nterms
         rval$cv.nvars           <- cv$nvars
         rval$cv.groups          <- cv$groups
         rval$cv.rsq.tab         <- cv$rsq.tab
         rval$cv.oof.rsq.tab     <- cv$oof.rsq.tab
         rval$cv.infold.rsq.tab  <- cv$infold.rsq.tab
         rval$cv.class.rate.tab  <- cv$class.rate.tab
         rval$cv.maxerr.tab      <- cv$maxerr.tab
         rval$cv.auc.tab         <- cv$auc.tab
         rval$cv.cor.tab         <- cv$cor.tab
         rval$cv.deviance.tab    <- cv$deviance.tab
         rval$cv.calib.int.tab   <- cv$calib.int.tab
         rval$cv.calib.slope.tab <- cv$calib.slope.tab
         rval$cv.oof.fit.tab     <- cv$oof.fit.tab # null unless varmod specified
         rval$varmod             <- cv$varmod      # null unless varmod specified
    }
    rval
}
earth.formula <- function(
    formula         = stop("no 'formula' arg"), # intercept will be ignored
    data            = NULL,
    weights         = NULL,
    wp              = NULL,
    subset          = NULL,
    na.action       = na.fail,
    keepxy          = FALSE,
    trace           = 0,
    glm             = NULL,
    ncross          = 1,
    nfold           = 0,
    stratify        = TRUE,
    varmod.method   = "none",
    varmod.exponent = 1,
    varmod.conv     = 1,
    varmod.clamp    = .1,
    varmod.minspan  = -5,
    Scale.y         = (NCOL(y)==1),
    ...)
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
    #--- earth.formula starts here
    trace <- check.trace.arg(trace)
    env <- parent.frame() # the environment from which earth was called
    Call <- make.call.generic(match.call(expand.dots=TRUE), "earth")
    if(trace >= 4)
        my.print.call("Call: ", Call)
    if(!is.null(Call[["x"]]))
        stop0("\"x\" argument not allowed in earth.formula")
    if(!is.null(Call[["y"]]))
        stop0("\"y\" argument not allowed in earth.formula")

    Call2 <- match.call(expand.dots=FALSE)

    # subset and weights handled in earth.fit, so match only on formula, data, and na.action

    m <- match(c("formula", "data", "na.action"), names(Call2), 0)
    mf <- Call2[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    if(!is.null(mf$na.action))
        stop0("\"na.action\" argument is not allowed (it is set internally to na.fail)")
    mf$na.action <- na.fail
    mf <- eval.parent(mf)

    # expand factors in x, convert to double matrix, add colnames

    x <- model.matrix(attr(mf, "terms"), mf)
    rownames(x) <- possibly.delete.rownames(x)
    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]     # silently discard intercept
    else
        warning0("ignored -1 in formula (earth objects always have an intercept)")
    # strip white space for better reading earth formula for e.g. earth(y~I(X-3))
    # because model.matrix inserts spaces around the minus
    if(!is.null(colnames(x)))
        colnames(x) <- strip.white.space(colnames(x))

    # expand factors in y, convert to double matrix, add colnames

    y <- model.response(mf, "any")  # "any" means factors are allowed
    ylevels <- get.ylevels(y, glm)
    terms. <- attr(mf, "terms")
    yname <- NULL
    if(!is.factor(y))
        yname <- names(attr(terms., "dataClasses"))[[attr(terms., "response")[1]]]
    y <- expand.arg(y, env, is.y.arg=TRUE, xname=yname)
    rownames(y) <- possibly.delete.rownames(y)

    rval <- earth.fit(x=x, y=y, weights=weights, wp=wp, subset=subset,
                      na.action=na.action, keepxy=keepxy, trace=trace, glm=glm,
                      Scale.y=Scale.y, ...)

    rval$namesx.org <- get.namesx(mf)
    rval$namesx <- make.unique(rval$namesx.org)
    rval$levels <- ylevels
    rval$terms <- terms.
    rval$call <- Call
    rval$wp <- wp

    keepxy <- check.boolean(keepxy)
    if(keepxy) {
        if(!is.null(data))
            rval$data <- data
        else # OLD: if(trace >= 1)
            warning0("No \"data\" argument to earth so \"keepxy\" is limited\n")
        rval$y <- y
        rval$subset <- subset
    }
    # earth.cv will return null unless nfold>1
    # subset parameter was already checked in earth.fit so it's safe to use it here
    cv <- earth.cv(rval,
                   if(is.null(subset)) x else x[subset,,drop=FALSE],
                   if(is.null(subset)) y else y[subset,,drop=FALSE],
                   subset, weights, na.action, keepxy, trace, glm,
                   ncross, nfold, stratify,
                   varmod.method, varmod.exponent, varmod.conv, varmod.clamp, varmod.minspan,
                   Scale.y, env, ...)
    if(!is.null(cv)) {
         rval$cv.list            <- cv$list.
         rval$cv.nterms          <- cv$nterms
         rval$cv.nvars           <- cv$nvars
         rval$cv.groups          <- cv$groups
         rval$cv.rsq.tab         <- cv$rsq.tab
         rval$cv.oof.rsq.tab     <- cv$oof.rsq.tab
         rval$cv.infold.rsq.tab  <- cv$infold.rsq.tab
         rval$cv.class.rate.tab  <- cv$class.rate.tab
         rval$cv.maxerr.tab      <- cv$maxerr.tab
         rval$cv.auc.tab         <- cv$auc.tab
         rval$cv.cor.tab         <- cv$cor.tab
         rval$cv.deviance.tab    <- cv$deviance.tab
         rval$cv.calib.int.tab   <- cv$calib.int.tab
         rval$cv.calib.slope.tab <- cv$calib.slope.tab
         rval$cv.oof.fit.tab     <- cv$oof.fit.tab # null unless varmod specified
         rval$varmod             <- cv$varmod      # null unless varmod specified
    }
    rval
}
# This is called from earth.default or earth.formula, not directly
# because the x and y args must be expanded for factors first.

earth.fit <- function(
    x       = stop("no 'x' arg"), # x and y already processed by model.matrix
    y       = stop("no 'y' arg"), # NAs are not allowed in x or y, an error msg if so
    weights = NULL,         # case weights (row weights)
    wp      = NULL,         # response weights (column weights)
    subset  = NULL,         # which rows in x to use
    na.action = na.fail,    # only legal value is na.fail
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

    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep"),
    nprune  = NULL,         # max nbr of terms (including intercept) in prune subset

                            # Following will not usually be used by the user

    Object  = NULL,         # if null, recreate earth model from scratch with forward pass
                            # if not null: no forward pass, just pruning pass

    Get.crit           = get.gcv,            # criterion func for model select during pruning
    Eval.model.subsets = eval.model.subsets, # function used to evaluate model subsets
    Scale.y            = (NCOL(y)==1), # TRUE to scale y in the forward pass
    Force.xtx.prune    = FALSE, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    Force.weights      = FALSE, # TRUE to force use of weight code in earth.c
    Use.beta.cache     = TRUE,  # TRUE to use beta cache, for speed
    Exhaustive.tol     = 1e-10,
    ...)                        # unused
{
    check.no.family.arg.to.earth(...)
    stop.if.dots.used("earth.fit", ...)
    if(is.logical(trace) && trace) {
        warning0("earth: converted trace=TRUE to trace=4")
        trace <- 4
    }
    Force.xtx.prune <- check.boolean(Force.xtx.prune)
    Use.beta.cache  <- check.boolean(Use.beta.cache)
    env <- parent.frame()
    y.org <- y
    pmethod <- match.arg1(pmethod) # check pmethod is legal, and expand
    # For binomial glms we have to drop paired y cols before passing y to
    # to the C earth routines.  The logical vector glm.bpairs keeps track
    # of which cols are used.
    glm.bpairs <- NULL                       # NULL means all columns used
    if(!is.null(glm)) {
        glm <- get.glm.arg(glm)
        glm.bpairs <- get.glm.bpairs(y, glm) # returns NULL if all cols used
    }
    if(trace >= 1) {
        if(trace >= 4)
            cat("\n")
        print.matrix.info("x", x, NULL, NULL,       details=(trace >= 4),
                          all.names=(trace >= 2), all.rows=(trace >= 5))
        if(trace >= 4)
            cat("\n")
        print.matrix.info("y", y, NULL, glm.bpairs, details=(trace >= 4),
                          all.names=(trace >= 2), all.rows=(trace >= 5))
        if(!is.null(weights)) {
            if(trace >= 4)
                cat("\n")
            print.matrix.info("weights", weights, NULL, NULL, details=(trace >= 4),
                              all.names=(trace >= 2), all.rows=(trace >= 5))
        }
        if(trace >= 3) # sic, gives good spacing later
            cat("\n")
    }
    # we do basic parameter checking here but much more in ForwardPass in earth.c
    check.integer.scalar(nk, min=1)
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop0("illegal \"na.action\", only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop0("illegal \"na.action\", only na.action=na.fail is allowed")
    na.action <- na.fail
    if(trace != -1 &&
       !is.null(weights) &&
       !isTRUE(all.equal(weights, rep(weights[1], length(weights))))) {
        warning0("support of weights is provisional in this version of earth")
    }
    if(trace == -1)
        trace <- 0
    n.allowedfunc.args <- check.allowed.arg(allowed)
    stopif(is.vector(x)) # should have been converted to matrix in earth.default or earth.formula
    stopif(is.vector(y)) # ditto
    if(nrow(x) == 0)
        stop("\"x\" has no rows")
    if(ncol(x) == 0)    # this happens for example for earth(Volume~Volume,data=trees)
        stop("\"x\" has no columns")
    if(nrow(x) != nrow(y))
        stop("nrow(x) ", nrow(x), " != nrow(y) ", nrow(y))
    # x and y must be double for calls to C functions
    if(!all(is.double(x)))
        stop("non double entries in \"x\" argument")
    if(!all(is.double(y)))
        stop("non double entries in \"y\" argument")
    weights.or.null <- weights
    if(is.null(weights))
        weights <- repl(1, nrow(x))
    weights <- check.weights(weights, "weights", nrow(x))
    if(!is.double(weights)) # weights must be double for calls to C functions
        weights <- as.double(weights)
    if(!is.null(subset)) {
        # duplicates are allowed in subsets so user can specify a bootstrap sample
        subset <- plotmo::check.index(subset, "subset", x, allow.dups=TRUE, allow.zeroes=TRUE)
        x <- x[subset, , drop=FALSE]
        y <- y[subset, , drop=FALSE]
        weights <- weights[subset]
    }
    yw <- NULL # weighted version of y
    if(Force.weights ||
       (!is.null(weights.or.null) && !all(abs(weights - weights[1]) < 1e-8))) {
        yw <- sqrt(weights) * y
    }
    if(!is.null(wp)) { # column weights
        wp <- check.weights(wp, "wp", ncol(y))
        wp <- sqrt(wp / mean(wp)) # normalize
        # multiply each column of y by its normalized weight
        y <- y  * outer(repl(1, nrow(y)), wp)
        if(!is.null(yw))
            yw <- yw * outer(repl(1, nrow(y)), wp)
    }
    if(!is.null(glm.bpairs)) {
        y <- y [, glm.bpairs, drop=FALSE]
        if(!is.null(yw))
            yw <- yw[, glm.bpairs, drop=FALSE]
    }
    if(is.null(Object)) {
        rval <- forward.pass(x, y, yw, weights,
                             trace, degree, penalty, nk, thresh,
                             minspan, endspan, newvar.penalty, fast.k, fast.beta,
                             linpreds, allowed, Use.beta.cache, Scale.y, Force.weights,
                             n.allowedfunc.args, env)
          bx     <- rval$bx
          dirs   <- rval$dirs
          cuts   <- rval$cuts
          reason <- rval$reason
    } else {
        # no forward pass: get here if update() called me with no forward pass params
        if(trace >= 1)
            cat("Skipped forward pass\n")
        check.classname(Object, deparse(substitute(Object)), "earth")
        bx      <- NULL
        dirs    <- Object$dirs
        cuts    <- Object$cuts
        reason  <- Object$reason
    }
    rval <- pruning.pass(if(is.null(yw)) x else sqrt(weights) * x,
                         if(is.null(yw)) y else yw,
                         subset,
                         trace, penalty, pmethod, nprune, bx, dirs, cuts, Get.crit,
                         Eval.model.subsets, Force.xtx.prune, Exhaustive.tol)
      bx             <- rval$bx.prune
      rss.per.subset <- rval$rss.per.subset # vector of RSSs for each model (index on subset size)
      gcv.per.subset <- rval$gcv.per.subset # vector of GCVs for each model
      prune.terms    <- rval$prune.terms    # triang mat: each row is a vector of term indices
      selected.terms <- rval$selected.terms # vec of term indices of selected model

    bx <- bx[, selected.terms, drop=FALSE]
    bx <- bx / sqrt(weights) # unweight bx

    # add names for returned values

    pred.names <- colnames(x)
    term.names <- get.earth.term.name(1:nrow(dirs), dirs, cuts, pred.names, x)
    colnames(bx) <- term.names[selected.terms]
    dimnames(dirs) <- list(term.names, pred.names)
    dimnames(cuts) <- list(term.names, pred.names)

    nresp <- ncol(y)  # number of responses
    nselected <- length(selected.terms)

    # Regress y on bx to get fitted.values etc.
    # The as.matrix calls after the call to lm are needed if y is
    # a vector so the fitted.values etc. are always arrays.

    if(!is.null(yw))
        lfit <- lm.wfit(bx, y, w=weights, singular.ok=FALSE)
    else
        lfit <- lm.fit(bx, y, singular.ok=FALSE)

    fitted.values <- as.matrix(lfit$fitted.values)
    residuals     <- as.matrix(lfit$residuals)
    coefficients  <- as.matrix(lfit$coefficients)
    response.names <- colnames(y)
    colnames(fitted.values) <- response.names
    colnames(residuals)     <- response.names
    colnames(coefficients)  <- response.names

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
    if(!is.null(glm)) {
        y.glm <- y.org
        if(!is.null(subset))
            y.glm <- y.glm[subset, , drop=FALSE]
        glm.list <- earth.glm(bx, y.glm, weights, na.action, glm,
                              trace, glm.bpairs, response.names[1], env)
        glm.coefs <- get.glm.coefs(glm.list, ncol(coefficients),
                                   selected.terms, term.names, response.names)
    }
    # prepare returned summary statistics

    rss <- rss.per.subset[nselected]
    rsq  <- get.rsq(rss, rss.per.subset[1])
    gcv  <- gcv.per.subset[nselected]
    grsq <- get.rsq(gcv, gcv.per.subset[1])
    rss.per.response  <- vector(mode="numeric", length=nresp)
    rsq.per.response  <- vector(mode="numeric", length=nresp)
    gcv.per.response  <- vector(mode="numeric", length=nresp)
    grsq.per.response <- vector(mode="numeric", length=nresp)

    for(iresp in 1:nresp) {
        rss.per.response[iresp]  <- ss(residuals[,iresp])
        tss                      <- ss(y[,iresp] - mean(y[,iresp]))
        rsq.per.response[iresp]  <- get.rsq(rss.per.response[iresp], tss)
        gcv.null                 <- Get.crit(tss, 1, penalty, nrow(x))
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

        rss.per.subset = rss.per.subset,# nprune x 1, RSS of each model, across all resp
        gcv.per.subset = gcv.per.subset,# nprune x 1, GCV of each model, across all resp

        fitted.values  = fitted.values, # ncases (after subset) x nresp
        residuals      = residuals,     # ncases (after subset) x nresp
        coefficients   = coefficients,  # selected terms only: nselected x nresp

        penalty        = penalty,          # copy of penalty argument
        nk             = nk,               # copy of nk argument
        thresh         = thresh,           # copy of thresh argument
        reason         = reason,           # reason we terminated the forward pass
        weights        = weights.or.null), # copy of weights argument, may be NULL
    class = "earth")

    if(!is.null(glm.list)) {
        rval$glm.list         <- glm.list   # list of glm models, NULL if none
        rval$glm.coefficients <- glm.coefs  # matrix of glm coefs, nselected x nresp
        rval$glm.bpairs       <- glm.bpairs # null unless paired binomial cols
    }
    rval
}
effective.nbr.of.params <- function(ntermsVec, nknotsVec, penalty)  # for GCV calculation
{
    if(penalty < 0) # special case: term and knots are free so GCV == RSS/ncases
        repl(0, length(ntermsVec))
    else
        ntermsVec + (penalty * nknotsVec)
}
# this returns the RSS and selected terms for each subset of size 1:nprune

eval.model.subsets <- function(     # this is the default Eval.model.subsets
    bx,      # weighted basis matrix
    y,       # weighted model response
    pmethod,
    nprune,  # max nbr of terms (including intercept) in prune subset, in range 1..nterms
    Force.xtx.prune, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    trace)
{
    stopifnot(nprune >= 1 && nprune <= nrow(bx))

    if(Force.xtx.prune)         # user explicitly asked for xtx subset evaluation
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune)

    else if(ncol(y) > 1)  {     # leaps cannot deal with multiple responses
        if(pmethod != "none" && pmethod != "backward")
            stop0("pmethod == \"", pmethod,
                  "\" is not allowed with multiple response models\n",
                  "       (y has ", ncol(y), " columns, use trace=4 to see y)")
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune)

    } else if(ncol(bx) <= 2) {  # leaps code gives an error for small number of cols
        if(pmethod != "none" && pmethod != "backward")
            pmethod <- "backward"
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune)

    } else
        eval.model.subsets.with.leaps(bx, y, pmethod, nprune, trace)
}
eval.model.subsets.with.leaps <- function(
    bx,
    y,
    pmethod,
    nprune,
    trace)
{
    convert.lopt <- function(lopt, nprune) # convert lopt format to prune.terms format
    {
        # Assignment fills matrix column wise. We want row wise, so
        # take upper triangle and then transpose.

        prune.terms <- matrix(0, nrow=nprune, ncol=nprune)
        prune.terms[upper.tri(prune.terms, diag=TRUE)] <- lopt
        t(prune.terms)
    }
    #--- eval.model.subsets.with.leaps starts here ---
    # Make warnings in the leaps routines behave as errors.
    # We have seen the leaps routines return bad data (?) after issuing
    # the warning "XHAUST returned error code -999",
    # and when that occurs we don't want to continue running.
    old.warn <- getOption("warn")
    on.exit(options(warn=old.warn))
    options(warn=2)

    rprune <- leaps.setup(x=bx, y=y,
        force.in=1,        # make sure intercept is in model
        force.out=NULL,
        intercept=FALSE,   # we have an intercept so leaps.setup must not add one
        nvmax=nprune, nbest=1, warn.dep=TRUE)

    rprune <- switch(pmethod,
        backward   = leaps.backward(rprune),
        none       = leaps.backward(rprune), # for stats, won't actually prune
        exhaustive = leaps.exhaustive(rprune, really.big=TRUE),
        forward    = leaps.forward(rprune),
        seqrep     = leaps.seqrep(rprune))

    rss.per.subset <- as.vector(rprune$ress) # convert from n x 1 mat to vector

    list(rss.per.subset = rss.per.subset, # vec of RSSs for each model (index on subset size)
         prune.terms = convert.lopt(rprune$lopt, nprune)) # each row is a vec of term indices
}
# This calls the earth.c routine EvalSubsetsUsingXtxR.
# Unlike the leaps code, it can deal with multiple responses (i.e. multiple y columns)

eval.subsets.xtx <- function(
    bx,
    y,
    pmethod,
    nprune,
    Force.xtx.prune)
{
    bad.pmethod <- function()
    {
        stop0("pmethod=\"", pmethod,
              "\" is not allowed with \"eval.subsets.xtx\"")
    }
    backward <- function(bx, y)
    {
        ncases <- nrow(bx)
        nterms <- ncol(bx)
        nresp <- ncol(y)
        stopifnot(is.double(bx))
        stopifnot(is.double(y))
        # TODO replace .C call with alternative interface that doesn't require DUP=TRUE
        rval <- .C("EvalSubsetsUsingXtxR",
            prune.terms = matrix(0, nrow=nterms, ncol=nterms), # double PruneTerms[]
            rss.per.subset = vector(mode="numeric", length=nterms),
            as.integer(ncases),     # const int *pnCases
            as.integer(nresp),      # const int *pnResp
            as.integer(nterms),     # const int *pnMaxTerms
            bx,                     # const double bx[]
            y,                      # const double y[]
            PACKAGE="earth")

        # above returns all subsets, so trim back to nprune below

        list(rss.per.subset = rval$rss.per.subset[1:nprune],
             prune.terms    = rval$prune.terms[1:nprune, 1:nprune, drop=FALSE])
    }
    #--- eval.subsets.xtx starts here ---

    rprune <- switch(pmethod,
        backward   = backward(bx, y),
        none       = backward(bx, y), # for stats, won't actually prune
        exhaustive = bad.pmethod(),
        forward    = bad.pmethod(),
        seqrep     = bad.pmethod())
}
do.scale.y <- function(y, yname)
{
    y.scaled <- scale(y)
    # make sure that scaling was ok
    i <- which(attr(y.scaled,"scaled:scale") == 0)
    if(length(i)) {
        if(ncol(y) > 1)
            stop0("cannot scale column ", i[1],
                  " of ", yname,
                  " (values are all equal to ", y[1,i], ")")
        else
            stop0("cannot scale ", yname,
                  " (values are all equal to ", y[1,1], ")\n",
                  "Try Scale.y=FALSE?")
    }
    y
}
forward.pass <- function(x, y, yw, weights, # must be double, but yw can be NULL
                         trace, degree, penalty, nk, thresh,
                         minspan, endspan, newvar.penalty, fast.k, fast.beta,
                         linpreds, allowed, Use.beta.cache, Scale.y, Force.weights,
                         n.allowedfunc.args, env)
{
    stopifnot(nrow(x) == nrow(y))
    if(!is.null(yw)) {
        stopifnot(nrow(y) == nrow(yw))
        stopifnot(ncol(y) == ncol(yw))
    }
    npreds <- ncol(x)
    fullset <- repl(0, nk) # element will be set TRUE if corresponding term used
    linpreds <- convert.linpreds.to.logical(linpreds, npreds, x)
    if(trace >= 2)
        print.linpreds(linpreds, x)

    if(Scale.y) {
        # scale y for better stability in the forward pass
        y.scaled <- do.scale.y(y, "y")
        if(!is.null(yw))
            yw.scaled <- do.scale.y(yw, "weighted y")
    } else {
        y.scaled <- y
        if(!is.null(yw))
            yw.scaled <- yw
    }
    # this calculation of nbytes if not accurate, it doesn't matter
    nbytes <- 8 * (nk^2 * ncol(x) + (nrow(x) * (3 + 2*nk + ncol(x)/2)))
    if(nbytes > 1e7)  # need more than 10 MBytes to build model?
        gc()

    # Mar 2012: R cmd check now complains if you pass R NULL via .C, so
    # we gyp this by passing a special value in my.null.
    # TODO is this reliable and safe?

    my.null = 999999L

    if(is.null(allowed))
        allowed <- my.null

    Force.weights <- check.boolean(Force.weights)

    pred.names <- if(trace >= 2 || n.allowedfunc.args >= 4) colnames(x) else my.null

    reason <- integer(length=1) # reason we terminated the forward pass

    on.exit(.C("FreeR",PACKAGE="earth")) # frees memory if user interrupts

    rval <- .C("ForwardPassR",
        fullset = as.integer(fullset),           # out: int FullSet[]
        bx   = matrix(0, nrow=nrow(x), ncol=nk), # out: double bx[]
        dirs = matrix(0, nrow=nk, ncol=npreds),  # out: double Dirs[]
        cuts = matrix(0, nrow=nk, ncol=npreds),  # out: double Cuts[]
        reason = reason,                         # out: int*
        x,                                       # in: const double x[]
        y.scaled,                                # in: const double y[]
        if(is.null(yw)) my.null else yw.scaled,  # in: const double yw[]
        weights,                        # in: const double WeightsArg[]
        as.integer(nrow(x)),            # in: const int *pnCases
        as.integer(ncol(y)),            # in: const int *pnResp
        as.integer(npreds),             # in: const int *pnPreds
        as.integer(degree),             # in: const int *pnMaxDegree
        as.double(penalty),             # in: const double *pPenalty
        as.integer(nk),                 # in: const int *pnMaxTerms
        as.double(thresh),              # in: double *pThresh
        as.integer(minspan),            # in: const int *pnMinSpan
        as.integer(endspan),            # in: const int *pnEndSpan
        as.integer(fast.k),             # in: const int *pnFastK
        as.double(fast.beta),           # in: const double *pFastBeta
        as.double(newvar.penalty),      # in: const double *pNewVarPenalty
        as.integer(linpreds),           # in: const int LinPreds[]
        allowed,                        # in: const SEXP Allowed
        as.integer(n.allowedfunc.args), # in: const int *pnAllowedFuncArgs
        env,                            # in: const SEXP Env for Allowed
        as.integer(Use.beta.cache),     # in: const int *pnUseBetaCache
        as.double(trace),               # in: const double *pTrace
        pred.names,                     # in: const char *sPredNames[], can be MyNull
        my.null,                        # in: const SEXP MyNull
        as.integer(Force.weights),      # in
        NAOK = TRUE, # we check for NAs etc. internally in C ForwardPass
        PACKAGE="earth")

    fullset <- as.logical(rval$fullset)

    list(reason = rval$reason,
         bx     = rval$bx[, fullset, drop=FALSE],
         dirs   = rval$dirs[fullset, , drop=FALSE],
         cuts   = rval$cuts[fullset, , drop=FALSE])
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
                paste0("x[,", ipred, "]")
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
                    s <- paste0(s, "*")
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
    #--- --- get.earth.term.name starts here --- ---
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
        warning0("duplicate term name \"", term.names[which(duplicated)[1]], "\"\n",
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

    # Removed in earth version 3.1-2:
    #     if(max(nparams, na.rm=TRUE) >= ncases)
    #         warning0("effective number ", max(nparams, na.rm=TRUE),
    #                  " of GCV parameters >= number ", ncases, " of cases")

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
get.ylevels <- function(y, glm)
{
    if(!is.null(levels(y)))
        return(levels(y))
    # following needed for predict.earth(type="class")
    if(is.logical(y))
        return(c(FALSE, TRUE))
    if(is.numeric(y)) {
        range <- range(y)
        if(range[2] - range[1] == 1)
            return(c(range[1], range[2]))
    }
    NULL
}
# Remove useless(?) "1" "2" "3" ... rownames for x (added by
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
print.pruning.pass <- function(trace, pmethod, penalty, nprune, selected.terms,
                               prune.terms, rss.per.subset, gcv.per.subset, dirs)
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
            cat0(if(iterm==nselected) "chosen " else "       ",
                 format(iterm, width=4),
                 sprintf("%12.4f ", grsq),
                 sprintf("%7.4f",   get.rsq(rss.per.subset[iterm], rss.per.subset[1])),
                 sprintf("%11.4f ", delta.grsq),
                 sprintf("%6d",     get.nused.preds.per.subset(dirs, selected)),
                 "  ")
            if(trace >= 4)
                cat(selected)
        cat("\n")
        }
    cat("\n")
    }
    if(trace >= 1) {
        cat0("Prune method \"", pmethod[1], "\" penalty ", penalty,
             " nprune ", nprune, ": selected ", nselected, " of ")
        selected <- prune.terms[nselected,]
        selected <- selected[selected != 0]
        cat(nrow(dirs), "terms, and",
            get.nused.preds.per.subset(dirs, selected),
            "of", ncol(dirs), "predictors\n")
        cat("After backward pass GRSq",
            format(get.rsq(gcv.per.subset[nselected], gcv.per.subset[1]), digits=3),
            "RSq",
            format(get.rsq(rss.per.subset[nselected], rss.per.subset[1]), digits=3),
            "\n")
    }
}
pruning.pass <- function(x, y, # x and y may be weighted
    subset, trace, penalty, pmethod, nprune, bx, dirs, cuts, Get.crit,
    Eval.model.subsets, Force.xtx.prune, Exhaustive.tol)
{
    possibly.print.pacifier <- function(bx, nprune) # print a reminder if pruning will be slow
    {
        if(pmethod == "exhaustive" && nprune > 1) {
            nsubsets <- 0 # approx, assumes brute force exhaustive search
            for(subset.size in 1:nprune)
                nsubsets <- nsubsets + choose(ncol(bx), subset.size)
            if(trace >= 1 || nsubsets > 1e9) {
                cat0("Exhaustive pruning: number of subsets ",
                     format(nsubsets, digits=2),
                     if(trace >= 1) " " else "\n")
            }
        }
    }
    # If pmethod is exhaustive and bx is ill conditioned, change pmethod to
    # backward. This prevents leaps.exhaustive returning error code -999.
    #
    # Note that bx should never be ill-conditioned (RegressAndFix should
    # take care of that).  However it seems that dqrdc2 (called by
    # RegressAndFix) does not detect certain types of ill conditioning (with
    # any tol). This is probably because we are near the numerical noise floor
    # and the column norms in dqrdc are not monotically decreasing.
    # This change was made in Apr 2011.
    #
    # TODO This would be better handled by simply removing collinear cols in bx?

    possibly.change.exhaustive.to.backward <- function(pmethod, bx)
    {
        if(pmethod == "exhaustive") {
            if(!is.numeric(Exhaustive.tol) || length(Exhaustive.tol) != 1) {
                cat("\n")
                stop0("illegal Exhaustive.tol, try something like Exhaustive.tol=1e-8")
            }
            if(Exhaustive.tol < 0 || Exhaustive.tol > .1) {
                cat("\n")
                stop0("illegal Exhaustive.tol ", Exhaustive.tol,
                      ", try something like Exhaustive.tol=1e-8")
            }
            sing.vals <- svd(bx)$d  # expensive
            cond <- sing.vals[length(sing.vals)] / sing.vals[1]
            if(is.na(cond) || cond < Exhaustive.tol) {
                if(trace >= 1)
                    cat("\n")
                warning0("forced pmethod=\"backward\" ",
                    "(bx is ill conditioned, sing val ratio ",
                    format(cond, digits=2), ")")
                pmethod <- "backward"
            } else if(trace >= 1)
                cat0("(bx sing val ratio ", format(cond, digits=2), ")\n")
        }
        pmethod
    }
    get.nprune <- function()  # convert user's nprune argument to a valid value
    {
        nterms <- nrow(dirs)
        if(is.null(nprune))
            nprune <- nterms
        if(nprune > nterms) {
            # Commented out Apr 2011
            # warning0("specified \"nprune\" ", nprune,
            #         " is greater than the number ", nterms, " of available model terms ",
            #         "\nForcing \"nprune\" to ", nterms)
            nprune <- nterms
        }
        if(nprune < 1)
            stop0("\"nprune\" is less than 1")
        nprune
    }
    #--- pruning.pass starts here ---
    if(!is.null(bx) && is.null(subset))
        bx.prune <- bx   # use bx created in forward pass
    else                 # else recreate bx with all terms (update invoked earth)
        bx.prune <- get.bx(x, 1:nrow(dirs), dirs, cuts) # all terms
    stopifnot(nrow(bx.prune) == nrow(y))
    nprune <- get.nprune()
    possibly.print.pacifier(bx.prune, nprune)
    pmethod <- possibly.change.exhaustive.to.backward(pmethod, bx.prune)
    flush.console() # make sure previous messages get seen, pruning make take a while
    rval <- Eval.model.subsets(bx.prune, y,
                               pmethod, nprune, Force.xtx.prune, trace)
      rss.per.subset <- rval$rss.per.subset # RSS for each subset (across all responses)
      prune.terms    <- rval$prune.terms    # each row is a vec of term indices

    # sanity checks here because Eval.model.subsets may be user written
    # rss.per.subset will be less than nprune if lin dep discovered by leaps
    stopifnot(length(rss.per.subset) <= nprune)
    nprune <- length(rss.per.subset)
    stopifnot(NROW(prune.terms) == nprune)
    stopifnot(NCOL(prune.terms) == nprune)
    prune.terms <- prune.terms[1:nprune, 1:nprune, drop=FALSE]
    stopif(any(prune.terms[,1] != 1))   # check intercept column
    gcv.per.subset <- Get.crit(rss.per.subset, 1:nprune, penalty, nrow(bx.prune))
    if(!all(is.finite(rss.per.subset)))
        warning0("earth: non finite RSS in model subsets ",
                 "(see the rss.per.subset returned by earth)")
    check.vec(rss.per.subset, "rss.per.subset")
    # else if(!all(is.finite(gcv.per.subset)))
    #     warning0("earth: non finite GCV in model subsets ",
    #              "(see the gcv.per.subset returned by earth)")
    selected.terms <- 1:nprune  # all terms
    if(pmethod != "none") {
        # choose the subset which has the lowest GCV in the vector of GCVS
        selected.terms <- prune.terms[which.min(gcv.per.subset),]
        selected.terms <- selected.terms[selected.terms != 0]
    }
    print.pruning.pass(trace, pmethod, penalty, nprune, selected.terms,
                       prune.terms, rss.per.subset, gcv.per.subset, dirs)

    list(bx.prune      = bx.prune,
        rss.per.subset = rss.per.subset, # vector of RSSs for each model (index on subset size)
        gcv.per.subset = gcv.per.subset, # vector of GCVs for each model (index on subset size)
        prune.terms    = prune.terms,    # triang mat: each row is a vector of term indices
        selected.terms = selected.terms) # vec of model terms in best model
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
    key.pair <- double(nterms)                  # put h(5-x1) before h(x1-5)
    if(nterms > 1)
        for(i in 2:nterms) {
            key.linpreds[i] <- -sum(dirs[i, ] == 2)
            key.cuts[i] <- cuts[i, first.fac.order[i]]
            ifirst.non.zero <- which(dirs[i, ] != 0)[1]
            stopifnot(length(ifirst.non.zero) == 1)
            key.pair[i] <- dirs[i, ifirst.non.zero] == 1
    }
    order(key.degrees, key.linpreds, key.x, key.cuts, key.pair)
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
        stop0("degree ", degree, " < 0")
    if(min.degree < 0)
        stop0("min.degree ", min.degree, " < 0")
    if(degree < min.degree)
        stop0("degree ", degree, " < min.degree ", min.degree)
    check.which.terms(x$dirs, which.terms)
    dirs <- x$dirs[which.terms, , drop=FALSE]
    new.order <- switch(match.arg1(decomp),
                   anova = reorder.terms.anova(
                                dirs, x$cuts[which.terms,,drop=FALSE]),
                   none  = 1:length(which.terms))
    degrees <- get.degrees.per.term(dirs[new.order, , drop=FALSE])
    new.order[degrees >= min.degree & degrees <= degree]
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
#    This default decision to do a forward pass or not can be overridden
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
    #--- update.earth starts here ---

    check.classname(object, deparse(substitute(object)), "earth")
    Call <- object$call
    stopif(is.null(Call))
    do.forward.pass <- FALSE
    if(!is.null(formula.)) {
        if(is.null(Call$formula))
            stop0("\"formula.\" argument is not allowed on ",
                  "objects created without a formula")
        Call$formula <- update.formula(formula(object), formula.)
        do.forward.pass <- TRUE
    }
    # figure out what trace should be

    this.call <- match.call(expand.dots=TRUE)
    trace <- get.update.arg(this.call$trace, "trace", object,
                            trace1=NULL, "update.earth", print.trace=FALSE)
    trace <- eval.parent(trace)
    if(is.name(trace))      # TODO needed when called from earth.cv with glm=NULL, why?
        trace <- eval.parent(trace)
    if(is.null(trace))
        trace <- 0
    if(is.name(Call$glm))   # TODO needed when called from earth.cv with glm=NULL, why?
        Call$glm <- eval.parent(Call$glm)
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
    if(check.boolean(ponly))
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
        warning0("deviance.earth: returning earth (not GLM) deviance")
    object$rss
}
effects.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning0("effects.earth: returning NULL")
    NULL
}
family.earth <- function(object, ...)
{
    family(object$glm.list[[1]])
}
anova.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning0("anova.earth: returning NULL")
    NULL
}
# use.names can have the following values:
#   TRUE:  return name if possible, else return x[,i] or x[i-1].
#   FALSE: return x[,i]
#   -1:    return x[i] with 0 based indexing (treat x as a C array)

variable.names.earth <- function(object, ..., use.names=TRUE)
{
    warn.if.dots.used("variable.names.earth", ...)
    ipred <- 1:ncol(object$dirs)
    if(length(use.names) != 1)
        stop0("illegal value for use.names")
    if(use.names == TRUE) {
        varname <- colnames(object$dirs)[ipred]
        if(!is.null(varname) && !is.na(varname))
            varname
        else
            paste0("x[,", ipred, "]")
    } else if(use.names == FALSE)
        paste0("x[,", ipred, "]")
    else if(use.names == -1)
        paste0("x[", ipred-1, "]")
    else
        stop0("illegal value for use.names \"", use.names, "\"")
}
case.names.earth <- function(object, ...)
{
    if(is.null(row.names(object$residuals)))
        paste(1:nrow(object$residuals))
    else
        row.names(object$residuals)
}
residuals.earth <- function(object = stop("no 'object' arg"), type = NULL, warn=TRUE, ...)
{
    glm.resids <- function(object, type)
    {
        g <- object$glm.list
        if(is.null(g))
            stop0("residuals.earth: type \"", type, "\" can be used ",
                  "only on earth-glm models")
        colnames <- ""
        for(imodel in seq_along(g)) {
            rval1 <- residuals(g[[imodel]], type)
            if(imodel == 1)
                rval <- rval1
            if(NROW(rval1) != NROW(rval)) # should never happen
                stop0("residuals.earth: glm.list[[", imodel, "]] does ",
                      "not conform to glm.list[[", 1, "]] ",
                      "(residuals have a different length)")
            if(imodel > 1) {
                colnames <- c(colnames)
                rval <- cbind(rval, rval1)
            }
        }
        rval
    }
    #--- residuals.earth starts here ---

    warn.if.dots.used("residuals.earth", ...)
    if(warn && is.null(type) && !is.null(object$glm.list))
        warning0("residuals.earth: returning earth (not glm) residuals")
    if(is.null(type))
        type <- "earth"
    types <- c("earth", "pearson", "deviance",
               "glm.pearson", "glm.working", "glm.response", "glm.partial")
    rval <- switch(match.choices(type, types),
        earth    = object$residuals,
        deviance = if(is.null(object$glm.list))
                       object$residuals
                   else
                       glm.resids(object, "deviance"),
        pearson      = object$residuals / get.se(object),
        glm.pearson  = glm.resids(object, "pearson"),
        glm.working  = glm.resids(object, "working"),
        glm.response = glm.resids(object, "response"),
        glm.partial  = glm.resids(object, "partial"))

    if(!is.matrix(rval))
        rval <- matrix(rval, ncol = 1)
    if(type != "glm.partial")
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
        warning0("extractAIC.earth: returning GCV instead of AIC")
    if(scale != 0)
        warning0("extractAIC.earth: ignored scale parameter ", scale)
    if(k != 2)
        warning0("extractAIC.earth: ignored k parameter ", k)
    warn.if.dots.used("extractAIC.earth", ...)
    nterms <- length(fit$selected.terms)
    c(effective.nbr.of.params(nterms, get.nknots(nterms), fit$penalty), fit$gcv)
}
hatvalues.earth <- function(model, ...)
{
    stop.if.dots.used("hatvalues.earth", ...)
    if(is.null(model$y))
        stop0("hatvalues.earth: earth object has no $y field.\n",
              "       Use keepxy=TRUE in the call to earth.")
    stopifnot(!is.null(ncol(model$y)))
    stopifnot(!is.null(model$bx))
    if(ncol(model$y) > 1)
        warning0("multiple response earth model: \n",
                 "getting hat values for only the first response")
    hatvalues(lm(model$y[,1] ~ model$bx))
}
coef.earth <- function(object, decomp="none", ...)
{
    warn.if.dots.used("coef.earth", ...)
    coef <- object$coefficients
    if(NCOL(coef) > 1)
        stop0("coef.earth: multiple response models not supported")
    new.order <- reorder.earth(object, decomp=decomp)
    names <- spaceout(rownames(coef))
    coef <- coef[new.order,]
    names(coef) <- names
    coef
}
weights.earth <- function(object, ...)
{
    warn.if.dots.used("weights.earth", ...)
    if(is.null(object$weights)) # weights arg to earth was NULL?
        repl(1, length(object$fitted.values))
    else
        object$weights
}

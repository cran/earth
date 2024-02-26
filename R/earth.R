# earth.R: an implementation of Friedman's Multivariate Adaptive
#          Regression Splines, commonly known as MARS.
#
# This code is derived from code in mda.R by Hastie and Tibshirani.
# Functions are in alphabetical order after earth.default and earth.formula.
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
#-----------------------------------------------------------------------------

# This is a list of those formal arguments of earth.fit that can be changed without
# requiring a new forward pass.
# NOTE: if you change the pruning formal arguments in earth.fit(), update this too!

prune.only.args <- c("glm", "trace", "nprune", "pmethod",
    "Use.beta.cache", "Force.xtx.prune", "Get.leverages", "Exhaustive.tol")

earth <- function(...) { UseMethod("earth") }

earth.default <- function( # user called earth with x y args (not formula)
    x               = stop("no 'x' argument"),
    y               = stop("no 'y' argument"),
    weights         = NULL,         # case weights
    wp              = NULL,         # response column weights
    subset          = NULL,         # which rows in x to use
    na.action       = na.fail,      # only legal value is na.fail
    pmethod         = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
    keepxy          = FALSE,        # true to retain x, y, etc in returned value
    trace           = 0,
    glm             = NULL,
    degree          = 1, # max degree of interaction (1=additive model) (Friedman's mi)
    nprune          = NULL,         # max nbr of terms (including intercept) in pruned subset
    nfold           = 0,            # number of folds per cross-validation
    ncross          = 1,            # number of cross-validations, ignored unless nfold>0
    stratify        = TRUE,         # stratify levels in cross-validation folds
    varmod.method   = "none",       # estimate cross-validation pred intervals
    varmod.exponent = 1,            # power transform applied to fitted response
    varmod.conv     = 1,            # max mean percent coef change for IRLS iterations
    varmod.clamp    = .1,           # min.sd predicted by varmod is varmod.clamp * mean(sd)
    varmod.minspan  = -3,           # minspan for varmod call to earth
    Scale.y         = NULL,         # TRUE to scale y in the forward pass
    ...)                            # passed on to earth.fit
{
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    env <- parent.frame() # the environment from which earth was called
    call <- make.call.generic(match.call(), "earth")
    if(trace >= 4)
        printcall("Call: ", call)
    if(!is.null(call$da)) # matches anything beginning with "da", doesn't alias with other args
        stop0("'data' argument not allowed in earth.default")
    if(is.character(na.action)) {
        if(is.na(pmatch(na.action, "na.fail")))
            stop0("illegal 'na.action', only na.action=na.fail is allowed")
    } else if(!identical(na.action, na.fail))
        stop0("illegal 'na.action', only na.action=na.fail is allowed")
    keepxy <- check.numeric.scalar(keepxy, min=0, max=2, logical.ok=TRUE)
    pmethod <- match.arg1(pmethod, "pmethod")
    if(pmethod == "cv")
        keepxy <- min(1, keepxy)
    if(keepxy) {
        xkeep <- x
        ykeep <- y
        # TODO this should be done in one place instead of here and also in earth.fit
        # TODO does lm process the weights before or after subset (I'm assuming before)
        if(!is.null(subset)) {
            # duplicates are allowed in subset so user can specify a bootstrap sample
            subset <- check.index(subset, "subset", xkeep, allow.dups=TRUE, allow.zeros=TRUE)
            xkeep <- if(is.null(dim(xkeep))) xkeep[subset] else xkeep[subset, , drop=FALSE]
            ykeep <- if(is.null(dim(ykeep))) ykeep[subset] else ykeep[subset, , drop=FALSE]
        }
    }
    xname <- trunc.deparse(substitute(x))
    if(is.matrix(x) && is.double(x[,1])) {
        # x is already a double matrix, save memory and time by skipping expand.arg
        colnames(x) <- gen.colnames(x, xname, "x")
        modvars <- get.identity.modvars(x)
    } else {
        # expand factors, convert to double matrix with column names
        ret <- expand.arg.modvars(x, env, trace=0, is.y.arg=FALSE, name=xname)
        x       <- ret$x
        modvars <- ret$modvars
        rownames(x) <- possibly.delete.rownames(x)
    }
    ylevels <- get.ylevels(y)
    y <- expand.arg(y, env, trace=0, is.y.arg=TRUE, name=trunc.deparse(substitute(y)))
    rownames(y) <- possibly.delete.rownames(y)
    # we need the diag of the hat matrix if varmod.method != "none"
    varmod.method <- match.choices(varmod.method,
                                   c("none", VARMOD.METHODS), "varmod.method")
    check.cv.args(nfold, ncross, pmethod, varmod.method)
    pmethod1 <- pmethod
    update.earth.called.me <- !is.null(dota("Object", DEF=NULL, ...))
    if(pmethod == "cv" && !update.earth.called.me) {
        trace1(trace,
               "=== pmethod=\"cv\": Preliminary model with pmethod=\"backward\" ===\n")
        if(nfold <= 1 || ncross < 1)
            stop0("the nfold argument must be specified when pmethod=\"cv\"")
        pmethod1 <- "backward"
    }
    rv <- earth.fit(x=x, y=y, weights=weights, wp=wp, subset=subset,
                    na.action=na.action, offset=NULL, pmethod=pmethod1,
                    keepxy=keepxy, trace=trace, glm=glm, degree=degree,
                    nprune=nprune, Scale.y=Scale.y, ...)

    rv$call <- call

    # namesx is no longer used by earth but provided for back compat with other apps
    # (modvars subsumes namesx)
    rv$namesx <- rownames(modvars)

    # The modvars field was added to earth models in Sep 2020 (earth 5.2.0).
    #
    # rownames(modvars) for earth.default differs from earth.formula when for example:
    #    earth(x, y) (where x has two columns, "num" and "sqrt(num)")
    #       then rownames(modvars) = c("num", "sqrt(num)")
    #    earth(y~num+sqrt(num), data=dat)  (where dat has two columns, "y" and "num")
    #      then rownames(modvars) = "num"

    rv$modvars <- modvars

    rv$levels <- ylevels
    rv$wp <- wp
    if(keepxy) {
        rv[["x"]] <- xkeep
        # following ruse is needed else vector y's lose their name
        # Jun 2015: but we don't use it for factors because as.matrix
        #           converts factors to character variables (!)
        if(is.null(dim(ykeep)) && !is.factor(ykeep)) {
            ykeep <- as.matrix(ykeep)
            colnames(ykeep) <- colnames(y)[1]
        }
        rv[["y"]] <- ykeep
        rv$subset <- subset
        # TODO consider doing the following
        # rv$.Environment <- parent.frame()
    }
    if(nfold > 1 && ncross >= 1) {
        if(!is.null(subset))
            stop0("'subset' cannot be used with 'nfold' (implementation restriction)")
        glm.arg <- process.glm.arg(glm)
        cv <- earth_cv(object=rv,
                x=if(is.null(subset)) x else x[subset,,drop=FALSE],
                y=if(is.null(subset)) y else y[subset,,drop=FALSE],
                subset=subset, weights=weights, na.action=na.action,
                pmethod=pmethod, keepxy=keepxy,
                trace=if(trace >= 4.1) trace else if(trace) .5 else 0,
                trace.org=trace,
                glm.arg=glm.arg, degree=degree,
                nfold=nfold, ncross=ncross, stratify=stratify,
                get.oof.fit.tab = varmod.method != "none",
                get.oof.rsq.per.subset = keepxy || pmethod == "cv",
                Scale.y=rv$Scale.y, env=env, ...)
        rv$cv.list <- cv$cv.list
        rv$cv.nterms.selected.by.gcv  <- cv$nterms.selected.by.gcv
        rv$cv.nvars.selected.by.gcv   <- cv$nvars.selected.by.gcv
        rv$cv.groups                  <- cv$groups # groups used for cross validation
        rv$cv.rsq.tab                 <- cv$rsq.tab
        rv$cv.maxerr.tab              <- cv$maxerr.tab
        if(!is.null(cv$class.rate.tab))  rv$cv.class.rate.tab  <- cv$class.rate.tab
        if(!is.null(cv$auc.tab))         rv$cv.auc.tab         <- cv$auc.tab
        if(!is.null(cv$cor.tab))         rv$cv.cor.tab         <- cv$cor.tab
        if(!is.null(cv$deviance.tab))    rv$cv.deviance.tab    <- cv$deviance.tab
        if(!is.null(cv$calib.int.tab))   rv$cv.calib.int.tab   <- cv$calib.int.tab
        if(!is.null(cv$calib.slope.tab)) rv$cv.calib.slope.tab <- cv$calib.slope.tab
        if(!is.null(cv$oof.rsq.tab))     rv$cv.oof.rsq.tab     <- cv$oof.rsq.tab
        if(!is.null(cv$infold.rsq.tab))  rv$cv.infold.rsq.tab  <- cv$infold.rsq.tab
        if(!is.null(cv$oof.fit.tab))     rv$cv.oof.fit.tab     <- cv$oof.fit.tab
        if(pmethod == "cv") {
            rv.backward <- rv
            tab <- rv$cv.oof.rsq.tab
            stopifnot(nrow(tab) > 1, ncol(tab) > 1)
            mean.oof.rsq.per.subset <- tab[nrow(tab),]
            # Sep 2020: nprune1 added else get (for some rand seeds):
            #    evimp: Error in object$prune.terms[isubset, -1] : subscript out of bounds
            nprune1 <- if(is.specified(nprune)) nprune else length(mean.oof.rsq.per.subset)
            nterms.selected.by.cv <- which.max(mean.oof.rsq.per.subset[1:nprune1])
            trace1(trace,
"\n=== pmethod=\"cv\": Calling update.earth internally for nterms selected by cv %g ===\n",
                nterms.selected.by.cv)
            trace2(trace, "\n")
            # July 2017 TODO following necessary for eval.parent(call) in update.earth
            penalty         <- dota("penalty",         DEF=if(degree > 1) 3 else 2, ...)
            nk              <- dota("nk",              DEF=min(200, max(20, 2 * ncol(x))) + 1, ...)
            thresh          <- dota("thresh",          DEF=0.001, ...)
            minspan         <- dota("minspan",         DEF=0, ...)
            endspan         <- dota("endspan",         DEF=0, ...)
            newvar.penalty  <- dota("newvar.penalty",  DEF=0, ...)
            fast.k          <- dota("fast.k",          DEF=20, ...)
            fast.beta       <- dota("fast.beta",       DEF=1, ...)
            linpreds        <- dota("linpreds",        DEF=FALSE, ...)
            allowed         <- dota("allowed",         DEF=NULL, ...)
            Object          <- dota("Object ",         DEF=NULL, ...)
            Adjust.endspan  <- dota("Adjust.endspan",  DEF=2, ...)
            Auto.linpreds   <- dota("Auto.linpreds",   DEF=TRUE, ...)
            Force.weights   <- dota("Force.weights",   DEF=FALSE, ...)
            Use.beta.cache  <- dota("Use.beta.cache",  DEF=TRUE, ...)
            Force.xtx.prune <- dota("Force.xtx.prune", DEF=FALSE, ...)
            Get.leverages   <- dota("Get.leverages",   DEF=NROW(x) < 1e5, ...)
            Exhaustive.tol  <- dota("Exhaustive.tol",  DEF=1e-10, ...)

            rv <- update.earth(rv, ponly=TRUE, trace=trace,
                               pmethod="cv", nprune=nterms.selected.by.cv,
                               nfold=0, ncross=1,
                               glm=glm, varmod.method="none",
                               # July 2017 TODO following necessary for eval.parent(call) in update.earth
                               penalty=penalty,
                               nk=nk,
                               thresh=thresh,
                               minspan=minspan,
                               endspan=endspan,
                               newvar.penalty=newvar.penalty,
                               fast.k=fast.k,
                               fast.beta=fast.beta,
                               linpreds=linpreds,
                               allowed=allowed,
                               Object=Object,
                               Adjust.endspan=Adjust.endspan,
                               Auto.linpreds=Auto.linpreds,
                               Force.weights=Force.weights,
                               Use.beta.cache=Use.beta.cache,
                               Force.xtx.prune=Force.xtx.prune,
                               Get.leverages=Get.leverages,
                               Exhaustive.tol=Exhaustive.tol)

            if(trace == .5)
                printf("%sGRSq %.3f RSq %.3f nterms selected by cv %g",
                    "Final model with pmethod=\"cv\": ",
                    rv$grsq , rv$rsq, length(rv$selected.terms))

            rv$call    <- call
            rv$pmethod <- "cv"
            rv$nprune  <- nprune
            rv$dirs           <- rv.backward$dirs
            rv$cuts           <- rv.backward$cuts
            rv$prune.terms    <- rv.backward$prune.terms
            rv$rss.per.subset <- rv.backward$rss.per.subset
            rv$gcv.per.subset <- rv.backward$gcv.per.subset
            rv$cv.list        <- rv.backward$cv.list
            rv$cv             <- rv.backward$cv
            # the number of terms that would have been selected by pmethod="backward"
            rv$backward.selected.terms <- rv.backward$selected.terms

            rv$cv.oof.fit.tab     = rv.backward$cv.oof.fit.tab
            rv$cv.infold.rsq.tab  = rv.backward$cv.infold.rsq.tab
            rv$cv.oof.rsq.tab     = rv.backward$cv.oof.rsq.tab

            # The following were calculated using the best model selected at each
            # fold using the fold's GCVs.  To minimize confusion, we delete them.
            rv$cv.rsq.tab         <- NULL
            rv$cv.maxerr.tab      <- NULL
            rv$cv.class.rate.tab  <- NULL
            rv$cv.auc.tab         <- NULL
            rv$cv.cor.tab         <- NULL
            rv$cv.deviance.tab    <- NULL
            rv$cv.calib.int.tab   <- NULL
            rv$cv.calib.slope.tab <- NULL
        }
        if(trace >= .5)
            printf("\n")
    }
    # TODO only do varmod for final model if pmethod="cv", similarly for glm
    if(varmod.method != "none") {
        oof.fit.tab <- rv$cv.oof.fit.tab
        stopifnot(!is.null(oof.fit.tab))
        model.var <- apply(oof.fit.tab,  1, var) # var  of each row of oof.fit.tab
        model.var <- matrix(model.var, ncol=1)
        rv$varmod <- varmod(rv, varmod.method, varmod.exponent,
                            varmod.conv, varmod.clamp, varmod.minspan,
                            trace, x, y, model.var)
    }
    rv$Scale.y <- NULL
    rv
}
earth.formula <- function( # user called earth with formula arg (not x y args)
    formula         = stop("no 'formula' argument"), # intercept will be ignored
    data            = NULL,
    weights         = NULL,
    wp              = NULL,
    subset          = NULL,
    na.action       = na.fail,
    pmethod         = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
    keepxy          = FALSE,
    trace           = 0,
    glm             = NULL,
    degree          = 1,    # max degree of interaction (1=additive model) (Friedman's mi)
    nprune          = NULL, # max nbr of terms (including intercept) in pruned subset
    nfold           = 0,
    ncross          = 1,
    stratify        = TRUE,
    varmod.method   = "none",
    varmod.exponent = 1,
    varmod.conv     = 1,
    varmod.clamp    = .1,
    varmod.minspan  = -3,
    Scale.y         = NULL,
    ...)
{
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    env <- parent.frame() # the environment from which earth was called
    call <- make.call.generic(match.call(), "earth")
    if(trace >= 4)
        printcall("Call: ", call)
    if(!is.null(call[["x"]]))
        stop0("'x' argument not allowed in earth.formula")
    if(!is.null(call[["y"]]))
        stop0("'y' argument not allowed in earth.formula")

    call2 <- match.call(expand.dots=FALSE)

    # subset is handled in earth.fit so it isn't included here
    # we handle weights here in the same way as source code of lm (so
    # weights are first searched for in the data passed to the formula)

    m <- match(c("formula", "data", "weights", "na.action", "offset"), names(call2), 0)
    formdat <- get.data.from.formula(mf=call2[c(1, m)], formula, data, env, trace)
        x       <- formdat$x
        y       <- formdat$y
        weights <- formdat$weights
        offset  <- formdat$offset
        terms   <- formdat$terms
        xlevels <- formdat$xlevels
        ylevels <- formdat$ylevels

    keepxy <- check.numeric.scalar(keepxy, min=0, max=2, logical.ok=TRUE)
    pmethod <- match.arg1(pmethod, "pmethod")
    # we need the diag of the hat matrix if varmod.method != "none"
    varmod.method <- match.choices(varmod.method,
                                   c("none", VARMOD.METHODS), "varmod.method")
    check.cv.args(nfold, ncross, pmethod, varmod.method)
    pmethod1 <- pmethod
    update.earth.called.me <- !is.null(dota("Object", DEF=NULL, ...))
    if(pmethod == "cv" && !update.earth.called.me) {
        trace1(trace,
"=== pmethod=\"cv\": Preliminary model with pmethod=\"backward\" ===\n")
        if(nfold <= 1 || ncross < 1)
            stop0("the nfold argument must be specified when pmethod=\"cv\"")
        pmethod1 <- "backward"
    }
    rv <- earth.fit(x=x, y=y, weights=weights, wp=wp, subset=subset,
                    na.action=na.action, offset=offset, pmethod=pmethod1,
                    keepxy=keepxy, trace=trace, glm=glm, degree=degree,
                    nprune=nprune, Scale.y=Scale.y, ...)

    rv$call <- call

    # namesx is no longer used by earth but provided for back compat with other apps
    # (modvars has subsumed namesx)
    rv$namesx <- rownames(formdat$modvars)

    rv$modvars <- formdat$modvars  # modvars added to earth models in Sep 2020 (earth 5.2.0).

    rv$terms <- terms
    # March 2019: added xlevels to match what lm does (and so does linmod.R in the plotmo tests)
    rv$xlevels <- xlevels
    rv$levels <- ylevels
    rv$wp <- wp
    if(keepxy) {
        if(!is.null(data))
            rv$data <- data
        else if(trace >= 0)
            warning0("No 'data' argument to earth so 'keepxy' is limited\n")
        rv[["y"]] <- y
        if(!is.null(subset)) {
            # duplicates are allowed in subset so user can specify a bootstrap sample
            subset <- check.index(subset, "subset", data, allow.dups=TRUE, allow.zeros=TRUE)
            rv$data <- data[subset, , drop=FALSE]
            rv[["y"]] <- rv[["y"]][subset, , drop=FALSE]
        }
        rv$subset <- subset
    }
    # TODO make the following code a subroutine, it's identical to code in earth.default
    if(nfold > 1 && ncross >= 1) {
        if(!is.null(subset))
            stop0("'subset' cannot be used with 'nfold' (implementation restriction)")
        glm.arg <- process.glm.arg(glm)
        cv <- earth_cv(object=rv,
                x=if(is.null(subset)) x else x[subset,,drop=FALSE],
                y=if(is.null(subset)) y else y[subset,,drop=FALSE],
                subset=subset, weights=weights, na.action=na.action,
                pmethod=pmethod, keepxy=keepxy,
                trace=if(trace >= 4.1) trace else if(trace) .5 else 0,
                trace.org=trace,
                glm.arg=glm.arg, degree=degree,
                nfold=nfold, ncross=ncross, stratify=stratify,
                get.oof.fit.tab = varmod.method != "none",
                get.oof.rsq.per.subset = keepxy || pmethod == "cv",
                Scale.y=rv$Scale.y, env=env, ...)
        rv$cv.list <- cv$cv.list
        rv$cv.nterms.selected.by.gcv  <- cv$nterms.selected.by.gcv
        rv$cv.nvars.selected.by.gcv   <- cv$nvars.selected.by.gcv
        rv$cv.groups                  <- cv$groups # groups used for cross validation
        rv$cv.rsq.tab                 <- cv$rsq.tab
        rv$cv.maxerr.tab              <- cv$maxerr.tab
        if(!is.null(cv$class.rate.tab))  rv$cv.class.rate.tab  <- cv$class.rate.tab
        if(!is.null(cv$auc.tab))         rv$cv.auc.tab         <- cv$auc.tab
        if(!is.null(cv$cor.tab))         rv$cv.cor.tab         <- cv$cor.tab
        if(!is.null(cv$deviance.tab))    rv$cv.deviance.tab    <- cv$deviance.tab
        if(!is.null(cv$calib.int.tab))   rv$cv.calib.int.tab   <- cv$calib.int.tab
        if(!is.null(cv$calib.slope.tab)) rv$cv.calib.slope.tab <- cv$calib.slope.tab
        if(!is.null(cv$oof.rsq.tab))     rv$cv.oof.rsq.tab     <- cv$oof.rsq.tab
        if(!is.null(cv$infold.rsq.tab))  rv$cv.infold.rsq.tab  <- cv$infold.rsq.tab
        if(!is.null(cv$oof.fit.tab))     rv$cv.oof.fit.tab     <- cv$oof.fit.tab
        if(pmethod == "cv") {
            rv.backward <- rv
            tab <- rv$cv.oof.rsq.tab
            stopifnot(nrow(tab) > 1, ncol(tab) > 1)
            mean.oof.rsq.per.subset <- tab[nrow(tab),]
            # Sep 2020: nprune1 added (necessary only for some rand seeds):
            #    evimp: Error in object$prune.terms[isubset, -1] : subscript out of bounds
            nprune1 <- if(is.specified(nprune)) nprune else length(mean.oof.rsq.per.subset)
            nterms.selected.by.cv <- which.max(mean.oof.rsq.per.subset[1:nprune1])
            trace1(trace,
"\n=== pmethod=\"cv\": Calling update.earth internally for nterms selected by cv %g ===\n",
                nterms.selected.by.cv)
            trace2(trace, "\n")
            # July 2017 TODO following necessary for eval.parent(call) in update.earth
            penalty         <- dota("penalty",         DEF=if(degree > 1) 3 else 2, ...)
            nk              <- dota("nk",              DEF=min(200, max(20, 2 * ncol(x))) + 1, ...)
            thresh          <- dota("thresh",          DEF=0.001, ...)
            minspan         <- dota("minspan",         DEF=0, ...)
            endspan         <- dota("endspan",         DEF=0, ...)
            newvar.penalty  <- dota("newvar.penalty",  DEF=0, ...)
            fast.k          <- dota("fast.k",          DEF=20, ...)
            fast.beta       <- dota("fast.beta",       DEF=1, ...)
            linpreds        <- dota("linpreds",        DEF=FALSE, ...)
            allowed         <- dota("allowed",         DEF=NULL, ...)
            Object          <- dota("Object ",         DEF=NULL, ...)
            Adjust.endspan  <- dota("Adjust.endspan",  DEF=2, ...)
            Auto.linpreds   <- dota("Auto.linpreds",   DEF=TRUE, ...)
            Force.weights   <- dota("Force.weights",   DEF=FALSE, ...)
            Use.beta.cache  <- dota("Use.beta.cache",  DEF=TRUE, ...)
            Force.xtx.prune <- dota("Force.xtx.prune", DEF=FALSE, ...)
            Get.leverages   <- dota("Get.leverages",   DEF=NROW(x) < 1e5, ...)
            Exhaustive.tol  <- dota("Exhaustive.tol",  DEF=1e-10, ...)

            # July 2017 TODO necessary when form is a local var in function calling earth.formula
            rv$call$formula <- formula
            rv$call$data <- data

            rv <- update.earth(rv, ponly=TRUE, trace=trace,
                               pmethod="cv", nprune=nterms.selected.by.cv,
                               nfold=0, ncross=1,
                               glm=glm, varmod.method="none",
                               # July 2017 TODO following necessary for eval.parent(call) in update.earth
                               penalty=penalty,
                               nk=nk,
                               thresh=thresh,
                               minspan=minspan,
                               endspan=endspan,
                               newvar.penalty=newvar.penalty,
                               fast.k=fast.k,
                               fast.beta=fast.beta,
                               linpreds=linpreds,
                               allowed=allowed,
                               Object=Object,
                               Adjust.endspan=Adjust.endspan,
                               Auto.linpreds=Auto.linpreds,
                               Force.weights=Force.weights,
                               Use.beta.cache=Use.beta.cache,
                               Force.xtx.prune=Force.xtx.prune,
                               Get.leverages=Get.leverages,
                               Exhaustive.tol=Exhaustive.tol)

            if(trace == .5)
                printf("%sGRSq %.3f RSq %.3f nterms selected by cv %g",
                    "Final model with pmethod=\"cv\": ",
                    rv$grsq , rv$rsq, length(rv$selected.terms))

            rv$call    <- call
            rv$pmethod <- "cv"
            rv$nprune  <- nprune
            rv$dirs           <- rv.backward$dirs
            rv$cuts           <- rv.backward$cuts
            rv$prune.terms    <- rv.backward$prune.terms
            rv$rss.per.subset <- rv.backward$rss.per.subset
            rv$gcv.per.subset <- rv.backward$gcv.per.subset
            rv$cv.list        <- rv.backward$cv.list
            rv$cv             <- rv.backward$cv
            # the number of terms that would have been selected by pmethod="backward"
            rv$backward.selected.terms <- rv.backward$selected.terms

            rv$cv.oof.fit.tab     = rv.backward$cv.oof.fit.tab
            rv$cv.infold.rsq.tab  = rv.backward$cv.infold.rsq.tab
            rv$cv.oof.rsq.tab     = rv.backward$cv.oof.rsq.tab

            # The following were calculated using the best model selected at each
            # fold using the fold's GCVs.  To minimize confusion, we delete them.
            rv$cv.rsq.tab         <- NULL
            rv$cv.maxerr.tab      <- NULL
            rv$cv.class.rate.tab  <- NULL
            rv$cv.auc.tab         <- NULL
            rv$cv.cor.tab         <- NULL
            rv$cv.deviance.tab    <- NULL
            rv$cv.calib.int.tab   <- NULL
            rv$cv.calib.slope.tab <- NULL
        }
        if(trace >= .5)
            printf("\n")
    }
    # TODO only do varmod for final model if pmethod="cv", similarly for glm
    if(varmod.method != "none") {
        oof.fit.tab <- rv$cv.oof.fit.tab
        stopifnot(!is.null(oof.fit.tab))
        model.var <- apply(oof.fit.tab,  1, var) # var  of each row of oof.fit.tab
        model.var <- matrix(model.var, ncol=1)
        rv$varmod <- varmod(rv, varmod.method, varmod.exponent,
                            varmod.conv, varmod.clamp, varmod.minspan,
                            trace, x, y, model.var)
    }
    rv$Scale.y <- NULL
    rv
}
# Like stats::.getXlevels but also works for Terms for multiple-response
# model terms made with Formula and with a "Response" attribute.
# If we use .getXlevels and not .getXlevels2, model.frame.default issues
#     Warning: variable 'pclass' is not a factor
# for e.g. a<-earth(pclass+age~sibsp, data=etitanic); plotmo(a, nresponse=1)
# This function is based on .getXlevels R version 3.5.3 (March 2019).

.getXlevelsMulti <- function(Terms, m)
{
    deparse2 <- function(x) { # copy of stats:::deparse2
        paste(deparse(x, width.cutoff = 500L,
              backtick = !is.symbol(x) && is.language(x)), collapse = " ")
    }
    xvars <- vapply(attr(Terms, "variables"), deparse2, "")[-1L]
    yvars <- attr(Terms, "response")
    if(is.null(yvars) || yvars[1] == 0)
        yvars <- attr(Terms, "Response")
    if(any(yvars) > 0) xvars <- xvars[-yvars]
    if(length(xvars)) {
        xlev <- lapply(m[xvars],
                function(x)
                    if(is.factor(x)) levels(x)
                    else if (is.character(x)) levels(as.factor(x))
                    else NULL)
        xlev[!vapply(xlev, is.null, NA)]
    } else NULL
}
check.which.terms <- function(dirs, which.terms) # ensure which.terms is valid
{
    if(is.null(which.terms))
        stop0("'which.terms' is NULL")
    if(length(which.terms) == 0)
        stop0("length(which.terms) == 0")
    if(which.terms[1] != 1)
        stop0("first element of 'which.terms' must be 1, the intercept term")
    if(NCOL(which.terms) > 1) {
        for(i in seq_len(NCOL(which.terms)))
            check.index(which.terms[,i], "which.terms", dirs,
                        allow.zeros=TRUE, allow.dups=TRUE)
    } else
            check.index(which.terms, "which.terms", dirs,
                        allow.zeros=TRUE, allow.dups=TRUE)
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
# called only by earth.formula

get.data.from.formula <- function(mf, formula, data, env, trace)
{
    mf[[1]] <- as.name("model.frame")
    if(!is.null(mf$na.action))
        stop0("'na.action' argument is not allowed (it is set internally to na.fail)")
    mf$na.action <- na.fail

    # for backward compat, use class "Formula" -- not "formula" -- only when necessary
    if(must.use.Formula(formula)) {
        # use class "Formula" (allows multiple responses separated by +)

        trace1(trace, "Using class \"Formula\" because lhs of formula has terms separated by \"+\"\n")
        Formula <- Formula::Formula(formula)
        if(length(attr(Formula, "lhs")) > 1) # e.g. y1 | y2 ~ .
            stop0("multiple parts on left side of formula (because \"|\" was used)")
        if(length(attr(Formula, "rhs")) > 1) # e.g. y ~ x1 | x2
            stop0("multiple parts on right side of formula (because \"|\" was used)")
        mf$formula <- Formula
        mf <- eval(mf, envir=env)
        terms <- terms(mf)
        x <- model.matrix(Formula, data=mf, rhs=1)
        # TODO work around for model.matrix.Formula which incorrectly includes
        # `(weights)` in x (Formula package version 1.2-3 March 2019).
        # Happens only if dot is used on rhs of formula?
        # e.g. d<-data.frame(x=1:9,y=1:9,z=1:9);earth(y+z~.,data=d,weights=1:9,trace=2)
        if(any(colnames(x) == "`(weights)`")) {
            trace2(trace, "Deleting `(weights)` column from 'Formula' model.matrix\n")
            x <- x[, colnames(x) != "`(weights)`", drop=FALSE]
        }
        y <- model.part(Formula, data=mf, lhs=1)
        # TODO Sep 2020: work around for model.matrix.Formula which incorrectly includes
        # `log(O3)` in x if log(O3) is used in y (i.e. on the lhs of the formula)
        # e.g. earth(log(O3) + wind ~ ., data=ozone1) is isn't handled correctly
        which <- which(colnames(y) != naken(colnames(y)))
        if(any(which))
            stop0("terms like \'", colnames(y)[which[1]],
                  "\' are not allowed on the LHS of a multiple-response formula")
        # add extra attributes to terms for use by earth
        attr(terms, "Formula") <- Formula
        attr(terms, "Response") <- 1:NCOL(y) # TODO is 1:NCOL(y) reliable here?
        iresp <- attr(terms, "Response") # is a vector if multiple responses
                                           # (empirically -1 works even with mult responses)
    } else {
        # use class "formula" (lhs does not have two terms separated by + or |)
        # we use formula not Formula here for backwards compatibility
        # could still be a multiple response e.g. cbind(survived, died)~.
        mf <- eval(mf, envir=env)
        terms <- terms(mf)
        # expand factors in x, convert to double matrix, add colnames
        x <- model.matrix(terms, data=mf)
        y <- model.response(mf, "any")  # "any" means factors are allowed
        iresp <- attr(terms, "response") # is 1 even with cbind(survived, died)~.

    }
    # the "assign" attribute has an entry for each column in x
    # giving the term in the formula which gave rise to the column
    xassign <- attr(x, "assign")
    xassign <- xassign[-1] # delete response (-1 correct even with multiple responses)

    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]     # silently discard intercept
    else
        warning0("ignored -1 in formula (earth objects always have an intercept)")

    rownames(x) <- possibly.delete.rownames(x)

    ylevels <- get.ylevels(y) # TODO this always returns NULL if Formula was used

    # expand factors in y, convert to double matrix, add colnames
    if(length(iresp) > 1)                                      # multiple columns
        yname <- "y"                                           # generic name
    else
        yname <- names(attr(terms, "dataClasses"))[[iresp[1]]] # name of variable
    y <- expand.arg(y, env, trace=0, is.y.arg=TRUE, name=yname)
    rownames(y) <- possibly.delete.rownames(y)

    # as.vector to avoid problems with Nx1 weights (same as lm source code)
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if(!is.null(offset)) {
        check.offset.var.is.in.data(terms, data)
        check.vec(offset, "formula offset", expected.len=NROW(y), logical.ok=FALSE)
    }
    list(x=x, y=y, weights=weights, offset=offset,
         terms=terms,
         modvars=get.modvars(x, xassign, terms), # Sep 2020
         xlevels=.getXlevelsMulti(terms, mf), ylevels=ylevels)
}
# modvars is a matrix specifying which input variables
# are used in each column of the model matrix.
#
# Called by get.data.from.formula()
# and by expand.arg.modvars() when the arg has to expanded by calling model.matrix.
#
# Columns correspond to columns of the model matrix (same as cols of earth$dirs).
# Rows correspond to unique "naked" variables in the formula.
#
# Example (cf earth.object.Rd):
#
#   formula: survived ~ age + pclass + sqrt(age) + sex + sex:parch + offset(off)
#
#   attr(terms,"factors"),  same colnames as earth$dirs:
#
#                 age pclass sqrt(age) sex sex:parch
#     survived      0      0         0   0         0  # response will be dropped
#     age           1      0         0   0         0
#     pclass        0      1         0   0         0
#     sqrt(age)     0      0         1   0         0 # sqrt(age) will be merged with age
#     sex           0      0         0   1         2 # 2 is inherited from attr(terms,"factors")
#     parch         0      0         0   0         1
#     offset(off)   0      0         0   0         0 # offset will be marked as 9999 below
#
#   modvars:
#             age pclass2nd pclass3rd sqrt(age) sexmale sexfemale:parch sexmale:parch
#     age       1         0         0         1       0               0             0
#     pclass    0         1         1         0       0               0             0
#     sex       0         0         0         0       1               2             2
#     parch     0         0         0         0       0               1             1
#     off    9999      9999      9999      9999    9999            9999          9999

get.modvars <- function(x, xassign, terms)
{
    factors <- attr(terms,"factors")
    offset <- 0
    if(!is.null(attr(terms,"offset")))
        offset <- attr(terms,"offset")
    if(offset) # mark offset row as special
        factors[offset,] <- rep(9999, length.out=ncol(factors))

    # drop response

    # following relies on extra attr "Response" added in get.data.from.formula
    is.Formula <- !is.null(attr(terms, "Response"))

    if(is.Formula) {    # Formula interface: factors includes response rows and cols
        iresp <- attr(terms, "Response")
        if(all(iresp) > 0)
            factors <- factors[-iresp, -iresp, drop=FALSE]
    } else {            # formula interface: includes response row but not response col
        iresp <- attr(terms, "response")
        if(all(iresp) > 0) # iresp will be 0 if  no respose e.g. formula is ~x1+x2
            factors <- factors[-iresp, , drop=FALSE]
    }
    # keep rows only for used variables
    # (note: variables will be unused if there is a "-" in the formula)
    # e.g. rownames num,int,fac,sqrt(num),ord,bool,offset(off)
    #      becomes  num,    fac,sqrt(num),ord,bool,offset(off)
    factors <- factors[which(rowSums(factors) > 0), , drop=FALSE]
    modvars <- matrix(0, nrow=nrow(factors), ncol=ncol(x))
    colnames(modvars) <- colnames(x)
    rownames(modvars) <- naken(rownames(factors)) # e.g. num,fac,sqrt(num),ord,bool,offset(off)
                                                  # to   num,fac,num       ord,bool,off
    nrow <- nrow(modvars)
    unique <- rep(TRUE, length.out=nrow)
    if(nrow > 0) {
        for(irow in 1:nrow(modvars))
            for(icol in 1:ncol(modvars))
                modvars[irow, icol] <- factors[irow, xassign[icol]]
    }
    if(nrow > 1) {
        # merge rows with duplicate rownames into first row with that rowname
        rownames <- rownames(modvars)
        for(i in 1:(nrow-1))
            for(j in (i+1):nrow)
                if(rownames[j] == rownames[i]) {    # row j is a duplicate of row i?
                    unique[j] <- FALSE
                    modvars[i,] <- modvars[i,] + modvars[j,] # merge row j into row i
                }
    }
    # drop rows with duplicated rownames
    modvars <- modvars[unique, , drop=FALSE]      # to    num,fac,          ord,bool,off
    modvars
}
get.identity.modvars <- function(x) # uses ncol(x) and colnames(x)
{
    modvars <- diag(ncol(x))
    colnames(modvars) <- rownames(modvars) <- colnames(x)
    modvars
}
get.ylevels <- function(y)
{
    if(!is.null(levels(y)))
        return(levels(y))
    # following needed for predict.earth(type="class")
    if(is.logical(y))
        return(c(FALSE, TRUE))
    if(is.numeric(y)) {
        range <- range(y, na.rm=TRUE) # forward pass will check NAs later
        if(range[2] - range[1] == 1)
            return(c(range[1], range[2]))
    }
    NULL
}
good.colname <- function(name)
{
    # The nchar check prevents super long names (60 is arb)
    # that are actually contents of vectors e.g. c(1,2,3,etc.)
    # The grep ensures that there are no more than three commas,
    # also to prevent using the contents of vectors.

    !is.null(name) && nchar(name) <= 60 && !grepany(",.*,.*,", name)
}
good.colnames <- function(x)
{
    colnames <- colnames(x)
    if(is.null(colnames))
        return(FALSE)
    for(i in seq_along(colnames))
        if(!good.colname(colnames[i]))
            return(FALSE)
    return(TRUE)
}
# If lhs of formula has two terms separated by + or |, use
# class "Formula" (to support multiple responses)
#     e.g. y+y2~.
#     e.g. y|y2~. (issue an error message if this is used, to help user)
#
# Otherwise, use class "formula" (for backwards compatability)
#     e.g. y~.              standard case
#     e.g. y/y2~.           because + isn't used
#     e.g. cbind(y+y2)~.    because + is internal to cbind
#     e.g. (y+y2)~.         because + is in parentheses
#     e.g. I(y+y2)~.        because + is in parentheses
#
# TODO this function depends on the implementation of all.names (it depends
# on the order in which all.names returns the names in the formula)

must.use.Formula <- function(formula)
{
    all.names <- all.names(formula)
    length(all.names) > 3 && all.names[1] == "~" &&
           (all.names[2] == "+" || all.names[2] == "|")
}
# Remove useless(?) "1" "2" "3" ... rownames for x (added by
# model.matrix) so earth.formula x is same as earth.default x,
# and to save memory (although not as important now that R
# hashes strings internally).

possibly.delete.rownames <- function(x)
{
    # decide by looking at first few names
    n <- length(rownames(x))
    if((n >= 1 && (is.null(rownames(x)[1]) || rownames(x)[1] == "1")) &&
       (n >= 2 && (is.null(rownames(x)[2]) || rownames(x)[2] == "2")) &&
       (n >= 3 && (is.null(rownames(x)[3]) || rownames(x)[3] == "3")))
        NULL
    else
        rownames(x)
}
# print a reminder if exhaustive pruning will be slow
possibly.print.exhaustive.pruning.reminder <- function(nprune, trace, bx, bx.cond)
{
    nsubsets <- 0 # approx, assumes brute force exhaustive search
    for(subset.size in seq_len(nprune))
        nsubsets <- nsubsets + choose(ncol(bx), subset.size)
    if(trace >= 1 || nsubsets > 1e9) {
        cat0("Exhaustive pruning: number of subsets ",
             format(nsubsets, digits=2),
             "   bx sing val ratio ", format(bx.cond, digits=2), "\n")
    }
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
    x           = stop("no 'x' argument"),
    which.terms = x$selected.terms,
    decomp      = c("anova", "none"),
    degree      = 99,       # max degree, 0 returns just the intercept
    min.degree  = 0,
    ...)                    # unused
{
    warn.if.dots(...)
    if(degree < 0)
        stop0("degree ", degree, " < 0")
    if(min.degree < 0)
        stop0("min.degree ", min.degree, " < 0")
    if(degree < min.degree)
        stop0("degree ", degree, " < min.degree ", min.degree)
    check.which.terms(x$dirs, which.terms)
    dirs <- x$dirs[which.terms, , drop=FALSE]
    new.order <- switch(match.arg1(decomp, "decomp"),
                   anova = reorder_terms_anova(
                                dirs, x$cuts[which.terms,,drop=FALSE]),
                   none  = 1:length(which.terms))
    degrees <- get.degrees.per.term(dirs[new.order, , drop=FALSE])
    new.order[degrees >= min.degree & degrees <= degree]
}
# return a vector of term numbers, ordered as per the "anova" decomposition

reorder_terms_anova <- function(dirs, cuts)
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
# b) This function also handle appropriately objects that were or were
#    not created using a formula i.e. were created by earth.formula() or
#    by earth.default().
#
# c) This function retrieves x and y from object$x and object$y if need be
#    and also data, weights, wp, and subset.

update.earth <- function(
    object   = stop("no 'object' argument"),
    formula. = NULL,    # formula. is optional
    ponly    = FALSE,   # force prune only, no forward pass
    ...,                # dots passed on to earth()
    evaluate = TRUE)    # for compatibility with generic update
{
    check.classname(object, substitute(object), "earth")
    call <- object$call
    stopifnot(!is.null(call))
    do.forward.pass <- FALSE
    if(!is.null(formula.)) {
        if(is.null(call$formula))
            stop0("'formula.' argument is not allowed on ",
                  "objects created without a formula")
        call$formula <- update.formula(formula(object), formula.)
        do.forward.pass <- TRUE
    }
    env <- parent.frame() # TODO should use model.env(object) here?
    # figure out what trace should be
    this.call <- match.call()
    trace <- get.update.arg(this.call$trace, "trace", object, env,
                            trace1=NULL, "update.earth", print.trace=FALSE)
    trace <- eval.parent(trace)
    if(is.name(trace))    # TODO needed when called from earth_cv with glm=NULL, why?
        trace <- eval.parent(trace)
    if(is.null(trace))
        trace <- 0
    if(is.name(call$glm)) # TODO needed when called from earth_cv with glm=NULL, why?
        call$glm <- eval.parent(call$glm)
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots) > 0) {
        if(anyNA(pmatch(names(dots), prune.only.args)))
            do.forward.pass <- TRUE
        else if(!is.null(dots$nfold) || !is.null(call$nfold)) {
            trace1(trace,
                "update.earth: forcing forward pass because nfold argument used\n")
            do.forward.pass <- TRUE
        }
        existing <- !is.na(match(names(dots), names(call)))
        for(i in names(dots)[existing])     # replace existing args
            call[[i]] <- dots[[i]]
        if(any(!existing)) {                # append new args
            call <- c(as.list(call), dots[!existing])
            call <- as.call(call)
        }
    }
    if(is.null(call$formula)) {
        call[["x"]] <- get.update.arg(this.call[["x"]], "x", object, env, trace)
        call[["y"]] <- get.update.arg(this.call[["y"]], "y", object, env, trace)
    } else
        call$data <- get.update.arg(this.call$data, "data", object, env, trace)
    call$subset  <- get.update.arg(this.call$subset,  "subset",  object, env, trace)
    call$weights <- get.update.arg(this.call$weights, "weights", object, env, trace)
    call$wp      <- get.update.arg(this.call$wp,      "wp",      object, env, trace)
    if(check.boolean(ponly))
        do.forward.pass <- FALSE
    call$Object <- if(do.forward.pass) NULL else substitute(object)
    if(evaluate)
        eval.parent(call)
    else
        call
}

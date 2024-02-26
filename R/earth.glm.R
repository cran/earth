# earth.glm.R: Generalized Linear Model support for earth

try.something.like <-
    "Try something like earth(y~x, glm=list(family=binomial))"

check.glm.model <- function(g, resp.name) # g is a model created by calling glm()
{
    # TODO following check is pointless? df is only defined for summary.glm?
    df <- if("df" %in% names(g)) g[["df"]] else NULL # avoid partial matching
    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        stop0("earth glm response \"", resp.name, "\": ", nsingular,
              " coefficients not defined because of singularities")
    glm.coef <- coef(g)
    if(length(glm.coef) == 0)
        stop0("earth glm response \"", resp.name, "\": no glm coefficients")
    check.vec(glm.coef, "glm coef")
}
# Check for a common user error: specifying a family argument
# to earth that is not wrapped in glm=list(family=...))
# Jan 2019: also check that epsilon and maxit are properly enclosed in glm list

check.no.family.arg.to.earth <- function(..., is.null.glm.arg)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(!is.null(dots$fa)) # partial match
        stop0("illegal 'family' argument to earth\n", try.something.like)
    if(!is.null.glm.arg) {
        if(!is.null(dots$eps))
            stop0("illegal 'epsilon' argument to earth\n",
                  "Try something like earth(y~x, glm=list(family=binomial, control=list(epsilon=1e-9)))")
        if(!is.null(dots$maxi))
            stop0("illegal 'maxit' argument to earth\n",
                  "Try something like earth(y~x, glm=list(family=binomial, control=list(maxit=99)))")
    }
}
# Note that on entry process.glm.arg has already checked the glm argument
# Most args are direct copies of args to earth.fit.

earth_glm <- function(bx, y, weights, na.action, offset,
                      glm.arg, trace, is.bpairs, env)
{
    hack.intercept.only.glm.model <- function(g)
    {
        # Skullduggery for intercept-only glm models.
        # Get the model into a form usable by later functions like predict().
        # We need to remove all references to EarthIntercept else predict()
        # will try to find it and complain because it cannot.

        g$coefficients <- g$coefficients[1]
        g$R            <- g$R[1,1]
        g$qr$qr        <- g$qr$qr[,1,drop=FALSE]
        g$model        <- g$model[,1,drop=FALSE]
        g$x            <- g$x[,1,drop=FALSE]
        g$data         <- NULL

        # yarg ~ EarthIntercept becomes yarg ~ yarg
        # TODO this approach used because glm() won't allow just yarg~
        g$terms[[3]] <- g$terms[[2]]

        # list(yarg, EarthIntercept) becomes list(yarg)
        attr(g$terms, "variables") <- call("list", quote(yarg))

        #        EarthIntercept
        # yarg                0   becomes an empty matrix

        attr(g$terms, "factors") <- matrix(nrow=0, ncol=0)

        # "EarthIntercept" becomes an empty character vector
        attr(g$terms, "term.labels") <- character(0)

        # list(yarg, EarthIntercept) becomes list(yarg)
        attr(g$terms, "predvars") <- call("list", quote(yarg))

        g
    }
    #--- earth_glm starts here ---
    if(trace >= 4)
        cat("\n")
    ncases <- nrow(bx)
    intercept.only <- ncol(bx) == 1
    if(intercept.only) {
        # glm() requires something on the rhs of the formula.
        # But this is an intercept-only model, so actually nothing on the rhs.
        # To work around that, give glm() the earth intercept, which will have
        # no effect on the glm model but will cause an extra coefficient etc. in
        # the value returned by glm.  We remove that extra data later (in
        # hack.intercept.only.glm.model).
        # Actually the fake intercept does have a small effect on the
        # model: dof is off by one (which also affects vals derived from dof).
        trace1(trace, "earth_glm: intercept-only earth model\n")
        bx.data.frame <- as.data.frame(bx) # bx has a single column, the earth intercept
        colnames(bx.data.frame) <- "EarthIntercept" # for sanity checking
    } else {
        # default operation: drop intercept with -1
        bx.data.frame <- as.data.frame(eval.parent(bx[, -1, drop=FALSE]))
    }
    # Convert args to form expected by glm().
    # We need to convert glm() args whose default is not NULL.

    control <- glm.arg$control
    if(is.null(control))
       control <- glm.control()
    # FIXED (earth 2.3-5): get control params
    if(!is.null(glm.arg$epsilon))
       control$epsilon <- glm.arg$epsilon
    if(!is.null(glm.arg$maxit))
       control$maxit <- glm.arg$maxit
    if(!is.null(glm.arg$trace))
       control$trace <- glm.arg$trace
    family <- get.glm.family(glm.arg$family, env=env)
    stopifnot(is.null(glm.arg$weights))

    glm.list <- list()  # returned list of glm models
    non.converged <- NULL
    for(iycol in seq_len(ncol(y))) { # for each y column
        yarg <- if(is.bpairs) y[,1:2]       # two columns
                else y[, iycol, drop=FALSE] # single column
        stopifnot(!is.null(colnames(yarg)))
        if(trace >= 4) {
            print_summary(yarg, "glm y", trace=2)
            printf("\n")
            print_summary(weights, "glm weights", trace=2)
            printf("\n")
        }
        # FIXED (earth 2.3-4): removed offset etc. arguments because of
        # difficulties evaluating them later in the correct environment
        # (process.glm.arg has already checked if such args were supplied by the user).

        # TODO consider setting x=FALSE, y=FALSE if keepxy=-1
        g <- glm(yarg ~ ., family=family, data=bx.data.frame,
                 weights=weights, na.action=na.action, offset=offset,
                 control=control, model=TRUE, trace=(trace>=2),
                 method="glm.fit", x=TRUE, y=TRUE, contrasts=NULL)

        if(intercept.only)
            g <- hack.intercept.only.glm.model(g)

        # if !converged remember this response for warning later
        if(!g$converged)
            non.converged <- c(non.converged, colnames(yarg)[1])

        if(trace >= 1) {
            printf("GLM %s devratio %.2f dof %d/%d iters %d\n",
                colnames(yarg)[1],
                get.devratio(g$null.deviance, g$deviance),
                g$df.residual, g$df.null, g$iter)
        }
        check.glm.model(g, colnames(yarg)[1])
        glm.list[[iycol]] <- g
        if(is.bpairs) {
            stopifnot(NCOL(y) == 2) # paranoia
            break                   # note break
        }
    }
    if(!is.null(non.converged))
        warning0("the glm algorithm did not converge for response",
                 if(length(non.converged) == 1) " " else "s ",
                 quotify(non.converged))
    glm.list
}
# process family here instead of in glm() so can give relevant error message

get.glm.family <- function(family, env)
{
    if(is.null(family))
        family <- gaussian
    if(is.character(family))
        family <- get(family, mode="function", envir=env)
    if(is.function(family))
        family <- family()
    if(!inherits(family, "family") || is.null(family$family))
        stop0("earth: illegal 'family' in 'glm' argument\n",
              try.something.like)
    family
}
# This returns earth's glm argument but with abbreviated names
# expanded to their full name.  It also checks that the glm argument is
# valid. Called before calling glm().  We want to make sure that the
# user hasn't specified, say, subset as a glm argument. The subset
# should only be specified as an earth argument so the subset is the
# same for earth and glm.
# FIXED (earth 2.3-4): disallow offset etc. arguments because of
# difficulties evaluating them later in the correct environment.

process.glm.arg <- function(glm.arg) # glm.arg is earth's glm argument
{
    # return glm.arg but with abbreviated names expanded to their full name.
    match.glm.arg <- function(glm.arg)
    {
        allowed.glm.args <- c("formula", "family", "data", "weights", "subset",
            "na.action", "control", "model", "method", "x", "y",
            "contrasts", "epsilon", "maxit", "trace", "bpairs")
        for(i in seq_along(glm.arg)) {
            name <- names(glm.arg)[[i]]
            j <- pmatch(name, "bpairs", nomatch=0)
            if(j != 0)
                warning0("earth: the '", name, "' argument is no longer supported\n",
                      "         (binomial pairs are determined automatically)\n",
                      "         See comments for 'bpairs' in the earth NEWS file (earth version 5.0.0)")
            j <- pmatch(name, allowed.glm.args, nomatch=0)
            if(j == 0)
                stop0("earth: '", name, "' is not supported in glm argument to earth")
            if(allowed.glm.args[j] != "bpairs")
                names(glm.arg)[[i]] <- allowed.glm.args[j]
        }
        # expand family argument if it is a string

        if(is.character(glm.arg$family)) {
            family.strings <-
              c("binomial", "gaussian",  "Gamma",  "inverse.gaussian",
                "poisson",  "quasi",  "quasibinomial",  "quasipoisson")

            i <- pmatch(glm.arg$family, family.strings, nomatch=0)
            if(i == 0)
                stop0("earth: illegal family '", glm.arg$family, "' in glm argument\n",
                      try.something.like)
            glm.arg$family <- family.strings[i]
        }
        glm.arg
    }
    #--- process.glm.arg starts here ---

    if(is.null(glm.arg))
        return(NULL)
    if(!is.list(glm.arg))
        stop0("earth: 'glm' argument must be a list\n", try.something.like)
    if(length(glm.arg) == 0)
        stop0("earth: 'glm' argument list is empty\n", try.something.like)
    argnames <- names(glm.arg)
    if(length(argnames) == 0)
        stop0("earth: no argument names in 'glm' argument list\n", try.something.like)
    glm.arg <- match.glm.arg(glm.arg)  # expand argument names to their full name
    if(is.null(glm.arg$family))
        stop0("earth: 'glm' argument must have a 'family' parameter\n",
              try.something.like)

    always.true.args <- c("x", "y", "model")

    imatch <- pmatch(always.true.args, argnames, nomatch=0)
    if(any(imatch))
        stop0("earth: illegal '", argnames[imatch[1]], "' in 'glm' argument\n",
              "These are always effectively TRUE")

    earths.args <- c("subset", "weights") # illegal because these get passed on to glm internally
    imatch <- pmatch(earths.args, argnames, nomatch=0)
    if(any(imatch)) {
        imatch <- imatch[imatch != 0]
        stop0("earth: illegal '", argnames[imatch[1]], "' in 'glm' argument\n",
              "       Use earth's '", argnames[imatch[1]], "' argument instead ",
              "(which will be passed on to glm internally)")
    }
    earths.args <- c("formula") # this is plain illegal
    imatch <- pmatch(earths.args, argnames, nomatch=0)
    if(any(imatch)) {
        imatch <- imatch[imatch != 0]
        stop0("earth: illegal '", argnames[imatch[1]], "' in 'glm' argument\n",
              "       Use earth's '", argnames[imatch[1]], "' argument instead")
    }
    glm.arg
}
# get.glm.coefs returns a ncoeffs * nresponses matrix

get.glm.coefs <- function(glm.list, nresp, selected.terms, term.names, resp.names)
{
    coefs <- matrix(nrow=length(selected.terms), ncol=nresp)
    col.names <- character(length=nresp)
    for(iresp in seq_len(nresp)) {
        coefs[,iresp] <- glm.list[[iresp]]$coefficients
        col.names[iresp] <- resp.names[iresp]
    }
    colnames(coefs) <- col.names
    rownames(coefs) <- term.names[selected.terms]
    coefs
}
is.binomial <- function(family) # return true if family is binomial or quasibinomial
{
    (
        (is.character(family) &&            # e.g. "binomial"
            (substr(family, 1, 1) == "b" ||
             substr(family, 1, 6) == "quasib")) # "quasib" excludes "quasip" (quasipoisson)
        ||
        (class(family)[1] == "function" &&  # e.g. binomial
            (identical(body(family), body(binomial)) ||
             identical(body(family), body(quasibinomial))))
        ||
        (class(family)[1] == "family" &&    # e.g. binomial()
            (family$family == "binomial" ||
             family$family == "quasibinomial"))
    )
}
is.poisson <- function(family) # return true if family is poisson or quasipoisson
{
    (
        (is.character(family) &&            # e.g. "poisson"
            (substr(family, 1, 1) == "p" ||
             substr(family, 1, 6) == "quasip")) # "quasip" excludes "quasib" (quasibinomial)
        ||
        (class(family)[1] == "function" &&  # e.g. poisson
            (identical(body(family), body(poisson)) ||
             identical(body(family), body(quasipoisson))))
        ||
        (class(family)[1] == "family" &&    # e.g. poisson()
            (family$family == "poisson" ||
             family$family == "quasipoisson"))
    )
}
# called from print.summary.earth

print_earth_glm <- function(object, digits, fixed.point, prefix.space)
{
    glm.list <- object$glm.list
    nresp <- length(glm.list)
    # print the maxit if any glm model did not converge
    maxit.msg <- ""
    for(iresp in seq_len(nresp)) {
        g <- glm.list[[iresp]]
        if(!g$converged) {
            maxit.msg <- sprint(", maxit=%g", g$control$maxit)
            break # note break
        }
    }
    cat0(if(prefix.space) "   " else "",
         "GLM (family ", glm.list[[1]]$family$family, ", link ",
         glm.list[[1]]$family$link, maxit.msg, "):\n")
    stopifnot(!is.null(object$glm.stats))
    print(tweak.glm.stats.for.printing(object$glm.stats, digits,
                                       fixed.point, prefix.space),
          digits=max(3, digits-1))
    cat0("\n")
}
get.glm.stats <- function(glm.list, response.names)
{
    # not sure if all glm models have an aic field, so play safe with aic field
    aic <- glm.list[[1]]$aic[1]
    has.aic <- !is.null(aic) && !anyNA(aic)

    nresp <- length(glm.list)
    stats <- matrix(nrow=nresp, ncol=7+has.aic)

    colnames <- c("nulldev", "df", "dev", "df", "devratio")
    if(has.aic)
        colnames <- c(colnames, "AIC")
    colnames(stats) <- c(colnames, "iters", "converged")

    if(nresp == 1)
        rownames(stats) <- ""
    else
        rownames(stats) <- response.names

    for(iresp in seq_len(nresp)) {
        g <- glm.list[[iresp]]
        devratio <- get.devratio(g$null.deviance, g$deviance)
        if(has.aic)
            stats[iresp,] <- c(g$null.deviance, g$df.null,
                               g$deviance, g$df.residual, devratio,
                               g$aic, g$iter, g$converged)
        else
            stats[iresp,] <- c(g$null.deviance, g$df.null,
                               g$deviance, g$df.residual, devratio,
                               g$iter, g$converged)
    }
    stats
}
get.devratio <- function(null.deviance, deviance)
{
    devratio <- (null.deviance - deviance) / null.deviance
    # TODO not sure what best boundary cases handling is
    if(null.deviance < 0 || is.nan(devratio)) # division by zero
        return(NA)
    if(devratio < 1e-4)
        devratio <- 0 # prevent e.g. -0.0 and unsightly things like "3.45e-6"
    devratio
}
tweak.glm.stats.for.printing <- function(stats, digits, fixed.point, prefix.space)
{
    has.aic <- "AIC" %in% colnames(stats)
    stats[,"nulldev"] <- signif(stats[,"nulldev"], digits)
    stats[,"dev"] <- signif(stats[,"dev"], digits)
    stats[,"devratio"] <- signif(stats[,"devratio"], max(2, digits-4))
    if(has.aic) # signif() matches code in stats::print.glm
        stats[,"AIC"] <- signif(stats[,"AIC"], max(3, digits-3))
    # space out columns for readability
    colnames(stats) <-
        if(has.aic)
            c(paste0(if(prefix.space) "  " else "", "nulldev"),
              "df", "      dev", "df", "  devratio", "    AIC", "iters", "converged")
        else
            c(paste0(if(prefix.space) "  " else "", "nulldev"),
              "df", "      dev", "df", "  devratio", "   iters", "converged")
    stats
}
# Called from print.summary.earth
# g is a glm object
# Most of the following was lifted from print.summary.glm
# but tweaked to include response names (necessary for multiple
# response glm earth models).

print_glm_details <- function(g, nresp, digits, fixed.point, resp.name)
{
    if(nresp > 1)
        prefix <- paste("GLM", resp.name)
    else
        prefix <- paste("GLM")
    sumg <- summary(g)
    cat0(prefix, " deviance residuals:\n")
    if(sumg$df.residual > 5) {
        sumg$deviance.resid <- quantile(sumg$deviance.resid,na.rm=TRUE)
        names(sumg$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print.default(sumg$deviance.resid, digits=digits, na.print="", print.gap=2)
    df <- if("df" %in% names(sumg)) sumg[["df"]] else NULL
    cat0("\n", prefix)
    cat0(" coefficients (family ", g$family$family, ", link ", g$family$link, ")\n")
    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        cat0(nsingular,  # should never happen for earth glm
             " coefficients not defined because of singularities\n")
    aliased <- is.na(coef(g))
    stopifnot(length(aliased) > 0, all(!aliased)) # already checked in check.glm.model
    coefs <- sumg$coefficients
    rownames(coefs) <- spaceout(rownames(coefs))
    # TODO can't use fixed.point here, would like to
    printCoefmat(coefs, digits=digits, signif.stars=FALSE, na.print="NA")
    if(sumg$dispersion != 1) # only show dispersion if it is not 1
        cat0("\n", prefix, " dispersion parameter for ", sumg$family$family,
             " family taken to be ", sumg$dispersion, "\n")
    cat("\n")
    NULL
}

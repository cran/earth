# earth.glm.R: Generalized Linear Model support for earth

try.something.like <-
    "Try something like earth(y~x, glm=list(family=binomial))"

check.glm.model <- function(g, response.name) # g is a model created by calling glm()
{
    df <- if("df" %in% names(g)) g[["df"]] else NULL

    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        stop1("earth glm ", response.name, ": ", nsingular,
              " coefficients not defined because of singularities")

    aliased <- is.na(coef(g))
    if(length(aliased) == 0)
        stop1("earth glm ", response.name, ": no glm coefficients");
    if(any(aliased))
        stop1("earth glm ", response.name, ": NA in glm coefficients")
}

# Check for a common user error: specifying a family argument
# to earth that is not wrapped in glm=list(family=...))

check.no.family.arg.to.earth <- function(...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(!is.null(dots$fa)) # partial match
        stop1("illegal \"family\" argument to earth\n", try.something.like)
}

# This duplicates some tests in binomial in family.R for
# better error reporting in the earth context

check.yarg.for.binomial.glm <- function(yarg, weights, mustart, more.y.columns)
{
    if(ncol(yarg) == 1 && any(yarg < 0 | yarg > 1)) {
        cat("Error: binomial glm model with a vector y:",
            "y values must be between 0 and 1\n")
        if(more.y.columns)
            cat("Possible remedy: pair this column with",
                "the next column using the \"bpairs\" arg\n")
        cat("The first few rows of the y argument to glm are\n", sep="")
        print(head(yarg))
        stop1("glm with earth failed, see above messages")
    }
}

# Note that on entry get.glm.arg has already checked the glm argument
# All args except bx, glm, and glm.bpairs are direct copies of args to earth.fit.

earth.glm <- function(bx, y, weights, na.action, glm,
                      trace, glm.bpairs, response.name)
{
    if(trace >= 2)
        cat("\n")

    if(ncol(bx) == 1)
        stop1("earth cannot build a glm model because the earth model\n",
              "is an intercept only model.  This typically means that\n",
              "all terms except the intercept were deleted during pruning.\n")

    ncases <- nrow(bx)
    # -1 below to drop intercept, eval.parent assumes that this is called by earth.fit
    bx.data.frame <- as.data.frame(eval.parent(bx[,-1,drop=FALSE]))

    # Convert args to form expected by glm().
    # We need to convert glm() args whose default is not NULL.

    if(is.null(weights))
       weights <- rep(1, ncases)
    control <- glm$control
    if(is.null(control))
       control <- glm.control()
    # FIXED (earth 2.3-5): get control params
    if(!is.null(glm$epsilon))
       control$epsilon <- glm$epsilon
    if(!is.null(glm$maxit))
       control$maxit <- glm$maxit
    if(!is.null(glm$trace))
       control$trace <- glm$trace
    env <- parent.frame()
    family <- get.glm.family(glm$family, env=env)
    is.binomial <- is.binomial(family)

    # Fit a glm model for each y column.  Except that if there are
    # paired y columns then fit a single glm for each pair.
    # Note that we don't need to look at earth's wp argument here
    # because each glm model is built independently.

    iycol <- 1          # y column index
    imodel <- 1         # model index
    glm.list <- list()  # returned list of glm models

    while(iycol <= ncol(y)) {
        # get yarg, the response for call to glm()
        # it will be a single column or a a paired binomial column

        if(is.null(glm.bpairs))
            yarg <- y[, iycol, drop=FALSE]      # single y column
        else {
            if(!glm.bpairs[iycol])
                stop1(
"unmatched FALSE value in \"bpairs\" for y column ", iycol, "\n",
"       Each FALSE in \"bpairs\" should be preceded by a TRUE\n",
"       Your bpairs is ", paste.with.space(glm.bpairs))

            if(iycol + 1 <= ncol(y) && !glm.bpairs[iycol+1]) {
                yarg <- y[, c(iycol,iycol+1)]   # paired y columns
                iycol <- iycol + 1
            } else                              # single y column
                yarg <- y[, iycol, drop=FALSE]

        }
        if(is.binomial)
            check.yarg.for.binomial.glm(yarg, weights, glm$mustart, iycol < ncol(y))
        iycol <- iycol + 1
        stopif(is.null(colnames(yarg)))
        if(trace >= 4)
            print.matrix.info("y arg to glm()", yarg, all.names=TRUE)

        # FIXED (earth 2.3-4): removed offset etc. arguments because of
        # difficulties evaluating them later in the correct environment.

        g <- glm(yarg ~ ., family=family, data=bx.data.frame,
                weights=glm$weights, na.action=na.action,
                control=control, model=TRUE, trace=(trace>=2),
                method="glm.fit", x=TRUE, y=TRUE, contrasts=NULL)

        if(trace == 0 && !g$converged) # give a message specific to this reponse
            cat("earth glm ", response.name,
                ": did not converge after ", g$iter," iterations\n", sep="")
        if(trace >= 1) {
            cat("GLM ", colnames(yarg)[1], ": ", sep="")
            print.one.earth.glm(g, digits=getOption("digits"))
        }
        check.glm.model(g, colnames(yarg)[1])
        glm.list[[imodel]] <- g
        imodel <- imodel + 1
    }
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
    if(is.null(family$family))
        stop1("earth: illegal \"family\" in \"glm\" argument\n",
              try.something.like)
    family
}

# This returns the glm argument but with abbreviated names
# expanded to their full name.  It also checks that the glm argument is
# valid. Called before calling glm().  We want to make sure that the
# user hasn't specified, say, subset as a glm argument. The subset
# should only be specified as an earth argument so the subset is the
# same for earth and glm.
# FIXED (earth 2.3-4): disallow offset etc. arguments because of
# difficulties evaluating them later in the correct environment.

get.glm.arg <- function(glm)    # glm arg is earth's glm arg
{
   # return glm but with abbreviated names expanded to their full name.

    match.glm.arg <- function(glm)
    {
        glm.args <- c("formula", "family", "data", "weights", "subset",
            "na.action", "control", "model", "method", "x", "y",
            "contrasts", "epsilon", "maxit", "trace", "bpairs")

        for(i in 1:length(glm)) {
            j <- pmatch(names(glm)[[i]], glm.args, nomatch=0)
            if(j == 0)
                stop1("earth: \"",
                      names(glm)[[i]], "\" is not supported in glm argument to earth")
            names(glm)[[i]] <- glm.args[j]
        }
        # expand family argument if it is a string

        if(is.character(glm$family)) {
            family.strings <-
              c("binomial", "gaussian",  "Gamma",  "inverse.gaussian",
                "poisson",  "quasi",  "quasibinomial",  "quasipoisson")

            i <- pmatch(glm$family, family.strings, nomatch=0)
            if(i == 0)
                stop1("earth: illegal family \"", glm$family, "\" in glm argument\n",
                      try.something.like)
            glm$family <- family.strings[i]
        }
        glm
    }
    # get.glm.arg starts here

    if(!is.list(glm))
        stop1("earth: \"glm\" argument must be a list\n", try.something.like)
    if(length(glm) == 0)
        stop1("earth: \"glm\" argument list is empty\n", try.something.like)
    arg.names <- names(glm)
    if(length(arg.names) == 0)
        stop1("earth: no argument names in \"glm\" argument list\n", try.something.like)
    glm <- match.glm.arg(glm)  # expand argument names to their full name
    if(is.null(glm$family))
        stop1("earth: \"glm\" argument must have a \"family\" parameter\n",
              try.something.like)

    always.true.args <- c("x", "y", "model")

    imatch <- pmatch(always.true.args, arg.names)
    imatch <- imatch[!is.na(imatch)]
    if(any(imatch))
        stop1("earth: illegal \"", arg.names[imatch[1]],
              "\" in \"glm\" argument\n",
              "These are always effectively TRUE")

    # FIXED Oct 2008: removed "weights" from list below to
    # allow weights for glm list. Needed for example in
    # library(segmented); data(down); fit.e<-earth(cases/births~age,data=down,glm=list(family="binomial", weights=down$births))
    # If no weights, get warning: non-integer #successes in a binomial glm
    earths.args <- c("formula", "subset")
    imatch <- pmatch(earths.args, arg.names)
    imatch <- imatch[!is.na(imatch)]
    if(any(imatch)) {
        stop1("earth: illegal \"", arg.names[imatch[1]],
              "\" in \"glm\" argument\n",
              "Use earth's \"", arg.names[imatch[1]], "\" argument instead")
    }
    glm
}

# get.glm.coefs returns a ncoeffs * nresponses matrix

get.glm.coefs <- function(glm.list, nresp, selected.terms, term.names, response.names)
{
    coefs <- matrix(nrow=length(selected.terms), ncol=nresp)
    col.names <- character(length=nresp)
    for(iresp in 1:nresp) {
        coefs[,iresp] <- glm.list[[iresp]]$coefficients
        col.names[iresp] <- response.names[iresp]
    }
    colnames(coefs) <- col.names
    rownames(coefs) <- term.names[selected.terms]
    coefs
}

# Return a boolean vector saying which cols in y must be passed
# on to the C earth routine.
# If the user explicitly specified bpairs then we use that bpairs.
# Else we try to figure out bpairs automatically.
# Returns NULL if all y columns should be used
# (and returns NULL if family is not binomial).

get.glm.bpairs <- function(y, glm)
{
    bpairs <- glm$bpairs
    if(!is.null(bpairs)) {              # bpairs provided by user?
        if(ncol(y) == 1)
            stop1("\"bpairs\" argument is illegal because y has only one column")
        if(!is.binomial(glm$family))
            stop1("\"bpairs\" argument is illegal because the family ",
                  "is not binomial or quasibinomial")
        check.index.vec("bpairs", bpairs, y,
                        check.empty=TRUE, use.as.col.index=TRUE)
        bpairs <- to.logical(bpairs, NCOL(y))
    } else {
        bpairs <- rep(TRUE,ncol(y))
        if(is.binomial(glm$family)) {
            # If two adjacent columns both have values <0 or >1 then
            # assume that the columns are paired

            if(nrow(y) < 1)
                return(NULL) # later check will issue err msg for too short y
            i <- 1           # column number
            repeat {
               if(i + 1 > ncol(y))
                   break;
               if(any(y[,i] < 0 | y[,i] > 1) && any(y[,i+1] < 0 | y[,i+1] > 1)) {
                   bpairs[i+1] <- FALSE # columns are paired, so discard 2nd column
                   i <- i + 1
               }
               i <- i + 1
            }
        }
    }
    if(all(bpairs))
        NULL                                # all y columns used
    else
        bpairs
}

is.binomial <- function(family) # return true if family is binom or quasibinom
{
    (
        (is.character(family) &&            # e.g. "binomial"
            (substr(family, 1, 1) == "b" ||
             substr(family, 1, 6) == "quasib"))
        ||
        (class(family) == "function" &&     # e.g. binomial
            (identical(body(family), body(binomial)) ||
             identical(body(family), body(quasibinomial))))
        ||
        (class(family) == "family" &&       # e.g. binomial()
            (family$family == "binomial" ||
             family$family == "quasibinomial"))
    )
}

is.poisson <- function(family) # return true if family is poisson or quasipoisson
{
    (
        (is.character(family) &&            # e.g. "poisson"
            (substr(family, 1, 1) == "p" ||
             substr(family, 1, 6) == "quasip"))
        ||
        (class(family) == "function" &&     # e.g. poisson
            (identical(body(family), body(poisson)) ||
             identical(body(family), body(quasipoisson))))
        ||
        (class(family) == "family" &&       # e.g. poisson()
            (family$family == "poisson" ||
             family$family == "quasipoisson"))
    )
}

# called from print.summary.earth

print.earth.glm <- function(obj, digits, fixed.point)    # obj is an earth object
{
    glm.list <- obj$glm.list
    nresp <- length(glm.list)

    cat("\nGLM ")
    if(nresp == 1)
        print.one.earth.glm(glm.list[[1]], digits)
    else {                                  # create a matrix and print that
        cat("(family ", glm.list[[1]]$family$family, ", link ",
             glm.list[[1]]$family$link, ")\n", sep="")

        a <- matrix(nrow=nresp, ncol=6)
        colnames(a) <- c("null.deviance", "df", "deviance", "df", "iters", "converged")
        rownames(a) <- colnames(obj$fitted.values)
        for(iresp in 1:nresp) {
            g <- glm.list[[iresp]]
            a[iresp,] <- c(g$null.deviance, g$df.null,
                           g$deviance, g$df.residual,
                           g$iter, g$converged)
        }
        if(fixed.point)
            a <- my.fixed.point(a, digits)
        print(a, digits=digits)
    }
}

# Called from print.summary.earth
# g is a glm object
# Most of the following was lifted from print.summary.glm
# but tweaked to include response names (necessary for multiple
# response glm earth models).

print.glm.details <- function(g, nresp, digits, fixed.point, response.name)
{
    if(nresp > 1)
        prefix <- paste("GLM", response.name)
    else
        prefix <- paste("GLM")
    sumg <- summary(g)
    cat(prefix, " deviance residuals:\n", sep="")
    if(sumg$df.residual > 5) {
        sumg$deviance.resid <- quantile(sumg$deviance.resid,na.rm=TRUE)
        names(sumg$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print.default(sumg$deviance.resid, digits=digits, na.print="", print.gap=2)
    df <- if("df" %in% names(sumg)) sumg[["df"]] else NULL
    cat("\n", prefix, sep="")
    cat(" coefficients (family ", g$family$family, ", link ",
         g$family$link, ")\n", sep="")
    if(!is.null(df) && (nsingular <- df[3] - df[1]))
        cat(nsingular,  # should never happen for earth glm
            " coefficients not defined because of singularities\n", sep="")
    aliased <- is.na(coef(g))
    stopif(length(aliased) == 0 || any(aliased)) # already checked in check.glm.model
    coefs <- sumg$coefficients
    rownames(coefs) <- spaceout(rownames(coefs))
    # TODO can't use fixed.point here, would like to
    printCoefmat(coefs, digits=digits, signif.stars=FALSE, na.print="NA")
    if(sumg$dispersion != 1) # only show dispersion if it is not 1
        cat("\n", prefix, " dispersion parameter for ", sumg$family$family,
            " family taken to be ", format(sumg$dispersion), "\n", sep="")
    cat("\n")
    NULL
}

print.one.earth.glm <- function(g, digits)
{
    cat("null.deviance ",  format(g$null.deviance, digits=digits),
         " (",             g$df.null, " dof)",
         "   deviance ",   format(g$deviance, digits=digits),
         " (",             g$df.residual, " dof)",
         "   iters ",      g$iter,
         if(!g$converged) " did not converge" else "",
         "\n", sep="")
}

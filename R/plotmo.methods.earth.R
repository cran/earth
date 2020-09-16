# plotmo.rpart.R: plotmo methods for earth objects

plotmo.singles.earth <- function(object, x, nresponse, trace, all1, ...)
{
    get.earth.vars.for.plotmo(ndegree=1, object, x, nresponse, trace, all=all1, ...)
}
plotmo.pairs.earth <- function(object, x, nresponse=1, trace=0, all2=FALSE, ...)
{
    get.earth.vars.for.plotmo(ndegree=2, object, x, nresponse, trace, all=all2, ...)
}
get.earth.vars.for.plotmo <- function(ndegree, object, x, nresponse, trace, all, ...)
{
    modvars <- object$modvars
    stopifnot(!is.null(modvars))

    stopifnot(ndegree == 1 || ndegree == 2)
    stopifnot(nrow(modvars) == length(object$namesx))

    # default return is all singles or all pairs
    def.return <- get.earth.vars.def.return(ndegree, modvars, x)

    if(all) # user wants all used predictors, not just those in ndegree terms?
        return(def.return)

    if(ncol(x) < nrow(modvars)) {
        if(ndegree == 1) { # so issue only one warning per invocation of plotmo
            format <- paste0("Cannot determine which variables to plot (use all1=TRUE?)\n",
                             "             ncol(x) %d < nrow(modvars) %d\n",
                             "             colnames(x)=%s\n",
                             "             rownames(modvars)=%s")
            warnf(format,
                  ncol(x), nrow(modvars),
                  paste.c(colnames(x), maxlen=100),
                  paste.c(rownames(modvars), maxlen=100))
        }
        return(def.return)
    }
    dirs <- object$dirs[object$selected.terms, , drop=FALSE]
    degree <- get.degrees.per.term(dirs) == ndegree # rows in dirs for terms of ndegree
    if(all(degree == 0)) {
        return(NULL)                    # no terms of ndegree
    }

    # set intercept row to 0 (else we will plot the
    # intercept for degree1 and degree2 plots)
    # we don't plot the offset by default because doing so can
    # supersize the ylim (thus compressing the curves on the other plots)
    # TODO inconsistent with glm models where we do plot the offset
    if(any(modvars == 9999)) {
        modvars[modvars == 9999] <- 0
        if(ndegree == 1 && trace >= 0)
            cat0("Note: the offset in the formula is not plotted\n",
                 "      (use all1=TRUE to plot the offset, ",
                 "or use trace=-1 to silence this message)\n\n")
    }

    dirs <- dirs[degree, , drop=FALSE]  # rows in dirs for terms of ndegree
    stopifnot(ncol(modvars) == ncol(object$dirs))
    singles <- NULL

    for(irow in seq_len(nrow(dirs))) {
        vars <- which(dirs[irow, ] != 0)
        if(length(vars) != ndegree)
            warnf("get.earth.vars.for.plotmo ndegree %d: irow %d length(vars) %d\n",
                  ndegree, irow, length(vars))
        # ivar1 will be length 2 for earth terms like x1:x2, or x1:x2 * h(3-x3)
        # because the term uses two variables, x1 and x2
        # (there was a x1:x2 in the formula)
        ivar1 <- which(modvars[,vars[1]] != 0)
        if(ndegree == 2) {
            ivar2 <- which(modvars[,vars[2]] != 0)
            singles <- c(singles, generate.all.pairs(ivar1, ivar2))
        } else {
            singles <- c(singles, ivar1)
        }
    }
    single.names <- rownames(modvars)[singles]
    single.names <- gsub("\`", "", single.names) # remove backticks if any

    # hack for booleans which get expanded by model.matrix from "bool" to "boolTRUE"
    # this currently only affects caret models (see above comments for caret models)
    if(any(grepl(".+TRUE$", single.names)) && !any(grepl(".+TRUE$", colnames(x)))) {
        if(trace >= 2) {
            format <- paste0("get.earth.vars.for.plotmo: ",
                             "deleting \"TRUE\" in single.names=%s\n",
                             "                           to match colnames(x)=%s\n")
            printf(format,
                   paste.c(single.names, maxlen=100),
                   paste.c(colnames(x), maxlen=100))
        }
        single.names <- gsub("TRUE", "", single.names)
    }
    match <- match(single.names, colnames(x), nomatch=0)
    match <- matrix(match, ncol=ndegree, byrow=TRUE)

    # sanity checks
    if(any(match == 0) ||
            length(match) %% ndegree != 0 ||
            length(match) != length(single.names)) {
        if(ndegree == 1) { # so issue only one warning per invocation of plotmo
            if(trace >= 2) {
                cat("\nAdditional information for warning below:\n\n")
                printf("any(match == 0): %d\n", any(match == 0))
                printf("length(match) %% ndegree != 0: %d\n", length(match) %% ndegree != 0)
                printf("length(match) != length(single.names): %d\n",
                    length(match) != length(single.names))
                printf("(ndegree == 1 && length(unique(match)) != length(match)): %d\n",
                    (ndegree == 1 && length(unique(match)) != length(match)))
                cat("\nsingles:", singles, "\n")
                cat("\nmatch:\n")
                print(match)
                cat("\nhead(x):\n")
                print(head(x))
                cat("\nhead(object$dirs):\n")
                print(head(object$dirs))
                cat("\nhead(modvars):\n")
                print(head(modvars))
                cat("\n")
            }
            format <- paste0("Cannot determine which variables to plot (use all1=TRUE?)\n",
                             "             single.names=%s\n",
                             "             colnames(x)=%s\n")
            warnf(format,
                  paste.c(single.names, maxlen=100),
                  paste.c(colnames(x), maxlen=100))
       }
        match <- def.return
    } # end of sanity checks

    return(match) # plotmo will remove duplicate rows and "toggled dups" like c(2,5) and c(5,2)
}
# example 1: ivar=5      ivar2=8     return c(5,8)
# example 2: ivar=c(5,6) ivar2=8     return c(5,8, 6,8)
# example 3: ivar=c(5,6) ivar2=(8,9) return c(5,8, 5,9, 6,8, 6,9)

generate.all.pairs <- function(ivar1, ivar2)
{
    pairs <- NULL
    for(i in ivar1)
        for(j in ivar2)
            if(i != j) # necessary for terms like h(num-4)*sqrt(num) (both vars are "num")
                pairs <- c(pairs, c(i, j))
    pairs
}
# get the default return for get.earth.vars.for.plotmo (all singles or all pairs)

get.earth.vars.def.return <- function(ndegree, modvars, x)
{
    # The min (for all.singles) below is necessary if ncol(x) < nrow(modvars)
    # This need for min currently only affects caret models.
    # This can happen with caret models because plotmo gets the caret data
    # from the formula call to caret, but caret calls earth.default with x
    # from the model.matrix() already applied to the formula.
    # For example a formula like y~fac+num (two vars)
    # gets passed to earth x as fac1 fac2 fac3 num (four vars).
    # This only affects formulas that have factors (and logicals) in them, because
    # the columns names in x don't match the variable names in the formula.

    all.singles <- seq_len(min(nrow(modvars), ncol(x))) # all singles

    def.return <- all.singles # default return, used if something goes wrong

    if(ndegree == 2) {
        if(length(all.singles) <= 1)
            def.return <- matrix(0, nrow=0, ncol=ndegree)  # no pairs
        else {
            # plotmo_doubles will remove redundant "toggled duplicates" like c(2,5) and c(5,2)
            all.pairs <- generate.all.pairs(all.singles, all.singles)
            def.return <-  matrix(all.pairs, ncol=ndegree, byrow=TRUE) # all pairs
            # limit number of pair plots to 20 (arb)
            def.return <- def.return[seq_len(min(nrow(def.return), 20)), , drop=FALSE]
        }
    }
    def.return
}
plotmo.y.earth <- function(object, trace, naked, expected.len, ...)
{
    temp <- plotmo::plotmo.y.default(object, trace, naked, expected.len)

    # plotmo.y.default returns list(field=y, do.subset=do.subset)
    # do the same processing on y as earth does, e.g. if y is a two
    # level factor, convert it to an indicator column of 0s and 1s

    colnames <- colnames(temp$field)

    temp$field <- expand.arg(temp$field, model.env(object), trace, is.y.arg=TRUE,
                             name=if(!is.null(colnames)) colnames else "y")

    temp
}
plotmo.pairs.bagEarth <- function(object, x, ...) # caret package
{
    pairs <- matrix(0, nrow=0, ncol=2)
    for(i in seq_along(object$fit))
        pairs <- rbind(pairs, plotmo.pairs.earth(object$fit[[i]], x))
    pairs[order(pairs[,1], pairs[,2]),]
}
plotmo.y.bagEarth <- function(object, trace, naked, expected.len, ...)
{
    plotmo.y.earth(object, trace, naked, expected.len)
}
# back compatibility
get.plotmo.pairs.bagEarth <- function(object, env, x, trace, ...)
{
    plotmo.pairs.bagEarth(object, x, ...)
}
get.plotmo.y.bagEarth <- function(object, env, y.column, expected.len, trace, ...)
{
    plotmo.y.bagEarth(object, trace)
}

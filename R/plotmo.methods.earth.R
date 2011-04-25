# plotmo.rpart.R: plotmo methods for earth objects
# See the descriptions of the methods in plotmo:::plotmo.methods.R.

get.plotmo.singles.earth <- function(object, env, x, trace, all1)
{
    if(length(coef(object)) == 1)
        warning0("the earth model appears to be an intercept only model")
    singles <- NULL
    max.degree <- 1
    if(all1) # user wants all used predictors, not just those in degree1 terms?
        max.degree <- 99
    selected <- object$selected.terms[
                    reorder.earth(object, degree=max.degree, min.degree=1)]
    if(length(selected) > 0) {
        degree1.dirs <- object$dirs[selected, , drop=FALSE]
        # column numbers of dirs that have predictors in degree1 terms
        icol <- which(degree1.dirs != 0, arr.ind=TRUE)[,2]
        if(!any(sapply(x, is.factor)))  # no factors in x?
            singles <- icol
        else {                          # factors in x
            colnames <- colnames(object$dirs)[icol]
            prednames <- object$namesx.org
            for(ipred in seq_along(prednames)) {
                if(is.factor(x[,ipred])) {
                    # This knows how to deal with expanded factor names because
                    # it e.g. looks for "^pclass" in "pclass3rd"
                    # TODO this can give extra predictors if variable names alias
                    #      e.g. "x" and "x1" are both variable names
                    if(length(grep(paste0("^", prednames[ipred]), colnames)) > 0)
                        singles <- c(singles, ipred)
                } else if(prednames[ipred] %in% colnames)
                    singles <- c(singles, ipred)
            }
        }
    }
    singles
}
get.plotmo.pairs.earth <- function(object, env, x, trace, ...)
{
    pairs <- matrix(0, nrow=0, ncol=2)      # no pairs
    selected <- object$selected.terms[      # selected is all degree 2 terms
                    reorder.earth(object, degree=2, min.degree=2)]
    pairs <- vector(mode="numeric")
    for(i in selected)                      # append indices of the two preds in term i
        pairs <- c(pairs, which(object$dirs[i,] != 0))
    pairs <- unique(matrix(pairs, ncol=2, byrow=TRUE))
    if(nrow(pairs) > 0 && any(sapply(x, is.factor))) { # any columns in x are factors?
        # pairs works off expanded factor names, so replace each name
        # with index of original variable name
        # TODO this can give wrong results if variable names alias
        #      e.g. if "x" and "x1" are both variable names this takes the LAST
        #      of the matching names so correct with "x" "x1" but not "x1" "x"
        dir.colnames <- colnames(object$dirs)
        prednames <- object$namesx.org
        prednames.hat <- paste0("^", prednames)
        for(i in 1:nrow(pairs))
            for(j in 1:2) {
                ipred1 <- 0
                for(ipred in seq_along(prednames.hat))
                    if(length(grep(prednames.hat[ipred], dir.colnames[pairs[i, j]])) > 0)
                        ipred1 <- ipred
                if(ipred1 == 0)
                    stop0("internal error: illegal ipred1 in get.plotmo.pairs.earth")
                pairs[i, j] <- ipred1
            }
    }
    pairs
}
get.plotmo.y.earth <- function(object, env, y.column, expected.len, trace)
{
    y <- plotmo:::get.plotmo.y.default(object, env, y.column, expected.len, trace)

    # do the same processing on y as earth does, e.g. if y is a two
    # level factor, convert it to an indicator column of 0s and 1s

    y <- expand.arg(y, env, is.y.arg=TRUE, colnames(y))
    if(length(colnames(y)) == 1 && colnames(y) == "y")
        colnames(y) <- NULL # remove artificial colname added by expand.arg
    # TODO revisit y.column handling here
    if(!is.null(object$glm.list[[1]])) # if an earth.glm model, use y.column 1
        y.column <- 1

    list(y=y, y.column=y.column)
}
get.plotmo.pairs.bagEarth <- function(object, env, x, trace, ...)
{
    pairs <- matrix(0, nrow=0, ncol=2)
    for(i in 1:length(object$fit))
        pairs <- rbind(pairs,
                       get.plotmo.pairs.earth(object$fit[[i]], env, x, trace))
    pairs[order(pairs[,1], pairs[,2]),]
}
get.plotmo.y.bagEarth <- function(object, env, y.column, expected.len, trace)
{
    get.plotmo.y.earth(object, env, y.column, expected.len, trace)
}

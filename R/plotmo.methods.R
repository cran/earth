# plotmo.methods.R:  method functions for plotmo, like get.singles.xxx
#                    and get.pairs.xxx

#------------------------------------------------------------------------------
# plotmo.prolog gets called at the start of plotmo

plotmo.prolog <- function(object, object.name) UseMethod("plotmo.prolog")

plotmo.prolog.default <- function(object, object.name)
{
    # Here we just establish with some sort of plausibility that object is a model obj.
    # The general idea is to let the user know what is going on if plotmo fails later.

    if(!is.list(object))
        stop1("'", object.name, "' is not a model object")
    else if(length(coef(object)) == 1)
        warning1("'", object.name, "' appears to be an intercept only model")

    NULL
}
#------------------------------------------------------------------------------
# Return a vector of indices of predictors for degree1 plots
# The indices are col numbers in the x matrix
# The default method simply returns the indices of all predictors

get.singles <- function(object, x, degree1, pred.names, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.singles\n\n")
    UseMethod("get.singles")
}

get.singles.default <- function(object, x, degree1, pred.names, trace)
{
    ifirst <- if(pred.names[1]=="(Intercept)") 2 else 1 # delete intercept, if any
    if(is.character(degree1) && !is.na(pmatch(degree1, "all")))
        ifirst:length(pred.names) # degree1 = "all"
    else {
        check.index.vec("degree1", degree1, pred.names)
        (ifirst:length(pred.names))[degree1]
    }
}
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
get.singles.rpart <- function(object, x, degree1, pred.names, trace)
{
    # get all the variables that are actually used in the tree
    irow <- as.integer(row.names(object$frame))
    var.names <- character(length=max(irow))
    var.names[irow] <- as.character(object$frame$var)
    ivar <- charmatch(var.names, pred.names)
    is.split <- !is.na(ivar) & ivar > 0 # same as var.names != "<leaf>" & var.names !=""
    if(sum(is.split) == 0)
        stop1("the rpart tree has no splits")
    ivar <- unique(ivar[is.split])
    if(is.character(degree1) && !is.na(pmatch(degree1, "all")))
        degree1 <- 1:length(ivar)
    check.index.vec("degree1", degree1, ivar)
    ivar[degree1]
}
#------------------------------------------------------------------------------
# Each row of the returned Pairs matrix is the indices of two predictors
# for a degree2 plot
# The indices are col numbers in the x matrix
# See also get.pairs.earth

get.pairs <- function(object, x, degree2, pred.names, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.pairs\n\n")
    UseMethod("get.pairs")
}
# Predictors x1 and x2 are considered paired if they appear in the formula
# in forms such as x1:x2 or I(x1*x2) or s(x1,x2)

get.pairs.default <- function(object, x, degree2, pred.names, trace=FALSE)
{
    if(is.character(degree2) && !is.na(pmatch(degree2, "all"))) {
        # special case: degree2=="all"
        used.vars <- get.singles.default(object, x, "all", pred.names, trace)
        if(length(used.vars) == 0)
            return(matrix(0, nrow=0, ncol=2)) # no pairs
        col1 <- rep(used.vars, times=length(used.vars))
        col2 <- rep(used.vars, each=length(used.vars))
        Pairs <- cbind(col1, col2)
        Pairs <- Pairs[col1 != col2, , drop=FALSE]
        return(unique(t(apply(Pairs, 1, sort)))) # remove duplicate pairs
    }
    Pairs <- matrix(0, nrow=0, ncol=2)   # no pairs
    term.labels <- NULL
    if(!is.null(object$call$formula)) {
        form <- object$call$formula
        if(typeof(form) != "language")
            form <- eval.parent(form, n=2)
        form <- as.formula(form)
        data <- get.formula.data(object, object$call$data, FALSE, trace)
        terms <- terms(form, data=data)
    }
    if(!is.null(terms))
        term.labels <- attr(terms, "term.labels")
    if(!is.null(term.labels))
        Pairs <- get.pairs.from.term.labels(term.labels, pred.names, trace)
    else {
        if(trace > 0)
            cat("no degree2 plots because no $call$formula$term.labels\n")
        if(is.specified(degree2))
            warning1("'degree2' specified but no degree2 plots ",
                     "(because no $call$formula$term.labels)")
    }
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("'degree2' specified but no degree2 plots")
    if(nrow(Pairs) > 0) {
        check.index.vec("degree2", degree2, Pairs)
        Pairs <- Pairs[degree2, , drop=FALSE]
    }
    Pairs
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
            # Pairs works off expanded factor names, so replace each name
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
# We consider rpart variables paired if one is the direct parent of the
# other in the tree.

get.pairs.rpart <- function(object, x, degree2, pred.names, trace)
{
    if(is.character(degree2) && !is.na(pmatch(degree2, "all"))) {
        # TODO this may be redundant for rpart models
        used.vars <- get.singles.rpart(object, x, "all", pred.names, trace)
        if(length(used.vars) == 0)
            return(matrix(0, nrow=0, ncol=2)) # no pairs
        col1 <- rep(used.vars, times=length(used.vars))
        col2 <- rep(used.vars, each=length(used.vars))
        Pairs <- cbind(col1, col2)
        Pairs <- Pairs[col1 != col2, , drop=FALSE]
        return(unique(t(apply(Pairs, 1, sort)))) # remove duplicate pairs
    }
    irow <- as.integer(row.names(object$frame))
    var.names <- character(length=max(irow))
    var.names[irow] <- as.character(object$frame$var)
    ivar <- charmatch(var.names, pred.names)
    # following is the same as var.names != "<leaf>" & var.names !=""
    is.split <- !is.na(ivar) & ivar > 0
    if(sum(is.split) == 0)
        stop1("the rpart tree has no splits")
    Pairs <- NULL
    for(i in 1:length(ivar)) {
        if(is.split[i]) {
            left <- 2 * i
            if(is.split[left] && ivar[i] != ivar[left])
                Pairs <- c(Pairs, ivar[i], ivar[left])
            right <- left + 1
            if(is.split[right] && ivar[i] != ivar[right])
                Pairs <- c(Pairs, ivar[i], ivar[right])
        }
    }
    if(is.null(Pairs))
        Pairs <- matrix(0, nrow=0, ncol=2)
    else {
        Pairs <- matrix(Pairs, ncol=2, byrow=TRUE)
        Pairs <- unique(t(apply(Pairs, 1, sort))) # remove duplicate pairs
        if(nrow(Pairs) > 0) {
            check.index.vec("degree2", degree2, Pairs)
            Pairs <- Pairs[degree2, , drop=FALSE]
        }
    }
    if(nrow(Pairs) == 0 && is.specified(degree2))
        warning1("'degree2' specified but no degree2 plots")
    Pairs
}
#------------------------------------------------------------------------------
# define method function because don't want bogus warnings from plotmo.prolog.default
plotmo.prolog.bagEarth <- function(object) NULL

get.pairs.bagEarth <- function(object, x, degree2, pred.names, trace)
{
    pairs <- matrix(0, nrow=0, ncol=2)
    for(i in 1:length(object$fit))     # TODO could probably be vectorized
        pairs <- rbind(pairs, get.pairs(object$fit[[i]], x, degree2, pred.names, trace))
    pairs <- unique(pairs)
    pairs[order(pairs[,1], pairs[,2]),]
}
#------------------------------------------------------------------------------
# Given the term.labels, return an npairs x 2 matrix specifying
# which predictors are pairs. The elements in the returned matrix are col indices of x.
#
# It works like this: extract substrings from each term.label that look
# like predictor pairs and qualify them as a valid pair if both predictors
# in the pair are in pred.names.
#
# The following combos of x1 and x2 are considered pairs: "x1*x2" "x1:x2" "s(x1,x2)"
#
# This routine is not infallible but works for the commonly used formulas.

get.pairs.from.term.labels <- function(term.labels, pred.names, trace=TRUE)
{
    if(trace > 0)
        cat("term.labels:", term.labels, "\n")
    Pairs <- matrix(0, nrow=0, ncol=2)          # no pairs
    for(i in 1:length(term.labels)) {
        s <- strip.white.space(term.labels[i])
        s <- gsub("[+*/,]", ":", s)             # replace + * / , with :
        s <- gsub("=[^,)]+", "", s)             # delete "=any"

        # get the indices of expressions of the form "ident1:ident2"
        igrep <- gregexpr(
            "[a-zA-Z._][a-zA-Z._0-9$]*:[a-zA-Z._][a-zA-Z._0-9$]*", s)[[1]]

        if(trace > 0)
            cat("considering", s)

        if(igrep[1] > 0) for(i in seq_along(igrep)) {
            # extract the i'th "ident1:ident2" into Pair
            start <- igrep[i]
            stop <- start + attr(igrep, "match.length")[i] - 1
            Pair <- substr(s, start=start, stop=stop)
            Pair <- strsplit(Pair, ":")[[1]]    # Pair is now c("ident1","ident2")
            ipred1 <- which(pred.names == Pair[1])
            ipred2 <- which(pred.names == Pair[2])
            if(trace > 0)
                cat("->", Pair, "at", if(length(ipred1)) ipred1 else NA,
                    if(length(ipred2)) ipred2 else NA)
            if(length(ipred1) == 1 && length(ipred2) == 1 && Pair[1] != Pair[2])
                Pairs <- c(Pairs, ipred1, ipred2)
        }
        if(trace > 0)
            cat("\n")
    }
    unique(matrix(Pairs, ncol=2, byrow=TRUE))
}
#------------------------------------------------------------------------------
plotmo.predict <- function(object, newdata, type, se.fit,
                            pred.names, ipred1, ipred2, trace)
{
    if(trace > 0) {
        if(ipred2 == 0)
            cat("\nplotmo.predict for predictor \"", pred.names[ipred1], "\" ", sep="")
        else
            cat("\nplotmo.predict for predictors \"",
                pred.names[ipred1], "\" and \"", pred.names[ipred2], "\" ", sep="")
        if(se.fit)
            cat("se.fit=TRUE ")
        cat("with newdata[", NROW(newdata), ",", NCOL(newdata), "]:\n", sep="")
        print(head(newdata, 3))
    }
    UseMethod("plotmo.predict")
}
plotmo.predict.default <- function(object, newdata, type, se.fit,
                                   pred.names, ipred1, ipred2, trace)
{
    # We do our own error handling to give some context to the error messages.
    # TODO but this means that traceback() is not helpful
    try1 <- try(
        if(se.fit)
            predict(object, newdata=newdata, trace=trace, type=type, se.fit=TRUE)
        else
            predict(object, newdata=newdata, trace=trace, type=type))

    if(is.try.error(try1)) {
        cat("\n")
        stop1("Call to predict.", class(object)[1], " failed")
    }
    try1
}
plotmo.predict.lda <- function(object, newdata, type, se.fit,
                               pred.names, ipred1, ipred2, trace)
{
    y <- plotmo.predict.default(object, newdata, type, se.fit,
                                pred.names, ipred1, ipred2, trace)
    get.lda.yhat(y, type, trace)
}
plotmo.predict.qda <- function(object, newdata, type, se.fit,
                               pred.names, ipred1, ipred2, trace)
{
    y <- plotmo.predict.default(object, newdata, type, se.fit,
                                pred.names, ipred1, ipred2, trace)
    get.lda.yhat(y, type, trace)
}
#------------------------------------------------------------------------------
# Return the data matrix for the given object with the response deleted.
#
# If the model has a call$formula, the columns of the returned matrix are in the same
# order as the predictors in the formula.
#
# The default function tries hard to get x regardless of the model.
# Note that the alternative approach of simply calling the standard
# model.matrix wouldn't get us what we want here because it can return
# columns with headings like "ns(x3,4)" whereas we want the "naked" predictor x3.
#
# The n=2 and n=3 in the calls to eval.parent() take us to the caller of plotmo.
#
# TODO get.plotmo.x.default uses a nasty "formula string manipulation"
#      hack, must be a better way?

get.plotmo.x <- function(object, trace=FALSE)
{
    if(trace > 0)
        cat("\n--get.plotmo.x\n\n")
    UseMethod("get.plotmo.x")
}

get.plotmo.x.default <- function(
    object = stop("no 'object' arg"),
    trace  = FALSE)
{
    # get x by calling model.frame() with a stripped formula

    get.x.from.formula <- function(object, trace)
    {
        Call <- object$call
        if(is.null(Call))
            return(NULL)    # error will be reported later
        m <- match(c("formula", "data"), names(Call), 0)
        if(all(m == 0))
            return(NULL)
        Call <- Call[c(1, m)]
        Call[[1]] <- as.name("model.frame")
        # TODO it would be nice to use whatever na handling the original model function
        # used, but there seems to be no general way of knowing what that is.
        # In the meantime the following hack suffices for my purposes.
        Call$na.action <- if(inherits(object, "rpart")) na.pass else na.fail
        form <- Call$formula
        if(is.null(form))
            return(NULL)
        # following "if" is needed for: form <- Volume ~ .; earth(form, data=trees)
        # fixes bug reported by Martin Maechler and Michael Amrein
        if(typeof(form) != "language")
            form <- eval.parent(form, n=3)
        formula.as.string <- format(form)
        stripped.formula <- strip.formula.string(formula.as.string)
        if(trace > 0)
            cat("formula ", formula.as.string, "\n",
                "stripped formula ", stripped.formula, "\n", sep="")
        Call$formula <- parse(text=stripped.formula)[[1]]
        Call$data <- get.formula.data(object, Call$data, FALSE, trace)
        if(trace > 0)
            my.print.call("about to call ", Call)
        if(length(grep(".+~", stripped.formula))) { # has response?
            x <- try(eval.parent(Call, n=3)[,-1])   # then eval without the response
        } else {
            warning1("formula has no response variable, formula is ", stripped.formula)
            x <- try(eval.parent(Call, n=3))
        }
        if(is.try.error(x)) {
            if(length(grep("missing", x)))
                stop1("could not evaluate formula")
            else
                stop1("could not evaluate formula (variables were deleted?)")
        }
        if(NCOL(x) == 1) {
            # if one predictor, model.matrix returns a vec with no colname, so fix it
            x <- data.frame(x)
            colnames(x) <- strip.formula.string(attr(object$terms, "term.labels")[1])
        }
        x
    }
    badx <- function(x, check.colnames)
    {
        is.null(x) || is.try.error(x) || NROW(x) == 0 ||
            (check.colnames && is.null(colnames(x)))
    }
    # get.plotmo.x.default starts here

    if(!is.list(object))
        stop1("get.plotmo.x.default cannot get x matrix --- object is not a model object")

    # look for an x with column names

    try.error.message <- NULL
    # use of brackets rather than $ below prevents incorrect partial matching
    x <- object[["x"]]
    if(!badx(x, TRUE) && trace > 0)
        cat("got x with colnames from object$x\n")
    if(badx(x, TRUE)) {
        x <- get.x.from.formula(object, trace > 0)
        if(!badx(x, TRUE) && trace > 0)
            cat("got x with colnames from object$call$formula\n")
    }
    if(badx(x, TRUE)) {
        x <- try(eval.parent(object$call[["x"]], n=2), silent=TRUE)
        if(!badx(x, TRUE) && trace > 0)
            cat("got x with colnames from object$call$x\n")
        if(is.try.error(x))
            try.error.message <- x
    }
    # if don't have an x with colnames look for one without colnames
    # the call to as.data.frame below will add V1, V2, ... colnames if necessary

    if(badx(x, TRUE)) {
        x <- object[["x"]]
        if(!badx(x, FALSE) && trace > 0)
            cat("got x without colnames from object$x\n")
        if(badx(x, FALSE)) {
            x <- get.x.from.formula(object, trace)
            if(!badx(x, FALSE) && trace > 0)
                cat("got x without colnames from object$call$formula\n")
        }
        if(badx(x, FALSE)) {
            x <- try(eval.parent(object$call[["x"]], n=2), silent=TRUE)
            if(!badx(x, FALSE) && trace > 0)
                cat("got x without colnames from object$call$x\n")
            if(is.try.error(x))
                try.error.message <- x
        }
    }
    if(badx(x, FALSE)) {
        if(trace > 0) {
            cat("Looked unsuccessfully for an x in the following places:\n")
            cat("\nobject$x:\n")
            print(head(object$x, 3))
            cat("\nobject$call$formula:\n")
            dput(object$call$formula)
            cat("\nobject$call$x:\n")
            if(is.null(try.error.message))
                print(head(eval.parent(object$call$x, n=2), 3))
            else
                cat(gsub("Error in ", "", try.error.message[1]))
            cat("\n")
        }
        stop1("get.plotmo.x.default cannot get x matrix --- ",
              "tried object$x, object$call$formula, and object$call$x")
    }
    x <- as.data.frame(x)
    weights <- weights(object)
    if(!is.null(weights) && any(abs(weights - weights[1]) > 1e-8))
        warning1("'weights' are not supported by 'plotmo', ignoring them")

    subset. <- get.subset(object, trace > 0)
    if(!is.null(subset.)) {
        check.index.vec("subset", subset., x, check.empty=TRUE, allow.duplicates=TRUE)
        x <- x[subset., , drop=FALSE]
    }
    x
}
#------------------------------------------------------------------------------
# get.plotmo.y is similar to model.response but can deal with models
# created without a formula
# The n=2 and n=3 in the calls to eval.parent() take us to the caller of plotmo.

get.plotmo.y <- function(
    object = stop("no 'object' arg"),
    ycolumn,            # which column of response to use if response has multiple cols
    expected.len,
    trace)
{
    if(trace > 0)
        cat("\n--get.plotmo.y\n\n")
    UseMethod("get.plotmo.y")
}
get.plotmo.y.default <- function(
    object = stop("no 'object' arg"),
    ycolumn,            # which column of response to use if response has multiple cols
    expected.len,
    trace)
{
    get.y.from.formula <- function(object)
    {
        Call <- object$call
        if(is.null(Call))
            return(NULL)    # error will be reported later
        m <- match(c("formula", "data"), names(Call), 0)
        if(all(m == 0))
            return(NULL)
        Call <- Call[c(1, m)]

        # replace the formula with the stripped formula
        form <- Call$formula
        if(is.null(form))
            return(NULL)
        if(typeof(form) != "language")
            form <- eval.parent(form, n=3)
        formula.as.string <- paste(format(form), collapse=" ")
        stripped.formula <- strip.formula.string(formula.as.string)
        Call$formula <- parse(text=stripped.formula)[[1]]
        if(trace > 0)
            cat("formula ", formula.as.string, "\n",
                "stripped formula ", stripped.formula, "\n", sep="")

        Call$data <- get.formula.data(object, Call$data, TRUE, trace)
        Call[[1]] <- as.name("model.frame")
        Call$na.action <- if(inherits(object, "rpart")) na.pass else na.fail # TODO hack
        stripped.formula <- strip.formula.string(formula.as.string)
        if(trace > 0)
            my.print.call("about to call ", Call)
        Call <- try(eval.parent(Call, n=3))
        if(!is.try.error(Call))
            model.response(Call, type="any")
    }
    bady <- function(y)
    {
        is.null(y) || is.try.error(y)
    }
    # get.plotmo.y.default starts here
    try.error.message <- NULL
    y <- object[["y"]]
    if(!bady(y) && trace > 0)
        cat("got y from object$y\n")
    if(bady(y)) {
        y <- get.y.from.formula(object)
        if(!bady(y) && trace > 0)
            cat("got y from object$call$formula\n")
    }
    if(bady(y)) {
        y <- try(eval.parent(object$call[["y"]], n=2), silent=TRUE)
        if(!bady(y) && trace > 0)
            cat("got y from object$call$y\n")
        if(is.try.error(y))
            try.error.message <- y
    }
    if(bady(y)) {
        if(trace > 0) {
            cat("Looked unsuccessfully for y in the following places:\n")
            cat("\nobject$y:\n")
            print(head(object$y, 3))
            cat("\nobject$call$formula:\n")
            dput(object$call$formula)
            cat("\nobject$call$y:\n")
            if(is.null(try.error.message))
                print(head(eval.parent(object$call$y, n=2), 3))
            else
                cat(gsub("Error in ", "", try.error.message[1]))
            cat("\n")
        }
        stop1("get.plotmo.y.default cannot get y --- ",
              "tried object$call$formula, object$call$y, and object$y")
    }
    subset. <- get.subset(object, trace)
    y <- check.and.print.y(y, "get.plotmo.y", ycolumn, expected.len, trace, subset.)
    if(!is.null(subset.)) {
        check.index.vec("subset", subset., y, check.empty=TRUE, allow.duplicates=TRUE)
        y <- y[subset.]
    }
    y
}
get.formula.data <- function(object, data.name, get.y, trace)
{
    data.is.good <- function(...)
    {
        # the length test is necessary for lm which saves x as an
        # empty list if its x arg is FALSE, don't know why
        good <- !is.null(x) && length(x)
        if(good && trace > 0)
            cat("get.formula.data: got",
                if(get.y) "y" else "x", "from", paste(..., sep=""), "\n")
        good
    }
    xname <- if(get.y) "y" else "x"
    # use of brackets rather than $ below prevents incorrect partial matching
    x <- object[["data"]]
    if(!data.is.good("object$data")) {
        x <- object[[xname]]
        if(!data.is.good("object$", xname)) {
            if(!is.null(data.name)) {
                x <- eval.parent(data.name, n=4)
                if(!data.is.good("data passed to original call to ", class(object)))
                    stop1("the original data \"", data.name,
                          "\" is no longer available",
                          if(inherits(object, "earth"))
                                " (use keepxy=TRUE)"
                          else if(inherits(object, "lm"))
                                paste(" (use ", xname, "=TRUE)", sep="")
                          else  "")
            }
        }
    }
    x
}
get.subset <- function(object, trace) # called by get.plotmo.x.default and get.plotmo.y.default
{
    subset. <- object$subset
    if(is.null(subset.)) {
        # the n=3 takes us to the caller of plotmo
        subset. <- try(eval.parent(object$call$subset, n=3), silent=TRUE)
        if(is.try.error(subset.))
            subset. <- NULL
        #TODO revisit the following, it converts function (x, ...) UseMethod("subset")
        else if(typeof(subset.) == "closure")
            subset. <- NULL
    }
    if(!is.null(subset.) && trace > 0) {
        cat("subset length " , length(subset.), sep="")
        try(cat(" min", min(subset.), "max", max(subset.)), silent=TRUE)
        cat(" values ")
        for(i in 1:min(10, length(subset.)))
            cat(subset.[i], "")
        cat("...\n")
    }
    subset.
}
# Given a formula (as string), return a string with the "naked" predictors.
#
# Example: y ~ x9 + ns(x2,4) + s(x3,x4,df=4) + x5:sqrt(x6)
# becomes: y ~ x9 + x2 + x3 + x4 + x5 + x6
# which will later result in a model.matrix with columns x9 x2 x3 x4 x5 x6.
#
# This routine is not infallible but works for the commonly used formulas.

strip.formula.string <- function(form)
{
    gsubi <- function(pat, rep, x) gsub(pat, rep, x, ignore.case=TRUE)

    igrep <- grep("[a-zA-Z._][a-zA-Z._0-9]\\$", form)   # check for "ident$"
    if(length(igrep) > 0) {
        # TODO formula has vars with $, this confuses predict() later, why?
        # they cause "Warning: after calling plotmo.predict, y has the wrong length"
        stop1("plotmo: names with \"$\" are not yet supported\n",
            "The unacceptable formula is ", form)
    }
    form <- strip.white.space(form)
    args <- gsubi(".*~", "", form)                  # extract everything after ~

    # We don't want to mess with anything in [square brackets]
    # Doing that with regular expressions is tricky, so we adopt this approach:
    # change "+" "-" "," in square brackets to #PLUS# #MINUS# #COMMA# to protect them

    args <- gsubi("(\\[.*)\\+(.*\\])", "\\1#PLUS#\\2", args)
    args <- gsubi("(\\[.*)\\-(.*\\])", "\\1#MINUS#\\2", args)
    args <- gsubi("(\\[.*)\\,(.*\\])", "\\1#COMMA#\\2", args)

    args <- gsubi("[-*/:]", "+", args)              # replace - / * : with +

    # next two gsubs allow us to retain "x=x1" but drop "df=99" from "bs(x=x1, df=99)"

    args <- gsubi("\\([a-z._0-9$]+=", "(", args)    # replace "(ident=" with "("
    args <- gsubi("[a-z._0-9$]+=[^,)]+", "", args)  # delete "ident=any"

    # replace ",ident" with ")+f(ident", thus "s(x0,x1)" becomes "s(x0)f(x1)"

    args <- gsubi(",([a-z._])", ")+s(\\1", args)

    args <- gsubi("[a-z._0-9$]*[(]", "", args)      # remove "ident("
    args <- gsubi("[,)][^+-]*", "", args)           # remove remaining ",arg1,arg2)"

    # change #PLUS# etc. back to what they where
    args <- gsubi("#MINUS#", "-", args)
    args <- gsubi("#PLUS#", "+", args)
    args <- gsubi("#COMMA#", ",", args)

    # workaround for "error: invalid type (list) for variable 'trees[, -3]'"
    # for a <- earth(trees[,3] ~ as.matrix(trees[,-3])); plotmo(a)
    #
    # TODO removed because although it fixes that problem we still get
    # Warning: 'newdata' had 10 rows but variable(s) found have 31 rows
    #
    # if(is.list(eval.parent(parse(text=args, n=4))))
    #   args<-paste("as.matrix(", args, ")")

    response <- ""
    if(length(grep("~", form)))                     # if there is a response
        response <- gsubi("~.*", "~", form)         # then extract all before ~

    # FIXED 7 Dec 2007 reported by Joe Retzer
    # collapse possible multiple element response and args into a single string

    strip.white.space(paste(response[1], paste(args, collapse=" "), collapse=" "))
}

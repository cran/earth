# model.matrix.earth.R: functions for diddling around with model matrices
#
# The main functions are:
#
# expand.arg(x, env, is.y.arg) expand factors in x and convert to double mat with col names
#      Called by earth.formula, earth.default, get.earth.x
#
# model.matrix(terms, data) standard R function to expand factors
#     Called by expand.arg earth.formula, get.earth.x, predict.earth(type="terms")
#
# model.matrix.earth(object, x, ...) x arg must not be expanded, returns bx
#     Called by predict.earth
#
# get.earth.x(object, data) returns x expanded for factors and all double
#      Called by model.matrix.earth
#
# get.bx(x, which.terms, dirs, cuts) x arg must be already expanded
#      Called by model.matrix.earth, pruning.pass
#-----------------------------------------------------------------------------

# If model.frame can't interpret the data passed to it it silently
# returns the fitted values.  This routine makes that not silent.

check.nrows <- function(expected.nrows, actual.nrows, fitted.nrows, Callers.name)
{
    if(actual.nrows != expected.nrows) {
        if(actual.nrows == fitted.nrows)
            stop1("model.frame.default could not interpret the data passed to ",
                  Callers.name)
        else  # can probably never get here
            warning1(Callers.name, " returned a number ", actual.nrows,
                    " of rows that was different to the number ",
                    expected.nrows,
                    " of rows in the data")
    }
}

# Generate a column name for each column in x.
# x could be a vector.
# if xname is specified and x is a vector then use xname
# else use the existing colnumn names where possible.

generate.colnames <- function(x, is.y.arg=FALSE, xname=NULL)
{
    good.name <- function(x)
    {
        # The nchar check prevents super long names
        # that are actually contents of vectors e.g. c(1,2,3,etc.)
        # The grep ensures that there are no more than three commas,
        # also to prevent using the contents of vectors.

        !is.null(x) &&
        nchar(x) <= 60 &&
        length(grep(",.*,.*,", x)) == 0
    }
    names. <- rep("", length=NCOL(x))

    if(NCOL(x) == 1 && good.name(xname[1]))
        names. <- xname
    else {
        # copy valid names in colnames(x) to names

        col.names <- colnames(x)
        for (i in 1:NCOL(x))
            if (good.name(col.names[i]))
                names.[i] <- col.names[i]
    }
    # if any name is still "", convert it to an "xN" style name

    basename <- if(is.y.arg) "y" else "x"
    which. <- which(names. == "")
    if(any(which.))
        if(length(names.) == 1)
            names. <- basename
        else
            names.[which.] <- paste(basename, (1:ncol(x))[which.], sep="")

    make.unique(strip.white.space(names.))
}

# Called from earth.fit just before doing the pruning pass
# Also called by model.matrix.earth (which returns bx)
# The x arg must be already expanded

get.bx <- function(x, which.terms, dirs, cuts)
{
    stopifnot(all(dirs[1,] == 0))   # intercept term dirs must all be 0
    check.which.terms(dirs, which.terms)
    stopifnot(NCOL(x) > 0)
    bx <- matrix(0, nrow=nrow(x), ncol=length(which.terms))
    ibx <- 1
    for(iterm in which.terms) {
        temp1 <- 1
        for(ipred in 1:ncol(x)) {
            if(dirs[iterm, ipred] == 2)  # predictor enters linearly?
                temp1 <- temp1 * x[, ipred]
            else if(dirs[iterm, ipred] == -1 || dirs[iterm, ipred] == 1) {
                temp2 <- dirs[iterm, ipred] * (x[, ipred] - cuts[iterm, ipred])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            } else if(dirs[iterm, ipred] != 0)
                stop1("illegal direction ", dirs[iterm, ipred], " in \"dirs\"")
        }
        bx[, ibx] <- temp1
        ibx <- ibx + 1
    }
    colnames(bx) <- rownames(dirs[which.terms,])
    bx
}

# Called only by model.matrix.earth

get.earth.x <- function(    # returns x expanded for factors
    object  = stop("no 'object' arg"),
    data    = NULL,         # can be a dataframe, matrix, or vector
    env,                    # environment for evaluation
    trace   = 0,
    Callers.name)           # caller's name for trace messages
{
    # Given an x mat, return an x mat with column names equal to
    # expected.colnames and with the columns in their correct order
    # So the user can hand us an x mat without column names,
    # or with named columns but in the wrong order, or an xmat containing
    # only a needed subset of all the columns, etc.
    # This code is a mess and doesn't handle all cases, just the common ones.

    fix.x.columns <- function(x, expected.colnames)
    {
        col.names <- colnames(x)
        ncol.names <- length(col.names)
        nexpected <- length(expected.colnames)
        if(is.null(col.names)) {
            if(trace)
                cat(Callers.name, ": x has no column names, ",
                    "adding column names: ", paste.with.space(expected.colnames),
                    "\n", sep="")
            col.names <- expected.colnames
         } else if(ncol.names < nexpected) {
            # CHANGED Oct 2008: allow user to specify less than the expected
            # nbr of columns -- which is ok if he specifies all predictors
            # actually used by the model.

            imatch <- pmatch(col.names, expected.colnames, nomatch=0)
            if (any(imatch == 0)) {
                # can't repair the error because there are colnames in x that aren't
                # in expected.names (tends to happen with expanded factor names)
                stop1(Callers.name, ": x has ", ncol.names,
                      " columns, expected ", length(expected.colnames),
                      " to match: ", paste.with.space(expected.colnames))
            }
            # Create a new x, putting the existing cols into their correct positions.
            # Cols that aren't in the original x will end up as all 999s in the
            # the recreated x; that doesn't matter if they are for predictors
            # that are unused in the earth model.
            # We use 999 instead of NA else the model matrix routines fail incorrectly
            # because they don't like NAs when doing predictions.
            if(trace)
                cat(Callers.name, ": x has missing columns, ",
                    "creating a new x with all cols\n", sep="")
            imatch <- pmatch(expected.colnames, col.names, nomatch=0)
            x.original <- x
            x <- matrix(data=999, nrow=nrow(x), ncol=nexpected)
            for (i in 1:nexpected)
                if (imatch[i])
                    x[,i] <- x.original[,imatch[i]]
            col.names <- expected.colnames
        } else if(ncol.names > nexpected) {
            # TODO not sure what to do here (do nothing so old regression tests pass)
        } else {
            imatch <- pmatch(col.names, expected.colnames, nomatch=0)
            if(all(imatch == 0)) {
                   if(trace)
                        cat(Callers.name,
                            ": unexpected x column names, renaming columns\n",
                            "    Old names: ", paste.with.space(col.names), "\n",
                            "    New names: ", paste.with.space(expected.colnames),
                            "\n", sep="")

                    col.names <- expected.colnames
            } else {
                 # replace indices for non-found predictor names with their value
                 # i.e. assume columns with unknown names are in their right position
                 for(i in 1:length(imatch))
                     if(imatch[i] == 0)
                         imatch[i] = i

                 # if any columns are in the wrong order then fix their order
                 # (imatch will be 1,2,3,... if columns are in the right order)

                 if(!all(imatch == seq_along(imatch))) {
                   s <- paste(Callers.name, ": x columns are in the wrong order%s\n",
                                "    Old columns: ", paste.with.space(col.names), "\n",
                                "    New columns: ", paste.with.space(expected.colnames),
                                "\n", sep="")
                   if(length(imatch) == ncol(x)) {
                       if(trace)
                           cat(sprintf(s, ", correcting the column order"), sep="")
                       x <- x[,imatch]
                       col.names <- col.names[imatch]
                   } else
                       warning1(sprintf(s, ""))
                }
            }
        }
        colnames(x) <- col.names
        x
    }
    # Return x with matrix dimensions.
    # If x is already a matrix this does nothing.
    # Else allow x to be a vector if its length is an integer
    # multiple of the number of columns in the original x

    my.as.matrix <- function(x)
    {
        if(is.null(ncol(x))) {
            nrows <- length(x) / length(object$namesx)
            if(floor(nrows) == nrows)
                dim(x) <- c(nrow=nrows, ncol=length(x) / nrows)
            else
                stop1(Callers.name,
                    ": could not convert vector x to matrix because\n",
                    "       length(x) ", length(x),
                    " is not a multiple of the number ", ncol(object$namesx),
                    " of predictors ",
                    "\n       Expected predictors: ",
                    paste.with.space(colnames(object$namesx)))
        }
        x
    }
    check.expanded.ncols <- function(x, object)
    {
        if(ncol(x) != ncol(object$dirs)) {
            stop1(Callers.name,
                     ": the number ", NCOL(x),
                     " of columns of x after factor expansion\n",
                     "does not match the number ", NCOL(object$dirs),
                     " of columns of the earth object",
                     "\n  expanded x:  ", paste.with.space(colnames(x)),
                     "\n  object$dirs: ", paste.with.space(colnames(object$dirs)),
                     "\nPossible remedy: check factors in the input data")
        }
    }
    # get.earth.x starts here

    trace <- get.update.arg(trace, "trace", object,
                            trace1=NULL, Callers.name, print.trace=FALSE)
    if(is.null(trace))
        trace <- 0
    this.call <- match.call(expand.dots=TRUE)
    if(is.null(object$terms)) {
        # object was created with earth.default, no formula

        x <- get.update.arg(data, "x", object, trace, Callers.name)
        x <- my.as.matrix(data)
        x <- fix.x.columns(x, object$namesx)
        if(trace)
            print.matrix.info("x", x, Callers.name)
        x <- expand.arg(x, env)
    } else {
        # object was created with earth.formula

        Terms <- delete.response(object$terms)
        data <- get.update.arg(data, "data", object, trace, Callers.name)
        data <- my.as.matrix(data)
        data <- fix.x.columns(data, object$namesx)
        data <- as.data.frame(data)
        expected.nrows <- nrow(data)
        if(trace)
            print.matrix.info("x", data, Callers.name)
        data <- model.frame(Terms, data)
        classes <- attr(Terms, "dataClasses")
        if(!is.null(classes)) {
             # use "try" to be lenient, allow numeric to be used for factors etc.
             z <- try(.checkMFClasses(classes, data), silent=FALSE)
             if(class(z) == "try-error") {
                 # error msg already printed by .checkMFClasses
                 cat("Forging on regardless, first few rows of x are\n")
                 print(head(data))
             }
        }
        x <- model.matrix(Terms, data)
        check.nrows(expected.nrows, nrow(x), nrow(object$fitted.values), Callers.name)
        intercept <- match("(Intercept)", colnames(x), nomatch=0)
        if(intercept)
            x <- x[, -intercept, drop=FALSE]    # silently discard intercept
    }
    if(nrow(x) == 0)
        stop1("empty model matrix")
    check.expanded.ncols(x, object)
    x
}

# Called by update.earth and get.earth.x
#
# Which x should we use? The precedence is [1] the x parameter, if any,
# in this call to update [2] the $x in the earth object (which exists
# if keepxy=TRUE was used the original call to earth) [3] the x found
# in the original call to earth.
# Same applies for y, subset, weights, and wp.
# The "arg" argument is from the current call to update or predict

get.update.arg <- function(arg, argname, object,
                           trace1, Callers.name="update.earth", print.trace=TRUE)
{
    if(!print.trace) # print.trace arg prevents recursion issues with trace
        trace1 = FALSE
    if(is.null(arg)) {
        temp. <- try(eval.parent(object[[argname, exact=TRUE]], n=2), silent=TRUE)
        if(!is.null(temp.) && class(temp.) != "try-error") {
            arg <- object[[argname, exact=TRUE]]
            if(trace1)
                cat(Callers.name, ": using ",
                    NROW(temp.), " by ", NCOL(temp.), " ", argname,
                    " saved by keepxy in original call to earth\n", sep="")
        } else {
            temp. <- try(eval.parent(object$call[[argname, exact=TRUE]], n=2), silent=TRUE)
            if(!is.null(temp.) && class(temp.) != "try-error") {
                arg <- object$call[[argname, exact=TRUE]]
                if(trace1)
                    cat(Callers.name, ": using ",
                        NROW(temp.), " by ", NCOL(temp.), " ", argname,
                        " argument from original call to earth\n", sep="")
             }
        }
    }
    arg
}

# Called by predict.earth and can also be called by users directly.
# Called object$bx if all x, subset, which.terms equal NULL.

model.matrix.earth <- function(     # returns bx
    object       = stop("no 'object' arg"),
    x            = NULL,            # x arg must not yet be expanded
    subset       = NULL,
    which.terms  = NULL,
    ...,                                 # unused, for generic method comparibility
    env          = parent.frame(),
    trace        = FALSE,
    Callers.name = "model.matrix.earth") # caller's name for trace messages
{
    warn.if.dots.used("model.matrix.earth", ...)
    check.classname(object, deparse(substitute(object)), "earth")
    trace <- check.trace.arg(trace)
    if(is.null(x) && is.null(subset) && is.null(which.terms)) {
        if(trace)
            cat(Callers.name, ": returning object$bx\n", sep="")
        return(object$bx)
    }
    x <- get.earth.x(object, data=x, env, trace, Callers.name)
    if(is.null(which.terms))
        which.terms <- object$selected.terms
    if(!is.null(subset)) {
        check.index.vec("subset", subset, x,
                        check.empty=TRUE, allow.duplicates=TRUE, allow.zeroes=TRUE)
        x <- x[subset, , drop=FALSE]
    }
    get.bx(x, which.terms, object$dirs, object$cuts)
}

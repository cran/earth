# model.matrix.earth.R: Functions for manipulating earth model matrices
#
# The main functions are:
#
# expand.arg(x, env, is.y.arg, name)    in expand.arg.R (not in this module)
#
#     Expand factors in x and convert to double mat with col names
#     Called by earth.formula, earth.default, get.earth.x
#
#
# stats::model.matrix(terms, data)      standard R function to expand factors
#
#     Called by expand.arg earth.formula, get.earth.x, predict.earth(type="terms")
#
#
# model.matrix.earth(object, x, ...)    x arg must not be expanded, returns bx
#
#     Called by predict.earth
#
#
# get.earth.x(object, data) returns     returns x expanded for factors and all double
#
#      Called by model.matrix.earth
#
#
# get.bx(x, which.terms, dirs, cuts)    x arg must be already expanded
#
#      Called by model.matrix.earth, pruning.pass
#
#-----------------------------------------------------------------------------

# Called from earth.fit just before doing the pruning pass
# Also called by model.matrix.earth (which returns bx)
# The x arg must be already expanded

get.bx <- function(x, which.terms, dirs, cuts)
{
    stopifnot(all(dirs[1,] == 0))   # intercept term dirs must all be 0
    check.which.terms(dirs, which.terms)
    stopifnot(NCOL(x) > 0)
    colnames <- rownames(dirs[which.terms,,drop=FALSE])
    bx <- matrix(0, nrow=nrow(x), ncol=length(which.terms),
                 dimnames=list(NULL, colnames))
    ibx <- 1
    for(iterm in which.terms) {
        temp1 <- 1
        for(ipred in seq_len(ncol(x))) {
            dir <- dirs[iterm, ipred]
            if(dir == 2)  # predictor enters linearly?
                temp1 <- temp1 * x[, ipred]
            else if(dir == -1 || dir == 1) {
                temp2 <- dir * (x[, ipred] - cuts[iterm, ipred])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            } else if(dir != 0)
                stop0("illegal direction ", dir, " in 'dirs'")
        }
        bx[, ibx] <- temp1
        ibx <- ibx + 1
    }
    bx
}
# called only by model.matrix.earth (used to generate a bx matrix)
# returns x expanded for factors
# data can be a dataframe, matrix, or vector

get.earth.x <- function(object, data=NULL, env, trace=0, Callers.name)
{
    trace <- get.update.arg(trace, "trace", object, env,
                            trace1=NULL, Callers.name, print.trace=FALSE)
    if(is.null(trace))
        trace <- 0

    this.call <- match.call()

    if(is.null(object$terms)) # model was created with earth.default, no formula?
        x <- get.earth.x.default(object, data, env, trace, Callers.name)
    else                      # model was created with earth.formula
        x <- get.earth.x.formula(object, data, env, trace, Callers.name)

    if(NROW(x) == 0)
        stop0("empty model matrix")

    # Fix: April 2010, allow earth to play nicely with fda with factors in x
    if(ncol(x) > ncol(object$dirs)) # too many columns?
        x <- x[, colnames(x) %in% colnames(object$dirs), drop=FALSE] # select only the columns in dirs

    check.expanded.ncols(x, object)

    x
}
# object was created with earth.default, no formula
# called only by get.earth.x

get.earth.x.default <- function(object, data, env, trace, Callers.name)
{
    x <- get.update.arg(data, "x", object, env, trace, Callers.name)
    namesx <- rownames(object$modvars)
    x <- possibly.convert.vector.to.matrix(x, namesx)
    # following allows data to be a list e.g. newdata=etitanic[1,,drop=TRUE]
    x <- possibly.convert.list.to.data.frame(x)
    x <- fix.newdata.cols(x, namesx, is.xy.model=TRUE, trace, Callers.name)
    if(trace >= 1) {
        print_summary(x, sprint("%s: x", Callers.name), trace=2)
        trace2(trace, "\n")
    }
    expand.arg(x, env, trace, is.y.arg=FALSE)
}
# object was created with earth.formula
# called only by get.earth.x

get.earth.x.formula <- function(object, data, env, trace, Callers.name)
{
    terms.without.response <- delete.Response(object$terms)
    data <- get.update.arg(data, "data", object, env, trace, Callers.name)
    namesx <- rownames(object$modvars)
    data <- possibly.convert.vector.to.matrix(data, namesx)
    # following allows data to be a list e.g. newdata=etitanic[1,,drop=TRUE]
    data <- possibly.convert.list.to.data.frame(data)
    data <- fix.newdata.cols(data, namesx, is.xy.model=FALSE, trace, Callers.name)
    data <- as.data.frame(data)
    expected.nrows <- nrow(data)
    if(trace >= 1) {
        print_summary(data, sprint("%s: x", Callers.name), trace=2)
        trace2(trace, "\n")
    }
    if(!is.null(attr(terms.without.response, "offset")))
        check.offset.var.is.in.data(terms.without.response, data)

    colnames(data) <- gsub("\`", "", colnames(data)) # remove backticks if any

    # March 2019: added xlev to match what lm does (and also linmod.R in the plotmo tests)
    # necessary for: mod <- earth(Sepal.Length~Species, data=iris);
    #                predict(mod, newdata=data.frame(Species="setosa")) # used to fail
    mf <- model.frame(terms.without.response, data=data, na.action=na.pass, xlev=object$xlevels)
    if(trace >= 1) {
        print_summary(mf,
                      sprint("%s: after call to model.frame: mf", Callers.name),
                      trace=2)
        trace2(trace, "\n")
    }
    classes <- attr(terms.without.response, "dataClasses")
    if(!is.null(classes)) {
        # Use "try" for leniency, to allow numeric to be used for factors etc.
        # There is special treatment for the following message because it seems to be benign:
        #   variable 'foo' was fitted with type "nmatrix.1" but type "numeric" was supplied
        try <- try(.checkMFClasses(classes, mf), silent=TRUE)
        if(is.try.err(try) && !grepl("\"nmatrix.1\" .* \"numeric\"", try[1])) {
            cat(try)
            cat("Continuing anyway, first few rows of modelframe are\n")
            print(head(mf))
        }
    }
    x <- model.matrix(terms.without.response, mf)
    check.nrows(expected.nrows, nrow(x), nrow(object$fitted.values), Callers.name)
    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]    # silently discard intercept
    x
}
# Like stats::delete.response but can handle multiple-response
# Formula objects with a "Response" attr.
# Can also handle conventional formula objects with a "response" attr.

delete.Response <- function (termobj, issue.warning=TRUE)
{
    a <- attributes(termobj)
    y <- a$response # response index
    if(is.null(y) || y[1] == 0)
        y <- a$Response # multiple resp termobj built with Formula
    if(is.null(y) || y[1] == 0) {
        if(issue.warning) {
            formula <- termobj
            attributes(formula) <- NULL
            warning0("formula has no response: ",
                     quotify(paste(formula, collapse=" ")))
        }
        return(termobj)
    }
    # following copied from stats::delete.response, R version 3.5.3 (March 2019)
    #
    # TODO Sep 2020: The R source code has changed for delete.response
    # Comment in new R source code (4.0.0 dev): "do this `by hand' as previous approach was vulnerable to re-ordering"
    # Therefore we need to update the code below to match the new R source code.
    a$response <- 0
    a$variables <- a$variables[-(1+y)]
    a$predvars <- a$predvars[-(1+y)]
    if(length(a$factors))
       a$factors <- a$factors[-y, , drop = FALSE]
    if(length(a$offset))
        a$offset <- ifelse(a$offset > y, a$offset-1, a$offset)
    if(length(a$specials)) {
        for(i in seq_along(a$specials)) {
            b <- a$specials[[i]]
            a$specials[[i]] <- ifelse(b > y, b-1, b)
        }
    }
    if(length(y) == 1)       # conventional formula object?
        termobj[[2]] <- NULL # termobj is list(~, response, rhs)
    else {                   # multiple response Formula object
        check.ymax <- function(len) {
            if(len < ymax) {
                attributes(termobj) <- NULL # for paste in error message
                stop0("Cannot delete response from ",
                      quotify(paste(termobj, collapse=" ")),
                      "\n       because ", deparse(substitute(len)),
                      " is ", len, ", expected length at least ", ymax)
            }
        }
        termobj <- strip_multiple_response_from_Formula(termobj)
        ymax <- max(y) # for error checking
        if(length(a$factors)) {
            check.ymax(NCOL(a$factors))
            a$factors <- a$factors[, -y, drop = FALSE]
        }
        check.ymax(length(a$term.labels))
        a$term.labels <- a$term.labels[-y]
        # TODO do we need the following?
        # check.ymax(length(a$order))
        # a$term.order <- a$order[-y]
        a$Response <- 0
    }
    attributes(termobj) <- a
    termobj
}
# Returns modified formula (the modified termobj) without attributes.
# TODO Is there are simpler way of doing this?

strip_multiple_response_from_Formula <- function(termobj) # termobj created by Formula
{
    check.class <- function(element, classes) {
        if(!class(element)[1] %in% classes)
            stop0("Cannot delete response from ",
                  quotify(paste(Formula, collapse=" ")),
                  "\n       because class(", deparse(substitute(element)),
                  ") is ", quotify(class(element)), " which is not in ", quotify(classes))
    }
    check.index <- function(element, i) {
        if(length(element) < i) {
            stop0("Cannot delete response from ",
                  quotify(paste(Formula, collapse=" ")),
                  "\n       because length(", deparse(substitute(element)),
                  ") is ", length(element), ", expected length at least ", i)
        }
    }
    Formula <- termobj # termobj is list(~,  list(     +, response,    rhs))
                       # index:          [1]  [2]  [2][1]    [2][2]  [2][3]

    attributes(Formula) <- NULL # converts class c("terms","formula") to "call"
    check.class(Formula, "call")
    check.index(Formula, 2)
    check.class(Formula[[2]], c("name", "call"))
    check.index(Formula[[2]], 3)
    check.class(Formula[[2]][[3]], c("name", "call"))
    Formula[[2]] <- Formula[[2]][[3]]  # extract rhs into 2nd elem of Formula
    Formula
}
check.expanded.ncols <- function(x, object) # called only by model.matrix.earth
{
    if(NCOL(x) != NCOL(object$dirs)) {
        format <- paste0("model.matrix.earth could not interpret the data\n",
                         "           model.matrix returned %d column%s %s\n",
                         "           need %d column%s %s")
        stopf(format,
              NCOL(x), if(NCOL(x) == 1) ":" else "s:",
              if(NCOL(x) == 0) "" else paste.with.quotes(colnames(x), maxlen=50),

              NCOL(object$dirs), if(NCOL(object$dirs) == 1) ":" else "s:",
              paste.with.quotes(colnames(object$dirs), maxlen=50))
    }
}
# Called by predict.earth and can also be called by users directly.
# Return object$bx if all x, subset, which.terms equal NULL.

model.matrix.earth <- function(     # returns bx
    object       = stop("no 'object' argument"),
    x            = NULL,            # x arg (not yet expanded)
    subset       = NULL,
    which.terms  = NULL,
    trace        = 0,
    ...,                            # unused, for generic method comparibility
    Env          = parent.frame(),
    Callers.name = "model.matrix.earth") # caller's name for trace messages
{
    warn.if.dots(...)
    check.classname(object, substitute(object), "earth")
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    if(is.null(x) && is.null(subset) && is.null(which.terms)) {
        if(trace >= 1)
            cat0(Callers.name, ": returning object$bx\n")
        return(object$bx)
    }
    x <- get.earth.x(object, data=x, Env, trace, paste("get.earth.x from", Callers.name))
    if(is.null(which.terms))
        which.terms <- object$selected.terms
    if(!is.null(subset)) {
        # duplicates are allowed in subsets so user can specify a bootstrap sample
        check.index("subset", subset, x, allow.dups=TRUE, allow.zeros=TRUE)
        x <- x[subset, , drop=FALSE]
    }
    get.bx(x, which.terms, object$dirs, object$cuts)
}
# Called by update.earth and get.earth.x
#
# Which x should we use? The precedence is [1] the x parameter, if any,
# in this call to update [2] the $x in the earth object (which exists
# if keepxy=TRUE was used the original call to earth) [3] the x found
# in the original call to earth.
# Same applies for y, subset, weights, and wp.
# The "arg" argument is from the current call to update or predict

get.update.arg <- function(arg, argname, object, env,
                           trace1, Callers.name="update.earth", print.trace=TRUE,
                           reeval=TRUE) # TODO hack to re-evaluate
{
    if(!print.trace) # print.trace arg prevents recursion issues with trace
        trace1 = FALSE
    if(is.null(arg)) {
        temp <- try(eval(object[[argname, exact=TRUE]], envir=env), silent=TRUE)
        if(!is.null(temp) && !is.try.err(temp)) {
            if(reeval)
                arg <- object[[argname, exact=TRUE]]
            else
                arg <- temp
            if(trace1 >= 1)
                cat0(Callers.name, ": using ",
                     NROW(temp), " by ", NCOL(temp), " ", argname,
                     " saved by keepxy in original call to earth\n")
        } else {
            temp <- try(eval(object$call[[argname, exact=TRUE]], envir=env), silent=TRUE)
            if(!is.null(temp) && !is.try.err(temp)) {
                if(reeval)
                    arg <- object$call[[argname, exact=TRUE]]
                else
                    arg <- temp
                if(trace1 >= 1)
                    cat0(Callers.name, ": using ",
                         NROW(temp), " by ", NCOL(temp), " ", argname,
                         " argument from original call to earth\n")
             }
        }
    }
    arg
}
# If stats::model.frame can't interpret the data passed to it it silently
# returns the fitted values.  This routine makes that not silent.
# Note that this won't work if where model.frame returns the wrong results
# but coincidentally returns actual.nrows.expected.nrows.

check.nrows <- function(expected.nrows, actual.nrows, fitted.nrows, Callers.name)
{
    if(actual.nrows != expected.nrows) {
        if(actual.nrows == fitted.nrows)
            stop0("model.frame.default could not interpret the data passed to ",
                  Callers.name,
                  "\n        (actual.nrows=", actual.nrows,
                  " expected.nrows=", expected.nrows,
                  " fitted.nrows=", fitted.nrows, ")")
        else  # can probably never get here
            warning0(Callers.name, " returned a number ", actual.nrows,
                     " of rows that was different from the number ",
                     expected.nrows,
                     " of rows in the data")
    }
}
# If x is already a matrix or data.frame this does nothing.
# Else if x is a vector, return a matrix with length(colnames) columns.

possibly.convert.vector.to.matrix <- function(x, colnames)
{
    if(is.null(ncol(x)) && !is.list(x)) {
        nrows <- length(x) / length(colnames)
        if(floor(nrows) == nrows)
            dim(x) <- c(nrow=nrows, ncol=length(colnames))
        else
            stop0("Could not convert vector x to matrix because ",
                  "length(x) ", length(x), "\n",
                  "       is not a multiple of the number ",
                  length(colnames), " of predictors ",
                  "\n       Expected predictors: ",
                  paste.with.quotes(colnames, maxlen=50))
    }
    x
}
possibly.convert.list.to.data.frame <- function(x)
{
    if(is.list(x) && !is.data.frame(x))
        x <- as.data.frame(x)
    x
}
# Given an x matrix or data.frame, return an x with column names equal to
# expected.colnames and with the columns in their correct order
# So the user can hand us an x without column names,
# or with named columns but in the wrong order, or an x containing
# only a needed subset of all the columns, etc.
# This code is a mess and doesn't handle all cases, just the common ones.

fix.newdata.cols <- function(x, namesx, is.xy.model, trace, Callers.name)
{
    colnames <- colnames(x)
    ncolnames <- length(colnames)
    nexpected <- length(namesx)

    if(is.null(colnames)) {
        if(trace >= 1)
            cat0(Callers.name, ": x has no column names, ",
                 "adding column names: ", paste.collapse(namesx),
                 "\n")
        ncol <- min(ncol(x), length(namesx))
        colnames(x)[1:ncol] <- namesx[1:ncol]
        colnames <- colnames(x)

     } else if(ncolnames < nexpected) {
        ret <- add.missing.newdata.cols(x, namesx, trace, Callers.name)
            x <- ret$x
            colnames <- ret$colnames

    } else if(ncolnames > nexpected) {
        NULL # TODO not sure what to do here (do nothing so old regression tests pass)

    } else {
        ret <- fix.newdata.colnames(x, namesx, trace, Callers.name)
            x <- ret$x
            colnames <- ret$colnames
    }

    colnames(x) <- colnames
    x
}
# Allow user to specify less than the expected
# nbr of columns -- which is ok if they specify all predictors
# actually used by the model.
#
# Called only by fix.newdata.cols
# which in turn is called only get.earth.x
# (via get.earth.x.default and get.earth.x.formula)

add.missing.newdata.cols <- function(x, namesx, trace, Callers.name)
{
    colnames <- colnames(x)
    nexpected <- length(namesx)

    imatch <- pmatch(colnames, namesx, nomatch=0)
    if(any(imatch == 0)) {
        # can't repair the error because there are colnames in x that aren't
        # in expected.names (tends to happen with expanded factor names)
        format <- paste0("could not interpret newdata\n",
                         "           model.matrix returned %d column%s %s\n",
                         "           need %d column%s %s")
        stopf(format,
              NCOL(x), if(NCOL(x) == 1) ":" else "s:",
              if(NCOL(x) == 0) "" else paste.with.quotes(colnames(x), maxlen=50),

              length(namesx), if(length(namesx) == 1) ":" else "s:",
              paste.with.quotes(namesx, maxlen=50))

    }
    # Create a new x, putting the existing cols into their correct positions.
    # Cols that aren't in the original x will end up as all NAs in the
    # the recreated x; that doesn't matter for predict.earth if they are
    # for predictors that are unused in the earth model.
    if(trace >= 1)
        cat0("newdata has missing columns, adding missing cols with all NAs\n")
    imatch <- pmatch(namesx, colnames, nomatch=0)
    x.original <- x
    x <- matrix(data=NA_real_, nrow=nrow(x), ncol=nexpected)
    for(i in seq_len(nexpected))
        if(imatch[i])
            x[,i] <- x.original[,imatch[i]]
    colnames <- namesx
    list(x=x, colnames=namesx)
}
# called only fix.newdata.cols (called when ncolnames == nexpected)

fix.newdata.colnames <- function(x, namesx, trace, Callers.name)
{
    colnames <- colnames(x)
    imatch <- pmatch(colnames, namesx, nomatch=0)
    if(all(imatch == 0)) {
        if(trace >= 1)
            cat0(Callers.name,
                 ": unexpected x column names, renaming columns\n",
                 "    Old names: ", paste.collapse(colnames), "\n",
                 "    New names: ", paste.collapse(namesx),
                 "\n")

        colnames <- namesx
    } else {
         # replace indices for non-found predictor names with their value
         # i.e. assume columns with unknown names are in their right position
         for(i in seq_along(imatch))
             if(imatch[i] == 0)
                 imatch[i] = i

         # if any columns are in the wrong order then fix their order
         # (imatch will be 1,2,3,... if columns are in the right order)

         if(!all(imatch == seq_along(imatch))) {
           s <- paste0(Callers.name, ": x columns are in the wrong order%s\n",
                       "    Old columns: ", paste.collapse(colnames), "\n",
                       "    New columns: ", paste.collapse(namesx),
                       "\n")
           if(length(imatch) == ncol(x)) {
               trace1(trace, s, ", correcting the column order")
               x <- x[, imatch, drop=FALSE]
               colnames <- colnames[imatch]
           } else
               warnf(s, "")
        }
    }
    list(x=x, colnames=colnames)
}
strip.func.call <- function(colnames) # e.g. "as.numeric(x3,99)" becomes "x3"
{
    regex <- ".+\\("                                 # matches foo(, does not match (Intercept)
    if(any(grepl(regex, colnames))) {
        colnames <- gsub(regex, "", colnames)        # replace foo(
        colnames <- gsub("[,)][^+-]*", "", colnames) # remove remaining ",arg1,arg2)"
    }
    colnames
}

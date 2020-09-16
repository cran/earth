# expand.arg.R:
#
# This module provides, amongst others, the following function:
#
# expand.arg(x, env, is.y.arg, name)
#       Expand factors in x and convert to double mat with col names
#       Called by earth.formula, earth.default, get.earth.x
#
#-----------------------------------------------------------------------------

# Return x with all values converted to double, and with factors expanded to ind cols.
#
# Always returns a matrix (never a vector) and always with column names.
#
# Factors in earth's y argument (is.y.arg==TRUE) are alays expanded
# using contr.earth.response(), even if the factors are ordered.
# i.e. one indicator column for each factor level.
#
# The internal strategy here is essentially:
#         if(x is already double)
#             return x unchanged (but add colnames if necessary)
#         else
#             convert x to a data.frame and invoke model.frame and model.matrix

expand.arg <- function(
    x,              # "x" is x or y arg to earth
    env,            # evironment for evaluation
    trace,          # passed to gen.colnames
    is.y.arg=FALSE, # is.y.arg is TRUE if y arg to earth
    name=NULL)      # used for colnames when x has no name
{
    expand.arg.modvars(x, env, trace, is.y.arg, name)$x
}

# like expand.arg but also return the modvars matrix
# (modvars translates from the expanded x colnames to original var names)

expand.arg.modvars <- function(
    x,              # "x" is x or y arg to earth
    env,            # evironment for evaluation
    trace,          # passed to gen.colnames
    is.y.arg=FALSE, # is.y.arg is TRUE if y arg to earth
    name=NULL)      # used for colnames when x has no name
{
    if(is.null(ncol(x))) # ensure x is a matrix, not a vector
        dim(x) <- c(nrow=length(x), ncol=1)
    if(is.y.arg) {
        # We must do this here else the call to model.matrix later generates
        # two columns for each logical or two-level factor column.
        x <- convert.two.level.resp.to.numeric(x)
    }
    if(is.double(x)) {
        # Already double so no need to convert.  Note that is.double() returns
        # TRUE for a matrix of doubles but always FALSE for data.frames.
        colnames(x) <- gen.colnames(x, name, if(is.y.arg) "y" else "x", trace)
        return(list(x=x, modvars=get.identity.modvars(x)))
    }
    if(is.y.arg) {
        # we always use contr.earth.response for the y argument (left side of formula)
        old.contrasts <- getOption("contrasts")
        on.exit(options(contrasts=old.contrasts))
        options(contrasts=c("contr.earth.response", "contr.earth.response"))
    }
    colnames(x) <- gen.colnames(x, name, if(is.y.arg) "y" else "x", trace)
    x <- as.data.frame(x)
    ncol.org <- ncol(x)
    mf <- call("model.frame", formula = ~., data=x, na.action=na.pass)
    mf <- eval(mf, envir=env)
    terms <- terms(mf)
    x <- model.matrix(object=attr(mf, "terms"), data=mf)
    # the "assign" attribute has an entry for each column in x
    # giving the term in the formula which gave rise to the column
    xassign <- attr(x, "assign")
    xassign <- xassign[-1] # delete response (-1 correct even with multiple responses)

    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]    # discard intercept

    # If x had only one col, model.matrix sometimes prepends "x" to the col names.
    # This seems to always happen if the column of x was a factor (which gets
    # expanded to multiple columns), but sometimes at other times. I don't know why.
    # If so, remove the "x" prefix.
    colnamesx <- colnames(x)
    if(ncol.org == 1 &&
            all(substr(colnamesx, 1, 1) == "x") &&
            # don't strip prefix if colnames are c("x","y") or c("x","x3")
            all(nchar(colnamesx[1]) > 1) &&
            # don't strip x. prefix from orded factors (which get expanded to x.L x.Q x.C ..)
            any(substr(colnamesx, 1, 2) != "x.")) {
        # remove the "x" prefix.
        colnames(x) <- substr(colnamesx, 2, 61) # strip 1st char of each colname
    }
    list(x=x, # all cols of x are now of type double, with column names
         modvars=get.modvars(x, xassign, terms))
}
# contr.earth.response returns an nlevels by nlevels diag matrix e.g.
#
#       A B C
#     A 1 0 0
#     B 0 1 0
#     C 0 0 1
#
# The base and contrasts arguments are ignored

contr.earth.response <- function(x, base, contrasts)
{
     contr <- array(0, c(length(x), length(x)), list(x, x))
     diag(contr) <- 1
     contr
}
# Here "two.level" means logical or two-level factor.
# These get converted to a numeric column of 0s and 1s.
# This code doesn't touch y if no changes are needed.

convert.two.level.resp.to.numeric <- function(y)
{
    stopifnot(!is.null(dim(y)))
    if(is.data.frame(y)) {
        # Dataframe, so handle each column independently.
        # Get here if y is a dataframe in call to earth.default.
        for(icol in seq_len(ncol(y))) {
            ycol <- y[,icol]
            if(is.logical(ycol))
                y[,icol] <- as.numeric(ycol)
            else if(is.factor(ycol) && nlevels(ycol) <= 2) # two-level factor?
                y[,icol] <- as.numeric(ycol) - 1
        }
    } else {
        # Not dataframe, must be a matrix.  All columns are of the same class.
        # Can't use above code because for example y[,icol] <- as.numeric(ycol)
        # generates NAs if y is a matrix of factors.
        convert.logical <- convert.factor <- FALSE
        colnames <- colnames(y)
        for(icol in seq_len(ncol(y))) {
            ycol <- y[,icol]
            if(is.character(ycol)) {
                # model.matrix in expand.arg would later convert the column to N columns
                # if nbr unique strings is N, which is incorrect, so block that here
                stop0("y is a character variable: ", paste.with.quotes(y, maxlen=40))
            } else if(is.logical(ycol)) {
                convert.logical <- TRUE
            } else if(is.factor(ycol) && nlevels(ycol) <= 2) { # two-level factor?
                convert.factor <- TRUE
                if(!is.null(colnames) || ncol(y) == 1)
                    colnames[icol] <- levels(ycol)[2]
            }
        }
        nrow <- nrow(y)
        ncol <- ncol(y)
        if(convert.logical) {
            y <- as.numeric(y)      # convert to 0s and 1s
            dim(y) <- c(nrow, ncol)
            colnames(y) <- colnames
        } else if(convert.factor) {
            y <- as.numeric(y) - 1  # minus 1 to convert to 0s and 1s
            dim(y) <- c(nrow, ncol)
            colnames(y) <- colnames
        }
    }
    y
}

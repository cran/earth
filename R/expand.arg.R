# expand.arg.R:

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
     len <- length(x)
     contr <- array(0, c(len, len), list(x, x))
     diag(contr) <- 1
     contr
}

# Return x but with factors expanded into dummy variables.
# and with all values converted to double.
# Always returns a matrix (never a vector) and always with column names.
# Factors in the y argument are treated in a non standard
# way --- see the earth man page for details.

expand.arg <- function(x,               # "x" is x or y arg to earth
                    env,                # evironment for evaluation
                    is.y.arg=FALSE,     # is.y.arg is TRUE if y arg to earth
                    xname=NULL)         # used for colname when x has no name
{
    # "two.level" here means logical or two level factor
    # They get converted to a single numeric column of 0s and 1s.
    # Must do this here else the call to model.matrix later generates
    # two cols for each logical or two level factor column.

    convert.two.level.to.numeric <- function(y)
    {
        stopif(is.null(dim(y)))
        if(class(y)[1] == "data.frame") {
            # dataframe, handle each column independently
            # get here if y is a dataframe in call to earth.default

            for(icol in 1:ncol(y)) {
                if(is.factor(y[,icol]) && nlevels(y[,icol]) <= 2)
                    y[,icol] <- as.numeric(y[,icol]) - 1
                else if(is.logical(y[,icol]))
                    y[,icol] <- as.numeric(y[,icol])
            }
        } else {
             # not dataframe, so all columns must be the same class

            nrows <- nrow(y)
            ncols <- ncol(y)
            convert.from.two.levels <- FALSE
            convert.from.logical <- FALSE
            colnames. <- colnames(y)
            for(icol in 1:ncol(y)) {
                ycol <- y[,icol]
                if(is.factor(ycol) && nlevels(ycol) <= 2) {
                    convert.from.two.levels <- TRUE
                    colnames.[icol] <- format(levels(ycol)[2])
                    break;
                }
                else if(is.logical(ycol)) {
                    convert.from.logical <- TRUE
                    break;
                }
            }
            if(convert.from.two.levels) {
                y <- as.numeric(y) - 1  # minus 1 to convert to 0s and 1s
                dim(y) <- c(nrows, ncols)
            } else if(convert.from.logical) {
                y <- as.numeric(y)      # convert to 0s and 1s
                dim(y) <- c(nrows, ncols)
            }
            colnames(y) <- colnames.
        }
        y
    }
    #--- expand.arg starts here ---

    if(is.null(ncol(x)))        # make sure x is a matrix, not a vector
        dim(x) <- c(nrow=length(x), ncol=1)
    if(is.y.arg)
        x <- convert.two.level.to.numeric(x)
    if(is.double(x)) {          # already double? then no need to convert
        colnames(x) <- generate.colnames(x, is.y.arg, xname)
        return(x)
    }
    if(is.y.arg) {
        old.contrasts <- getOption("contrasts")
        on.exit(options(contrasts=old.contrasts))
        options(contrasts=c("contr.earth.response", "contr.earth.response"))
    }
    is.data.frame <- class(x)[1] == "data.frame"
    if(is.data.frame)
        mf <- call("model.frame", formula = ~., data=x, na.action=na.pass)
    else
        mf <- call("model.frame", formula = ~x, na.action=na.pass)
    mf <- eval(mf, env)
    mf.has.colnames <- !is.null(colnames(mf))
    x <- model.matrix(object=attr(mf, "terms"), data=mf)
    intercept <- match("(Intercept)", colnames(x), nomatch=0)
    if(intercept)
        x <- x[, -intercept, drop=FALSE]     # discard intercept

    # If !is.data.frame, model.matrix prepends "x" to the column names,
    # so remove the "x", but only if x had column names to begin with.

    if(!is.data.frame && mf.has.colnames) {
        # strip 1st char of each column name
        colnames(x) <- sapply(colnames(x), substr, 2, 99)
    }
    colnames(x) <- generate.colnames(x, is.y.arg, xname)

    x   # all columns are now double with column names
}

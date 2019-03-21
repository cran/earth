# bpairs.R: code to support "binomial pairs" in glm binomial models

is.bpairs <- function(y, family, trace, is.earth) # true if response y is a binomial pair
{
    msg <- NULL
    is.bpairs <- FALSE
    is.nonneg <- TRUE
    is.binomial <- is.binomial(family)
    if(is.binomial)
        check.no.na.in.mat(y) # prevent confusing errors later from comparisons
    if(is.binomial && NROW(y) > 1 && NCOL(y) == 2) {
        y1 <- y[, 1, drop=TRUE]
        y2 <- y[, 2, drop=TRUE]
        if(!(is.numeric(y1) || is.logical(y1)) ||
           !(is.numeric(y2) || is.logical(y2))) {
            return(FALSE) # note return
        }
        rowsums <- rowSums(y)
        is.int <- all((round(y1) == y1) & (round(y2) == y2))
        is.nonneg <- all((y1 >= 0) & (y2 >= 0))
        is.rowsums.greater.than.1 <- any(rowsums > 1)
        is.bpairs <- is.int && is.nonneg && is.rowsums.greater.than.1
        # TODO this and following stop sometimes cause duplicate messages, even with trace=0
        msg <- bpairs.msg(is.bpairs, y, rowsums, is.earth, is.nonneg,
                          is.int, is.rowsums.greater.than.1, trace)
    }
    if(is.binomial && !is.bpairs &&
            (is.numeric(y) || is.logical(y)) && (!all(y >= 0 & y <= 1) || !is.nonneg)) {
        cat0("\nprint(head(y)):\n")
        print(head(y))
        cat0("\n")
        # This preempts the following error from within glm() later:
        #   Error in eval(family$initialize) : y values must be 0 <= y <= 1
        stop0("Binomial response (see above): all values should be between 0 and 1, or a binomial pair",
              if(is.null(msg)) "" else paste0("\n       ", msg))
    }
    is.bpairs
}
# is.binomial is true and y has two columns when this is called
bpairs.msg <- function(is.bpairs, y, rowsums, is.earth, is.nonneg,
                        is.int, is.rowsums.greater.than.1, trace)
{
    msg <- NULL
    if(!is.bpairs) {
         earth.msg <-
            if(!is.earth) ""
            else sprint("\nEarth will build two GLM models with responses \"%s\" and \"%s\"",
                        colname(y, 1), colname(y, 2))
        if(!is.int) { # glm will give a warning later (Warning: non-integer #successes in a binomial glm)
            cat0("\nprint(head(y)):\n")
            print(head(y))
            cat0("\n")
            msg <- "Response has two columns but is not a binomial pair because not all values are integers"
            printf("%s%s\n\n", msg, earth.msg)
        } else if(!is.nonneg) { # will see error message below
            msg <- "Response has two columns but is not a binomial pair because some values are negative"
            trace1(trace, "%s\n", msg)
        } else if(!is.rowsums.greater.than.1) { # no warning from glm later
            msg <- "Response has two columns but is not a binomial pair because no rows sum to greater than 1"
           trace1(trace, "%s%s\n\n", msg, earth.msg)
        }
    } else if(any(rowsums == 0))
        trace1(trace,
            "Note: Both entries in row %d %sof the %s and %s response are zero\n",
            which(rowsums==0)[1],
            if(length(which(rowsums==0)) == 1) "" else "(and others) ",
            colname(y, 1), colname(y, 2))
    msg
}
# When expanding the binomial pair, the first column of the
# short y is considered to be "true", the second "false".
#
# Example short data:
#   dose temp survived died
# 1    5   20        1    3
# 2    2   20        2    3
# 3    2   30        0    1  # note 0 survived (to test)
# 4    9   20        2    0  # note 0 died
# 5    5   20        2    1  # note that predictors same as row 1 (dose=5 temp=20)
# 6    9   30        0    0  # both rows 0
#
# Equivalent long data:
#    dose temp survived
# 1     5   20        1
# 2     5   20        0
# 3     5   20        0
# 4     5   20        0
#
# 5     2   20        1
# 6     2   20        1
# 7     2   20        0
# 8     2   20        0
# 9     2   20        0
#
# 10    2   30        0
#
# 11    9   20        1
# 12    9   20        1
#
# 13    5   20        1
# 14    5   20        1
# 15    5   20        0
#
# 16    9   30        0  # both rows zero in short data, so treat as a "false"

expand.bpairs <- function(...) { UseMethod("expand.bpairs") }

expand.bpairs.formula <- function(formula=stop("no 'formula' argument"), data=NULL, sort=FALSE, ...)
{
    stop.if.dots(...)
    call <- match.call(expand.dots=FALSE)
    imatch <- match(c("formula", "data"), names(call), 0)
    mf <- call[c(1, imatch)]
    # we use the Formula package to allow multiple responses
    # TODO this is not exactly consistent with earth, which uses
    #      Formula only for formulas with + (else earth uses formula)
    Formula <- Formula::Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- Formula
    mf$na.action <- na.pass # NAs are allowed, they get propogated as is
    mf <- eval.parent(mf)

    terms <- terms(Formula, data=data)
    varnames <- attr(terms, "term.labels")
    if(length(attr(Formula, "rhs")) > 1)
        stop0("invalid formula: too many terms on the right hand side")
    x <- model.part(Formula, data=mf, rhs=1)
    if(NCOL(x) == 0)
        stopf("expand.bpairs: the right side of the formula does not have any variables")
    x <- as.data.frame(x)

    if(length(attr(Formula, "lhs")) > 1)
        stop0("invalid formula: too many terms on the left hand side")
    y <- model.part(Formula, data=mf, lhs=1)
    # following handles when lhs of formula is a matrix e.g. cbind(success, fail) ~ .
    if(NCOL(y[[1]]) > 1)
        y <- y[[1]]
    respname <- colnames(y)[1]
    if(is.null(respname) || length(respname) != 1 || nchar(respname) == 0)
        stop0("expand.bpairs: cannot get response from formula") # paranoia
    if(NCOL(y) != 2)
        stopf("expand.bpairs: '%s' does not have two columns", respname)
    y <- as.data.frame(y)
    # trace=1 below to give extra info if we invoke stop()
    if(!is.bpairs(y, family="binomial", trace=1, is.earth=FALSE))
        stopf("expand.bpairs: the left side of the formula is not a two-column matrix of binomial pairs")

    expand.bpairs.aux(x, y, sort)
}
expand.bpairs.default <- function(data = stop("no 'data' argument"), y = NULL, sort=FALSE, ...)
{
    help.msg <- paste0(
"The y argument should be one of:\n",
"\n",
"     o Two column matrix or dataframe of binomial pairs.\n",
"\n",
"     o Two-element numeric vector specifying the response columns in 'data'.\n",
"\n",
"     o Two-element character vector specifying the response column names in 'data'.\n",
"       The full names must be used (partial matching isn't supported).")

    stop.if.dots(...)
    if(is.null(y))
        stop0("expand.bpairs: no y argument\n\n", help.msg)
    stopifnot(is.vector(data) || is.matrix(data) || is.data.frame(data))
    data.name <- trunc.deparse(substitute(data)) # for possible error message
    if(length(y) == 2 && is.numeric(y) &&
            round(y[1]) == y[1] && round(y[2]) == y[2]) {
        # y is a two element numeric vector specifying two columns in data
        ycolumns <- y
        check.index(ycolumns, "ycolumns", data, is.col.index=1, allow.negatives=FALSE)
        y <- data[,  ycolumns, drop=FALSE]
        yarg.name <- sprint("%s[,c(%g,%g)]", data.name, ycolumns[1], ycolumns[2])
        check.index(-ycolumns, "ycolumns", data, is.col.index=1, allow.negatives=TRUE)
        data <- data[, -ycolumns, drop=FALSE]
        if(ncol(data) < 1)
            stop0("expand.bpairs: x is empty after removing response columns")
    } else if(length(y) == 2 && is.character(y)) {
        # y is a two element character vector specifying two columns in data
        ycolumns <- y
        check.index(ycolumns, "ycolumns", data, is.col.index=2, allow.negatives=FALSE)
        y <- data[,  ycolumns, drop=FALSE]
        colnames <- colnames(data)
        i1 <- match(ycolumns[1], colnames)
        i2 <- match(ycolumns[2], colnames)
        # following check is probably unnecessary after above call to check.index
        if(length(i1) != 1 || length(i2) != 1 || anyNA(i1) || anyNA(i2))
            stopf("expand.bpairs: cannot find '%s' or '%s' in colnames(data)", ycolumns[1], ycolumns[2])
        yarg.name <- sprint("%s[,c(\"%s\",\"%s\")]", data.name, ycolumns[1], ycolumns[2])
        ycolumns <- c(i1,i2)
        check.index(-ycolumns, "ycolumns", data, is.col.index=1, allow.negatives=TRUE)
        data <- data[, -ycolumns, drop=FALSE]
        if(ncol(data) < 1)
            stop0("expand.bpairs: x is empty after removing response columns")
    } else { # y is the response
        yarg.name <- trunc.deparse(substitute(y))
        # basic error checking to preempt confusing messages from is.bpairs() below
        if(NCOL(y) != 2 || NROW(y) < 2)
            stop0("expand.bpairs: bad y argument '", unquote(yarg.name), "'\n\n", help.msg)
    }
    data <- as.data.frame(data)
    y <- as.data.frame(y)
    # trace=1 below to give extra info if we invoke stop()
    if(!is.bpairs(y, family="binomial", trace=1, is.earth=FALSE))
        stopf("expand.bpairs: %s is not a two-column matrix of binomial pairs", yarg.name)

    expand.bpairs.aux(data, y, sort)
}
expand.bpairs.aux <- function(x, y, sort) # returns a data.frame with attributes
{
    sort <- check.boolean(sort)
    stopifnot(ncol(x) >= 1)
    if(nrow(x) != nrow(y))
        stopf("expand.bpairs: x has %d row%s, but the response y has %d row%s",
             nrow(x), if(nrow(x)==1) "" else "s",
             nrow(y), if(nrow(y)==1) "" else "s")

    # Remove columns in x that match colnames in y.
    # This allows expand.bpairs(formula=y~.,data=x)
    # when colnames(y) is c("survived", "died") and
    # "survived" and "died" are also columns in x.
    x <- possibly.delete.column(x, y, 1)
    x <- possibly.delete.column(x, y, 2)
    stopifnot(ncol(x) >= 1)
    colname.y1 <- colnames(y)[1]
    if(is.null(colname.y1) || nchar(colname.y1) == 0)
        colname.y1 <- "true"
    rowsums <- rowSums(y)
    # For simplicity, if both values in a y row are zero, we treat
    # this as a "false".  Properly we should ignore the entry, but that gets
    # very complicated, because we would need to work with a subset of the data.
    rowsums[rowsums == 0] <- 1 # include row even if both values in row are 0

    n <- sum(rowsums)               # length of long data
    ylong <- logical(length=n)      # long form of y
    xlong <- x[FALSE, , drop=FALSE] # data.frame with all variables, but zero rows
    rownames <- character(length=n) # rownames in long data
    bpairs.index <- repl(0L, nrow(y)) # for recompacting bpairs later
    i <- 1 # index into long data
    for(ishort in 1:nrow(y)) {
        bpairs.index[ishort] <- as.integer(i)
        ntrue <- y[ishort, 1, drop=TRUE] # drop is needed if y is a data.frame
        nfalse <- y[ishort, 2, drop=TRUE]
        if(ntrue + nfalse == 0) # both values zero?
            nfalse <- 1         # treat as false
        if(nfalse > 0) {
            i2 <- i + nfalse - 1
            ylong[i:i2] <- FALSE
            xlong[i:i2,] <- x[ishort,]
            rownames[i:i2] <- sprint("row%d.%d", ishort, 1:(i2-i+1))
            i <- i + nfalse
        }
        if(ntrue > 0) {
            i2 <- i + ntrue - 1
            ylong[i:i2] <- TRUE
            xlong[i:i2,] <- x[ishort,]
            rownames[i:i2] <- sprint("row%d.%d", ishort, (nfalse+1):(nfalse+i2-i+1))
            i <- i + ntrue
        }
    }
    ylong <- as.matrix(ylong)
    colnames(ylong) <- colname.y1
    df <- data.frame(ylong, xlong)
    rownames(df) <- rownames
    if(sort) {
        stopifnot(ncol(df) >= 2)
        icol <- c(2:ncol(df), 1) # want y column to be last in sort order
        # this sorts on variable values (variables on left take precedence)
        # see example on help page for "order"
        df <- df[ do.call(order, df[,icol]), ]
    } else {
        attr(df, "bpairs.index") <- bpairs.index
    }
    attr(df, "ynames") <- colnames(y)
    df
}
possibly.delete.column <- function(x, y, icol)
{
    stopifnot(is.data.frame(x) && is.data.frame(y))
    stopifnot(icol >= 1 && icol <= 2 && ncol(y) == 2)
    colname.y <- colnames(y)[icol]
    if(!is.null(colname.y) && nchar(colname.y) > 0) {
        imatch <- match(colname.y, colnames(x), 0)
        if(length(imatch) > 1) # paranoia
            stopf("multiple columns match '%s'", colname.y)
        if(imatch) {
            # warnf("dropping column '%s' from x because it matches a column name in y", colnames(x)[imatch])
            x[[imatch]] <- NULL # delete the column in x
        }
    }
    x
}

# offset.R: misc functions for handling an offset term in the earth formula

# If an offset variable is specified in the formula, the variable must be in the data
# (it can't be passed as a global variable, independently of the data).
# Without this check, the earth model will build ok, but predict.earth later fails,
# because stat::model.frame fails (with a confusing message) if offset is not in the data.
# (Models built with lm have this problem, but with earth we instead help the user here.)

check.offset.var.is.in.data <- function(terms, data)
{
    if(is.null(data))
        stopf("if an offset is specified in the formula, the 'data' argument must be used")
    offset.index <- attr(terms,"offset")
                            stopifnot(!is.null(offset.index))
    if(length(offset.index) > 1)
        stop0("only one offset is allowed")
    varnames <- rownames(attr(terms, "factors"))
                            stopifnot(!is.null(varnames))
                            stopifnot(offset.index >= 1)
                            stopifnot(offset.index <= length(varnames))
    offset.term <- varnames[offset.index]
                            stopifnot(grepl("^offset\\(", offset.term))
    # convert "offset(foo)" to "foo"
    offset.varname <- substring(offset.term, 8, nchar(offset.term)-1)
    offset.varname <- naken.collapse(offset.varname) # convert e.g. "log(Holders)" to "Holders"
    if(!(offset.varname %in% colnames(data)))
        stopf("the offset variable '%s' in '%s' must be in the data",
              offset.varname, offset.term)
}
# get offset specified in model formula, if any
get.predict.offset <- function(object, newdata, trace)
{
    terms <- object$terms
    if(is.null(terms))
        return(NULL)
    ioffset <- attr(terms, "offset")
    if(is.null(ioffset))
        return(NULL)
    # following should have already been caught by stop.if.dots in earth.fit
    stopifnot(is.null(object$call$offset))
    # following should have already been caught in earth.formula
    stopifnot(length(ioffset) == 1) # only one offset is allowed

    offset <- eval(attr(terms, "variables")[[ioffset+1]], envir=newdata)

    if(trace >= 1)
        cat0("predict.earth: offset: ", as.char(offset), "\n")
    offset
}

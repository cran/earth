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
    offset.varname <- naken(offset.varname) # convert e.g. "log(Holders)" to "Holders"
    if(!(offset.varname %in% colnames(data)))
        stopf("the offset variable '%s' in '%s' must be in the data",
              offset.varname, offset.term)
}
# This converts an offset of say "log(Holders)" to "Holders".
# It was copied from the plotmo code.
# It's overkill for nudifying offset specs, but works for that purpose.
# TODO Could use base::all.vars here instead?

naken <- function(s) # e.g. "s(x3,x4,df=4)" becomes "x3+x4"
{
    s <- paste.collapse(strip.space(s))

    # We don't want to mess with anything in [square brackets].
    # So we replace the bracketed expression with "#BRACKETS#",
    # and then replace that back again at the end.
    # Needed for e.g. lm(trees[,3]~trees[,1:2])

    brackets <- replace.brackets("\\[.*\\]", "#BRACKETS#", s)
    s <- brackets$s

    s <- gsub("[-*/:]", "+", s)                 # replace - / * : with +

    # next two gsubs allow us to retain "x=x1" but drop "df=99" from "bs(x=x1, df=99)"

    s <- gsub("\\(._$[[:alnum:]]+=", "(", s)    # replace "(ident=" with "("
    s <- gsub("[._$[:alnum:]]+=[^,)]+", "", s)  # delete "ident=any"

    # replace ",ident" with ")+f(ident", thus "s(x0,x1)" becomes "s(x0)f(x1)"

    s <- gsub(",([._$[:alpha:]])", ")+f(\\1", s)

    if(grepl("[._$[:alnum:]]*[(]", s)) {
        s <- gsub("[._$[:alnum:]]*[(]", "", s)  # replace ident(
        s <- gsub("[,)][^+-]*", "", s)          # remove remaining ",arg1,arg2)"
    }
    # s is now something like x1+x2, split it on "+" for further processing

    s <- strsplit(s, "+", fixed=TRUE)[[1]]
    s <- unique(s) # remove duplicates
    # remove numbers e.g. sin(x1*x2/12) is nakened to x1+x1+12, we don't want the 12
    is.num <- sapply(s, function(x) grepl("^([0-9]|\\.[0-9])", x))
    # but keep the intercept if there is one
    which1 <- which(s == "1")
    is.num[which1] <- FALSE
    s <- paste0(s[!is.num], collapse="+")

    sub("#BRACKETS#",  brackets$brackets, s) # change #BRACKETS# back to what it was
}
replace.brackets <- function(pattern, place.holder, s) # utility for naken
{
    brackets <- ""
    i <- regexpr(pattern, s)
    if(i > 0) {
        last <- i + attr(i,"match.length") - 1
        stopifnot(last > i)
        brackets <- substr(s, i, last)          # remember the bracketed expression
        s <- paste0(substr(s, 1, i-1), place.holder, substring(s, last+1))
                                                # replace [.*] with #BRACKETS#
    }
    return(list(s=s, brackets=brackets))
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

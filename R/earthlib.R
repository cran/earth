# earthlib.R: general purpose routines that are needed by the earth package
#
# Comments containing "$$" mark known issues.
# Stephen Milborrow Mar 2007 Petaluma
#
# $$ the make.space functions should take into account char height and other par settings

#--------------------------------------------------------------------------------------------
# Miscellaneous utilities

# The function stopif is intended for catching programmer errors.
# For user errors, we try to give a more informative message.

stopif <- function(...)
    stopifnot(!(...))

stop1 <- function(...)          # call.=FALSE so use traceback() to see call
    stop(..., call.=FALSE)

warning1 <- function(...)       # set options(warn=2) and traceback() to see the call
    warning(..., call.=FALSE)

paste.with.space <- function(s) paste(s, collapse = " ")

paste.quoted.names <- function(names) # add quotes and comma seperators
    paste("\"", paste(names, collapse="\" \""), "\"", sep="")

strip.white.space <- function(s) gsub("[ \t\n]", "", s)

pastef <- function(s, format, ...)          # paste the c printf style args to s
    paste(s, sprintf(format, ...), sep="")

warn.if.dots.used <- function(func.name, ...)
{
    dots <- match.call(expand.dots = FALSE)$...
    if(length(dots) == 1)
        warning1(func.name, " ignored argument \"", names(dots), "\"")
    else if(length(dots) > 1)
        warning1(func.name, " ignored arguments ", paste.quoted.names(names(dots)))
}

check.classname <- function(object, object.name, class.names)
{
    if(!inherits(object, class.names))
        stop1("the class of \'", object.name, "\' is not \"",
            paste(class.names, collapse="\" or \""), "\"")

}

# Remove the "..." arguments from call.
#
# This is used in model functions that don't use "..." and issue a
# warning if a "..." argument is used.
# If we don't remove the "..." argument from the call, subsequent
# calls to the model function via update continue will re-issue the
# warning, and there is no easy way for the user to get rid of
# the unwanted argument.

strip.dots.from.call <- function(call)
{
    idots <- which(names(call)=="...")
    if(length(idots) > 0)
        call[-idots]
    else
        call
}

# Modify call to refer to the generic e.g. "foo.default" becomes "foo".
# This means that functions like update() call foo() and not foo.default().
# An advantage is that we don't have to export foo.default().

make.call.generic <- function(call, func.name)
{
    call[[1]] <- as.name(func.name)
    call
}

# warn.if.not.all.finite is intended to help clarify possibly confusing
# messages from within the bowels of calls to other functions later

warn.if.not.all.finite <- function(x, text="unknown")
{
    if(any(is.factors <- sapply(x, is.factor)))
        x <- x[, !is.factors]               # remove factor columns
    if(!all(sapply(x, is.finite))) {
        warning1("non finite value in ", text)
        return(TRUE)                        # return TRUE if warning issued
    }
    FALSE
}

# Check that an index vector specified by the user is ok to index an object.
# We want to preclude confusing R messages or behaviour later.
# An example is when max(indexVec) > len(object) which quietly returns NA
# and can cause confusing downstream behaviour.

check.index.vec <- function(index.name, indexVec, object, check.empty = FALSE)
{
    if(is.null(indexVec)) {
        if(check.empty)
            stop1("'", index.name, "' is NULL")
        return(NULL)
    }
    if(any(is.na(indexVec)))
        stop1("NA in '", index.name, "'")
    if(!(NROW(indexVec) == 1 || NCOL(indexVec) == 1))
        stop1("'", index.name, "' must be a vector not a matrix (it has dimensions ",
            NROW(indexVec), " x ", NCOL(indexVec), ")")

    # assume that if object is an array then subset chooses rows (not cols)
    len <- if(is.vector(object)) length(object) else NROW(object)

    if(is.logical(indexVec)) {
        if(check.empty) {
            if(length(indexVec) == 0)
                stop1("length(", index.name, ") == 0")
            if(length(indexVec[indexVec == TRUE]) == 0)
                stop1("'", index.name, "' is all FALSE")
        }
        if(length(indexVec) > len)
            stop1("logical index vector '", index.name,
                "' is too long (its length is ", length(indexVec),
                "and the max allowed length is ", len, ")")
    } else if(is.numeric(indexVec)) {
        if(check.empty && length(indexVec) == 0)
            stop1("length(", index.name, ") == 0")
        if(any(indexVec < 0) && any(indexVec > 0))
            stop1("mixed negative and positive values in '", index.name, "'")
        if(any(indexVec == 0) && length(indexVec) != 1)
            warning1("ignored zero in '", index.name, "'")
        if(any(duplicated(indexVec)))
            warning1("duplicates in '", index.name, "'")
        if(any(abs(indexVec) > len))
            stop1("out of range value in '", index.name,
                "' (allowed range is 1:",  len, ", or negative)")
    } else
        warning1("index vector '", index.name,
            "' has an unusual class \"", class(indexVec), "\"")
    indexVec
}

exists.and.not.null <- function(object, mode="any", argname="")
{
    # following "if" is like is.null(object) but no error msg if object doesn't exist

    if(paste("'", object, "'", sep="") == "'NULL'")
        return(FALSE)

    if(paste("'", object, "'", sep="") == "'NA'")
        if(length(argname))
            stop1(argname, "=NA")
        else
            stop1(object, "illegal NA")

# $$ removed until I can get this to work reliably
#
#   if(!exists(paste("'", object, "'", sep=""), where=parent.frame(), mode=mode))
#       if(length(argname))
#           stop1("you specified ", argname, "=", object,
#               " but there is no such ", if(mode=="any") "object" else mode)
#       else
#           stop1(object, ": no such", if(mode=="any") "object" else mode)

    return(TRUE)
}

# Return a vector of n clearly distinguishable colors.  Here distinguishability
# comes before aesthetics and ideas of "matching colors" or of "ordered" sequences.
# The first three are also distinguishable on (my) monochrome printer.

discrete.plot.cols <- function(ncolors=5)
{
    cols <- c(1, "grey60", "brown", "lightblue", "pink", "green")
    if(ncolors > length(cols))   # won't really be distinguishable
        cols <- c(cols, heat.colors(ncolors - length(cols)))
    cols[1:ncolors]
}

#--------------------------------------------------------------------------------------------

match.arg1 <- function(arg)     # match.arg1 returns an integer
{
    formal.args <- formals(sys.function(sys.parent()))
    arg.name=deparse(substitute(arg))
    match.choices(arg, choices=eval(formal.args[[arg.name]]),  arg.name=arg.name)
}

match.choices <- function(arg, choices, arg.name)   # choices is a vector of strings
{
    stopif(is.null(choices))    # indicates a programming error
    if(all(arg == choices))
        return(1)
    i <- pmatch(arg, choices)
    if(any(is.na(i)))
        stop1(paste("bad ", arg.name, " argument \"", arg, "\"\n",
            "Choose one of: ", paste.with.space(choices), sep=""))
    if(length(i) > 1)
        stop1("'", arg.name, "' is ambiguous")
    i
}

#--------------------------------------------------------------------------------------------
get.sub.caption.from.call <- function(sub.caption, object)
{
    if(is.null(sub.caption))
        sub.caption <- strip.white.space(paste(deparse(object$call), sep="", collapse=""))
    sub.caption
}

# Call this only after a plot is on the screen to avoid
# an error message "plot.new has not been called yet"
#$$ the trimming code overtrims

show.sub.caption <- function(sub.caption, trim=FALSE)
{
    len.sub.caption <- nchar(sub.caption)
    if(len.sub.caption > 0) {
        if(trim) {
            # trim sub.caption to fit
            len <- len.sub.caption / strwidth(sub.caption, "figure")
            sub.caption <- substr(sub.caption, 1, len)
            # append ellipsis if chars deleted
            if(len < len.sub.caption)
                sub.caption <- paste(sub.caption, "...", sep="")
        }
        mtext(sub.caption, outer=TRUE, font=2, line=1.5, cex=1)
    }
    NULL
}

# the make.space functions should only be called if do.par is FALSE

make.space.for.sub.caption <- function(sub.caption)
{
    oma <- par("oma")
    if(nchar(sub.caption) > 0 && oma[3] < 3) {
        oma[3] <- 3
        par(oma=oma)
    }
    NULL
}

make.space.for.bottom.axis <- function()
{
    mar <- par("mar")
    if(mar[1] < 3) {
        mar[1] <- 3
        par(mar=mar)
    }
    NULL
}

make.space.for.left.axis <- function()
{
    mar <- par("mar")
    if(mar[2] < 3) {
        mar[2] <- 3
        par(mar=mar)
    }
    NULL
}

make.space.for.right.axis <- function()
{
    mar <- par("mar")
    if(mar[4] < 3.5) {
        mar[4] <- 3.5
        par(mar=mar)
    }
    NULL
}

#--------------------------------------------------------------------------------------------
# Given a list of objects, return a vector of strings.  Each string shows where
# the $call argument of the object differs from the call of the first object.
# (Thus the first object is used as a reference).

get.arg.strings <- function(
        objects,    # list of objects with $call arguments
        maxchars=16)
{
    # the gsub discards white space and the final ")"
    get.call <- function(iobj)
        gsub("[ \t\n)]",  "",  paste(format(objects[[iobj]]$call), collapse=" "))

    stopif(length(objects) == 0)
    call1 <- get.call(1)
    if(length(objects) == 1)
        return(substr(call1, 1, maxchars))
    call2 <- get.call(2)
    i <- first.non.matching.arg(call1, call2)
    if(i == 0)
        rval <- c("", "")
    else
        rval <- c(substr(call1, i, i+maxchars), substr(call2, i, i+maxchars))
    if(length(objects) > 2)
        for(iobj in 3:(length(objects))) {
            call2 <- get.call(iobj)
            i <- first.non.matching.arg(call1, call2)
            rval <- c(rval, if(i==0) "" else substr(call2, i, i+maxchars))
        }
    rval
}

# Return the position of the first non matching arg between two function call strings.
#
# More precisely, find the first non matching characters in s1 and s2.
# When it is found, step back until a comma or "(" is found.
# Return the index of the character after the comma.
#
# Example: s1 = lm(formula=O3~.,data=ozone
#          s2 = lm(formula=O3~.-wind,data=ozone
#
#          return index of "formula=O3~.-wind,data=ozone"
#          because formula is the first argument with a differing argument

first.non.matching.arg <- function(s1, s2)
{
    len <- min(nchar(s1), nchar(s2))
    if(len == 0)
        return(0)
    for(i in 1:len)
        if(substr(s1, i, i) != substr(s2, i, i))
            break;
    if(i == len || i == 1)  # no difference or all different?
        return(1)
    while (i >= 1 && substr(s2, i, i) != "," && substr(s2, i, i) != "(")
        i <- i - 1  # move backwards to character following comma or "("
    return(i+1)
}

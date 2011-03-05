# earthlib.R: general purpose routines that are needed by the earth package
#
# Comments containing "TODO" mark known issues.
# Stephen Milborrow Mar 2007 Petaluma
#
# TODO the make.space functions should take into account char height and other par settings

#------------------------------------------------------------------------------
# Miscellaneous utilities

any1 <- function(x) any(x != 0) # like any but no warning if x not logical

# The function stopif is intended for catching programmer errors.
# For user errors, we try to give a more informative message.

stopif <- function(...) stopifnot(!(...))

stop1 <- function(...)          # call.=FALSE so use traceback() to see call
    stop(..., call.=FALSE)

warning1 <- function(...)       # set options(warn=2) and traceback() to see the call
    warning(..., call.=FALSE)

paste.with.space <- function(s)
    paste(s, collapse=" ")

paste.with.comma <- function(s)
    paste(s, collapse=", ")

paste.quoted.names <- function(names) # add quotes and comma seperators
    paste("\"", paste(names, collapse="\" \""), "\"", sep="")

strip.white.space <- function(s)
    gsub("[ \t\n]", "", s)

printf <- function(format, ...) # like c printf
    cat(sprintf(format, ...))

pastef <- function(s, format, ...) # paste the c printf style args to s
    paste(s, sprintf(format, ...), sep="")

is.try.error <- function(obj)
   class(obj)[1] == "try-error"

paste1 <- function(...) # paste with no added spaces
    paste(..., sep="")

# is.na.or.zero's main purpose is to see if a plot component should be
# drawn, i.e., to see if the component "has a color"
# We use identical() and not is.na() below because is.na(x) gives warnings
# for certain x's, e.g if x is a function, and x == 0 gives warnings if x
# is a vector or a function etc.

is.na.or.zero <- function(x) identical(x, NA) || identical(x, 0)

warn.if.dots.used <- function(func.name, ...)
{
    dots <- match.call(expand.dots = FALSE)$...
    if(length(dots) == 1)
        warning1(func.name, " ignored unrecognized argument \"", names(dots), "\"")
    else if(length(dots) > 1)
        warning1(func.name, " ignored unrecognized arguments ",
                 paste.quoted.names(names(dots)))
}

check.classname <- function(object, object.name, class.names)
{
    if(!inherits(object, class.names))
        stop1("the class of \"", object.name, "\" is not \"",
            paste(class.names, collapse="\" or \""), "\"")
}

to.logical <- function(x, len)
{
    xlogical <- rep(FALSE, len)
    xlogical[x] <- TRUE
    xlogical
}

check.trace.arg <- function(trace) # make sure trace is a one element vector
{
    if(!is.vector(trace))
        warning1("bad \"trace\" argument")
    if(length(trace) != 1)
        warning1("\"trace\" argument has more than one element")
    as.numeric(trace[1])
}

# Modify call to refer to the generic e.g. "foo.default" becomes "foo".
# This means that functions like update() call foo() and not foo.default().
# An advantage is that we don't have to export foo.default().

make.call.generic <- function(Call, func.name)
{
    Call[[1]] <- as.name(func.name)
    Call
}

# warn.if.not.all.finite is intended to help clarify possibly confusing
# messages from within the bowels of calls to other functions later
# Return TRUE if warning issued.

warn.if.not.all.finite <- function(x, text="unknown")
{
    is.factors <- sapply(x, is.factor)
    if(any(is.factors)) {
        if(NCOL(x) == 1 || all(is.factors)) #TODO suspect
            return(FALSE)
        x <- x[, !is.factors]               # remove factor columns before is.finite check
    }
    if(any(sapply(x, is.na))) {
        warning1("NA in ", text)
        return(TRUE)
    }
    if(!all(sapply(x, is.finite))) {
        warning1("non finite value in ", text)
        return(TRUE)
    }
    FALSE
}

my.fixed.point <- function(x, digits)
{
    if(NROW(x) > 2) # if only intercept term and one other then don't use fixed point
        x <- apply(x, 2, zapsmall, digits+1)
    x
}

# Check that an index vector specified by the user is ok to index an object.
# We want to preclude confusing R messages or behaviour later.
# An example is when max(indexVec) > len(object) which quietly returns NA
# and can cause confusing downstream behaviour.

check.index.vec <- function(index.name, indexVec, object,
                        check.empty = FALSE, use.as.col.index=FALSE,
                        allow.negative.indices = TRUE,
                        allow.duplicates = FALSE,
                        allow.zeroes = FALSE)
{
    if(is.null(indexVec)) {
        if(check.empty)
            stop1("\"", index.name, "\" is NULL and cannot be used as an index vector")
        return(NULL)
    }
    if(any(is.na(indexVec)))
        stop1("NA in \"", index.name, "\"")
    if(!(NROW(indexVec) == 1 || NCOL(indexVec) == 1))
        stop1("\"", index.name, "\" must be a vector not a matrix ",
            "(\"", index.name, "\" has dimensions ",
            NROW(indexVec), " x ", NCOL(indexVec), ")")

    if(use.as.col.index)
        len <- NCOL(object)         # index is for cols of object
    else if(is.vector(object))
        len <- length(object)
    else
        len <- NROW(object)         # index is for rows of object

    if(is.logical(indexVec)) {
        if(check.empty) {
            if(length(indexVec) == 0)
                stop1("length(", index.name, ") == 0")
            if(length(indexVec[indexVec == TRUE]) == 0)
                stop1("\"", index.name, "\" is all FALSE")
        }
        if(length(indexVec) > len)
            stop1("logical index vector \"", index.name, "\" is too long\n",
                "       Its length is ", length(indexVec),
                " and the max allowed length is ", len)
    } else if(is.numeric(indexVec)) {
        if(check.empty) {
            if(length(indexVec) == 0)
                stop1("length(", index.name, ") == 0")
            else if(all(indexVec == 0))
                if(length(indexVec) == 1)
                    stop1("\"", index.name, "\" is 0")
                else
                    stop1("\"", index.name, "\" is all zeroes")
        }
        if(any(indexVec < 0) && any(indexVec > 0))
            stop1("mixed negative and positive values in \"", index.name, "\"")
        if(!allow.zeroes && any(indexVec == 0) && length(indexVec) != 1)
            warning1("zero in \"", index.name, "\"")
        if(!allow.duplicates && any(duplicated(indexVec)))
            warning1("duplicates in \"", index.name, "\"")
        if(!allow.negative.indices && any(indexVec < 0))
            stop1("negative value in \"", index.name, "\"")
        if(any(abs(indexVec) > len))
            if(len == 1)
                stop1("out of range value in \"", index.name,
                    "\" (the only legal value is 1)")
            else
                stop1("out of range value in \"", index.name,
                    "\" (allowed index range is 1:",  len, ")")
    } else
        warning1("index vector \"", index.name,
            "\" has an unusual class \"", class(indexVec), "\"")
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

    # TODO removed until I can get this to work reliably
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
# comes before aesthetics and "matching colors" or "ordered" sequences.
# The first three are also distinguishable on (my) monochrome printer.

discrete.plot.cols <- function(ncolors=5)
{
    cols <- c(1, "grey60", "brown", "lightblue", "pink", "green")
    if(ncolors > length(cols))   # won't really be distinguishable
        cols <- c(cols, heat.colors(ncolors - length(cols)))
    cols[1:ncolors]
}


# Example of my.print.call:
#
# Call: earth(formula=O3~., data=ozone1, trace=4, linpreds=c(3,
#       8), degree=2)
#
# Note that the 2nd line is horizontally aligned with the first.

my.print.call <- function(msg, Call)
{
    # don't print x or y if they are too long
    # TODO there must be a better way of checking length
    if(!is.null(Call$x)) {
        x. <- Call$x
        if(length(paste(substitute(x.))) > 100)
            Call$x <- paste("[", NROW(Call$x), ",", NCOL(Call$x),
                            "]-too-long-to-display", sep="")
    }
    if(!is.null(Call$y)) {
        y. <- Call$y
        if(length(paste(substitute(y.))) > 100)
            Call$y <- paste("[", NROW(Call$y), ",", NCOL(Call$y),
                            "]-too-long-to-display", sep="")
    }
    Call$na.action <- NULL # don't want to print the na.action
    s <- format(Call)
    if(length(s) > 8) {
        s <- s[1:8]
        s[8] <- paste(s[8], "\netc.")
    }
    s <- gsub("[ \t\n]", "", s)                 # remove white space

    # add newlines and prefix (spaces prefix all lines except the first)
    spaces. <- sprintf("%*s", nchar(msg), " ")   # nchar spaces

    s <- gsub(",", ", ", s)                     # replace comma with comma space
    s <- paste(s, collapse=paste("\n", spaces., sep=""), sep="")
    cat(msg, s, "\n", sep="")
}

#------------------------------------------------------------------------------

match.arg1 <- function(arg)     # match.arg1 returns an integer
{
    formal.args <- formals(sys.function(sys.parent()))
    arg.name=deparse(substitute(arg))
    match.choices(arg, choices=eval(formal.args[[arg.name]]),  arg.name=arg.name)
}

match.choices <- function(arg, choices, arg.name)   # choices is a vector of strings
{
    stopifnot(is.character(arg))
    stopif(is.null(choices))    # indicates a programming error
    if(all(arg == choices))
        return(1)
    i <- pmatch(arg, choices)
    if(any(is.na(i)))
        stop1(paste("bad \"", arg.name, "\" argument \"", arg, "\"\n",
            "Choose one of: ", paste.with.space(choices), sep=""))
    if(i == 0)
        stop1(paste("the \"", arg.name, "\" argument is ambiguous\n",
              "Choose one of: ", paste.with.space(choices), sep=""))
    i
}

#------------------------------------------------------------------------------
get.caption.from.call <- function(caption, object)
{
    if(is.null(caption))
        caption <- strip.white.space(paste(deparse(object$call), sep="", collapse=""))
    caption
}

# Call this only after a plot is on the screen to avoid
# an error message "plot.new has not been called yet"
# TODO the trimming code overtrims

show.caption <- function(caption, trim=0, show=TRUE)
{
    len.caption <- nchar(caption)
    if(len.caption > 0) {
        if(trim) {
            if(is.logical(trim))
                trim <- 1
            # trim caption to fit
            len <- len.caption * trim / strwidth(caption, "figure")
            caption <- substr(caption, 1, len)
            # append ellipsis if chars deleted
            if(len < len.caption)
                caption <- paste(caption, "...", sep="")
        }
        if(show)
            mtext(caption, outer=TRUE, font=2, line=1.5, cex=1)
    }
    caption
}

# the make.space functions should only be called if do.par is TRUE
# (otherwise par is not restored correctly)

make.space.for.caption <- function(caption)
{
    oma <- par("oma")
    if(nchar(caption) > 0 && oma[3] < 3) {
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

#------------------------------------------------------------------------------
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
    Call <- get.call(1)
    if(length(objects) == 1)
        return(substr(Call, 1, maxchars))
    Call2 <- get.call(2)
    i <- first.non.matching.arg(Call, Call2)
    if(i == 0)
        rval <- c("", "")
    else
        rval <- c(substr(Call, i, i+maxchars), substr(Call2, i, i+maxchars))
    if(length(objects) > 2)
        for(iobj in 3:(length(objects))) {
            Call2 <- get.call(iobj)
            i <- first.non.matching.arg(Call, Call2)
            rval <- c(rval, if(i==0) "" else substr(Call2, i, i+maxchars))
        }
    rval
}

# Return the position of the first non matching arg between two function call strings.
#
# More precisely, find the first non matching characters in s1 and s2.
# When it is found, step back until a comma or "(" is found.
# Return the index of the character after the comma or "(".
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
            break
    if(i == len || i == 1)  # no difference or all different?
        return(1)
    while(i >= 1 && substr(s2, i, i) != "," && substr(s2, i, i) != "(")
        i <- i - 1  # move backwards to character following comma or "("
    return(i+1)
}

# summarize a matrix
# TODO don't print all colnames if too many

print.matrix.info <- function(xname, x, Callers.name=NULL, bpairs=NULL,
                              details=TRUE, all.rows=FALSE, all.names=FALSE)
{
    if(!is.null(Callers.name))
        cat(Callers.name, ": ", sep="")
    cat(xname, " is a ", NROW(x), " by ", NCOL(x), " matrix: ", sep="")
    if(is.null(bpairs))
        bpairs <- rep(TRUE, NCOL(x))
    colnames <- colnames(x)
    if(is.null(colnames))
        colnames <- rep("", NCOL(x))
    stopifnot(length(bpairs) == length(colnames))
    icol <- 0
    n.names <- length(colnames)
    if(!all.names && !details)
        n.names <- min(5, n.names)
    for(i in 1:n.names) {
        if(bpairs[i]) {
            icol <- icol+1
            cat(icol, "=", sep="")
        } else
            cat(" (paired with ")
        if(nchar(colnames[i]))
            cat(colnames[i])
        else
            cat("UNNAMED")
        if(NCOL(x) > 1)
            xi <- x[,i]
        else
            xi <- x
        if(!(class(xi) == "numeric" || class(xi) == "matrix"))
            cat(" (", class(xi), ")", sep="")
        else if(typeof(xi) != "double")
            cat(" (", mode(xi), ")", sep="")
        if(!bpairs[i])
            cat(")")
        if(i < length(colnames) && (i == length(colnames) || bpairs[i+1]))
            cat(", ")
    }
    if(i < length(colnames))
        cat("...")
    cat("\n")
    if(details) {
        rownames(x) <- NULL
        if(all.rows || NROW(x) <= 6) { # head prints 6 rows
            cat("Contents of ", xname, " are\n", sep="")
            print(x)
        } else {
            rowstring <- if(class(x) == "numeric" || class(x) == "factor")
                             "elements"
                         else
                             "rows"
            cat("First few ", rowstring, " of ", xname, " are\n", sep="")
            print(head(x))
        }
    }
}

# earthlib.R: general purpose routines for the earth package
#
# TODO the make.space functions should take into account char height and other par settings

stop0 <- function(...) stop(..., call.=FALSE)

# like stop0, but pass args in printf format
stopf <- function(...) stop(sprintf(...), call.=FALSE)

warning0 <- function(...) warning(..., call.=FALSE)

# like warning0, but pass args in printf format
warnf <- function(...) warning(sprintf(...), call.=FALSE)

any1 <- function(x) any(x != 0) # like any but no warning if x not logical

# The function stopif is intended for catching programmer errors.
# For user errors, we try to give a more informative message.

stopif <- function(...) stopifnot(!(...))

cat0 <- function(...) cat(..., sep="")      # cat with no added spaces

paste0 <- function(...) paste(..., sep="")  # paste with no added spaces

paste.with.space <- function(s) paste(s, collapse=" ")

paste.with.comma <- function(s) paste(s, collapse=", ")

paste.quoted.names <- function(names) # add quotes and comma seperators
    paste0("\"", paste(names, collapse="\" \""), "\"")

quote.with.c <- function(names) # return "x" or c("x1", "x2")
{
    if(length(names) == 1)
        sprintf("\"%s\"", names)
    else
        sprintf("c(%s)", paste0("\"", paste(names, collapse="\", \""), "\""))
}
paste.with.c <- function(x) # return x or c(x1, x2)
{
    if(length(x) == 1)
        sprintf("%d", x)
    else
        sprintf("c(%s)", paste(x, collapse=", "))
}
printf <- function(format, ...) cat(sprintf(format, ...)) # like c printf

pastef <- function(s, format, ...) # paste the c printf style args to s
    paste0(s, sprintf(format, ...))

strip.white.space <- function(s) gsub("[ \t\n]", "", s)

# Note: there is no string line type corresponding to 1, so this
# converts 1 to "1" which is an illegal lty, so must be specially
# handled in functions which use the lty string.

lty.as.char <- function(lty)
{
    char <- lty
    if(is.numeric(lty)) {
        char <- NULL
        tab <- c("1", "44", "13", "1343", "73", "2262") # from par man page
        stopifnot(length(lty) > 0)
        for(i in seq_along(lty)) {
            stopifnot(lty[i] >= 1 && lty[i] <= length(tab))
            char <- c(char, tab[lty[i]])
        }
    }
    char
}
get.mean.rsq <- function(rss, tss, wp)
{
    if(is.null(wp))
        wp <- repl(1, length(rss))
    stopifnot(length(rss) == length(tss) && length(wp) == length(tss))
    total.rsq <- 0
    for(iresp in seq_along(rss))
        total.rsq <- total.rsq + wp[iresp] * get.rsq(rss[iresp], tss[iresp])
    sum(total.rsq) / sum(wp)
}
get.rsq <- function(rss, tss)
{
    rsq <- 1 - rss / tss
    # following makes testing easier across machines in presence of numerical error
    rsq[rsq > -1e-5 & rsq < 1e-5] <- 0
    rsq
}
get.weighted.rsq <- function(y, yhat, w=NULL)
{
    stopifnot(length(y) == length(yhat))
    if(is.null(w)) {
        rss <- sum((y - yhat)^2)
        tss <- sum((y - mean(y))^2)
    } else {
        stopifnot(length(w) == length(yhat))
        rss <- sum(w * (y - yhat)^2)
        tss <- sum(w * (y - weighted.mean(y, w))^2)
    }
    get.rsq(rss, tss)
}
weighted.mean <- function(x, w) sum(w * x) / sum(w)

ss <- function(x) sum(as.vector(x^2)) # sum of squares

# is.specified's main purpose is to see if a plot component should be
# drawn, i.e., to see if the component "has a color"
# We use identical() and not is.na() below because is.na(x) gives warnings
# for certain x's, e.g if x is a function, and x == 0 gives warnings if x
# is a vector or a function etc.

is.specified <- function(x) !identical(x, NA) && !identical(x, 0) && !is.null(x)

is.try.error <- function(obj) class(obj)[1] == "try-error"

warn.if.dots.used <- function(func.name, ...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots))
        warning0(func.name, " ignored unrecognized argument \"", names(dots)[1], "\"")
}
stop.if.dots.used <- function(func.name, ...)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots))
        stop0(func.name, ": unrecognized argument \"", names(dots)[1], "\"")
}
check.classname <- function(object, object.name, class.names)
{
    if(is.null(object))
        stop0("object \"", substr(object.name, 1, 30), "\" is NULL\n",
              "Expected an object of class \"",
              paste(class.names, collapse="\" or \""), "\"")
    if(!inherits(object, class.names)) {
        stop0("object \"", substr(object.name, 1, 30),
              "\" has an unexpected class\n",
              "Expected an object of class \"",
              paste(class.names, collapse="\" or \""),
              "\"\nbut got an object of class \"",
              paste(class(object), collapse="\", \""), "\"")
    }
}
to.logical <- function(x, len)
{
    xlogical <- repl(FALSE, len)
    xlogical[x] <- TRUE
    xlogical
}
repl <- function(x, length.out)
{
    check.numeric.scalar(length.out)
    stopifnot(floor(length.out) == length.out)
    stopifnot(length.out > 0)
    rep(x, length.out=length.out)
}
check.boolean <- function(b) # b==0 or b==1 is also ok
{
    if(length(b) != 1)
        stop0("the ", deparse(substitute(b)),
              " argument is not FALSE or TRUE or 0 or 1")
    if(!(is.logical(b) || is.numeric(b)) || is.na(b) || !(b == 0 || b == 1))
        stop0("the ", deparse(substitute(b)),
            " argument is not FALSE or TRUE or 0 or 1")
    b != 0 # convert to logical
}
check.numeric.scalar <- function(x, null.ok=FALSE, na.ok=FALSE, xname=NULL)
{
    if(is.null(xname))
        xname <- deparse(substitute(x))
    if(is.null(x)) {
        if(!null.ok)
            stop0(xname, "=NULL is not allowed")
    } else if(any(is.na(x))) {
        if(!na.ok)
            stop0(xname, "=NA is not allowed")
    } else if(!is.numeric(x))
        stop0("the ", xname, " argument must be numeric")
    else if(length(x) != 1)
        stop0("the ", xname, " argument must be scalar")
    x
}
check.integer.scalar <- function(x, min=NULL, null.ok=FALSE, na.ok=FALSE, xname=NULL)
{
    if(is.null(xname))
        xname <- deparse(substitute(x))
    check.numeric.scalar(x, null.ok, na.ok, xname)
    if(!is.null(min) && x < min)
        stop0("the ", xname, " argument must be at least ", min,
              " (you have ", xname, " = ", x, ")")
}
check.trace.arg <- function(trace) # make sure trace is a one element vector
{
    if(!is.vector(trace))
        warning0("illegal \"trace\" argument")
    if(length(trace) != 1)
        warning0("\"trace\" argument has more than one element")
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
my.fixed.point <- function(x, digits)
{
    if(NROW(x) > 2) # if only intercept term and one other then don't use fixed point
        x <- apply(x, 2, zapsmall, digits+1)
    x
}
# Return a vector of n clearly distinguishable colors.  Here distinguishability
# comes before aesthetics and "matching colors" or "ordered" sequences.
# The first three are also distinguishable on (my) monochrome printer.

discrete.plot.cols <- function(ncolors=5)
{
    cols <- c(1, "gray60", "brown", "lightblue", "pink", "green")
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
            Call$x <- paste0("[", NROW(Call$x), ",", NCOL(Call$x),
                             "]-too-long-to-display")
    }
    if(!is.null(Call$y)) {
        y. <- Call$y
        if(length(paste(substitute(y.))) > 100)
            Call$y <- paste0("[", NROW(Call$y), ",", NCOL(Call$y),
                             "]-too-long-to-display")
    }
    Call$na.action <- NULL # don't want to print the na.action
    s <- format(Call)
    if(length(s) > 8) {
        s <- s[1:8]
        s[8] <- paste(s[8], "\n...")
    }
    s <- gsub("[ \t\n]", "", s)                 # remove white space

    # add newlines and prefix (spaces prefix all lines except the first)
    spaces. <- sprintf("%*s", nchar(msg), " ")   # nchar spaces

    s <- gsub(",", ", ", s)                     # replace comma with comma space
    s <- paste(s, collapse=paste("\n", spaces., sep=""), sep="")
    cat0(msg, s, "\n")
}
match.choices1 <- function(arg, arg.name, choices) # returns the expanded arg
{
    if(!is.character(arg) || length(arg) != 1)
         stop0("illegal \"", arg.name, "\" argument.\n",
               "Choose one of: ", paste.quoted.names(choices))
    imatch <- pmatch(arg, choices)
    if(any(is.na(imatch))) {
        imatch <- NULL
        for(i in seq_along(choices))
            if(pmatch(arg, choices[i], nomatch=0))
                imatch <- c(i, imatch)
        if(length(imatch) == 0)
            stop0(arg.name, "=\"", arg, "\" is not allowed.\n",
                  "Choose one of: ", paste.quoted.names(choices))
        if(length(imatch) > 1)
            stop0(paste0(arg.name, "=\"", arg, "\" is ambiguous.\n",
                         "Choose one of: ", paste.quoted.names(choices)))
    }
    choices[imatch]
}
match.choices <- function(arg, choices) # choices is a vector of strings
{
    arg.name <- deparse(substitute(arg))
    match.choices1(arg, arg.name, choices)
}

match.arg1 <- function(arg) # returns the expanded arg
{
    formal.args <- formals(sys.function(sys.parent()))
    arg.name <- deparse(substitute(arg))
    match.choices1(arg[1], arg.name, eval(formal.args[[arg.name]]))
}
get.caption.from.call <- function(caption, object)
{
    if(is.null(caption) && !is.null(object$call))
        caption <- strip.white.space(paste0(deparse(object$call), collapse=""))
    caption
}
# the make.space functions should only be called if do.par is TRUE
# (otherwise par is not restored correctly)

make.space.for.caption <- function(caption="CAPTION")
{
    oma <- par("oma")
    needed <- 3
    # adjust for newlines in caption
    newlines <- grep("\n", caption)
    if(length(newlines) > 0)
        needed <- needed + .5 * newlines # .5 seems enough although 1 in theory
    if(!is.null(caption) && any(nchar(caption)) && oma[3] <= needed) {
        oma[3] <- needed
        par(oma=oma)
    }
}
# Call this only after a plot is on the screen to avoid
# an error message "plot.new has not been called yet"

show.caption <- function(caption, trim=TRUE, show=TRUE, cex=1)
{
    if(!is.null(caption) && any(nchar(caption))) {
        if(trim) {
            # trim each line of caption to fit
            caption <- strsplit(caption, "\n")[[1]]
            for(i in seq_along(caption)) {
                nchar.org <- nchar(caption[i])
                caption[i] <- substr(caption[i], 1, 60)
                if(nchar(caption[i]) < nchar.org) # append ellipsis if chars deleted
                    caption[i] <- paste0(caption[i], " ...")
            }
            caption <- paste(caption, collapse="\n")
        }
        if(show)
            mtext(caption, outer=TRUE, font=2, line=1, cex=cex)
    }
    caption
}
# like text, but with a white background
# TODO sign of adj is backwards?

text.on.white <- function(x, y, label, cex=1, adj=.5, font=1, ...)
{
    stopifnot(length(label) == 1)
    width       <- strwidth(label, cex=cex, font=font)
    char.width  <- strwidth("X", cex=cex, font=font)
    height      <- strheight(label, cex=cex, font=font)
    char.height <- strheight("X", cex=cex, font=font)
    if(length(adj) == 1)
        adj <- c(adj, .5)
    rect(x - adj[1]         * width  - .2 * char.width,
         y - 1 * adj[2]     * height - .5 * char.height,
         x + (1-adj[1])     * width  + .2 * char.width,
         y + 1 * (1-adj[2]) * height + .1 * char.height,
         col="white", border=NA)
    text(x=x, y=y, labels=label, cex=cex, adj=adj, font=font, ...)
}
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

print.matrix.info <- function(xname, x, Callers.name=NULL, bpairs=NULL,
                              details=TRUE, all.rows=FALSE, all.names=FALSE,
                              prefix=TRUE)
{
    if(prefix) {
        if(!is.null(Callers.name))
            cat0(Callers.name, ": ")
        cat0(xname, " is a ", NROW(x), " by ", NCOL(x), " matrix: ")
        if(is.null(bpairs))
            bpairs <- repl(TRUE, NCOL(x))
        colnames <- colnames(x)
        if(is.null(colnames))
            colnames <- repl("", NCOL(x))
        stopifnot(length(bpairs) == length(colnames))
        icol <- 0
        n.names <- length(colnames)
        if(!all.names && !details)
            n.names <- min(5, n.names)
        else # never print more than 20 columns
            n.names <- min(20, n.names)
        for(i in 1:n.names) {
            if(bpairs[i]) {
                icol <- icol+1
                cat0(icol, "=")
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
            if(!(class(xi)[1] == "numeric" || class(xi)[1] == "matrix"))
                cat0(" (", class(xi)[1], ")")
            else if(typeof(xi) != "double")
                cat0(" (", mode(xi), ")")
            if(!bpairs[i])
                cat(")")
            if(i < length(colnames) && (i == length(colnames) || bpairs[i+1]))
                cat(", ")
        }
        if(i < length(colnames))
            cat("...")
        cat("\n")
    }
    if(details) {
        rownames(x) <- NULL
        ncol.save <- 0
        if(NCOL(x) > 20) { # never print more than 20 columns
            ncol.save <- ncol(x)
            x <- x[, 1:20]
        }
        if(all.rows || NROW(x) <= 6) { # head prints 6 rows
            cat0("Contents of ", xname, " are\n")
            print(x)
        } else {
            rowstring <- if(class(x)[1] == "numeric" || class(x)[1] == "factor")
                             "elements"
                         else
                             "rows"
            cat0("First few ", rowstring, " of ", xname, " are\n")
            print(head(x))
        }
        if(ncol.save)
            printf("Not all %d columns were printed\n", ncol.save)
    }
}
# If xlim[1] == xlim[2], then plot() issues a message.  We don't want that,
# so use this function to make sure xlim[2] is different to xlim[1].
# We also use this function for ylim.

fix.lim <- function(xlim)
{
    if(!is.null(xlim)) {
        # constants below are arbitrary
        min <- -.001
        max <- if(is.null(xlim)) .001 else xlim[2]
        small <- max(1e-6, .001 * abs(xlim[1]), .001 * abs(xlim[2]))
        if(abs(xlim[2] - xlim[1]) < small) # illegal xlim?
            xlim <- c(xlim[1] - small, xlim[2] + small)
    }
    xlim
}
# We allow 20% of x to be nonpositive, useful if the response is essentially
# positive, but the predicted response has a few nonpositive values at the extremes.
# Needed for example if we will later take log(x) or sqrt(x).

check.that.most.are.positive <- function(x, xname, user.arg, non.positive.msg, frac.allowed=.2)
{
    check.numeric.scalar(frac.allowed)
    stopifnot(frac.allowed >= 0 && frac.allowed <= 1)
    nonpos <- x <= 0
    if(sum(nonpos, na.rm=TRUE) > frac.allowed * length(x)) { # more than frac.allowed nonpos?
        ifirst <- which(nonpos)[1]
        stop0(sprintf(
                "%s is not allowed because too many %ss are %s\n",
                user.arg, xname, non.positive.msg),
              sprintf(
                "       %.2g%% are %s (%g%% is allowed)\n",
                100 * sum(nonpos) / length(x), non.positive.msg, 100 * frac.allowed),
               sprintf("       e.g. %s[%d] is %g", xname, ifirst, x[ifirst]))
    }
}
get.cex.points <- function(npoints, len)
{
    n <- if(npoints <= 0) len else min(npoints, len)

    cex.points <-
        if     (n >= 5000) .2
        else if(n >= 3000) .4
        else if(n >= 1000) .6
        else if(n >= 300)  .8
        else               1

    cex.points / par("cex")
}
check.deprecated <-function(new, new.is.missing, old)
{
    if(length(old) != 1 || !is.na(old)) {
        new.name <- deparse(substitute(new))
        old.name <- deparse(substitute(old))
        warnf("%s is deprecated.  Please use %s instead.", old.name, new.name)
        if(!new.is.missing)
            warnf("%s and %s both specified.  Please use just %s.",
                  old.name, new.name, new.name)
        new <- old
    }
    new
}
check <- function(x, xname, check.name, check.func, allow.na=FALSE)
{
    any <- check.func(x)
    if(allow.na)
        any <- any[!is.na(any)]
    else {
        which.na <- which(is.na(any))
        if(length(which.na)) {
            stopf("NA in %s\n       %s[%d] is %g",
                  xname, xname, which.na[1], x[which.na[1]])
        }
    }
    if(any(any)) {
        which <- which(check.func(x))
        stopifnot(length(which) > 0)
        stopf("%s in %s\n       %s[%d] is %g",
              check.name, xname, xname, which[1], x[which[1]])
    }
}
check.vec <- function(x, xname, expected.len=NA, allow.logical=TRUE, allow.na=FALSE)
{
    if(!(NROW(x) == 1 || NCOL(x) == 1))
        stop0(xname, " is not a vector\n       ",
              xname, " has dimensions ", NROW(x), " by ", NCOL(x))
    if(!((allow.logical && is.logical(x)) || is.numeric(x)))
        stop0(xname, " is not numeric")
    if(!is.na(expected.len) && length(x) != expected.len)
        stop0(xname, " has the wrong length ",
              length(x), ", expected ", expected.len)
    if(allow.na)
        x[is.na(x)] <- 1 # prevent check is.finite from complaining
    else
        check(x, xname, "NA", is.na)
    check(x, xname, "non-finite value", function(x) {!is.finite(x)})
}
par.for.plot <- function(do.par, nfigs, cex.main, caption="CAPTION", right.axis=FALSE)
{
    stopifnot(length(do.par) == 1 &&
              (do.par == 0 || do.par == 1 || do.par == 2))
    stopifnot(nfigs > 0)
    old.par <- par(no.readonly=TRUE)
    if(do.par) {
        nrows <- ceiling(sqrt(nfigs))
        par(mfrow=c(nrows, nrows))
        par(mar=c(4, 4, 2, if(right.axis) 3 else 1)) # small margins to pack figs in
        par(mgp=c(1.6, .6, 0))                       # flatten axis elements
        par(cex.main=cex.main)
        make.space.for.caption(caption)
    }
    old.par
}
get.response.given.formula <- function(formula, data, subset)
{
    call. <- match.call(expand.dots=FALSE)
    mf <- call.[c(1, match(c("formula", "data"), names(call.), 0))]
    mf[[1]] <- as.name("model.frame")
    y <- try(model.response(eval.parent(mf), "any"))
    if(is.try.error(y))
        stop0("cannot get the model response from newdata")
    if(is.null(y)) # probably not needed
        stop0("cannot get the model response from newdata")
    as.matrix(y)
}
# the model was created with the x,y interface (no formula)

get.response.given.xy.model <- function(object, newdata)
{
    stopifnot(is.null(object$terms)) # shouldn't be here if model has formula
    colnames.newdata <- colnames(newdata)
    if(is.null(colnames.newdata))
        stop0("cannot get response from newdata because newdata has no column names\n",
              "Possible remedy: use a formula when building the model")
    fitted <- fitted(object)
    response.name <- colnames(fitted)
    if(is.null(response.name)) {
        y <- plotmo::get.plotmo.y.wrapper(object, parent.env, 1, NROW(fitted), trace=01)
        response.name <- colnames(y)
    }
    if(is.null(response.name))
        stopf("%s\n%s\n%s",
              "cannot get response from newdata",
              "Remedy 1: use a formula when building the model",
              "Remedy 2: use y with a column name when building the model")
    if(length(response.name) != 1)
        stop0("multiple response models are not supported here")
    which <- which(colnames.newdata == response.name)
    if(length(which) == 0)
        stop0("No column names in newdata match the original response name\n",
              sprintf("       Response name: %s\n", response.name),
              "       Column names in newdata: ", paste.with.space(colnames.newdata))
    if(length(which) > 1)
        stopf("multiple column names in newdata match the original response name %s",
              response.name)
    y <- as.matrix(newdata[, colnames.newdata[which]])
    stopifnot(!is.null(y))
    y <- as.matrix(y)
    colnames(y) <- response.name
    y
}
# extract the response column from newdata, using the model object
# to figure out which column is the reponse column

get.response <- function(object, newdata)
{
    stopifnot(!is.null(newdata))
    if(is.null(object$terms))
        get.response.given.xy.model(object, newdata)
    else
        get.response.given.formula(formula(object$terms), newdata)
}

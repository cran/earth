# print.earth.R: functions for summarizing and printing earth objects

# print.earth's first arg is actually a model object but called x for consistency with generic

print.earth <- function(x, print.glm=TRUE, digits=getOption("digits"), fixed.point=TRUE, ...)
{
    form <- function(x, pad)
    {
        sprint("%-*s", digits+pad,
                format(if(abs(x) < 1e-20) 0 else x, digits=digits))
    }
    #--- print.earth starts here
    check.classname(x, substitute(x), "earth")
    warn.if.dots(...)
    nresp <- NCOL(x$coefficients)
    is.cv <- !is.null(x$cv.list)
    nselected <- length(x$selected.terms)
    if(is.null(x$glm.list))     # glm.list is a list of glm models, null if none
        cat("Selected ")
    else
        cat("Earth selected ")  # remind user that these are for the earth not glm model

    cat(length(x$selected.terms), "of", nrow(x$dirs), "terms, and",
        get.nused.preds.per.subset(x$dirs, x$selected.terms),
        "of", ncol(x$dirs), "predictors")

    if(x$pmethod == "cv")
        printf(" using pmethod=\"cv\"\n")
    else if(x$pmethod == "none")
        printf(" (pmethod=\"none\")\n")
    else
        cat("\n")

    print.termcond(x)             # print reason we terminated the forward pass
    print.one.line.evimp(x)       # print estimated var importances on a single line
    # "try" below is paranoia so don't completely blow out if problem in print.offset
    try(print.offset(x, digits))  # print offset term in formula, if any
    try(print.weights(x, digits)) # print case weights, if any

    nterms.per.degree <- get.nterms.per.degree(x, x$selected.terms)
    cat("Number of terms at each degree of interaction:", nterms.per.degree)

    cat0(switch(length(nterms.per.degree),
           " (intercept only model)",
           " (additive model)"),
           "\n")

    if(nresp > 1)
        print.earth.multiple.response(x, digits=digits, fixed.point=fixed.point, ...)
    else
        print.earth.single.response(x, digits=digits, fixed.point=fixed.point, ...)

    if(print.glm) {
        if(!is.null(x$glm.list))
            print.earth.glm(x, digits, fixed.point)
        if(x$pmethod == "cv") {
            # The digits +1 below is an attempt to be compatible with the above prints.
            # The number of digits won't always match exactly, it's not critical.
            print.would.have(x, if(nresp > 1) digits+1 else digits)
        }
    }
    invisible(x)
}
print.earth.single.response <- function(x, digits, fixed.point, ...)
{
    is.cv <- !is.null(x$cv.list)
    spacer <- if(is.cv) "  " else "    "
    nselected <- length(x$selected.terms)
    if(!is.null(x$glm.list))
        cat("Earth ")   # remind user
    if(x$pmethod == "cv") {
        ilast <- nrow(x$cv.oof.rsq.tab)
        cat0("GRSq ", format(x$grsq, digits=digits),
             spacer, "RSq ",  format(x$rsq,  digits=digits),
             spacer, "mean.oof.RSq ",
             format(x$cv.oof.rsq.tab[ilast,nselected], digits=digits),
             " (sd ",
             format(sd(x$cv.oof.rsq.tab[-ilast,nselected], na.rm=TRUE), digits=3),
             ")\n")
    } else {
        ilast <- nrow(x$cv.rsq.tab)
        cat0("GCV ",   format(x$gcv,  digits=digits),
             spacer, "RSS ",  format(x$rss,  digits=digits),
             spacer, "GRSq ", format(x$grsq, digits=digits),
             spacer, "RSq ",  format(x$rsq,  digits=digits))
        if(is.cv)
            cat0(spacer, "CVRSq ",
                 format(x$cv.rsq.tab[ilast,ncol(x$cv.rsq.tab)], digits=digits))
        cat("\n")
    }
}
print.earth.multiple.response <- function(x, digits, fixed.point, ...)
{
    nresp <- NCOL(x$coefficients)
    stopifnot(nresp > 1)
    is.cv <- !is.null(x$cv.list)
    nselected <- length(x$selected.terms)
    # create a data.frame and print that
    mat <- matrix(nrow=nresp+1, ncol=4 + is.cv + (x$pmethod == "cv"))
    rownames(mat) <- c(colnames(x$fitted.values), "All")
    colnames <- c("GCV", "RSS", "GRSq", "RSq")
    if(is.cv)
        colnames <- c(colnames,
            if(x$pmethod=="cv") c("  mean.oof.RSq", "sd(mean.oof.RSq)")
            else "CVRSq")
    colnames(mat) <- colnames
    for(iresp in seq_len(nresp)) {
        mat[iresp,1:4] <- c(x$gcv.per.response[iresp],
                            x$rss.per.response[iresp],
                            x$grsq.per.response[iresp],
                            x$rsq.per.response[iresp])
        if(is.cv) {
            if(x$pmethod == "cv") {
                ilast <- nrow(x$cv.oof.rsq.tab)
                mat[iresp,5] <- NA
                mat[iresp,6] <- NA
            } else {
                ilast <- nrow(x$cv.rsq.tab)
                mat[iresp,5] <- x$cv.rsq.tab[ilast,iresp]
            }
        }
    }
    # final row for "All"
    mat[nresp+1,1:4] <- c(x$gcv, x$rss, x$grsq, x$rsq)
    if(is.cv) {
        if(x$pmethod == "cv") {
            mat[nresp+1,5] <- signif(x$cv.oof.rsq.tab[ilast,nselected],
                                   digits=digits)
            mat[nresp+1,6] <-
                signif(sd(x$cv.oof.rsq.tab[-ilast,nselected], na.rm=TRUE),
                       digits=digits)
        } else
            mat[nresp+1,5] <- x$cv.rsq.tab[ilast,ncol(x$cv.rsq.tab)]
    }
    cat("\n")
    if(!is.null(x$glm.list))
        cat("Earth\n") # remind user
    if(fixed.point)
       mat <- my.fixed.point(mat, digits)
    df <- as.data.frame(mat)
    # the following converts the matrix from numeric to character
    df[is.na(df)] <- "" # print NAs as blanks (in mean.oof.RSq column)
    print(df, digits=digits)
}
print.termcond <- function(object) # print reason we terminated the forward pass
{
    printf("Termination condition: ")
    if(is.null(object$termcond)) {
        printf("Unknown\n") # model was created by mars.to.earth
        return()
    }
    termcond <- object$termcond
    check.numeric.scalar(termcond)
    nk <- object$nk
    check.numeric.scalar(nk)
    nterms.before.pruning <- nrow(object$dirs)
    check.numeric.scalar(nterms.before.pruning)
    thresh <- object$thresh
    check.numeric.scalar(thresh)
    terms.string = if(nterms.before.pruning == 1) "term" else "terms"
    if(termcond == 1)
        printf("Reached nk %d\n", nk)
    else if(termcond == 2)
        printf("GRSq -Inf at %d %s\n", nterms.before.pruning, terms.string)
    else if(termcond == 3)
        printf("GRSq -10 at %d %s\n", nterms.before.pruning, terms.string)
    else if(termcond == 4)
        printf("RSq changed by less than %g at %d %s\n",
            thresh, nterms.before.pruning, terms.string)
    else if(termcond == 5)
        printf("Reached maximum RSq %.4f at %d %s\n",
               1-thresh, nterms.before.pruning, terms.string)
    else if(termcond == 6)
        printf("No new term increases RSq at %d %s\n",
               nterms.before.pruning, terms.string)
    else if(termcond == 7)
        printf("Reached nk %d\n", nk)
    else
        printf("Unknown (termcond %d)\n", termcond) # should never happen
}
# The first arg is actually an object but called x for consistency with generic

print.summary.earth <- function(
    x            = stop("no 'x' argument"),     # "summary.earth" object
    details      = x$details,
    decomp       = x$decomp,
    digits       = x$digits,
    fixed.point  = x$fixed.point,
    newdata      = x$newdata,
    ...)
{
    nresp <- NCOL(x$coefficients)
    warn.if.dots(...)
    if(!is.null(newdata)) {
        # print short summary on newdata
        printf("RSq %.3f on newdata (%d cases)\n", x$newrsq, NROW(newdata))
        if(!is.null(x$varmod)) {
            printf("\n")
            print.varmod(x$varmod, newdata=newdata, digits=digits)
        }
        return(invisible(x))
    }
    printcall("Call: ", x$call)
    cat("\n")
    is.glm <- !is.null(x$glm.list)   # TRUE if embedded GLM model(s)
    new.order <- reorder.earth(x, decomp=decomp)

    if(!is.glm || details)
        print.earth.coefficients(x, digits, fixed.point, new.order)
    if(is.glm)
        print.summary.earth.glm(x, details, digits, fixed.point, new.order)
    if(details)
        cat0("Number of cases: ", nrow(x$residuals), "\n")
    print.earth(x, digits, print.glm=FALSE)
    if(!is.null(x$glm.list))
        print.earth.glm(x, digits, fixed.point)
    if(!is.null(x$cv.list) && x$pmethod != "cv")
        print.cv(x)
    if(x$pmethod == "cv")
        print.would.have(x, if(nresp > 1) digits+1 else digits)
    if(!is.null(x$varmod)) {
        printf("\nvarmod: ")
        x$varmod <- print.varmod(x$varmod, digits=digits)
    }
    invisible(x)
}
print.offset <- function(object, digits)
{
    # currently only earth.formula objects can have an offset
    terms <- object$terms
    if(is.null(terms))
        return()
    offset.index <- attr(terms, "offset")
    if(is.null(offset.index) || is.null(object$offset))
        return()
    varnames <- rownames(attr(terms, "factors"))
    if(is.null(varnames) || offset.index < 1 || offset.index > length(varnames))
        return()
    offset.term <- varnames[offset.index]
    # convert "offset(foo)" to "foo"
    offset.term.name <- substring(offset.term, 8, nchar(offset.term)-1)
    print.values("Offset", offset.term.name, object$offset, digits)
}
print.weights <- function(object, digits)
{
    if(is.null(object$weights))
        return()
    print.values("Weights", NULL, object$weights, digits)
}
print.values <- function(name, term.name, values, digits)
{
    s <- if(is.null(term.name))
            paste0(name, ": ")
         else
            paste0(name, ": ", strip.space.collapse(term.name), " with values ")
    cat0(s)
    n <- min(30, length(values)) # save time by formatting a max of only 30 values
    stopifnot(n > 0)
    svalues <- character(n)
    maxlen <- max(25, getOption("width") - nchar(s) - 5) # -5 for ", ..."
    s <- if(!is.null(term.name) && substring(term.name, 1, 4) == "log(") { # common case
            # format each value individually (not aggregrated like format on a vector)
            # this is nice for example if just one weight is very small
            for(i in seq_along(svalues))
                svalues[i] <- strip.space(format(exp(values[i]), digits=digits))
            paste.trunc("log(", svalues, ")",
                        sep="", collapse=", ", maxlen=maxlen)
         } else {
            for(i in seq_along(svalues))
                svalues[i] <- strip.space(format(values[i], digits=digits))
            paste.trunc(svalues, sep="", collapse=", ", maxlen=maxlen)
        }
    cat0(s, "\n")
}
print.would.have <- function(x, digits)
{
    form <- function(x) format(x, digits=digits)
    nselected <- length(x$backward.selected.terms)
    ilast <- nrow(x$cv.oof.rsq.tab)
    cat0("\npmethod=\"backward\" would have selected",
        if(nselected == length(x$selected.terms)) " the same model" else "",
        ":\n    ",
        nselected,
        " terms ",
        get.nused.preds.per.subset(x$dirs, x$backward.selected.terms),
        " preds,  GRSq ",
        form(get.rsq(x$gcv.per.subset[nselected], x$gcv.per.subset[1])),
        "  RSq ",
        form(get.rsq(x$rss.per.subset[nselected], x$rss.per.subset[1])),
        "  mean.oof.RSq ",
        form(x$cv.oof.rsq.tab[ilast, nselected]),
        "\n")
    if(anyNA(x$cv.oof.rsq.tab[ilast, nselected]))
        printf(
"    (mean.oof.RSq is NA because most fold models have less than %g terms)\n",
               nselected)

}
print.earth.coefficients <- function(x, digits, fixed.point, new.order)
{
    nresp <- NCOL(x$coefficients)
    if(!is.null(x$strings)) {       # old style expression formatting?
        resp.names <- colnames(x$fitted.values)
        for(iresp in seq_len(nresp)) {
            cat0(resp.names[iresp], " =\n")
            cat(x$strings[iresp])
            cat("\n")
       }
    } else {
        rownames(x$coefficients) <- spaceout(rownames(x$coefficients))
        coef <- x$coefficients[new.order, , drop=FALSE]
        if(fixed.point)
            coef <- my.fixed.point(coef, digits)
        if(!is.null(x$glm.list))        # embedded GLM model(s)?
            cat("Earth coefficients\n") # remind user what these are
        else if(nresp == 1)
            colnames(coef) = "coefficients"
        print(coef, digits=digits)
        cat("\n")
    }
}
print.summary.earth.glm <- function(x, details, digits, fixed.point, new.order)
{
    nresp <- NCOL(x$coefficients)
    resp.names <- colnames(x$fitted.values)
    if(!is.null(x$strings)) {      # old style expression formatting?
       for(iresp in seq_len(nresp)) {
           g <- x$glm.list[[iresp]]
           cat("GLM ")
           cat0(resp.names[iresp], " =\n")
           cat(x$strings[nresp+iresp]) # glm strings index is offset by nresp
           cat("\n")
       }
    } else {
        cat("GLM coefficients\n")
        rownames(x$glm.coefficients) <- spaceout(rownames(x$glm.coefficients))
        coef <- x$glm.coefficients[new.order, , drop=FALSE]
        if(fixed.point)
            coef <- my.fixed.point(coef, digits)
        print(coef, digits=digits)
        cat("\n")
    }
    if(details)
        for(iresp in seq_len(nresp))
            print.glm.details(x$glm.list[[iresp]], nresp, digits,
                              my.fixed.point, resp.names[iresp])
}
# put some spaces into term names for readability
#     convert h(x1-5860)*h(x2--15)
#     to      h(x1-5860) * h(x2- -15)

spaceout <- function(rownames.)
{
    rownames. <- gsub("\\*", " * ", rownames.)   # spaces around *
    rownames. <- gsub("--", "- -", rownames.)    # spaces between --
    gsub("`", "", rownames.)                     # remove backquotes
}
# TODO Add an inverse.func arg to summary.earth, similar to plotmo.

summary.earth <- function(   # returns a superset, not a summary in the strict sense
    object       = stop("no 'object' argument"),
    details      = FALSE,
    style        = c("h", "pmax", "max", "C", "bf"),
    decomp       = "anova",
    digits       = getOption("digits"),
    fixed.point  = TRUE,
    newdata      = NULL,
    ...) # unused
{
    check.classname(object, substitute(object), "earth")
    details     <- check.boolean(details)
    fixed.point <- check.boolean(fixed.point)
    rval <- object
    rval$strings <- switch(match.arg1(style, "style"),
        "h"    = { stop.if.dots(...) },
        "pmax" = format.earth(x=object, style=style, decomp=decomp, digits=digits, ...),
        "max"  = format.earth(x=object, style=style, decomp=decomp, digits=digits, ...),
        "C"    = format.earth(x=object, style=style, decomp=decomp, digits=digits, ...),
        "bf"   = format.earth(x=object, style=style, decomp=decomp, digits=digits, ...))
    rval$details     <- details  # pass details arg on to print.summary.earth
    rval$decomp      <- decomp
    rval$digits      <- digits
    rval$fixed.point <- fixed.point
    if(!is.null(newdata)) {
        rval$newdata <- newdata
        rval$newrsq  <- plotmo::plotmo_rsq(object, newdata, ...)
    }
    class(rval) <- c("summary.earth", "earth")
    rval
}

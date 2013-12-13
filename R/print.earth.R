# print.earth.R: functions for summarizing and printing earth objects

# print.earth's 1st arg is actually an obj but called x for consistency with generic

print.earth <- function(x, print.glm=TRUE, digits=getOption("digits"), fixed.point=TRUE, ...)
{
    form <- function(x, pad)
    {
        sprintf("%-*s", digits+pad,
                format(if(abs(x) < 1e-20) 0 else x, digits=digits))
    }
    #--- print.earth starts here
    warn.if.dots.used("print.earth", ...)
    check.classname(x, deparse(substitute(x)), "earth")
    if(is.null(x$glm.list))     # glm.list is a list of glm models, null if none
        cat("Selected ")
    else
        cat("Earth selected ")  # remind user that these are for the earth not glm model

    cat(length(x$selected.terms), "of", nrow(x$dirs), "terms, and",
        get.nused.preds.per.subset(x$dirs, x$selected.terms),
        "of", ncol(x$dirs), "predictors ")

    nlinpreds <- sum(x$dirs[x$selected.terms,] == 2)
    if(nlinpreds == 1)
        cat0(nlinpreds, " linear predictor")
    else if(nlinpreds > 1)
        cat0(nlinpreds, " linear predictors")
    cat("\n")

    print.one.line.evimp(x) # print estimated var importances on a single line

    nterms.per.degree <- get.nterms.per.degree(x, x$selected.terms)
    cat("Number of terms at each degree of interaction:", nterms.per.degree)

    cat0(switch(length(nterms.per.degree),
           " (intercept only model)",
           " (additive model)"),
         "\n")

    nresp <- NCOL(x$coefficients)
    is.cv <- !is.null(x$cv.list)
    ilast <- nrow(x$cv.rsq.tab)

    if(nresp > 1) {
        # create a matrix and print that

        a <- matrix(nrow=nresp+1, ncol=4 + is.cv)
        rownames(a) <- c(colnames(x$fitted.values), "All")
        colnames. <- c("GCV", "RSS", "GRSq", "RSq")
        if(is.cv)
            colnames. <- c(colnames., "cv.rsq")
        colnames(a) <- colnames.
        for(iresp in 1:nresp) {
            a[iresp,1:4] <- c(x$gcv.per.response[iresp],
                              x$rss.per.response[iresp],
                              x$grsq.per.response[iresp],
                              x$rsq.per.response[iresp])
            if(is.cv)
                a[iresp,5] <- x$cv.rsq.tab[ilast,iresp]
        }
        # final row for "All"
        a[nresp+1,1:4] <- c(x$gcv, x$rss, x$grsq, x$rsq)
        if(is.cv)
            a[nresp+1,5] <- x$cv.rsq.tab[ilast,ncol(x$cv.rsq.tab)]
        cat("\n")
        if(!is.null(x$glm.list))
            cat("Earth\n") # remind user
        if(fixed.point)
           a <- my.fixed.point(a, digits)
        print(a, digits=digits)
    } else {
        if(!is.null(x$glm.list))
            cat("Earth ")   # remind user
        spacer <- if(is.cv) "  " else "    "
        cat0("GCV ",   format(x$gcv,  digits=digits),
             spacer, "RSS ",  format(x$rss,  digits=digits),
             spacer, "GRSq ", format(x$grsq, digits=digits),
             spacer, "RSq ",  format(x$rsq,  digits=digits))
        if(is.cv)
            cat0(spacer, "cv.rsq ",
                 format(x$cv.rsq.tab[ilast,1], digits=digits))
        cat("\n")
    }
    if(print.glm && !is.null(x$glm.list))
        print.earth.glm(x, digits, fixed.point)

    invisible(x)
}

# The first arg is actually an object but called x for consistency with generic

print.summary.earth <- function(
    x           = stop("no 'x' arg"),     # "summary.earth" object
    details     = x$details,
    decomp      = x$decomp,
    digits      = x$digits,
    fixed.point = x$fixed.point,
    ...)
{
    my.print.call("Call: ", x$call)
    nresp <- NCOL(x$coefficients)
    cat("\n")
    warn.if.dots.used("print.summary.earth", ...)
    is.glm <- !is.null(x$glm.list)   # TRUE if embedded GLM model(s)
    new.order <- reorder.earth(x, decomp=decomp)
    response.names <- colnames(x$fitted.values)

    # print coefficients
    if(!is.glm || details) {
        if(!is.null(x$strings)) {      # old style expression formatting?
            for(iresp in 1:nresp) {
                cat0(response.names[iresp], " =\n")
                cat(x$strings[iresp])
                cat("\n")
           }
        } else {
            rownames(x$coefficients) <- spaceout(rownames(x$coefficients))
            coef <- x$coefficients[new.order, , drop=FALSE]
            if(fixed.point)
                coef <- my.fixed.point(coef, digits)
            if(is.glm)
                cat("Earth coefficients\n") # remind user what these are
            else if(nresp == 1)
                colnames(coef) = "coefficients"
            print(coef, digits=digits)
            cat("\n")
        }
    }
    if(is.glm) {
        if(!is.null(x$strings)) {      # old style expression formatting?
           for(iresp in 1:nresp) {
               g <- x$glm.list[[iresp]]
               cat("GLM ")
               cat0(response.names[iresp], " =\n")
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
        if(details) for(iresp in 1:nresp)
           print.glm.details(x$glm.list[[iresp]], nresp, digits,
                             my.fixed.point, response.names[iresp])
    }
    if(details)
        cat0("Number of cases: ", nrow(x$residuals), "\n")
    print.earth(x, digits, print.glm=FALSE)
    if(!is.null(x$glm.list))
        print.earth.glm(x, digits, fixed.point)
    if(!is.null(x$cv.list))
        print.cv(x)
    invisible(x)
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

summary.earth <- function(   # returns a superset, not a summary in the strict sense
    object  = stop("no 'object' arg"),
    details = FALSE,
    style   = c("h", "pmax", "max", "C", "bf"),
    decomp  = "anova",        # see reorder.earth for legal decomp values
                              # use "pmax" for old earth style earth expression formatting
    digits  = getOption("digits"),
    fixed.point = TRUE,       # see help page
    ...)                      # extra args passed on to format.earth
{
    rval <- object
    rval$strings <- switch(match.arg1(style),
        { warn.if.dots.used("summary.earth", ...); NULL },                      # "h"
        format.earth(x=object, style=style, decomp=decomp, digits=digits, ...), # "pmax"
        format.earth(x=object, style=style, decomp=decomp, digits=digits, ...), # "max"
        format.earth(x=object, style=style, decomp=decomp, digits=digits, ...), # "C"
        format.earth(x=object, style=style, decomp=decomp, digits=digits, ...)) # "bf"
    rval$details     <- details  # ditto
    rval$decomp      <- decomp
    rval$digits      <- digits   # allows us to pass digits arg on to print.summary.earth
    rval$fixed.point <- fixed.point
    class(rval)      <- c("summary.earth", "earth")
    rval
}

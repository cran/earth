# print.earth.R: functions for summarizing and printing earth objects

# print.earth's 1st arg is actually an obj but called x for consistency with generic

print.earth <- function(x, print.glm=TRUE, digits=getOption("digits"), fixed.point=TRUE, ...)
{

    form <- function(x, pad)
    {
        sprintf("%-*s", digits+pad,
                format(if(abs(x) < 1e-20) 0 else x, digits=digits))
    }

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
        cat(nlinpreds, " linear predictor", sep="")
    else if(nlinpreds > 1)
        cat(nlinpreds, " linear predictors", sep="")
    cat("\n")

    print.estimated.predictor.importance(x)

    nterms.per.degree <- get.nterms.per.degree(x, x$selected.terms)
    cat("Number of terms at each degree of interaction:", nterms.per.degree)

    cat(switch(length(nterms.per.degree),
            " (intercept only model)",
            " (additive model)"),
        "\n", sep="")

    nresp <- NCOL(x$coefficients)

    if(nresp > 1) {
        # create a matrix and print that

        a <- matrix(nrow=nresp+1, ncol=4)
        rownames(a) <- c(colnames(x$fitted.values), "All")
        colnames(a) <- c("GCV", "RSS", "GRSq", "RSq")
        for(iresp in 1:nresp)
             a[iresp,] <- c(x$gcv.per.response[iresp],
                            x$rss.per.response[iresp],
                            x$grsq.per.response[iresp],
                            x$rsq.per.response[iresp])

        a[nresp+1,] <- c(x$gcv, x$rss, x$grsq, x$rsq)  # final row for "All"
        cat("\n")
        if(!is.null(x$glm.list))
            cat("earth\n") # remind user
        if(fixed.point)
           a <- my.fixed.point(a, digits)
        print(a, digits=digits)
    } else {
        if(!is.null(x$glm.list))
            cat("earth ")   # remind user
        cat("GCV ",   format(x$gcv,  digits=digits),
            "    RSS ",  format(x$rss,  digits=digits),
            "    GRSq ", format(x$grsq, digits=digits),
            "    RSq ",  format(x$rsq,  digits=digits),
            "\n", sep="")
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
    cat("\n")
    warn.if.dots.used("print.summary.earth", ...)
    nresp <- NCOL(x$coefficients)
    is.glm.model <- !is.null(x$glm.list)   # TRUE if embedded GLM model(s)
    new.order <- reorder.earth(x, decomp=decomp)
    response.names <- colnames(x$coeff)

    # print coefficients

    if(!is.glm.model || details) {
        if(!is.null(x$strings)) {      # old style expression formatting?
           for(iresp in 1:nresp) {
               cat(response.names[iresp], " =\n", sep="")
               cat(x$strings[iresp])
               cat("\n")
           }
        } else {
            if(is.glm.model)
                cat("earth coefficients\n") # remind user what these are
            rownames(x$coefficients) <- spaceout(rownames(x$coefficients))
            coef <- x$coefficients[new.order, , drop=FALSE]
            if(fixed.point)
                coef <- my.fixed.point(coef, digits)
            print(coef, digits=digits)
            cat("\n")
        }
    }
    if(is.glm.model) {
        if(!is.null(x$strings)) {      # old style expression formatting?
           for(iresp in 1:nresp) {
               g <- x$glm.list[[iresp]]
               cat("GLM ")
               cat(response.names[iresp], " =\n", sep="")
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
           print.glm.details(x$glm.list[[iresp]], nresp, digits, my.fixed.point, response.names[iresp])
    }
    if(details)
        cat("Number of cases: ", nrow(x$residuals), "\n", sep="")
    print.earth(x, digits, print.glm=FALSE)
    if(!is.null(x$glm.list))
        print.earth.glm(x, digits, fixed.point)
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
    decomp  = "anova",        # see reorder.earth for legal decomp values
    style   = c("h", "pmax"), # use "pmax" for old earth style earth expression formatting
    digits  = getOption("digits"),
    fixed.point = TRUE,       # see help page
    ...)                      # extra args passed on to format.earth
{
    rval <- object
    style = switch(match.arg1(style),
        NULL,
        rval$strings <- format.earth(x=object, digits=digits, decomp=decomp, style=style))
    rval$details     <- details  # ditto
    rval$decomp      <- decomp
    rval$digits      <- digits   # allows us to pass digits arg on to print.summary.earth
    rval$fixed.point <- fixed.point
    class(rval)      <- c("summary.earth", "earth")
    rval
}

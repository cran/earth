# format.earth.R

# Return a vector s of strings, length of vector is nresponses.
# But if there are embedded GLM model(s) then length of s is 2 * nresponses
# and strings for the GLM model(s) start at s[nresponses].
#
# For each model string, there is one term per line.  Each term (except
# the intercept) is made up of a coefficent which multiplies one or more
# hockey stick funcs.
#
# For the default style="h" the result looks like this:
#
#         23
#         +  5.7 * h(Girth-12.9)
#         -  2.9 * h(12.9-Girth)
#         + 0.72 * h(Height-76)
#
# For style="pmax" the result looks like this:
#
#         23.208244
#         +  5.7459616 * pmax(0,  Girth -   12.9)
#         -  2.8664516 * pmax(0,   12.9 -  Girth)
#         + 0.71833643 * pmax(0, Height -     76)
#
# decomp argument: see reorder.earth()
#
# The first arg is actually an object but called x for consistency with generic
#
# TODO would be nice to add an option to first list bx terms and then combine them
# TODO would be nice to add an option to print in term importance order

format.earth <- function(
    x           = stop("no 'x' arg"),   # "earth" object
    digits      = getOption("digits"),
    use.names   = TRUE,     # use predictor names, else "x[,1]" etc
    decomp      = "anova",  # see reorder.earth for legal decomp values
    style       = "h",      # see get.term.strings
    valid.names = FALSE,    # if TRUE, convert to valid C names
    ...)                    # unused, for consistency with generic
{
    check.classname(x, deparse(substitute(a)), "earth")
    warn.if.dots.used("format.earth", ...)
    nresp <- NCOL(x$coefficients)
    s <- vector(mode = "character", length=nresp)
    for(iresp in 1:nresp)
        s[iresp] <- format.one.response(iresp, x, digits, use.names,
                                        decomp, style=style,
                                        valid.names=valid.names, coefs=NULL)

    if(!is.null(x$glm.list))   # embedded GLM model(s)?
        for(iresp in 1:nresp)
            s[nresp + iresp] <- format.one.response(iresp, x, digits, use.names,
                                        decomp, style=style,
                                        valid.names=valid.names,
                                        coefs=x$glm.list[[iresp]]$coefficients)
    s
}

format.one.response <- function( # called by format.earth
    iresp,          # response index i.e. column in y matrix
    obj,            # "earth" obj
    digits,
    use.names,      # use predictor names, else "x[,1]" etc
    decomp,         # see reorder.earth for legal decomp values
    style,          # see get.term.strings
    valid.names,    # if TRUE, convert to valid C names
    coefs=NULL)     # if not NULL use these instead of obj$coefficients
{
    new.order <- reorder.earth(obj, decomp=decomp)
    if(is.null(coefs))
        coefs <- obj$coefficients[, iresp]
    coefs <- coefs[new.order]
    which.terms <- obj$selected.terms[new.order]
    dirs <- obj$dirs
    check.which.terms(dirs, which.terms)
    term.names <- get.term.strings(obj, digits, use.names, style, new.order)
    if(valid.names) # convert invalid chars to underscore
        term.names <- make.unique(gsub(":", "_", term.names), sep="_")
    coef.width <- get.coef.width(coefs[-1], digits)
    s <- ""         # result goes into this string
    s <- pastef(s, "  %.*g\n", digits=digits, coefs[1])
    iterm <- 2
    while(iterm <= length(which.terms)) {
        coef <- coefs[iterm]
        if(coef < 0)
            s <- pastef(s, "  - %s ",
            format(-coef, justify="left",w=coef.width,digits=digits,format="%g"))
        else
            s <- pastef(s, "  + %s ",
                        format(coef, justify="left",
                               w=coef.width,digits=digits,format="%g"))
        s <- pastef(s, "* %s", term.names[iterm])
        s <- pastef(s, "\n")
        iterm <- iterm + 1
    }
    s
}

get.coef.width <- function(coefs, digits)   # get print width for earth coefs
{
    if(length(coefs) > 0)
        max(nchar(format(abs(coefs), digits=digits)))
    else
        10 # arbitrary width if no coefs
}

# style argument:
#   "h"    gives "h(survived-0) * h(16-age)"
#   "pmax" gives "pmax(0, survived - 0) * pmax(0, 16 - age)"

get.term.strings <- function(obj, digits, use.names, style = c("h", "pmax"), neworder)
{
    style = switch(match.arg1(style),
        get.term.strings.h(obj, digits, use.names, neworder),    # "h"
        get.term.strings.pmax(obj, digits, use.names, neworder)) # "pmax"
}

get.term.strings.h <- function(obj, digits, use.names, new.order)
{
    # digits is unused

    if(!use.names)
        warning("use.names=FALSE ignored because style=\"h\"")

    s <- colnames(obj$bx)[new.order]
}

# TODO need to add factor simplification to this routine

get.term.strings.pmax <- function(obj, digits, use.names, new.order)
{
    # get.width returns the width for printing elements of the earth expression.
    # This is used to keep things lined up without too much white space.
    # This returns the widest of all possible printed elements.

    get.width <- function(which.terms, dirs, var.names, cuts, digits)
    {
        if(length(which.terms) == 1)
            return(10)  # return arbitrary width for intercept only model
        used.dirs <- dirs[which.terms, , drop=FALSE]
        # used.preds is a logical index vector which selects used x predictors
        used.preds <- apply(used.dirs, 2, any1)
        # as.list is needed so format treats each cut independently
        max(nchar(var.names[used.preds]),
            nchar(format(as.list(cuts[which.terms, used.preds]), digits=digits)))
    }
    which.terms <- obj$selected.terms[new.order]
    cuts <- obj$cuts
    var.names <- variable.names.earth(obj, use.names=use.names)
    which.terms <- obj$selected.terms[new.order]
    dirs <- obj$dirs
    width <- get.width(which.terms, dirs, var.names, cuts, digits)
    nterms <- length(which.terms)
    s <- character(nterms)
    s[1] = "(Intercept)"
    iterm <- 2
    while(iterm <= nterms) {
        isel.term <- which.terms[iterm]
        dir <- dirs[isel.term, , drop=FALSE]
        cut <- cuts[isel.term, , drop=FALSE]
        npreds <- ncol(cuts)
        prefix <- ""
        for(ipred in 1:npreds) {
            if(dir[ipred]) {
                if(dir[ipred] == 2)     # linear predictor?
                    s[iterm] <- pastef(s[iterm], "%s%-*s %*s            ",
                                    prefix, width=width,
                                    var.names[ipred], width=width, "")
                else if(dir[ipred] == -1)
                    s[iterm] <- pastef(s[iterm], "%spmax(0, %s - %*s) ",
                                     prefix, format(cut[ipred],
                                     width=width, digits=digits),
                                     width, var.names[ipred])
                else if(dir[ipred] == 1)
                    s[iterm] <- pastef(s[iterm], "%spmax(0, %*s - %s) ",
                                     prefix, width=width, var.names[ipred],
                                     format(cut[ipred], width=width,
                                     digits=digits))
                else
                    stop1("illegal direction ", dir[ipred], " in 'dirs'")

                prefix <- "* "
            }
        }
        iterm <- iterm + 1
    }
    s
}

# Return a string representing the linear model.
# Example: a <- lm(Volume ~ ., data = trees); cat(format(a))
# which yields:
#
#   -58
#   +  4.71 * Girth
#   + 0.339 * Height
#
# The first arg is actually an object but called x for consistency with generic
#
# TODO this function doesn't really belong in the earth package

format.lm <- function(
  x           = stop("no 'x' arg"),    # "lm" object, also works for "glm" objects
  digits      = getOption("digits"),
  use.names   = TRUE,
  valid.names = FALSE,                 # if TRUE, convert to valid C names
  ...)                                 # unused, for consistency with generic
{
  get.name <- function(ipred)
  {                                    # return "name" if possible, else "x[,i]"
      pred.name <- pred.names[ipred]
      if(is.null(pred.name) || is.na(pred.name) || !use.names)
          pred.name <- paste("x[,", ipred, "]", sep="")
      if(valid.names) # convert invalid chars to underscore
          pred.name <- make.unique(gsub(":", "_", pred.name), sep="_")
      pred.name
  }
  format1 <- function(coef)
  {
      format(coef, justify="left", w=coef.width, digits=digits, format="%g")
  }
  check.classname(x, deparse(substitute(x)), "lm")
  dataClasses <- attr(x$terms, "dataClasses")
  if(any((dataClasses == "factor") | (dataClasses == "ordered")))
      stop("a predictor has class 'factor' and format.lm cannot handle that")
  coefs <- coef(x)
  if(!is.vector(coefs) || NCOL(coefs) > 1)
      stop("format.lm can only handle single response models")
  intercept <- 0;
  pred.names <- variable.names(x)
  intercept.index <- match("(Intercept)", names(coefs), nomatch=0)
  if(intercept.index) {
      stopifnot(intercept.index == 1)
      intercept <- coefs[1]
      pred.names <- variable.names(x)[-1] # drop intercept
      coefs <- coefs[-1]
  }
  s <- sprintf("  %.*g\n", digits=digits, intercept)
  coef.width <- get.coef.width(coefs, digits)
  for(ipred in seq_along(coefs)) {
      coef <- coefs[ipred]
      if(coef < 0)
          s <- pastef(s, "  - %s ", format1(-coef))
      else
          s <- pastef(s, "  + %s ", format1(coef))
      s <- pastef(s, "* %s", get.name(ipred))
      s <- pastef(s, "\n")
  }
  s
}

# earth.methods.R: miscellaneous earth methods

anova.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning0("anova.earth: returning NULL")
    NULL
}
case.names.earth <- function(object, ...)
{
    if(is.null(row.names(object$residuals)))
        paste(seq_len(nrow(object$residuals)))
    else
        row.names(object$residuals)
}
coef.earth <- function(object, type  = c("response", "earth", "glm"), decomp="none", ...)
{
    warn.if.dots(...)
    type <- match.arg1(type, "type")
    coef <- object$glm.coefficients
    if(is.null(coef)) { # not a glm model?
        if(type == "glm")
            stop0("coef.earth: type == \"glm\" is not allowed because this is not an earth-glm model")
        coef <- object$coefficients
    } else if(type == "earth")
        coef <- object$coefficients
    if(NCOL(coef) > 1) {
        warning0("coef.earth: multiple response model: returning coefficients for just the first response")
        coef <- coef[,1,drop=FALSE]
    }
    new.order <- reorder.earth(object, decomp=decomp)
    names <- spaceout(rownames(coef))
    coef <- coef[new.order,]
    names(coef) <- names
    coef
}
deviance.earth <- function(object, warn=TRUE, ...)
{
    if(warn && !is.null(object$glm.list))
        warning0("deviance.earth: returning earth (not GLM) deviance")
    object$rss
}
effects.earth <- function(object, warn=TRUE, ...)
{
    if(warn)
        warning0("effects.earth: returning NULL")
    NULL
}
# Fake the AIC by returning the GCV.  This is enough for step() to work.
extractAIC.earth <- function(fit, scale = 0, k = 2, warn=TRUE, ...)
{
    if(warn)
        warning0("extractAIC.earth: returning GCV instead of AIC")
    if(scale != 0)
        warning0("extractAIC.earth: ignored scale parameter ", scale)
    if(k != 2)
        warning0("extractAIC.earth: ignored k parameter ", k)
    warn.if.dots(...)
    nterms <- length(fit$selected.terms)
    c(effective.nbr.of.params(nterms, get.nknots(nterms), fit$penalty), fit$gcv)
}
family.earth <- function(object, ...)
{
    stopifnot(!is.null(object$glm.list))
    family(object$glm.list[[1]])
}
hatvalues.earth <- function(model, ...)
{
    stop.if.dots(...)
    if(is.null(model$leverages))
        stop0("this earth model does not have leverages")
    model$leverages
}
fitted.earth <- function(object, type="response", ...)
{
    predict.earth(object, newdata=NULL, type=type, ...)
}
fitted.values.earth <- function(object, type="response", ...)
{
    predict.earth(object, newdata=NULL, type=type, ...)
}
# use.names can have the following values:
#   TRUE:  return name if possible, else return x[,i] or x[i-1].
#   FALSE: return x[,i]
#   -1:    return x[i] with 0 based indexing (treat x as a C array)

variable.names.earth <- function(object, ..., use.names=TRUE)
{
    warn.if.dots(...)
    ipred <- seq_len(ncol(object$dirs))
    if(length(use.names) != 1)
        stop0("illegal value for use.names")
    if(use.names == TRUE) {
        varname <- colnames(object$dirs)[ipred]
        if(!is.null(varname) && !anyNA(varname))
            varname
        else
            paste0("x[,", ipred, "]")
    } else if(use.names == FALSE)
        paste0("x[,", ipred, "]")
    else if(use.names == -1)
        paste0("x[", ipred-1, "]")
    else
        stop0("illegal value for use.names \"", use.names, "\"")
}
weights.earth <- function(object, ...)
{
    warn.if.dots(...)
    if(is.null(object$weights)) # weights arg to earth was NULL?
        repl(1, length(object$fitted.values[,1]))
    else
        object$weights
}

# predict.earth.R

predict.earth <- function(
    object   = stop("no 'object' argument"),
    newdata  = NULL,
    type     = c("link", "response", "earth", "class", "terms"),
    interval = "none",
    level    = .95,     # only used if interval != none
    thresh   = .5,      # only used if type="class"
    trace    = FALSE,
    ...)                # unused, for compatibility with generic predict
{
    check.classname(object, substitute(object), "earth")
    warn.if.dots(...)
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    type <- match.arg1(type, "type")
    env <- parent.frame() # the environment from which predict.earth was called
    if(type == "terms")
        fit <- predict_earth_terms(object, newdata, env, trace)
    else
        fit <- predict_earth_aux(object, newdata, env, type, thresh, trace)
    interval <- match.choices(interval,
                    c("none", "pint", "cint", "se", "abs.residual"), "interval")
    if(interval == "none") {
        if(!missing(level))
            stop0("predict.earth: level=", as.char(level),
                  " was specified but interval=\"none\"")
        return(fit) # note return
    }
    # the interval argument was used
    if(is.null(object$varmod))
        stop0("no prediction intervals because ",
              "the earth model was not built with varmod.method")
    if(type == "class" || type == "terms")
        stop0("predict.earth: the interval argument is not allowed ",
              "with type=\"", type, "\"")
    if(!is.null(object$glm.list) && type != "earth")
        stop0("predict.earth: with earth-glm models, use type=\"earth\" ",
              "when using the interval argument")
    if(NCOL(fit) != 1)
        stop0("predict.earth: the interval argument is not supported ",
              "for multiple response models")
    predict.varmod(object$varmod, newdata=newdata, type=interval, level=level)
}
predict_earth_aux <- function(object, newdata, env, type, thresh, trace)
{
    type.is.class <- type=="class"
    if(type.is.class) {
        type <- "response" # we want predicted probabilities
        ylevels <- object$levels
        if(is.null(ylevels))
            ylevels <- c(FALSE, TRUE)
    }
    if(is.null(newdata)) # no newdata?
        fit <- predict_earth_without_newdata(object, type, trace)
    else # user supplied newdata
        fit <- predict_earth_with_newdata(object, newdata, env, type, trace)

    if(type.is.class)
        fit <- convert.predicted.response.to.class(fit, ylevels,
                            colnames(object$coefficients)[1], thresh)
    fit
}
predict_earth_without_newdata <- function(object, type, trace)
{
    if(is.null(object$glm.list) || type=="earth") {
        print_returning_earth(object, trace, "fitted.values")
        fit <- object$fitted.values
    } else {    # glm predictions
        trace1(trace, "predict.earth: returning glm fitted.values\n")
        fit <- matrix(0, nrow=nrow(object$fitted.values),
                         ncol=ncol(object$fitted.values))
        colnames(fit) <- colnames(object$fitted.values)
        for(i in seq_along(object$glm.list))
            fit[,i] = predict.glm(object$glm.list[[i]], type=type)
    }
    fit
}
predict_earth_with_newdata <- function(object, newdata, env, type, trace)
{
    bx <- model.matrix.earth(object, newdata, trace=trace, Env=env,
                             Callers.name="model.matrix.earth from predict.earth")
    if(trace >= 1) {
        print_summary(bx, "predict.earth with newdata: bx", trace=2)
        trace2(trace, "\n")
    }
    offset <- get.predict.offset(object, newdata, trace)
    if(is.null(object$glm.list) || type=="earth") {
        print_returning_earth(object, trace, "predictions")
        fit <- bx %*% object$coefficients
        if(!is.null(offset)) {
            stopifnot(NROW(fit) == NROW(offset))
            fit <- fit + offset
        }
    } else { # glm predictions
        if(trace >= 1)
            cat("predict.earth: returning glm", type, "predictions\n")
        fit <- matrix(0, nrow=nrow(bx),
                         ncol=ncol(object$fitted.values))
        colnames(fit) <- colnames(object$fitted.values)
        bx <- eval(bx[,-1, drop=FALSE], envir=env) # -1 to drop intercept
        bx.data.frame <- as.data.frame(bx)
        if(!is.null(offset)) {
            stopifnot(NROW(bx.data.frame) == NROW(offset))
            bx.data.frame$offset <- offset
        }
        for(i in seq_along(object$glm.list)) {
            fit[,i] = predict.glm(object$glm.list[[i]],
                                  newdata=bx.data.frame, type=type)
            check.nrows(nrow(bx), nrow(fit),
                        nrow(object$fitted.values), "predict.earth")
        }
    }
    fit
}
print_returning_earth <- function(object, trace, msg)
{
    if(trace >= 1) {
        if(is.null(object$glm.list))
            cat("predict.earth: returning earth", msg, "\n")
        else
            cat("predict.earth: returning earth (not glm)", msg, "\n")
    }
}
convert.predicted.response.to.class <- function(resp, ylevels, resp.name, thresh=.5)
{
    which1 <- function(row, thresh) # row is a scalar or a vector
    {
        if(length(row) > 1)
            which. <- which.max(row)
        else
            which. <- if(row > thresh) 2 else 1
        which.
    }
    if(is.null(ylevels)) # should never happen
        stop0("predict.earth: cannot use type=\"class\" with this model")
    check.numeric.scalar(thresh)
    resp <- ylevels[apply(resp, 1, which1, thresh)]
    if(is.character(ylevels))
        resp <- factor(resp, levels = ylevels)
    fit <- as.matrix(resp, ncol=1)
    colnames(fit) <- resp.name
    fit
}
# type="terms" was passed to predict.earth, return just enough for termplot to work
predict_earth_terms <- function(object, newdata, env, trace)
{
    if(!is.null(object$glm.list))
        warning0("predict.earth: returning the earth (not glm) terms")
    bx <- model.matrix.earth(object, x=newdata, trace=trace,
                             Env=env, Callers.name="predict.earth")
    dirs <- object$dirs[object$selected.terms, , drop=FALSE]
    # retain only additive terms
    additive.terms <- get.degrees.per.term(dirs) == 1
    bx <- bx[, additive.terms, drop=FALSE]
    dirs <- dirs[additive.terms, , drop=FALSE]
    coefs <- object$coefficients[additive.terms, 1, drop=FALSE]
    additive.preds <- colSums(abs(dirs)) != 0
    dirs <- dirs[, additive.preds, drop=FALSE]
    var.names <- variable.names(object, use.names=TRUE)[additive.preds]
    termMat <- matrix(0, nrow=nrow(bx), ncol=ncol(dirs))
    colnames(termMat) <- var.names
    if(ncol(bx) >= 1)
        for(ipred in seq_len(ncol(dirs)))
            for(iterm in seq_len(ncol(bx)))
                if(dirs[iterm, ipred])
                    termMat[, ipred] =
                        termMat[, ipred] + coefs[iterm] * bx[, iterm]
    termMat
}

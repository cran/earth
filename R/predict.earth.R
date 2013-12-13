# predict.earth.R

# predict.earth returns multiple columns for multiple response models

predict.earth <- function(
    object  = stop("no 'object' arg"),
    newdata = NULL,
    type    = c("link", "response", "earth", "class", "terms"),
                        # "terms" always returns the earth not glm terms
                        # "terms" returns just the additive terms!
                        # and just the first response if more than one
    thresh  = .5,       # only used if type="class"
    trace   = FALSE,
    ...)                # unused, for compatibility with generic predict
{
    print.is.earth <- function(trace, object, msg)
    {
        if(trace >= 1)
            if(is.null(object$glm.list))
                cat("predict.earth: returning earth", msg, "\n")
            else
                cat("predict.earth: returning earth (not glm)", msg, "\n")
    }
    get.predicted.response <- function(object, newdata, type)
    {
        is.type.class <- FALSE
        ylevels <- object$levels
        if(type=="class") {
            is.type.class <- TRUE
            type <- "response" # we want predicted probabilities
            if(is.null(ylevels))
                ylevels <- c(FALSE, TRUE)
        }
        if(is.null(newdata)) {
            if(is.null(object$glm.list) || type=="earth") {
                print.is.earth(trace, object, "fitted.values")
                rval <- object$fitted.values
            } else {    # glm predictions
                if(trace >= 1)
                    cat("predict.earth: returning glm fitted.values\n")
                rval <- matrix(0, nrow=nrow(object$fitted.values),
                                  ncol=ncol(object$fitted.values))
                colnames(rval) <- colnames(object$fitted.values)
                for(i in 1:length(object$glm.list))
                    rval[,i] = predict.glm(object$glm.list[[i]], type=type)
            }
        } else {        # user supplied newdata
            bx <- model.matrix.earth(object, newdata, env=env,
                                     trace=trace,
                                     Callers.name="model.matrix.earth from predict.earth")
            if(trace >= 1)
                 print.matrix.info("bx", bx, "predict.earth", all.names=trace>=2)
            if(is.null(object$glm.list) || type=="earth") {
                print.is.earth(trace, object, "predictions")
                rval <- bx %*% object$coefficients
            } else {    # glm predictions
                if(trace >= 1)
                    cat("predict.earth: returning glm", type, "predictions\n")
                rval <- matrix(0, nrow=nrow(bx),
                                  ncol=ncol(object$fitted.values))
                colnames(rval) <- colnames(object$fitted.values)
                bx <- eval(bx[,-1, drop=FALSE], env) # -1 to drop intercept
                bx.data.frame <- as.data.frame(bx)
                for(i in 1:length(object$glm.list)) {
                    rval[,i] = predict.glm(object$glm.list[[i]],
                                           newdata=bx.data.frame, type=type)
                    check.nrows(nrow(bx), nrow(rval),
                                nrow(object$fitted.values), "predict.earth")
                }
            }
        }
        if(is.type.class)
            rval <- convert.predicted.response.to.class(rval, ylevels, thresh)
        rval
    }
    # returns just enough for termplot to work
    get.terms <- function(object, newdata)
    {
        if(!is.null(object$glm.list))
            warning0("predict.earth: returning the earth (not glm) terms")
        bx <- model.matrix.earth(object, x=newdata, env=env,
                                 trace=trace, Callers.name="predict.earth")
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
            for(ipred in 1:ncol(dirs))
                for(iterm in 1:ncol(bx))
                    if(dirs[iterm, ipred])
                        termMat[, ipred] =
                            termMat[, ipred] + coefs[iterm] * bx[, iterm]
        termMat
    }
    #--- predict.earth starts here ---

    check.classname(object, deparse(substitute(object)), "earth")
    warn.if.dots.used("predict.earth", ...)
    trace <- check.trace.arg(trace)
    env <- parent.frame() # the environment from which predict.earth was called
    switch(match.arg1(type),
        get.predicted.response(object, newdata, "link"),
        get.predicted.response(object, newdata, "response"),
        get.predicted.response(object, newdata, "earth"),
        get.predicted.response(object, newdata, "class"),
        get.terms(object, newdata))           # "terms"
}
# resp is a matrix, return a vector
convert.predicted.response.to.class <- function(resp, ylevels, thresh=.5)
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
        stop0("cannot use type=\"class\" with this model")
    resp <- ylevels[apply(resp, 1, which1, thresh)]
    if(is.character(ylevels))
        resp <- factor(resp, levels = ylevels)
    resp
}

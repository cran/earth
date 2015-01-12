# residuals.earth.R:

residuals.earth <- function(object = stop("no 'object' arg"), type = NULL, warn=TRUE, ...)
{
    glm.resids <- function(object, type)
    {
        g <- object$glm.list
        if(is.null(g))
            stop0("residuals.earth: type \"", type, "\" can be used ",
                  "only on earth-glm models")
        colnames <- ""
        for(imodel in seq_along(g)) {
            rval1 <- residuals(g[[imodel]], type)
            if(imodel == 1)
                rval <- rval1
            if(NROW(rval1) != NROW(rval)) # should never happen
                stop0("residuals.earth: glm.list[[", imodel, "]] does ",
                      "not conform to glm.list[[", 1, "]] ",
                      "(residuals have a different length)")
            if(imodel > 1) {
                colnames <- c(colnames)
                rval <- cbind(rval, rval1)
            }
        }
        rval
    }
    #--- residuals.earth starts here ---
    warn.if.dots.used("residuals.earth", ...)
    if(warn && is.null(type) && !is.null(object$glm.list))
        warning0("residuals.earth: returning earth (not glm) residuals")
    if(is.null(type))
        type <- "earth"
    types <- c("earth", "student", "delever", "deviance",
               "glm.pearson", "glm.working", "glm.response", "glm.partial")
    resids <- switch(match.choices(type, types),
        earth        = object$residuals,
        student      = { rinfo <- get.rinfo(object, 1, student=TRUE, delever=FALSE,
                                            "will be returned as NA")
                         rinfo$scale * rinfo$resids
                       },
        delever      = { rinfo <- get.rinfo(object, 1, student=FALSE, delever=TRUE,
                                            "will be returned as NA")
                         rinfo$scale * rinfo$resids
                       },
        deviance     = if(is.null(object$glm.list))
                           object$residuals
                       else
                           glm.resids(object, "deviance"),
        glm.pearson  = glm.resids(object, "pearson"),
        glm.working  = glm.resids(object, "working"),
        glm.response = glm.resids(object, "response"),
        glm.partial  = glm.resids(object, "partial"))

    if(!is.matrix(resids))
        resids <- matrix(resids, ncol = 1)
    if(type != "glm.partial")
        colnames(resids) <- colnames(object$residuals)
    rownames(resids) <- case.names(object)
    resids
}
resid.earth <- function(object = stop("no 'object' arg"), type = NULL, warn=TRUE, ...)
{
    residuals.earth(object, type, warn, ...)
}
get.student.scale <- function(object) # scale for studentization, inf if leverage is 1
{
    if(inherits(object, "earth")) {
        if(is.null(object$varmod))
            stop0("\"student\" is not allowed because\n",
                  "the model was not built with varmod.method")
        se <- predict(object, type="earth", interval="se")
    } else if(inherits(object, "lm")) {
        predict <- predict(object, se.versus=TRUE)
        se <- predict$residual.scale # TODO correct?
    } else
       stop0("\"student\" is not allowed because\n",
             "get.se doesn't know how to get the stderrs of a \"",
             class(object)[1], "\" object")

    1 / (se * sqrt(1 - hatvalues(object)))
}
possible.leverage.warning <- function(scale, object, leverage.msg)
{
    if(any(is.na(scale))) {
        which <- which(hatvalues(object) == 1)
        # following check is paranoia in case NA was caused by something else
        bad.hat <- which(hatvalues(object) == 1)
        if(length(bad.hat) > 0)
            warnf("response[%s] has a leverage of one and will be %s",
                  paste.with.c(bad.hat), leverage.msg)
        else # probably can never get here, but it doesn't matter if we do
            warnf("response[%s] has a NA scale and will be %s",
                  paste.with.c(which(is.na(scale))), leverage.msg)
    }
}
# get information on the residuals (the residuals, and their scale and name)

get.rinfo <- function(object, nresponse, student, delever, leverage.msg="ignored")
{
    resids <- residuals(object, warn=FALSE)
    if(is.null(resids))
        stop0("cannot get residuals for \"", class(object), "\" object")
    if(NCOL(resids) == 1)
        resids <- matrix(resids, ncol=1)
    nresponse <- plotmo::check.index(nresponse, "nresponse", resids, is.col.index=TRUE)
    resids <- resids[, nresponse, drop=FALSE]
    scale <- repl(1, length(resids))
    name <- "Residual"
    student <- check.boolean(student)
    if(student) {
        if(is.null(object$varmod))
            stop0("no studentized residuals because ",
                  "the earth model was not built with varmod.method")
        scale <- get.student.scale(object)
        name <- "Studentized Residual"
    }
    delever <- check.boolean(delever)
    if(delever) {
        if(student) # don't allow double denormalization
            stop0("the student and delever arguments cannot both be set")
        scale <- 1 / sqrt(1 - hatvalues(object))
        name <- "Delevered Residual"
    }
    # leverages of 1 cause an infinite scale, change that to NA for easier handling later
    scale[is.infinite(scale)] <- NA
    check.vec(scale, "scale", length(resids), allow.na=TRUE)
    check(scale, "scale", "non-positive value", function(x) { x <= 0 }, allow.na=TRUE)
    possible.leverage.warning(scale, object, leverage.msg)
    list(resids = resids, # raw resids (student and delever not applied)
         scale  = scale,  # will be 1 unless student or delever set
         name   = name)   # "Residual" or "Delevered Residual" etc.
}

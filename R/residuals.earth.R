# residuals.earth.R:

residuals.earth <- function(object=stop("no 'object' argument"), type=NULL, warn=TRUE, ...)
{
    warn.if.dots(...)
    warn <- check.boolean(warn)
    if(warn && is.null(type) && !is.null(object$glm.list))
        warning0("residuals.earth: returning earth (not glm) residuals")
    if(is.null(type))
        type <- "earth"
    types <- c("earth", "deviance", "response",
               "standardize", "delever",
               "pearson", "working", "partial",
               "glm.response", "glm.pearson", "glm.working", "glm.partial")
    if(is.null(object$residuals)) # I think this can only happen for cv models
        stop0("earth object has no residuals field.\n",
              "       Use keepxy=TRUE in the call to earth.")

    resids <- switch(match.choices(type, types, "type"),
        earth        = object$residuals,
        deviance     = if(is.null(object$glm.list))
                           object$residuals
                       else
                           glm.resids(object$glm.list, "deviance"),
        response     = if(is.null(object$glm.list))
                            object$residuals
                       else
                            glm.resids(object$glm.list, "response"),

        standardize  = plotmo::plotmo_standardizescale(object) * object$residuals,
        delever      = object$residuals / sqrt(1 - hatvalues(object)),

        pearson      = glm.resids(object$glm.list, "pearson"),
        working      = glm.resids(object$glm.list, "working"),
        partial      = glm.resids(object$glm.list, "partial"),

        glm.response = glm.resids(object$glm.list, "response"),
        glm.pearson  = glm.resids(object$glm.list, "pearson"),
        glm.working  = glm.resids(object$glm.list, "working"),
        glm.partial  = glm.resids(object$glm.list, "partial"))

    if(!is.matrix(resids))
        resids <- matrix(resids, ncol = 1)
    if(type != "partial" && type != "glm.partial")
        colnames(resids) <- colnames(object$residuals)
    rownames(resids) <- case.names(object)
    resids
}
glm.resids <- function(glm.list, type)
{
    if(is.null(glm.list))
        stop0("residuals.earth: type \"", type, "\" can be used ",
              "only on earth-glm models")
    colnames <- ""
    for(imodel in seq_along(glm.list)) {
        rval1 <- residuals(glm.list[[imodel]], type) # invokes residuals.glm
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
resid.earth <- function(object=stop("no 'object' argument"), type=NULL, warn=TRUE, ...)
{
    residuals.earth(object, type, warn, ...)
}

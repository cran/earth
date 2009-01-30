# earth.cv.R: functions for cross validation of earth models

# earth.cv returns null unless nfold > 1

earth.cv <- function(x, y, weights, wp, scale.y, subset, na.action,
                     glm, trace, keepxy, nfold, stratify, ...)
{
    if(!is.numeric(nfold))
        stop1("nfold argument must be numeric")
    if(floor(nfold) != nfold)
        stop1("nfold argument must be an integer, you have nfold=", nfold)
    if(nfold < 0 || nfold > 10000) # 10000 is arbitrary
        stop1("nfold argument must be between 0 and 10000, you have nfold=", nfold)
    # if(nfold == 1)
    #   warning1("nfold=1 treated as nfold=0 i.e. no cross validation")
    if(nfold <= 1)
        return (NULL)

    # get here if must do the cross validation

    wp.original <- wp
    if(is.null(wp))
        wp <- rep(1, ncol(y))   # all ones vector
    list. <- list()             # returned list of cross validated models
    ncases <- nrow(x)
    nresp <- ncol(y)            # number of responses
    # ndigits aligns prints without too much white space
    ndigits <- ceiling(log10(ncases - ncases/nfold + .1))
    groups <- sample(rep(1:nfold, length=ncases))
    if(stratify) {
        # Get (roughly) equal number of folds for each non-zero entry in each y column
        # If y was originally a factor before expansion to multiple columns, this is
        # equivalent to having the same numbers of each factor level in each fold.
        for(iresp in 1:ncol(y)) {
            yset = y[,iresp] != 0
            groups[yset] <- sample(rep(1:nfold, length=sum(yset)))
        }
    }
    if(any(table(groups) == 0))
        stop1("Not enough data to do ", nfold, 
              " fold cross validation (a CV test set is empty)")
    is.binomial <- is.poisson <- FALSE
    if(!is.null(glm)) {
        glm1 <- get.glm.arg(glm)
        family <- get.glm.family(glm1$family, env=parent.frame())
        is.binomial <- is.binomial(family)
        is.poisson <- is.poisson(family)
    }
    augmented.resp.names <- c(colnames(y), "mean") # response names plus "mean"
    format.string <- if(nfold < 10) "%d" else "%2d"
    augmented.fold.names <- c(paste("fold ", sprintf(format.string, 1:nfold), sep=""), "mean")
    nvars <- double(nfold+1)    # number of used predictors in each CV model
    nterms <- double(nfold+1)   # number of selected terms in each CV model
    names(nvars) <- names(nterms) <- augmented.fold.names

    rsq.tab <- matrix(0, nrow=nfold+1, ncol=1+nresp) # table of cv results, +1 for means
    colnames(rsq.tab) <- augmented.resp.names
    rownames(rsq.tab) <- augmented.fold.names

    maxerr.tab <- matrix(0, nrow=nfold+1, ncol=1+nresp) # table of cv results, +1 for max
    colnames(maxerr.tab) <- c(colnames(y), "max") # response names plus "max"
    rownames(maxerr.tab) <- c(paste("fold ", sprintf(format.string, 1:nfold), sep=""), "max")

    # the following remain NULL unless is.binomial or is.poisson
    deviance.tab <- calib.int.tab <- calib.slope.tab <- test.tab <- NULL 
    if(is.binomial || is.poisson) {
        deviance.tab    <- matrix(0, nrow=nfold+1, ncol=1+nresp) # mean deviance
        calib.int.tab   <- matrix(0, nrow=nfold+1, ncol=1+nresp)
        calib.slope.tab <- matrix(0, nrow=nfold+1, ncol=1+nresp)
        test.tab        <- matrix(0, nrow=nfold+1, ncol=1+nresp) # binomial auc, poisson cor
        colnames(deviance.tab) <- colnames(calib.int.tab) <- 
            colnames(calib.slope.tab) <- colnames(test.tab) <- augmented.resp.names
        rownames(deviance.tab) <- rownames(calib.int.tab) <- 
            rownames(calib.slope.tab) <- rownames(test.tab) <- augmented.fold.names
    }
    ilast <- ncol(rsq.tab) # index of final "mean" column in tables
    for(ifold in 1:nfold) {
        if(trace && trace < 1)
            printf("CV fold %2d: ", ifold)
        else if(trace >= 1) # newline etc. to distinguish this from other trace prints
            printf("\nCV fold %2d -----------------------------%s", ifold,
                   "---------------------------------------\n")
        test.subset <- (1:ncases)[which(groups == ifold)]
        train.subset <- (1:ncases)[-test.subset]
        fit <- earth.fit(x=x[train.subset,,drop=FALSE],
                      y=y[train.subset,,drop=FALSE],
                      weights=weights, wp=wp.original, scale.y=scale.y,
                      subset=subset, na.action=na.action, glm=glm,
                      trace=trace, ...)

        ytest <- y[test.subset,,drop=FALSE]
        yhat <- predict(fit, newdata=x[test.subset,,drop=FALSE], type="earth")
        if(is.binomial || is.poisson)
            response.yhat <- 
                predict(fit, newdata=x[test.subset,,drop=FALSE], type="response")

        # fill in this fold's row in summary tabs

        for(iresp in 1:nresp) {
            rsq.tab[ifold, iresp] <-
                get.rsq(var(yhat[,iresp] - ytest[,iresp]), var(ytest[,iresp]))
            if(is.binomial) {
                deviance.tab[ifold, iresp] <- 
                    get.binomial.deviance(response.yhat[,iresp], ytest[,iresp])
                calib <- get.binomial.calib(response.yhat[,iresp], ytest[,iresp])
                calib.int.tab[ifold, iresp]   <- calib[1]
                calib.slope.tab[ifold, iresp] <- calib[2]
                maxerr.tab[ifold, iresp] <- 
                    get.maxerr(ytest[,iresp] - response.yhat[,iresp])
                test.tab[ifold, iresp] <-
                    get.auc(response.yhat[,iresp], ytest[,iresp])
            }
            else if(is.poisson) {
                deviance.tab[ifold, iresp] <- 
                    get.poisson.deviance(response.yhat[,iresp], ytest[,iresp])
                calib <- get.poisson.calib(response.yhat[,iresp], ytest[,iresp])
                calib.int.tab[ifold, iresp]   <- calib[1]
                calib.slope.tab[ifold, iresp] <- calib[2]
                maxerr.tab[ifold, iresp] <-
                    get.maxerr(ytest[,iresp] - response.yhat[,iresp])
                test.tab[ifold, iresp] <- 
                    cor(response.yhat[,iresp], ytest[,iresp])
            }
            else {
                maxerr.tab[ifold, iresp] <- 
                    get.maxerr(ytest[,iresp] - yhat[,iresp])
	    	}
        }
        nvars[ifold] <- get.nused.preds.per.subset(fit$dirs, fit$selected.terms)
        nterms[ifold] <- length(fit$selected.terms)

        # init last column of summary tables
        # take weighted mean of RSq etc. across all responses (TODO suspect)

        rsq.tab[ifold, ilast]    <- sum(wp * rsq.tab[ifold, -ilast])    / sum(wp)
        maxerr.tab[ifold, ilast] <- get.maxerr(maxerr.tab[ifold, -ilast])
        if(is.binomial || is.poisson) {
            deviance.tab[ifold, ilast]    <- sum(wp * deviance.tab   [ifold, -ilast]) / sum(wp)
            calib.int.tab[ifold, ilast]   <- sum(wp * calib.int.tab  [ifold, -ilast]) / sum(wp)
            calib.slope.tab[ifold, ilast] <- sum(wp * calib.slope.tab[ifold, -ilast]) / sum(wp)
            test.tab[ifold, ilast]        <- sum(wp * test.tab       [ifold, -ilast]) / sum(wp)
        }
        if(!keepxy) {
            # zap big elements in fit so list of models does not take too much memory
            fit$bx <- NULL
            fit$residuals <- NULL
            fit$prune.terms <- NULL
        }
        list.[[ifold]] <- fit

        if(trace >= .5)
            trace.fold(ifold, trace, y, train.subset, test.subset, ndigits, rsq.tab[ifold,])
    }
    # init last row of summary tables

    ilast <- nfold+1 # index of final "mean" row in tables
    nvars[ilast] <- mean(nvars[-ilast])
    nterms[ilast] <- mean(nterms[-ilast])
    rsq.tab[ilast,]             <- colMeans(rsq.tab        [-ilast,])
    maxerr.tab[ilast,]          <- get.maxerr(maxerr.tab   [-ilast,])
    if(is.binomial || is.poisson) {
        deviance.tab[ilast,]    <- colMeans(deviance.tab   [-ilast,])
        calib.int.tab[ilast,]   <- colMeans(calib.int.tab  [-ilast,])
        calib.slope.tab[ilast,] <- colMeans(calib.slope.tab[-ilast,])
        test.tab[ilast,]        <- colMeans(test.tab       [-ilast,])
    }
    if(trace >= 1)
        cat("\n")
    if(trace >= .5) {
        trace.fold(-1, trace, y, TRUE, TRUE, ndigits, rsq.tab[ilast,])
        cat("\n")
    }
    list(rsq.tab=rsq.tab, 
         maxerr.tab=maxerr.tab, 
         deviance.tab=deviance.tab,
         calib.int.tab=calib.int.tab,
         calib.slope.tab=calib.slope.tab,
         auc.tab=if(is.binomial) test.tab else NULL,
         cor.tab=if(is.poisson) test.tab else NULL,
         nterms=nterms,
         nvars=nvars,
         groups=groups, # groups used for cross validation
         list.=list.) # list of earth objects built during cross validation
}
# print results for the current fold (ifold=-1 means "all")
# only called if trace >= .5

trace.fold <- function(ifold, trace, y, train.subset, test.subset, ndigits, rsq.row)
{
    if(ifold < 0)
        printf("CV all:     ")
    else if(trace >= 1)
        printf("CV fold %2d: ", ifold)
    printf("CV-RSq %-6.3f ", rsq.row[length(rsq.row)])
    nresp <- length(rsq.row) - 1 # -1 for "mean"
    if(nresp > 1) {
        cat("Per response CV-RSq ")
        for(iresp in 1:nresp)
            printf("%-6.3f ", rsq.row[iresp])
    }
    if(nresp > 1)
        cat("\n                          ");
    cat("ntrain-nz ")
    for(iresp in 1:nresp)
        printf("%*.0f ", ndigits, sum(y[train.subset, iresp] != 0))
    if(ifold >= 0) {
        cat("ntest-nz ")
        for(iresp in 1:nresp)
            printf("%*.0f ", ndigits, sum(y[test.subset, iresp] != 0))
    }
    cat("\n")
}

# Following functions lifted from Elith Leathwick code.
# In that code, the AUC calculation was adapted from Ferrier, Pearce and Watson.
# See Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance
# of habitat models developed using logistic regression.

get.binomial.deviance <- function(yhat, y) # yhat is predicted, y is observed
{
    deviance.contribs <- y * log(yhat) + (1 - y) * log(1 - yhat)
    mean(-2 * sum(deviance.contribs)) / length(y) # TODO length(y) ok, should use dof?
}

get.poisson.deviance <- function(yhat, y)
{
    deviance.contribs <- ifelse(y == 0, 0, y * log(y/yhat)) - (y - yhat)
    2 * sum(deviance.contribs) / length(y)
}

get.auc <- function(yhat, y) # area under ROC curve
{
    nx <- length(y[y == 0])
    ny <- length(y[y == 1])
    xy <- c(yhat[y == 0], yhat[y == 1])
    wilc <- nx * ny   +    (nx * (nx + 1)) / 2   -   sum(rank(xy)[1:nx])
    wilc / (nx * ny)
}

get.binomial.calib <- function(yhat, y) # returns c(intercept, slope)
{
    yhat <- yhat + 1e-005 # TODO why? (from Elith Leathwick code)
    yhat[yhat >= 1] <- 0.99999
    fit <- glm(y ~ log(yhat / (1 - yhat)), family = binomial)$coefficients
}

get.poisson.calib <- function(yhat, y) # returns c(intercept, slope)
{
    fit <- glm(y ~ log(yhat), family = poisson)$coefficients
}
get.maxerr <- function(errs) # get signed max absolute err; if matrix then of each col
{
    if(NCOL(errs) == 1)
        rv <- errs[which.max(abs(errs))]
    else {
		rv <- double(NCOL(errs))
		for(icol in 1:NCOL(errs)) {
			col <- errs[,icol]
			rv[icol] <- col[which.max(abs(col))]
		}
   	}
	rv
}
print.cv <- function(x) # called from print.earth for cross validated models
{
    stopif(is.null(x$cv.list))
    nresp <- NCOL(x$coefficients)
    ilast <- nrow(x$cv.rsq.tab) # index of "mean" row in summary tables
    cat("\n")
    printf("Cross validation: nterms %.2f sd %.2f  nvars %.2f sd %.2f\n\n", 
        x$cv.nterms[ilast], sd(x$cv.nterms[-ilast]), 
        x$cv.nvars[ilast], sd(x$cv.nvars[-ilast]))

    # is.binomial is true if x is a cross-validated binomial model
    is.binomial <- is.binomial(x$glm.list[[1]]$family$family)
    is.poisson <- is.poisson(x$glm.list[[1]]$family$family)

    # create a matrix and print that

    tab <- matrix(0, ncol=4 + (is.binomial || is.poisson) * 8, nrow=nresp+1)
    colnames. <- c("CV-RSq", "sd")
    if(is.binomial || is.poisson)
        colnames. <- c(colnames., "MeanDev", "sd", 
                                  "CalibInt", "sd", "CalibSlope", "sd")
    if(is.binomial)
        colnames. <- c(colnames., "AUC", "sd")
    else if(is.poisson)
        colnames. <- c(colnames., "cor", "sd")
    colnames. <- c(colnames., "MaxErr", "sd")
    colnames(tab) <- colnames.
    rownames(tab) <- c(colnames(x$fitted.values), "mean")

    tab[,1] <- x$cv.rsq.tab[ilast,]
    tab[,2] <- apply(x$cv.rsq.tab[-ilast,], 2, sd)
    if(is.binomial || is.poisson) {
        tab[,3] <- x$cv.deviance.tab[ilast,]
        tab[,4] <- apply(x$cv.deviance.tab[-ilast,], 2, sd)
        tab[,5] <- x$cv.calib.int.tab[ilast,]
        tab[,6] <- apply(x$cv.calib.int.tab[-ilast,], 2, sd)
        tab[,7] <- x$cv.calib.slope.tab[ilast,]
        tab[,8] <- apply(x$cv.calib.slope.tab[-ilast,], 2, sd)
    }
    if(is.binomial) {
        tab[,9]  <- x$cv.auc.tab[ilast,]
        tab[,10] <- apply(x$cv.auc.tab[-ilast,], 2, sd)
    } else if(is.poisson) {
        tab[,9]  <- x$cv.cor.tab[ilast,]
        tab[,10] <- apply(x$cv.cor.tab[-ilast,], 2, sd)
    }
    tab[,ncol(tab)-1] <- x$cv.maxerr.tab[ilast,]
    tab[,ncol(tab)]   <- apply(x$cv.maxerr.tab[-ilast,], 2, sd)
    if(nresp == 1)
        print(tab[1,,drop=FALSE]) # don't print "mean" row
    else
        print(tab)
    NULL
}

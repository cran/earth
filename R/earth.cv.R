# earth.cv.R: Functions for cross validation of earth models.
#             Note that earth.cv returns null unless nfold > 1.

earth.cv <- function(obj, x, y,
                     subset, weights, na.action, keepxy, trace, glm,
                     ncross, nfold, stratify,
                     varmod.method=c("none", VARMOD.METHODS),
                     varmod.exponent, varmod.conv, varmod.clamp, varmod.minspan,
                     Scale.y, env, degree=1, ...)
{
    get.fold.rsq.per.subset <- function(foldmod, oof.y, max.nterms, trace, must.print.dots)
    {
        wp.expanded <- wp.expanded / sum(wp.expanded)
        oof.rsq.per.subset <- infold.rsq.per.subset <- repl(0, max.nterms)
        # nrow(foldmod$dirs) is the number of terms in this fold's model before pruning
        for(nterms in 1:min(max.nterms, nrow(foldmod$dirs))) {
            trace.get.fold1(trace, must.print.dots, nterms)
            # penalty=-1 to enforce strict nprune
            # glm=NULL for speed, ok because we don't need the glm submodel
            # TODO 70% of cv time is spent in update.earth
            pruned.foldmod <- update(foldmod, nprune=nterms, penalty=-1, ponly=TRUE,
                                     glm=NULL, trace=max(0, trace-1))
            fit <- predict(pruned.foldmod, newdata=x, type="earth")
            oof.fit    <- fit[oof.subset,  , drop=FALSE]
            infold.fit <- fit[infold.subset, , drop=FALSE]
            for(i in 1:NCOL(fit)) { # for each response
                oof.rsq.per.subset[nterms] <- oof.rsq.per.subset[nterms] +
                    wp.expanded[i] * get.rsq(ss(oof.y[,i] - oof.fit[,i]),
                                             ss(oof.y[,i] - mean(oof.y[,i])))
                infold.rsq.per.subset[nterms] <- infold.rsq.per.subset[nterms] +
                    wp.expanded[i] * get.rsq(ss(infold.y[,i] - infold.fit[,i]),
                                             ss(infold.y[,i] - mean(infold.y[,i])))
            }
            if(nrow(foldmod$dirs) < max.nterms)
                for(nterms in (nrow(foldmod$dirs)+1): max.nterms)
                    oof.rsq.per.subset[nterms] <- infold.rsq.per.subset[nterms] <- NA
        }
        trace.get.fold2(trace, must.print.dots, nterms)
        list(oof.rsq.per.subset=oof.rsq.per.subset,
             infold.rsq.per.subset=infold.rsq.per.subset)
    }
    #--- earth.cv starts here ---
    stratify <- check.boolean(stratify)
    varmod.method <- match.arg1(varmod.method) # check varmod.method is legal and expand
    check.cv.args(ncross, nfold, varmod.method)
    if(ncross < 1 || nfold <= 1)
        return (NULL)
    if(nfold > nrow(x))
        nfold <- nrow(x)
    if(!is.null(obj$glm.bpairs))
        stop0("earth does not yet support cross validation of paired binomial responses")
    max.nterms <- nrow(obj$dirs)
    wp <- wp.expanded <- obj$wp
    if(is.null(wp))
        wp.expanded <- repl(1, ncol(y)) # all ones vector
    list. <- list()                     # returned list of cross validated models
    ncases <- nrow(x)
    nresp <- ncol(y)                    # number of responses
    # ndigits aligns trace prints without too much white space
    ndigits <- ceiling(log10(ncases - ncases/nfold + .1))
    # print pacifier dots if get.fold.rss.per.subset will be slow
    must.print.dots <- trace >= .5 && trace <= 1 && keepxy && nrow(x) * max.nterms > 50e3
    trace.cv.header(obj, nresp, trace, must.print.dots)
    n.oof.digits <- ceiling(log10(1.2 * ncases / nfold)) # 1.2 allows for diff sized subsets
    groups <- matrix(NA, nrow=ncross*ncases, ncol=2)
    colnames(groups) <- c("cross", "fold")
    for(icross in 1:ncross) {
        start <- ((icross-1) * ncases) + 1
        groups[start:(start+ncases-1), 1] <- icross
        groups[start:(start+ncases-1), 2] <- get.groups(y, nfold, stratify)
    }
    is.binomial <- is.poisson <- FALSE
    if(!is.null(glm)) {
        glm1 <- get.glm.arg(glm)
        family <- get.glm.family(glm1$family, env=env)
        is.binomial <- is.binomial(family)
        is.poisson <- is.poisson(family)
    }
    must.get.class.rate <- !is.null(obj$levels)
    # the final summary row of the tables is "mean", "all" or "max", depending on the statistic.
    if(ncross > 1)
        fold.names <- paste0(rep(paste0("fold", 1:ncross, "."), each=nfold),
                             rep(1:nfold, times=ncross))
    else
        fold.names <- sprintf("fold%d", 1:nfold)
    fold.names.plus.mean <- c(fold.names, "mean")
    fold.names.plus.all  <- c(fold.names, "all")
    fold.names.plus.max  <- c(fold.names, "max")
    resp.names.plus.mean <- c(colnames(y), "mean") # response names plus "mean"
    resp.names.plus.max  <- c(colnames(y), "max")
    ncross.fold <- ncross * nfold
    nvars <- double(ncross.fold+1)    # number of used predictors in each CV model
    nterms <- double(ncross.fold+1)   # number of selected terms in each CV model
    names(nvars) <- names(nterms) <- fold.names.plus.mean

    rsq.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # table of cv results, +1 for means
    colnames(rsq.tab) <- resp.names.plus.mean
    rownames(rsq.tab) <- fold.names.plus.mean

    maxerr.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # table of cv results, +1 for max
    colnames(maxerr.tab) <- resp.names.plus.max
    rownames(maxerr.tab) <- fold.names.plus.max

    deviance.tab <- calib.int.tab <- calib.slope.tab <- test.tab <- class.rate.tab <- NULL
    if(is.binomial || is.poisson) {
        deviance.tab    <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # mean deviance
        calib.int.tab   <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        calib.slope.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        test.tab        <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp) # binomial auc, poisson cor
        colnames(deviance.tab) <- colnames(calib.int.tab) <-
            colnames(calib.slope.tab) <- colnames(test.tab) <- resp.names.plus.mean
        rownames(deviance.tab) <- rownames(calib.int.tab) <-
            rownames(calib.slope.tab) <- rownames(test.tab) <- fold.names.plus.mean
    }
    if(must.get.class.rate) {
        class.rate.tab <- matrix(0, nrow=ncross.fold+1, ncol=1+nresp)
        colnames(class.rate.tab) <- resp.names.plus.mean
        rownames(maxerr.tab) <- fold.names.plus.all
    }
    oof.rsq.tab <- infold.rsq.tab <- oof.fit.tab <- NULL
    if(keepxy) {
        oof.rsq.tab <- infold.rsq.tab <- matrix(0, nrow=ncross.fold+1, ncol=max.nterms)
        colnames(oof.rsq.tab) <- colnames(infold.rsq.tab) <- paste0("nterms", 1:max.nterms)
        rownames(oof.rsq.tab) <- rownames(infold.rsq.tab) <- fold.names.plus.all
    }
    if(varmod.method != "none") {
        oof.fit.tab <- matrix(0, nrow=nrow(x), ncol=ncross) # preds on oof data
        colnames(oof.fit.tab) <- paste0("icross", 1:ncross)
    }
    for(icross in 1:ncross) {
        this.group <- get.this.group(icross, ifold, ncases, groups)
        for(ifold in 1:nfold) {
            icross.fold <- ((icross-1) * nfold) + ifold
            oof.subset <- (1:ncases)[which(this.group == ifold)]
            infold.subset <- (1:ncases)[-oof.subset]
            trace.fold.header(trace, ncross, icross, ifold)
            infold.x <- x[infold.subset,,drop=FALSE]
            infold.y <- y[infold.subset,,drop=FALSE]
            infold.weights <- weights[infold.subset]
            foldmod <- earth.default(x=infold.x, y=infold.y, weights=infold.weights,
                                     wp=wp, Scale.y=Scale.y, subset=subset,
                                     glm=glm, trace=trace, degree=degree, ...)
            foldmod$icross <- icross
            foldmod$ifold <- ifold
            oof.x <- x[oof.subset,,drop=FALSE]
            oof.y <- y[oof.subset,,drop=FALSE]
            oof.fit <- predict(foldmod, newdata=oof.x, type="earth")
            # fill in subset of entries in this icross column of oof.fit.tab
            # note that we use only the first response when there are multiple responses
            if(!is.null(oof.fit.tab))
                oof.fit.tab[oof.subset, icross] <- oof.fit[,1]
            oof.fit.resp <- NULL
            if(is.binomial || is.poisson)
                oof.fit.resp <- predict(foldmod, newdata=oof.x, type="response")
            else if(must.get.class.rate) # not glm but has binary response?
                oof.fit.resp <- oof.fit

            # fill in this fold's row in summary tabs

            for(iresp in 1:nresp) {
                rsq.tab[icross.fold, iresp] <-
                    get.rsq(sum((oof.fit[,iresp] - oof.y[,iresp])^2),
                            sum((oof.y[,iresp] - mean(oof.y[,iresp]))^2))
                if(is.binomial) {
                    deviance.tab[icross.fold, iresp] <-
                        get.binomial.deviance(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib <- get.binomial.calib(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib.int.tab[icross.fold, iresp]   <- calib[1]
                    calib.slope.tab[icross.fold, iresp] <- calib[2]
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit.resp[,iresp])
                    test.tab[icross.fold, iresp] <-
                        get.auc(oof.fit.resp[,iresp], oof.y[,iresp])
                }
                else if(is.poisson) {
                    deviance.tab[icross.fold, iresp] <-
                        get.poisson.deviance(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib <- get.poisson.calib(oof.fit.resp[,iresp], oof.y[,iresp])
                    calib.int.tab[icross.fold, iresp]   <- calib[1]
                    calib.slope.tab[icross.fold, iresp] <- calib[2]
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit.resp[,iresp])
                    test.tab[icross.fold, iresp] <-
                        cor(oof.fit.resp[,iresp], oof.y[,iresp])
                }
                else
                    maxerr.tab[icross.fold, iresp] <-
                        get.maxerr(oof.y[,iresp] - oof.fit[,iresp])
            } # end for iresp

            nvars[icross.fold] <- get.nused.preds.per.subset(foldmod$dirs, foldmod$selected.terms)
            nterms[icross.fold] <- length(foldmod$selected.terms)
            if(must.get.class.rate)
                class.rate.tab[icross.fold, ] <- get.class.rate(oof.fit.resp, oof.y, obj$levels)
            if(keepxy) {
                temp <-
                    get.fold.rsq.per.subset(foldmod, oof.y, max.nterms, trace, must.print.dots)
                oof.rsq.tab[icross.fold,]  <- temp$oof.rsq
                infold.rsq.tab[icross.fold,] <- temp$infold.rsq
            }
            # init last column of summary tables
            ilast.col <- nresp+1 # index of final (summary) column in tables
            rsq.tab[icross.fold, ilast.col]    <- weighted.mean(rsq.tab[icross.fold, -ilast.col], wp.expanded)
            maxerr.tab[icross.fold, ilast.col] <- get.maxerr(maxerr.tab[icross.fold, -ilast.col])
            if(is.binomial || is.poisson) {
                deviance.tab[icross.fold, ilast.col]    <- weighted.mean(deviance.tab   [icross.fold, -ilast.col], wp.expanded)
                calib.int.tab[icross.fold, ilast.col]   <- weighted.mean(calib.int.tab  [icross.fold, -ilast.col], wp.expanded)
                calib.slope.tab[icross.fold, ilast.col] <- weighted.mean(calib.slope.tab[icross.fold, -ilast.col], wp.expanded)
                test.tab[icross.fold, ilast.col]        <- weighted.mean(test.tab       [icross.fold, -ilast.col], wp.expanded)
            }
            if(!keepxy) # reduce memory by getting rid of big fields
                foldmod$bx <- foldmod$residuals <- foldmod$prune.terms <- NULL
            trace.fold(icross, ifold, trace, y, infold.subset, oof.subset, ncross, ndigits,
                       rsq.tab[icross.fold,], n.oof.digits, must.print.dots)
            list.[[icross.fold]] <- foldmod
        } # end for ifold
    } # end for icross

    # init last row of summary tables

    ilast <- ncross.fold+1 # index of last row in tables
    nvars [ilast]  <- mean(nvars[-ilast])
    nterms[ilast]  <- mean(nterms[-ilast])
    rsq.tab   [ilast,] <- colMeans(rsq.tab     [-ilast,])
    maxerr.tab[ilast,] <- get.maxerr(maxerr.tab[-ilast,])
    if(is.binomial || is.poisson) {
        deviance.tab   [ilast,] <- colMeans(deviance.tab   [-ilast,])
        calib.int.tab  [ilast,] <- colMeans(calib.int.tab  [-ilast,])
        calib.slope.tab[ilast,] <- colMeans(calib.slope.tab[-ilast,])
        test.tab       [ilast,] <- colMeans(test.tab       [-ilast,])
    }
    if(must.get.class.rate)
        class.rate.tab[ilast,]  <- colMeans(class.rate.tab [-ilast,])
    oof.rsq.per.subset <- NULL
    if(keepxy) {
        # there will be NAs in oof.rsq.tab if max terms in a fold is
        # less than max terms in full model
        oof.rsq.tab[ilast,]  <- colMeans(oof.rsq.tab[-ilast,], na.rm=TRUE)
        infold.rsq.tab[ilast,] <- colMeans(infold.rsq.tab[-ilast,], na.rm=TRUE)
    }
    varmod <- NULL
    if(varmod.method != "none") {
        meanfit   <- apply(oof.fit.tab,  1, mean) # mean of each row of oof.fit.tab
        meanfit   <- matrix(meanfit, ncol=1)      # convert to a n x 1 mat, for consistency
        model.var <- apply(oof.fit.tab,  1, var)  # var  of each row of oof.fit.tab
        model.var <- matrix(model.var, ncol=1)
        varmod <- varmod(obj,
                         varmod.method, varmod.exponent,
                         varmod.conv, varmod.clamp, varmod.minspan,
                         trace, x, y, meanfit, model.var, degree=degree)
    }
    if(trace >= 1)
        cat("\n")
    trace.fold(icross, -1, trace, y, TRUE, TRUE, ncross, ndigits,
               rsq.tab[ilast,], n.oof.digits, must.print.dots)
    if(trace >= .5)
        cat("\n")
    names(list.) <- fold.names
    list(list.           = list.,   # list of earth objects built during cross validation
         nterms          = nterms,
         nvars           = nvars,
         groups          = groups, # groups used for cross validation
         rsq.tab         = rsq.tab,
         oof.rsq.tab     = oof.rsq.tab,
         infold.rsq.tab  = infold.rsq.tab,
         class.rate.tab  = class.rate.tab,
         maxerr.tab      = maxerr.tab,
         auc.tab         = if(is.binomial) test.tab else NULL,
         cor.tab         = if(is.poisson) test.tab else NULL,
         deviance.tab    = deviance.tab,
         calib.int.tab   = calib.int.tab,
         calib.slope.tab = calib.slope.tab,
         varmod          = varmod,      # null unless varmod specified
         oof.fit.tab     = oof.fit.tab) # null unless varmod specified
}
check.cv.args <- function(ncross, nfold, varmod.method)
{
    if(!is.numeric(ncross))
        stop0("ncross must be numeric")
    if(length(ncross) != 1)
        stop0("ncross must be an integer")
    if(floor(ncross) != ncross)
        stop0("ncross must be an integer, you have ncross=", ncross)
    if(ncross < 0)
        stop0("ncross ", ncross, " is less than 0")
    if(ncross > 1000) # 1000 is arbitrary
        stop0("ncross ", ncross, " is too big")

    if(!is.numeric(nfold))
        stop0("nfold must be numeric")
    if(length(nfold) != 1)
        stop0("nfold must be an integer")
    if(floor(nfold) != nfold)
        stop0("nfold must be an integer, you have nfold=", nfold)
    if(nfold < 0)
        stop0("nfold ", nfold, " is less than 0")
    # if(nfold > 10000) # 10000 is arbitrary
    #     stop0("nfold ", nfold, " is too big")

    if(ncross > 1 && nfold < 2)
        stop0("ncross=", ncross, " yet nfold=", nfold)
    if(ncross < 1 && nfold > 1)
        stop0("ncross=", ncross, " yet nfold=", nfold)

    if(varmod.method != "none") {
        if(nfold <= 1)
            stop0("varmod.method=\"", varmod.method, "\" requires nfold greater than 1")
        if(ncross < 3)
            stop0("ncross=", ncross,
                  "\nWhen varmod.method is used, ncross should be larger\n",
                  "(suggest 100, 30 for a quick look, 3 for debugging code)")
    }
}
get.groups <- function(y, nfold, stratify)
{
    groups <- sample(repl(1:nfold, nrow(y)))
    if(stratify) {
        # Get (roughly) equal number of folds for each non-zero entry in each y column
        # If y was originally a factor before expansion to multiple columns, this is
        # equivalent to having the same numbers of each factor level in each fold.
        for(iresp in 1:ncol(y)) {
            yset <- y[,iresp] != 0
            groups[yset] <- sample(repl(1:nfold, sum(yset)))
        }
    }
    if(any(table(groups) == 0))
        stop0("Not enough data to do ", nfold,
              " fold cross validation (an out-of-fold set is empty)")
    groups
}
get.this.group <- function(icross, ifold, ncases, groups)
{
    start <- ((icross-1) * ncases) + 1
    groups[start:(start+ncases-1), 2]
}
trace.cv.header <- function(obj, nresp, trace, must.print.dots)
{
    if(trace == .5) {
        printf("Full model grsq %5.3f  rsq %5.3f\n", obj$grsq , obj$rsq)
        if(must.print.dots || nresp > 1)
            cat("\n")
    }
}
trace.fold.header <- function(trace, ncross, icross, ifold)
{
    if(trace >= .5 && trace < 1) {
        if(ncross > 1)
            printf("CV fold %2d.%-2d ", icross, ifold)
        else
            printf("CV fold %-2d ", ifold)
    } else if(trace >= 1) { # newline etc. to distinguish this from other trace prints
        if(ncross > 1)
            printf("\nCV fold %d.%d -----------------------------%s", icross, ifold,
                   "---------------------------------------\n")
        else
            printf("\nCV fold %d -----------------------------%s", ifold,
                   "---------------------------------------\n")
    }
}
# print results for the current fold (ifold=-1 means "all")

trace.fold <- function(icross, ifold, trace, y, infold.subset, oof.subset, ncross, ndigits,
                       rsq.row, n.oof.digits, must.print.dots)
{
    if(trace < .5)
        return()
    if(ifold < 0) {
        icross <- if(ncross > 1) "   " else ""
        printf("%s%s", "CV all     ", icross)
    } else if(trace >= 1) {
        if(ncross > 1)
            icross <- sprintf("%2d.", icross)
        else
            icross <- ""
        printf("CV fold %s%-2d ", icross, ifold)
    }
    printf("CVRSq %-6.3f ", rsq.row[length(rsq.row)])
    nresp <- length(rsq.row) - 1 # -1 for "all"
    if(nresp > 1) {
        cat("Per response CVRSq ")
        for(iresp in 1:nresp)
            printf("%-6.3f ", rsq.row[iresp])
    }
    if(nresp > 1) {
        if(ncross > 1)
            cat("\n                            ")
        else
            cat("\n                         ")
    }
    if(ifold < 0)
        printf("      %.*s       ", n.oof.digits, "                    ")
    else
        printf("n.oof %*.0f %2.0f%%  ", n.oof.digits,
               length(infold.subset),
               100 * (length(y[, 1]) - length(infold.subset)) / length(y[, 1]))
    cat("n.infold.nz ")
    if(nresp == 1)
        printf("%*.0f %2.0f%%  ", ndigits, sum(y[infold.subset, 1] != 0),
               100 * sum(y[infold.subset, 1] != 0) / length(y[infold.subset, 1]))
    else for(iresp in 1:nresp)
        printf("%*.0f ", ndigits, sum(y[infold.subset, iresp] != 0))

    if(ifold >= 0) {
        cat("n.oof.nz ")
        if(nresp == 1)
            printf("%*.0f %2.0f%%  ", ndigits, sum(y[oof.subset, 1] != 0),
                   100 * sum(y[oof.subset, 1] != 0) / length(y[oof.subset, 1]))
        else for(iresp in 1:nresp)
            printf("%*.0f  ", ndigits, sum(y[oof.subset, iresp] != 0))
        if(nresp > 1)
            cat("\n")
    }
    if(must.print.dots && trace == .5 && ifold > 0)
        cat("\n")
    cat("\n")
}
trace.get.fold1 <- function(trace, must.print.dots, nterms)
{
    if(trace >= 2)
        cat0("\nget.fold.rss.per.subset nterms=", nterms, "\n")
    else if(must.print.dots) {
        cat0(".")
        if(nterms %% 40 == 0) {
            cat0("\n")
            if(trace == .5)
                cat("            ")
        }
        flush.console()
    }
}
trace.get.fold2 <- function(trace, must.print.dots, nterms)
{
    if(must.print.dots && nterms %% 40) { # nterms %% 40 avoids double newline
        cat0("\n")
        if(trace == .5)
            cat("            ")
    }
}
print.cv <- function(x) # called from print.earth for cross validated models
{
    get.field <- function(field) round(x[[field]][ilast,], 2)
    get.sd    <- function(field) round(apply(x[[field]][-ilast,], 2, sd), 3)
    #--- print.cv starts here ---
    stopif(is.null(x$cv.list))
    ilast <- nrow(x$cv.rsq.tab) # index of "all" row in summary tables
    cat("\nNote: the cross-validation sd's below are standard deviations across folds\n\n")

    printf("Cross validation:   nterms %.2f sd %.2f    nvars %.2f sd %.2f\n\n",
        x$cv.nterms[ilast], sd(x$cv.nterms[-ilast]),
        x$cv.nvars[ilast], sd(x$cv.nvars[-ilast]))

    wide.spacing <- is.null(x$cv.deviance.tab) # if printing little then use wide spacing

    # create a data.frame and print that
    tab <- if(wide.spacing) data.frame("    CVRSq"=get.field("cv.rsq.tab"), check.names = FALSE)
           else             data.frame("CVRSq"  =get.field("cv.rsq.tab"))

    tab$sd.1=get.sd("cv.rsq.tab")

    if(!is.null(x$cv.class.rate.tab)) {
        if(wide.spacing) tab$"    ClassRate" <- get.field("cv.class.rate.tab")
        else             tab$"ClassRate"   <- get.field("cv.class.rate.tab")

        tab$sd.2 <- get.sd("cv.class.rate.tab")
    }
    if(wide.spacing)
        tab$"    MaxErr" <- x$cv.maxerr.tab[ilast,]
    else
        tab$"MaxErr"   <- x$cv.maxerr.tab[ilast,]

    tab$sd <- apply(x$cv.maxerr.tab[-ilast,], 2, sd)

    if(!is.null(x$cv.auc.tab)) {
        tab$"AUC" <- get.field("cv.auc.tab")
        tab$sd.3  <- get.sd("cv.auc.tab")
    }
    if(!is.null(x$cv.cor.tab)) {
        tab$"cv.cor.tab" <- get.field("cv.cor.tab")
        tab$sd.4         <- get.sd("cv.cor.tab")
    }
    if(!is.null(x$cv.deviance.tab)) {
        tab$"MeanDev"    <- x$cv.deviance.tab[ilast,]
        tab$sd.5         <- apply(x$cv.deviance.tab[-ilast,], 2, sd)
        tab$"CalibInt"   <- get.field("cv.calib.int.tab")
        tab$sd.6         <- get.sd("cv.calib.int.tab")
        tab$"CalibSlope" <- get.field("cv.calib.slope.tab")
        tab$sd.7         <- get.sd("cv.calib.slope.tab")
    }
    # change "sd.N" to plain "sd"
    names <- names(tab)
    names[grep("sd", names)] <- "sd"
    names(tab) <- names

    rownames(tab) <- c(colnames(x$fitted.values), "All")

    digits <- min(getOption("digits"), 2)

    if(NCOL(x$coefficients) == 1)   # single response model?
        print(tab[1,,drop=FALSE], digits=digits, row.names=FALSE) # don't print "All" row
    else                            # multiple response model
        print(tab, digits=digits)
}

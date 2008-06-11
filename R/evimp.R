# evimp.R: estimate variable importances in an earth object


# Return a vector of predictor indexes for predictors that are used
# in the final model

get.used.preds <- function(obj)   # obj is an earth object
{
    which(apply(obj$dirs[obj$selected.terms,,drop=FALSE],2,any1))
}

# Print predictors in order of decreasing estimated importance.
# A one line summary --- print up to 10 predictors.
# Called by print.summary.earth.

print.estimated.predictor.importance <- function(obj) # obj is an "earth" obj
{
    evimp <- row.names(evimp(obj, trim=FALSE))

    cat("Importance: ")
    nprint <- min(10, length(evimp))
    if(nprint == 0)
        cat("no predictors")
    else {
        cat(paste.with.comma(evimp[1:nprint]))
        if(nprint < length(evimp))
            cat(", ...")
    }
    cat("\n")
}

evimp <- function(obj, trim=TRUE) # see help page for description
{
    # convert col numbers in predtab to col numbers in importances
    as.icriti <- function(icrit) c(3,4,6)[icrit]

    check.classname(obj, deparse(substitute(obj)), "earth")
    nsubsets <- length(obj$selected.terms)
    dirs <- obj$dirs
    pred.names <- colnames(dirs)

    # tagged.pred.names is a copy of pred.names but with unused
    # predictors renamed by adding a "-unused" suffix.
    # By unused, we mean unused in the final model.

    used.preds <- to.logical(get.used.preds(obj), len=length(pred.names))
    tagged.pred.names <- pred.names
    tagged.pred.names[!used.preds] <-
            paste(tagged.pred.names[!used.preds], "-unused", sep="")

    # deltas[isubset, icrit] is the change in criterion value
    # for isubset using criterion icrit

    stopifnot(nsubsets >= 1)
    deltas <- matrix(nrow=nsubsets-1, ncol=3)
    colnames(deltas) <- c("nsubsets", "gcv", "rss")
    deltas[,"nsubsets"] <- rep(1, times=nsubsets-1)
    deltas[,"gcv"]      <- -diff(obj$gcv.per.subset[1:nsubsets])
    deltas[,"rss"]      <- -diff(obj$rss.per.subset[1:nsubsets])

    # preds.in.each.term[iterm] is the indices of predictors in term iterm

    preds.in.each.term <- apply(obj$dirs, 1, function(row) which(row != 0))

    # importances is the matrix we return

    importances <- matrix(0, nrow=length(pred.names), ncol=7)
    colnames(importances) <- c("col", "used", "nsubsets", "gcv", "", "rss", "")
    rownames(importances) <- tagged.pred.names
    importances[, "col"] <- 1:nrow(importances)
    importances[used.preds, "used"] <- 1

    if(nsubsets > 1) for(isubset in 2:nsubsets) {
        terms.in.this.subset <- obj$prune.terms[isubset,-1]  # -1 drops intercept
        preds.in.this.subset <-
            unique(unlist(preds.in.each.term[terms.in.this.subset]))

        for(icrit in 1:3) {
            icriti <- as.icriti(icrit)
            importances[preds.in.this.subset, icriti] <-
                importances[preds.in.this.subset, icriti] +
                deltas[isubset-1, icrit]
        }
    }
    # sort rows in "importances" by the nsubsets criteria
    # and with the "gcv" criterion as a secondary sort key

    order.nsubsets <- order(importances[,"nsubsets"], importances[,"gcv"], decreasing=TRUE)
    importances <- importances[order.nsubsets, , drop=FALSE]

    if(nrow(importances) > 1)
        for(icrit in 2:3) {
            # tag importances where gcv or rss ordering disagrees with nsubsets ordering

            icriti <- as.icriti(icrit)
            importances[, icriti+1] <- 1
            for(i in 2:nrow(importances))
                if(importances[i,icriti] > importances[i-1,icriti])
                    importances[i, icriti+1] <- 0

            # normalize importances

            max <- max(abs(importances[,icriti]))
            if(max != 0)
                importances[,icriti] <- 100 * importances[,icriti] / max
        }

    if(trim) {
        # keep only rows for predictors that are used in at least one subset

        in.at.least.one.subset <- importances[,"nsubsets"] != 0
        importances <- importances[in.at.least.one.subset, , drop=FALSE]
    }
    importances
}

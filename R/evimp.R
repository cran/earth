# evimp.R: estimate variable importances in an earth object

# Return a vector of column numbers for predictors that are used
# in the final model

get.used.preds <- function(obj)   # obj is an earth object
{
    which(apply(obj$dirs[obj$selected.terms,,drop=FALSE],2,any1))
}

# Print predictors in order of decreasing estimated importance.
# A one line summary --- print up to 10 predictors.
# Called by print.summary.earth.

print.one.line.evimp <- function(obj) # obj is an "earth" obj
{
    if(is.null(obj$prune.terms)) # this occurs on mars.to.earth objects
        return()
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

evimp <- function(obj, trim=TRUE, sqrt.=FALSE) # see help page for description
{
    # convert col numbers in predtab to col numbers in importances
    as.icriti <- function(icrit) c(3,4,6)[icrit]

    check.classname(obj, deparse(substitute(obj)), "earth")
    if(is.null(obj$prune.terms)) # this occurs on mars.to.earth objects
        stop1("object has no prune.terms, cannot get variable importances")
    nsubsets <- length(obj$selected.terms)
    dirs <- obj$dirs
    pred.names <- generate.colnames(dirs)

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
            if(max != 0) {
                if(sqrt.) {
                    temp <- sqrt(abs(importances[,icriti]) / max)
                    signs <- ifelse(importances[,icriti] < 0, -1, 1)
                    importances[,icriti] <- 100 * signs * temp
                } else
                    importances[,icriti] <- 100 * importances[,icriti] / max
            }
        }

    if(trim) {
        # keep only rows for predictors that are used in at least one subset

        in.at.least.one.subset <- importances[,"nsubsets"] != 0
        importances <- importances[in.at.least.one.subset, , drop=FALSE]
    }
    class(importances) <- "evimp"   # allows use of plot.evimp
    importances
}

# utility to print an evimp object without printing the class

print.evimp <- function(x = stop("no 'x' arg"), ...) # x is an "evimp" obj
{
    class(x) <- NULL
    print(x, ...)
}

# TODO this would be better if rotated clockwise 90 degrees so could easily read var names

plot.evimp <- function(
    x = stop("no 'x' arg"),     # an evimp object (called x for consistency with generic)
    cex.var = 1,                # cex for variable names, make smaller if have lots of vars

    type.nsubsets = "l",        # plot type for nsubsets graph, "b" is quite nice too, "n" for none
    col.nsubsets = "black",     # color of nsubsets line
    lty.nsubsets = 1,           # line type of nsubsets line

    type.gcv = "l",             # plot type for gcv graph,
    col.gcv = "lightblue",      # as above but for the gcv plot
    lty.gcv = 1,

    type.rss = "l",             # as above but for the rss plot
    col.rss = "gray60",
    lty.rss = 1,

    cex.legend = 1,             # cex for legend strings, use if want the legend to be smaller
    x.legend = nrow(x),         # x position of legend, use 0 for no legend
    y.legend = x[1,"nsubsets"], # y position of legend

    main = "Variable importance", # main title
    do.par = TRUE,              # call par() as appropriate
    ...)                        # extra args passed to plotting and method funcs
{
    # make sure that all evimp columns are present (extra columns are ok)
    if(any(pmatch(c("col", "used", "nsubsets", "gcv"), colnames(x), nomatch=0) == 0))
        stop("x is not an evimp matrix")

    max.subsets <- x[1, "nsubsets"]
    varlabs <- paste(rownames(x), sprintf("%3d", x[,"col"]))
    par <- par("mar", "cex")
    on.exit(par(par))
    cex.var <- par$cex * cex.var    # cex.var is relative to current cex
    if(do.par) {
        # TODO what is the best way of doing the bottom.margin calculation?
        # The .5 is a hack to convert nchars to line heights, as required by mar
        mar <- par$mar
        mar[1] <- cex.var * .5 * max(nchar(varlabs) + 6)    # bottom margin
        mar[4] <- mar[4] + 3                                # right margin
        par(mar=mar) # big bottom and right margins
    }
    plot(x[, "nsubsets"], ylim=c(0, max.subsets), type=type.nsubsets,
         xlab="", xaxt="n", ylab="nsubsets",
         main=main, lty=lty.nsubsets, col=col.nsubsets)
    lines(max.subsets * x[,"rss"] / 100, type=type.rss, lty=lty.rss, col=col.rss)
    # plot gcv second so it goes on top of rss (gcv arguably more important than rss)
    lines(max.subsets * x[,"gcv"] / 100, type=type.gcv, lty=lty.gcv, col=col.gcv)
    if(!is.null(x.legend) && !is.na.or.zero(x.legend))
        legend(x=x.legend, y = y.legend, xjust=1,   # top right corner by default
               legend=c("nsubsets", "gcv", "rss"),
               col=c(col.nsubsets, col.gcv, col.rss),
               lty=c(lty.nsubsets, lty.gcv, lty.rss),
               bg="white", cex=cex.legend)
    # right hand axis: normalized rss/gcv values, always 0...100
    # TODO how to get the x position in the call to text correct for all window sizes?
    axis(side=4,
         at=c(0,.2*max.subsets,.4*max.subsets,.6*max.subsets,.8*max.subsets,max.subsets),
         labels=c(0,20,40,60,80,100))
    text(x=nrow(x) + 1.8, y=max.subsets/2, "normalized gcv or rss",
         xpd=NA, # no clip to plot region
         srt=90) # rotate text
    # bottom axis: variable names
    # axis() ignores the cex parameter (a bug?), so set cex globally, on.exit will restore it
    par(cex=cex.var)
    axis(side=1, at=seq(1, nrow(x), by=1), labels=varlabs, las=3)
    invisible()
}

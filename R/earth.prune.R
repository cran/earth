# earth.prune.R:
#
# Functions are in alphabetical order.


# This exposes a problem, possibly in the leaps fortran code.
# It reports when more than one term changes in a single pruning pass step.
# To see run the call to earth below, and in the pruning pass note that at
# step 8 several terms are added and deleted (but really only one term
# should be added per step if pmethod="backward" or "forward").
#     earth(O3 ~ ., data = ozone1, degree=2, trace=5)

check.one.term.per.step <- function(prune.terms, trace)
{
    for(i in 2:nrow(prune.terms)) {
        which1 <- which2 <- repl(FALSE, ncol(prune.terms))
        which1[prune.terms[i-1,]] <- TRUE
        which2[prune.terms[i,]]   <- TRUE
        xor <- xor(which1, which2) # true for every term that changed
        if(sum(xor) > 1) { # more than one term changed?
            printf("%g terms changed between steps %g and %g of the pruning pass\n",
                   sum(xor), i-1, i)
            break
        }
    }
}
# Convert lopt format to prune.terms format

convert.lopt <- function(lopt, nprune)
{
    # Assignment fills matrix column wise. We want row wise, so
    # take upper triangle and then transpose.
    prune.terms <- matrix(0, nrow=nprune, ncol=nprune)
    prune.terms[upper.tri(prune.terms, diag=TRUE)] <- lopt
    t(prune.terms)
}
# This returns the RSS and selected terms for each subset of size 1:nprune

eval.model.subsets <- function(
    bx,      # weighted basis matrix
    y,       # weighted model response
    pmethod,
    nprune,  # max nbr of terms (including intercept) in prune subset, in range 1..nterms
    Force.xtx.prune, # TRUE to always call EvalSubsetsUsingXtx rather than leaps
    trace)
{
    stopifnot(nprune >= 1, nprune <= nrow(bx))

    if(ncol(y) > 1)  {           # leaps cannot handle multiple responses
        if(pmethod != "none" && pmethod != "backward")
            stop0("pmethod=\"", pmethod, "\" is not allowed with multiple response models\n",
                  "       (y has ", ncol(y), " columns, use trace=4 to see y)")
        trace2(trace, "Using EvalSubsetsUsingXtx (rather than leaps) because this is a multiple response model\n")
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune, trace)

    } else if(ncol(bx) <= 2) {   # leaps code gives an error for small number of cols
        if(pmethod != "none" && pmethod != "backward")
            pmethod <- "backward"
        trace2(trace, "Using EvalSubsetsUsingXtx (rather than leaps) because ncol(bx) <= 2\n")
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune, trace)

    } else if(Force.xtx.prune) { # user explicitly asked for xtx subset evaluation
        trace2(trace, "Using EvalSubsetsUsingXtx (rather than leaps) because Force.xtx.prune=TRUE\n")
        eval.subsets.xtx(bx, y, pmethod, nprune, Force.xtx.prune, trace)

    } else
        eval.model.subsets.with.leaps(bx, y, pmethod, nprune)
}
eval.model.subsets.with.leaps <- function(bx, y, pmethod, nprune)
{
    rprune <- leaps.setup(x=bx, y=y,
        force.in=1,        # make sure intercept is in model
        force.out=NULL,
        intercept=FALSE,   # we have an intercept so leaps.setup must not add one
        nvmax=nprune, nbest=1, warn.dep=TRUE)

    rprune <- switch(pmethod,
        backward   = leaps.backward(rprune),
        none       = leaps.backward(rprune), # for stats, won't actually prune
        exhaustive = leaps.exhaustive(rprune),
        forward    = leaps.forward(rprune),
        seqrep     = leaps.seqrep(rprune))

    list(rss.per.subset = as.vector(rprune$ress), # convert nx1 mat to vec
         # each row of prune.terms is a vec of term indices
         prune.terms    = convert.lopt(rprune$lopt, nprune))
}
# This calls the earth.c routine EvalSubsetsUsingXtxR.
# Unlike the leaps code, it can handle multiple responses (i.e. multiple y columns)

eval.subsets.xtx <- function(
    bx,
    y,
    pmethod,
    nprune,
    Force.xtx.prune,
    trace)
{
    bad.pmethod <- function()
    {
        stop0("pmethod=\"", pmethod, "\" is not allowed with 'eval.subsets.xtx'")
    }
    backward <- function(bx, y)
    {
        ncases <- nrow(bx)
        nterms <- ncol(bx)
        nresp <- ncol(y)
        stopifnot(is.double(bx))
        stopifnot(is.double(y))
        on.exit(.C("FreeEarth", PACKAGE="earth")) # if error or user interrupt, free mem
        # TODO replace .C call with alternative interface that doesn't require DUP=TRUE
        rv <- .C("EvalSubsetsUsingXtxR",
            prune.terms = matrix(0, nrow=nterms, ncol=nterms), # double PruneTerms[]
            rss.per.subset = vector(mode="numeric", length=nterms),
            as.integer(ncases),       # const int* pnCases
            as.integer(nresp),        # const int* pnResp
            as.integer(nterms),       # const int* pnMaxTerms
            bx,                       # const double bx[]
            y,                        # const double y[]
            as.double(max(trace, 0)), # in: const double* pTrace
            PACKAGE="earth")

        # above returns all subsets, so trim back to nprune below

        list(rss.per.subset = rv$rss.per.subset[seq_len(nprune)],
             prune.terms    = rv$prune.terms[seq_len(nprune), seq_len(nprune),
                                             drop=FALSE])
    }
    #--- eval.subsets.xtx starts here ---

    rprune <- switch(pmethod,
        backward   = backward(bx, y),
        none       = backward(bx, y), # for stats, won't actually prune
        exhaustive = bad.pmethod(),
        forward    = bad.pmethod(),
        seqrep     = bad.pmethod())
}
get.nused.preds.per.subset <- function(dirs, which.terms)
{
    # object was converted from mars? if so, ugly hack to allow plot routines to work
    if(is.null(which.terms))
        which.terms <- matrix(seq_len(ncol(dirs)), ncol(dirs), ncol(dirs))

    # allow which.terms to be a vector or matrix
    if(NROW(which.terms) == 1 || NCOL(which.terms) == 1)
        which.terms <- matrix(which.terms, nrow=1,
                              ncol=NROW(which.terms) * NCOL(which.terms == 1))

    nmodels <- NROW(which.terms)
    stopifnot(nmodels > 0)
    nused <- vector(mode="numeric", nmodels)
    for(i in seq_len(nmodels)) {
        check.which.terms(dirs, which.terms)
        nused[i] <- sum(0 != colSums(abs(
                             dirs[which.terms[i,,drop=FALSE], , drop=FALSE])))
    }
    nused
}
# If pmethod is exhaustive and bx is ill conditioned, change pmethod to
# backward. This prevents leaps.exhaustive returning error code -999.
#
# Note that bx should never be ill-conditioned (RegressAndFix should
# take care of that).  However it seems that dqrdc2 (called by
# RegressAndFix) does not detect certain types of ill conditioning (with
# any tol). This is probably because we are near the numerical noise floor
# and the column norms in dqrdc are not monotically decreasing.
# This change was made in Apr 2011.
#
# TODO This would be better handled by simply removing collinear cols in bx?

preprocess.exhaustive <- function(pmethod, nprune, Exhaustive.tol, trace, bx)
{
    pmethod <- "exhaustive"
    check.numeric.scalar(Exhaustive.tol)
    if(Exhaustive.tol < 0 || Exhaustive.tol > .1)
        stop0("illegal Exhaustive.tol ", Exhaustive.tol,
              ", try something like Exhaustive.tol=1e-8")
    sing.vals <- svd(bx)$d  # expensive
    bx.cond <- sing.vals[length(sing.vals)] / sing.vals[1]
    if(is.na(bx.cond) || bx.cond < Exhaustive.tol) {
        trace1(trace, "\n")
        warning0("forced pmethod=\"backward\" ",
            "(bx is ill conditioned, sing val ratio ",
            format(bx.cond, digits=2), ")")
        trace1(trace, "\n")
        pmethod <- "backward"
    }
    if(pmethod == "exhaustive")
        possibly.print.exhaustive.pruning.reminder(nprune, trace, bx, bx.cond)
    pmethod
}
print_pruning_pass <- function(trace, pmethod, penalty, nprune, selected.terms,
                               prune.terms, rss.per.subset, gcv.per.subset, dirs)
{
    nselected <- length(selected.terms)
    prev.grsq <- 0
    if(trace >= 3 && trace <= 7) {
        cat("Subset size        GRSq     RSq  DeltaGRSq nPreds")
        if(trace >= 4)
            cat("  Terms (col nbr in bx)")
        cat("\n")
        for(iterm in seq_along(rss.per.subset)) {
            grsq <- get.rsq(gcv.per.subset[iterm], gcv.per.subset[1])
            delta.grsq <- grsq - prev.grsq
            prev.grsq <- grsq
            selected <- prune.terms[iterm,]
            selected <- selected[selected != 0]
            cat0(if(iterm==nselected) "chosen "
                 else                 "       ",
                 format(iterm, width=4),
                 sprint("%12.4f ", grsq),
                 sprint("%7.4f",   get.rsq(rss.per.subset[iterm], rss.per.subset[1])),
                 sprint("%11.4f ", delta.grsq),
                 sprint("%6d",     get.nused.preds.per.subset(dirs, selected)),
                 "  ")
            if(trace >= 4)
                cat(selected)
            cat("\n")
        }
        cat("\n")
    }
    if(trace >= 5 && (pmethod == "backward" || pmethod == "forward"))
        check.one.term.per.step(prune.terms)
    if(trace >= 1) {
        cat0("Prune ", pmethod, " penalty ", penalty)
        if(pmethod != "cv")
            cat0(" nprune ", if(is.null(nprune)) "null" else nprune)
        cat0(": selected ", nselected, " of ")
        selected <- prune.terms[nselected,]
        selected <- selected[selected != 0]
        cat(nrow(dirs), "terms, and",
            get.nused.preds.per.subset(dirs, selected),
            "of", ncol(dirs), "preds\n")
        cat0("After pruning pass GRSq ",
            format(get.rsq(gcv.per.subset[nselected], gcv.per.subset[1]), digits=3),
            " RSq ",
            format(get.rsq(rss.per.subset[nselected], rss.per.subset[1]), digits=3),
            "\n")
    }
}
# This is called pruning.pass (not backward.pass) because pmethod may not be "backward".
#
# Note that pmethod="none" is equivalent to "backward" except that by
# default it retains all the terms created by the forward pass.  If nprune
# is specified, then it selects the terms that backward would have
# selected if backward selected nprune terms.
#
# If pmethod=="cv" we first do a normal pmethod="backward" pass with nprune=nterms,
# then select the subset using the nprune passed to this routine, rather than
# with the GCV.  The nprune passed to this routine will be machine generated
# as the nprune that gives the max mean oof rsq.

pruning.pass <- function(x, y, bx, # x, y, and bx are weighted if weights arg was used
                         pmethod, penalty, nprune,
                         trace, dirs, Force.xtx.prune, Exhaustive.tol)
{
    stopifnot(nrow(bx) == nrow(y))
    nterms <- nrow(dirs)
    check.integer.scalar(nprune, null.ok=TRUE, min=1)
    nprune.org <- nprune
    if(is.null(nprune) || pmethod=="cv")
        nprune <- nterms
    # else
    #     trace1(trace, "nprune=%g\n", nprune)
    nprune <- min(nprune, nterms)
    # Sep 2020: Keep best subset for all sizes up to nterms, not just nprune.
    #   This gives clearer graphs in earth_plotmodsel for nprune models.
    #   Don't do it for exhaustive because that's too slow.
    nprune.all <- if(pmethod != "exhaustive") nterms else nprune
    if(pmethod == "exhaustive")
        pmethod <- preprocess.exhaustive(pmethod, nprune, Exhaustive.tol, trace, bx)
    flush.console() # make sure previous messages get seen, pruning make take a while
    rv <- eval.model.subsets(bx, y,
                pmethod=if(pmethod=="cv") "backward" else pmethod,
                nprune.all, Force.xtx.prune, trace)
    rss.per.subset <- rv$rss.per.subset # RSS for each subset (across all responses)
    prune.terms    <- rv$prune.terms    # each row is a vec of term indices
    stopifnot(length(rss.per.subset) <= nprune.all)
    nprune     <- min(nprune,     length(rss.per.subset)) # probably unnecessary
    nprune.all <- min(nprune.all, length(rss.per.subset)) # probably unnecessary
    prune.terms <- prune.terms[seq_len(nprune.all), seq_len(nprune.all), drop=FALSE]
    stopifnot(all(prune.terms[,1] == 1)) # check intercept column
    gcv.per.subset <- get.gcv(rss.per.subset, seq_len(nprune.all), penalty, nrow(bx))
    if(!all(is.finite(rss.per.subset)))
        warning0("earth: non finite RSS in model subsets ",
                 "(see the rss.per.subset returned by earth)")
    check.vec(rss.per.subset, "rss.per.subset")
    selected.terms <- seq_len(nprune)
    if(pmethod == "cv") {
        # choose the subset at the nprune passed to this routine
        check.integer.scalar(nprune.org, min=1, max=nterms)
        selected.terms <- prune.terms[nprune.org,]
        selected.terms <- selected.terms[selected.terms != 0]
    } else if(pmethod != "none") {
        # choose the subset which has the lowest GCV in the vector of GCVS
        selected.terms <- prune.terms[which.min(gcv.per.subset[1:nprune]),]
        selected.terms <- selected.terms[selected.terms != 0]
    }
    print_pruning_pass(trace, pmethod, penalty, nprune.org, selected.terms,
                       prune.terms, rss.per.subset, gcv.per.subset, dirs)

    list(rss.per.subset = rss.per.subset, # vector of RSSs for each model (index on subset size)
         gcv.per.subset = gcv.per.subset, # vector of GCVs for each model (index on subset size)
         prune.terms    = prune.terms,    # triang mat: each row is a vector of term indices
         selected.terms = selected.terms) # vec of model terms in best model
}

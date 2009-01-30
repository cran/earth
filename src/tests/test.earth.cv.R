# test.earth.cv.R: test earth cross validation

library(earth)
if(!interactive())
    postscript()
data(ozone1)
data(trees)
data(etitanic)
options(warn=2)
options(width=200)

printh <- function(x, expect.warning=FALSE, max.print=0) # like print but with a header
{
    cat("===", deparse(substitute(x)))
    if(expect.warning)
        cat(" expect warning -->")
    else if (NROW(x) > 1)
        cat("\n")
    if (max.print > 0)
        print(head(x, n=max.print))
    else
        print(x)
}

cat("a0: trees\n\n")

set.seed(23)
a0 <- earth(Volume ~ ., data = trees, trace=0.5, nfold=3)
printh(a0$cv.rsq.tab)
printh(a0)
printh(summary(a0))

cat("a0a: trees with matrix interface\n\n")

set.seed(23)
a0a <- earth(trees[,-3], trees[,3], trace=0.5, nfold=3)
stopifnot(!identical(a0$cv.rsq.tab,a0a$cv.rsq.tab))
printh(a0a)
printh(summary(a0a))

cat("a1: trees with trace enabled\n\n")

set.seed(1)
a1 <- earth(Volume ~ ., data = trees, trace=1, nfold=3)
stopifnot(!identical(a0$cv.rsq.tab,a1$cv.rsq.tab))
printh(a1)
printh(summary(a1))

# test correct operation of update

cat("a2 <- update(a0) # should do cv\n")
set.seed(2)
a2 <- update(a0)
cat("a3 <- update(a0) # should do cv\n")
set.seed(3)
a3 <- update(a0, formula=Volume~.-Height)
printh(a3$cv.rsq.tab)
printh(a3)
printh(summary(a3))
cat("a4 <- update(a0, nfold=0, trace=.5)  # should not do cv\n")
set.seed(4)
a4 <- update(a0, nfold=0, trace=.5)
cat("a5 <- update(a4, trace=.5)           # should not do cv\n")
set.seed(5)
a5 <- update(a4)
cat("a5a <- update(a4, nfold=2, trace=.5) # should do cv\n")
set.seed(2)
a5a <- update(a4, nfold=2, trace=.5)

cat("a6: titanic data, one logical response\n\n")
survived. <- as.logical(etitanic$survived)
set.seed(6)
a6 <- earth(survived. ~ ., data=etitanic[,-2], degree=2, glm=list(family="binomial"), trace=0.5, nfold=3)
# TODO extend this list
printh(a6$cv.nterms)
printh(a6$cv.nvars)
printh(a6$cv.rsq.tab)
printh(a6$cv.maxerr.tab)
printh(a6$cv.deviance.tab)
printh(a6$cv.auc.tab)
printh(a6$cv.calib.int.tab)
printh(a6$cv.calib.slope.tab)
printh(a6$cv.auc.tab)
printh(a6$cv.cor.tab)
printh(a6$cv.groups)
printh(a6$cv.groups)
printh(a6$cv.list.)
printh(a6)
printh(summary(a6))
# expect a note: plotmo cannot use factors as axes for degree2 plots...
plotmo(a6)
printh(a6$cv.list[[2]])
printh(summary(a6$cv.list[[2]]))

cat("a6a: stratify=FALSE\n\n")
set.seed(6)
a6a <- earth(survived. ~ ., data=etitanic[,-2], degree=2, glm=list(family="binomial"), trace=0.5, nfold=3, stratify=FALSE)
printh(a6a$cv.rsq.tab)
printh(a6a)
printh(summary(a6a))

cat("a7: titanic data, multiple responses (i.e. 3 level factor)\n\n")
set.seed(3)
# keepxy is needed for summary and plotmo of submodels
a7 <- earth(pclass ~ ., data=etitanic, degree=2, glm=list(family="binomial"), trace=0.5, nfold=10, keepxy=TRUE)
printh(a7$cv.nterms)
printh(a7$cv.nvars)
printh(a7$cv.rsq.tab)
printh(a7$cv.maxerr.tab)
printh(a7$cv.deviance.tab)
printh(a7$cv.auc.tab)
printh(a7$cv.calib.int.tab)
printh(a7$cv.calib.slope.tab)
printh(a7$cv.auc.tab)
printh(a7$cv.cor.tab)
printh(a7$cv.groups)
printh(a7$cv.groups)
printh(a7$cv.list.)
printh(a7)
printh(summary(a7))
# expect a note: plotmo cannot use factors as axes for degree2 plots...
plotmo(a7, ycolumn=1)
printh(a7$cv.list[[3]])
printh(summary(a7$cv.list[[3]]))

cat("a7.wp: as above but with wp parameter\n\n")
set.seed(3)
a7.wp <- earth(pclass ~ ., data=etitanic, degree=2, glm=list(family="binomial"), trace=0.5, nfold=10, wp=c(1,3,1))
printh(a7.wp)
printh(summary(a7.wp))

# poisson models

counts <- c(18,  17,  15,  20,  10,  20,  25,  13,  12,
            18+2,17+2,15+2,20+2,10+2,20+2,25+2,13+2,12+2,
            18+3,17+3,15+3,20+3,10+3,20+3,25+3,13+3,12+3,
            18+4,17+4,15+4,20+4,10+4,20+4,25+4,13+4,12+4)
counts2 <- c(181,171,151,201,101,201,251,131,121,
             189,179,159,209,109,209,259,139,121,
             185,175,155,205,105,205,255,135,125,
             183,173,153,203,103,203,253,133,123)
outcome <- gl(3,1,4*9)
treatment <- gl(3,4*3)
d.AD <- data.frame(treatment, outcome, counts, counts2)

# one response poisson model
cat("a8p: one response poisson model\n\n")
a8p <- earth(counts ~ outcome + treatment, glm=list(family=poisson()), trace=0.5, pmethod="none", nfold=3)
printh(a8p)
printh(summary(a8p))
# two response poisson model
cat("a10: two response poisson model\n\n")
a10 <- earth(cbind(counts, counts2) ~ outcome + treatment, glm=list(fam="po"), trace=0.5, pmethod="none", nfold=3)
printh(a10)
printh(summary(a10))

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

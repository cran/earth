# test.weights.R

library(earth)
library(mda)
options(warn=1) # print warnings as they occur
check.equal <- function(x, y, msg="")
{
    diff <- x - y
    if (any(abs(diff) > 1e-8)) {
        cat(msg, "\n1st matrix:\n", sep="")
        print(x)
        cat("\n2nd matrix:\n")
        print(y)
        cat("\ndiff:\n")
        print(diff)
        stop("check.equal failed for ", msg, call.=FALSE)
    }
}
check.models.equal <- function(lm.mod, earth.mod)
{
    lm.mod.name <- deparse(substitute(lm.mod))
    earth.mod.name <- deparse(substitute(earth.mod))
    msg <- sprintf("%s vs %s", lm.mod.name, earth.mod.name)
    check.equal(lm.mod$coefficients,       earth.mod$coefficients,       msg=sprintf("%s coefficients", msg))
    check.equal(lm.mod$rss,                earth.mod$rss,                msg=sprintf("%s rss", msg))
    check.equal(lm.mod$residuals,          earth.mod$residuals,          msg=sprintf("%s residuals", msg))
    check.equal(summary(lm.mod)$r.squared, earth.mod$rsq,                msg=sprintf("%s rsq", msg))
    check.equal(summary(lm.mod)$r.squared, earth.mod$rsq.per.reponse[1], msg=sprintf("%s rsq.per.response", msg))
}
if(!interactive())
    postscript(paper="letter")

# artifical data
xxx <- 1:9
yyy <- 1:9
yyy[5] <- 9
data <- data.frame(x=xxx, y=yyy)
colnames(data) <- c("x", "y")

# Check against a linear model with weights, using linpreds.
# This also checks the backward pass's handling of weights.

lm1 <- lm(y~., data=data)
a1 <- earth(y~., data=data, linpreds=TRUE)
check.models.equal(lm1, a1)

weights <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
lm2 <- lm(y~., data=data, weights=weights)
a2  <- earth(y~., data=data, linpreds=TRUE, weights=weights)
check.models.equal(lm2, a2)

# check that we can get the weights from the data as per lm
lm2.a <- lm(y~xxx, data=data, weights=x) # weights from model frame
a2.a  <- earth(y~xxx, data=data, linpreds=TRUE, weights=x) # weights from model frame
a2.b  <- earth(y~xxx, data=data, linpreds=TRUE, weights=xxx) # weights from global env
check.models.equal(lm2.a, a2.a)
check.models.equal(a2.b, a2.a)

weights <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
lm3 <- lm(y~., data=data, weights=weights)
a3  <- earth(y~., data=data, linpreds=TRUE, weights=weights, trace=-1)
check.models.equal(lm3, a3)

lm4 <- lm(y~., data=data, weights=.1 * weights)
a4  <- earth(y~., data=data, linpreds=TRUE, weights=.1 * weights,
             minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
check.models.equal(lm4, a4)

# We want to see the effect only on the forward pass, so disable the
# backward pass with penalty=-1.  This also prevents "termination of the
# forward pass with a negative GRSq" with this artifical data.
#
# We can't use thresh=0, because then very small weights will still cause a usable
# reduction in RSq (remember that weights of zero are changed to very small weights
# in the current implementation).  So instead we use thresh=1e-8.
# This is a problem only with this very artifical data.  With real data, we
# want to use the standard thresh=.001, even with weights.

cat("=== a5.noweights ===\n")
par(mfrow = c(2, 2))
par(mar = c(3, 3, 3, 1))
par(mgp = c(1.5, 0.5, 0))
a5.noweights <- earth(y~., data=data,
                      minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a5.noweights, col.response=2, do.par=F, main="a5.noweights", grid.col="gray", jitter=0)
# TODO why does this model differ from the above model?
a5.noweights.force <- earth(y~., data=data, Force.weights=T,
                      minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a5.noweights.force, col.response=2, do.par=F, main="a5.noweights.force", grid.col="gray", jitter=0)

cat("=== a6.azeroweight ===\n")
a6.azeroweight  <- earth(y~., data=data, weights=c(1, 1, 1, 1, 0, 1, 1, 1, 1),
                         minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a6.azeroweight, col.response=2, do.par=F, main="a6.azeroweight", grid.col="gray", jitter=0)

cat("=== a7.asmallweight ===\n") # different set of weights (pick up notch in data, but with different forward pass RSq's)
a7.asmallweight  <- earth(y~., data=data, weights=c(1, 1, 1, 1, .5, 1, 1, 1, 1),
                          minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a7.asmallweight, col.response=2, do.par=F, main="a7.asmallweight", grid.col="gray", jitter=0)

cat("=== a7.xy.asmallweight ===\n") # x,y interface
a7.xy.asmallweight  <- earth(xxx, yyy, weights=c(1, 1, 1, 1, .5, 1, 1, 1, 1),
                          minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
check.models.equal(a7.xy.asmallweight, a7.xy.asmallweight)

cat("=== a8 ===\n")
par(mfrow = c(3, 2)) # new page
par(mar = c(3, 3, 3, 1))
par(mgp = c(1.5, 0.5, 0))
data$y <- c(0, 0, 0, 1, 0, 1, 1, 1, 1) != 0

# glm models first without weights
a8 <- earth(y~., data=data,
            linpreds=TRUE, glm=list(family=binomial),
            minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a8,
       col.response=2, do.par=F, main="a8 glm no weights\ntype=\"response\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
plotmo(a8, type="earth",
       col.response=2, do.par=F, main="a8 glm no weights\ntype=\"earth\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
glm <- glm(y~., data=data, family=binomial)
stopifnot(coefficients(a8$glm.list[[1]]) == coefficients(glm))

cat("=== a8.weights ===\n")
# now glm models with weights
glm.weights <- c(.8,1,1,.5,1,1,1,1,1)
a8.weights <- earth(y~., data=data,
                    linpreds=TRUE, glm=list(family=binomial),
                    weights=glm.weights,
                    minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a8.weights, type="response",
       col.response=2, do.par=F, main="a8.weights glm\ntype=\"response\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
plotmo(a8.weights, type="earth",
       col.response=2, do.par=F, main="a8.weights glm\ntype=\"earth\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
glm <- glm(y~., data=data, weights=glm.weights, family=binomial)
# TODO this fails if a weight is 0 in glm.weights
stopifnot(coefficients(a8.weights$glm.list[[1]]) == coefficients(glm))

cat("=== a8.weights ===\n")
# now glm models with weights
glm.weights <- c(.8,1,1,0,1,1,1,1,1)
a8.azeroweight <- earth(y~., data=data,
                    linpreds=TRUE, glm=list(family=binomial),
                    weights=glm.weights,
                    minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a8.azeroweight, type="response",
       col.response=2, do.par=F, main="a8.azeroweight glm\ntype=\"response\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
plotmo(a8.azeroweight, type="earth",
       col.response=2, do.par=F, main="a8.azeroweight glm\ntype=\"earth\"",
       grid.col="gray", ylim=c(-.2, 1.2), jitter=0)
# glm <- glm(y~., data=data, weights=glm.weights, family=binomial)
# print(coefficients(a8.azeroweight$glm.list[[1]]))
# print(coefficients(glm))
# # TODO this fails if a weight is 0 in glm.weights
# stopifnot(coefficients(a8.azeroweight$glm.list[[1]]) == coefficients(glm))

cat("=== plot.earth with weights ===\n")
# we also test id.n=TRUE and id.n=-1 here
old.par <- par(mfrow=c(2,2), mar=c(4, 3.2, 3, 3), mgp=c(1.6, 0.6, 0), oma=c(0,0,3,0), par(cex=1))
plot(a3, id.n=TRUE, SHOWCALL=TRUE, caption="compare a3 to to lm3", do.par=FALSE,
     which=c(3,4), caption.cex=1.5)
plot(lm3, id.n=9, which=c(1,2), sub.caption="")
par(old.par)

cat("=== plot.earth with earth-glm model and weights ===\n")
plot(a8, id.n=TRUE, caption="a8")
plot(a8.weights, id.n=TRUE, caption="a8.weights")
plot(a8.azeroweight, id.n=TRUE, caption="a8.azeroweight")
plot(a8.azeroweight, id.n=TRUE, delever=TRUE, caption="a8.azeroweight delever=TRUE")

# multivariate models

noise <- .01 * c(1,2,3,2,1,3,5,2,0)
data <- data.frame(x1=c(1,2,3,4,5,6,7,8,9), x2=c(1,2,3,3,3,6,7,8,9), y=(1:9)+noise)
data[5,] <- c(5, 5, 6)
colnames(data) <- c("x1", "x2", "y")

weights <- c(3, 2, 1, 1, 2, 3, 1, 2, 3)
lm20 <- lm(y~., data=data, weights=weights)
a20  <- earth(y~., data=data, linpreds=TRUE, weights=weights,
              minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
check.models.equal(lm20, a20)

a21.noweights <- earth(y~., data=data, # no weights for comparison
                       minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a21.noweights, col.resp=2, trace=-1, caption="a21.noweights", jitter=0)

weights <- c(1, 1, 1, 1, .5, 1, 1, 1, 1)
a10  <- earth(y~., data=data, weights=weights,
              minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a10, col.resp=2, caption="a10", jitter=0)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

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
plotmo(a5.noweights, col.response=2, do.par=F, main="a5.noweights", grid="gray")
# TODO why does this model differ from the above model?
a5.noweights.force <- earth(y~., data=data, Force.weights=T,
                      minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a5.noweights.force, col.response=2, do.par=F, main="a5.noweights.force", grid="gray")

cat("=== a6.azeroweight ===\n")
a6.azeroweight  <- earth(y~., data=data, weights=c(1, 1, 1, 1, 0, 1, 1, 1, 1),
                         minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a6.azeroweight, col.response=2, do.par=F, main="a6.azeroweight", grid="gray")

cat("=== a7.asmallweight ===\n") # different set of weights (pick up notch in data, but with different forward pass RSq's)
a7.asmallweight  <- earth(y~., data=data, weights=c(1, 1, 1, 1, .5, 1, 1, 1, 1),
                          minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=3)
plotmo(a7.asmallweight, col.response=2, do.par=F, main="a7.asmallweight", grid="gray")

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
       grid="gray", ylim=c(-.2, 1.2))
plotmo(a8, type="earth",
       col.response=2, do.par=F, main="a8 glm no weights\ntype=\"earth\"",
       grid="gray", ylim=c(-.2, 1.2))
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
       grid="gray", ylim=c(-.2, 1.2))
plotmo(a8.weights, type="earth",
       col.response=2, do.par=F, main="a8.weights glm\ntype=\"earth\"",
       grid="gray", ylim=c(-.2, 1.2))
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
       grid="gray", ylim=c(-.2, 1.2))
plotmo(a8.azeroweight, type="earth",
       col.response=2, do.par=F, main="a8.azeroweight glm\ntype=\"earth\"",
       grid="gray", ylim=c(-.2, 1.2))
# glm <- glm(y~., data=data, weights=glm.weights, family=binomial)
# print(coefficients(a8.azeroweight$glm.list[[1]]))
# print(coefficients(glm))
# # TODO this fails if a weight is 0 in glm.weights
# stopifnot(coefficients(a8.azeroweight$glm.list[[1]]) == coefficients(glm))

cat("=== plot.earth with weights ===\n")
plot(a8)
plot(a8.weights)
plot(a8.azeroweight)
plot(a8.azeroweight, delever=TRUE)

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
plotmo(a21.noweights, col.resp=2, trace=-1, caption="a21.noweights")

weights <- c(1, 1, 1, 1, .5, 1, 1, 1, 1)
a10  <- earth(y~., data=data, weights=weights,
              minspan=1, endspan=1, penalty=-1, thresh=1e-8, trace=-1)
plotmo(a10, col.resp=2, trace=-1, caption="a10")

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

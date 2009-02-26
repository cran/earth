# test.plotd: test earth with a biggish model

library(earth)

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

data(etitanic)
if(!interactive())
    postscript()
old.par <- par(no.readonly=TRUE)
par(mfrow=c(2,2))
par(cex = 0.7)
par(mar = c(4, 3, 1.7, 0.5))    # small margins and text to pack figs in
par(mgp = c(1.6, 0.6, 0))       # flatten axis elements

# test plotd

a1 <- earth(survived ~ ., data=etitanic, degree=2, glm=list(family=binomial))
printh(summary(a1))

plotd(a1) # default, no main
plotd(a1, main="2 lev fac, default")

plotd(a1, main="2 lev fac, kernel=epan adjust=.3",
      kernel="epan", adjust=.3, legend.names=c("leg1", "leg2"), 
      legend.pos=c(.3,4), legend.extra=TRUE)

plotd(a1, main="2 lev fac, type=earth", type="earth", xlab="my xlab", ylab="my ylab", 
      xlim=c(-.5, 1.5), zero.line=TRUE, vline.col="brown")

plotd(a1, main="2 lev fac, type=link", type="link", legend=FALSE, 
      col=c("pink", "lightblue"), lty=c(1,2), 
      vline.thresh=1, vline.col="gray", vline.lty=2)

# test plotd with histograms

plotd(a1, main="2 lev fac, default hist", hist=TRUE)

plotd(a1, main="2 lev fac, hist", 
      hist=TRUE, borders=c("green", "red"),
      xlab="my xlab", ylab="my ylab", xlim=c(-.2, 1.2),
      col=c("pink", "black"), lty=c(1,3), 
      vline.thresh=.65, vline.col="brown", vline.lty=2,
      legend=FALSE, breaks=5)

# test plotd with a numeric response and not glm

survivedn <- as.numeric(etitanic$survived)
a1a <- earth(survivedn ~ .-survived, data=etitanic, degree=2)
plotd(a1, main="numeric response")
plotd(a1, main="numeric response, type=class", type="class")
plotd(a1, main="hist, numeric response, type=class", hist=1, type="class")

a1b <- earth(survivedn - .5 ~ .-survived, data=etitanic, degree=2)
plotd(a1b, main="survivedn-.5, type=class, thresh=0", hist=1, type="class",thresh=0,vline.col="brown")
plotd(a1b, main="survivedn-.5, type=class, thresh=0.3", hist=1, type="class",thresh=0.3,vline.col="brown")

# examples from the man page

a <- earth(survived ~ ., data=etitanic, degree=2, glm=list(family=binomial))
printh(summary(a))

plotd(a, main="plotd default")
plotd(a, main="histogram", hist=TRUE)

plotd(a, main="histogram, type=class", 
      hist=TRUE, type="class", 
      col=c("brown", "black"), 
      legend.pos=c(.2,500), legend.cex=.7, legend.extra=TRUE, labels=TRUE)

# test plotd with a logical response

logical.survived <- as.logical(etitanic$survived) 
a2 <- earth(logical.survived ~ . - survived, data=etitanic, degree=2, glm=list(family=binomial))
printh(summary(a2))
plotd(a2, main="binary, default plotd")
plotd(a2, main="binary, default hist", hist=TRUE)

# test plotd with lm models

lm1 <- lm(survived ~ ., data=etitanic)
plotd(lm1)
plotd(lm1, main="lm1, survived")
plotd(lm1, hist=1, main="lm1, survived, hist=1, labels=1", labels=1)

lm2 <- lm(unclass(pclass)-1 ~ ., data=etitanic)
plotd(lm2)
plotd(lm2, main="lm2, pclass")

lm3 <- lm(cbind(survived, sin(age)) ~ ., data=etitanic) # nonsense model
plotd(lm3, xlim=c(-.5,1.5), thresh=.2, hist=1, main="lm3, NCOL(y)==2, thresh=.2")

lm4 <- lm(cbind(survived, sin(age), cos(age)) ~ ., data=etitanic) # nonsense model
plotd(lm4, xlim=c(-.5,1.5), thresh=0, hist=1, main="lm4, NCOL(y)==3, thresh=0")

# test plotd with glm models

glm.1 <- glm(survived ~ ., data=etitanic, family=binomial)
plotd(glm.1)
plotd(glm.1, main="glm.1, survived")

glm.2 <- glm(pclass ~ ., data=etitanic, family=binomial)
plotd(glm.2, main="glm.2, pclass")

# test plotd with a 3 level factor

a3 <- earth(pclass ~ ., data=etitanic, glm=list(family=binomial))
printh(summary(a3))

plotd(a3, main="3 lev fac, default")

plotd(a3, main="3 lev fac",
      xlab="my xlab", ylab="my ylab", xlim=c(-.2, 1.2),
      col=c("pink", "black", "green"), lty=c(1,3,1), 
      vline.thresh=.2, vline.col="blue", vline.lty=3,
      adjust=.3)

plotd(a3, main="3 lev fac, hist, default", hist=TRUE)

plotd(a3, main="3 lev fac, hist", 
      hist=TRUE, borders=c("green", "red", "blue"),
      xlab="my xlab", ylab="my ylab", xlim=c(-.2, 1.2),
      col=c("pink", "black", "green"), lty=c(1,2,3), 
      vline.thresh=.65, vline.col="gray", vline.lty=1,
      legend=FALSE, breaks=5)

plotd(a3, type="class", main="3 lev fac, type=class")
plotd(a3, type="class", main="3 lev fac, hist, type=class", hist=TRUE)

par(old.par)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

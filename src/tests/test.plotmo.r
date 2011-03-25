# test.plotmo.R: regression tests for plotmo
# Many of these tests are culled from man page examples and modified to try to confuse plotmo.
# Many of the plots are plotted twice so you can visually check by comparing
# plots in the same window, they should be substantially the same.
# Stephen Milborrow, Petaluma Jan 2007

print(R.version.string)

Trace <- 0

dopar <- function(nrows, ncols, caption = "")
{
    cat("                             ", caption, "\n")
    earth:::make.space.for.caption(caption)
    par(mfrow=c(nrows, ncols))
    par(mar = c(3, 3, 1.7, 0.5))
    par(mgp = c(1.6, 0.6, 0))
    par(cex = 0.7)
}
library(earth)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur
if(!interactive())
    postscript(paper="letter")

caption <- "basic earth test of plotmo"
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, degree1=2, degree2=4, caption=caption, trace=TRUE)

caption <- "test 5 x 5 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=51, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", caption=caption, trace=1)

caption <- "test 4 x 4 layout with ylab"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=30, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="ozone", caption=caption, trace=Trace)

caption <- "test 3 x 3 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=16, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", caption=caption, trace=Trace)

caption <- "test 2 x 2 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=9, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", caption=caption, trace=Trace)

caption <- "test 1 x 1 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=4, pmethod="n", degree=2)
plotmo(a, xlab="", ylab="", caption=caption, trace=Trace)

caption <- "test plotmo basic params"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(3,2,caption)
set.seed(1) # needed for reproducibility because of sample for rug in plotmo
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, caption=caption,
        main="test main", xlab="test xlab", ylab="test ylab", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=F, degree2=4, grid.func=mean, col.persp="white", ngrid2=10, phi=40, trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, lty.degree1=2, col.degree1=2, nrug=50, degree2=F, main="nrug=50", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, main="nrug=-1", trace=Trace)
set.seed(1)
plotmo(a, do.par=FALSE, degree1=1, nrug=20, ngrid1=50, degree2=F, main="ngrid1=50 nrug=20", trace=Trace)
plotmo(a, do.par=FALSE, degree1=NA, degree2=1, phi=60, box=F, r=100) # dots args

caption <- "test plotmo ylim"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(3,3,caption)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, caption=caption, xlab="ylim=default", trace=Trace)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=NA, xlab="ylim=NA", trace=Trace)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=c(0,20), xlab="ylim=c(0,20)", trace=Trace)

# term.plot calls predict.earth with an se parameter, even with termplot(se=FALSE)

caption <- "basic earth test against termplot"
dopar(4,4,caption)
earth:::make.space.for.caption("test caption1")
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, do.par=FALSE, ylim=NA, caption=caption, degree2=FALSE, trace=Trace)
cat("Ignore two warnings: predict.earth ignored argument \"se.fit\"\n")
termplot(a)

caption <- "test change order of earth predictors and cex"
dopar(4,4,caption)
a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
plotmo(a, do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,2), trace=Trace, cex=1)
termplot(a)

caption <- "test all1=TRUE"
a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
plotmo(a, caption=caption, all1=TRUE)
caption <- "test all2=TRUE"
print(summary(a))
plotmo(a, caption=caption, all2=TRUE)

a.global <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2, minspan=-1)
func1 <- function()
{
    caption <- "test call plotmo from within a function\nglobal dataframe (predict doy)"
    dopar(4,4,caption)
    a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2, minspan=-1)
    plotmo(a,        do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
           trace=Trace, type2="image", col.response=3, pch.response=20)
    plotmo(a.global, do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
          trace=Trace, type2="image", col.response=3, pch.response=20)

    caption <- "test call plotmo from within a function\nlocal dataframe, clip=FALSE"
    dopar(4,4,caption)
    ozone1.local <- ozone1[,c(1,2,3,4,5,6,7,8,10)]  # drop vis
    a <- earth(doy ~ humidity + temp + wind, data=ozone1.local, degree=2, minspan=-1)
    plotmo(a,        do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
           trace=Trace, type2="image", col.response=3, pch.response=20, clip=FALSE)
    plotmo(a.global, do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
           trace=Trace, type2="image", col.response=3, pch.response=20, clip=FALSE)

    caption <- "test call plotmo from within a function\nlocal dataframe using x,y interface"
    dopar(4,4,caption)
    x <- ozone1.local[,c(4,5,3)]    # humidty temp wind
    y <- ozone1.local[,9]           # doy
    a <- earth(x, y, degree=2, minspan=-1)
    plotmo(a,        do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
           trace=Trace, type2="image", col.response=3, pch.response=20)
    plotmo(a.global, do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,3),
           trace=Trace, type2="image", col.response=3, pch.response=20)
}
func1()

caption <- "test earth formula versus x,y model"
dopar(4,4,caption)
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, do.par=FALSE, caption=caption, trace=Trace)
a <- earth(ozone1[, -1], ozone1[,1], degree=2)
plotmo(a, do.par=FALSE, trace=Trace)

# single predictor
caption <- "test earth(O3~wind, data=ozone1, degree=2), single predictor"
dopar(2,2,caption)
a <- earth(O3~wind, data=ozone1, degree=2)
plotmo(a, trace=Trace)

caption = "se=2, earth(doy~humidity+temp+wind, data=ozone1) versus termplot (expect no se lines)"
dopar(3,2,caption)
a <- earth(doy~humidity + temp + wind, data=ozone1, degree=2)
cat("Ignore warning: predict.earth ignored argument \"se\"\n")
termplot(a)
cat("Ignore two warnings: predict.earth ignored argument \"se.fit\"\n")
plotmo(a, se=2, do.par=FALSE, ylim=NA, degree2=c(1:2), clip=FALSE, caption=caption, trace=Trace)

# test fix to bug reported by Joe Retzer, FIXED Dec 7, 2007
N <- 650
set.seed(2007)
q_4    <- runif(N, -1, 1)
q_2102 <- runif(N, -1, 1)
q_2104 <- runif(N, -1, 1)
q_3105 <- runif(N, -1, 1)
q_3106 <- runif(N, -1, 1)
q_4104 <- runif(N, -1, 1)
q_6101 <- runif(N, -1, 1)
q_6103 <- runif(N, -1, 1)
q_7104 <- runif(N, -1, 1)
q_3109 <- runif(N, -1, 1)
q_4103 <- runif(N, -1, 1)
q_2111 <- runif(N, -1, 1)
q_3107 <- runif(N, -1, 1)
q_3101 <- runif(N, -1, 1)
q_3104 <- runif(N, -1, 1)
q_7107 <- runif(N, -1, 1)
depIndex <- sin(1.0 * q_4 + rnorm(650, sd=.8)) + sin(1.8 * q_2102 + rnorm(650, sd=.8)) + sin(1.3 * q_2104 + rnorm(650, sd=.8)) + sin(1.4 * q_3105 + rnorm(650, sd=.8)) +
            sin(1.5 * q_3106 + rnorm(650, sd=.8)) + sin(1.6 * q_4104 + rnorm(650, sd=.8)) + sin(1.8 * q_6101 + rnorm(650, sd=.8)) + sin(1.8 * q_6103 + rnorm(650, sd=.8)) +
            sin(1.9 * q_7104 + rnorm(650, sd=.8)) + sin(2.0 * q_3109 + rnorm(650, sd=.8))

regDatCWD <- as.data.frame(cbind(depIndex, q_4, q_2102, q_2104, q_3105, q_3106, q_4104, q_6101, q_6103, q_7104, q_3109, q_4103, q_2111, q_3107, q_3101, q_3104, q_7107))
earthobj <- earth(depIndex ~  q_4+q_2102+q_2104+q_3105+q_3106+q_4104+q_6101+q_6103+q_7104+q_3109+q_4103+q_2111+q_3107+q_3101+q_3104+q_7107, data=regDatCWD)
print(summary(earthobj, digits = 2))
plotmo(earthobj)

# long predictor names

a.rather.long.in.fact.very.long.name.q_4 <- q_4
a.rather.long.in.fact.very.long.name.q_2102 <- q_2102
a.rather.long.in.fact.very.long.name.q_2104 <- q_2104
a.rather.long.in.fact.very.long.name.q_3105 <- q_3105
a.rather.long.in.fact.very.long.name.q_3106 <- q_3106
a.rather.long.in.fact.very.long.name.q_4104 <- q_4104
a.rather.long.in.fact.very.long.name.q_6101 <- q_6101
a.rather.long.in.fact.very.long.name.q_6103 <- q_6103
a.rather.long.in.fact.very.long.name.q_7104 <- q_7104
a.rather.long.in.fact.very.long.name.q_3109 <- q_3109
a.rather.long.in.fact.very.long.name.q_4103 <- q_4103
a.rather.long.in.fact.very.long.name.q_2111 <- q_2111
a.rather.long.in.fact.very.long.name.q_3107 <- q_3107
a.rather.long.in.fact.very.long.name.q_3101 <- q_3101
a.rather.long.in.fact.very.long.name.q_3104 <- q_3104
a.rather.long.in.fact.very.long.name.q_7107 <- q_7107
a.rather.long.in.fact.very.long.name.for.the.response <- depIndex
a.rather.long.in.fact.very.long.name.for.the.dataframe <-
        as.data.frame(cbind(
                a.rather.long.in.fact.very.long.name.for.the.response,
                a.rather.long.in.fact.very.long.name.q_4,
                a.rather.long.in.fact.very.long.name.q_2102,
                a.rather.long.in.fact.very.long.name.q_2104,
                a.rather.long.in.fact.very.long.name.q_3105,
                a.rather.long.in.fact.very.long.name.q_3106,
                a.rather.long.in.fact.very.long.name.q_4104,
                a.rather.long.in.fact.very.long.name.q_6101,
                a.rather.long.in.fact.very.long.name.q_6103,
                a.rather.long.in.fact.very.long.name.q_7104,
                a.rather.long.in.fact.very.long.name.q_3109,
                a.rather.long.in.fact.very.long.name.q_4103,
                a.rather.long.in.fact.very.long.name.q_2111,
                a.rather.long.in.fact.very.long.name.q_3107,
                a.rather.long.in.fact.very.long.name.q_3101,
                a.rather.long.in.fact.very.long.name.q_3104,
                a.rather.long.in.fact.very.long.name.q_7107))

a.rather.long.in.fact.very.long.name.for.the.modelA <-
        earth(a.rather.long.in.fact.very.long.name.for.the.response ~
                a.rather.long.in.fact.very.long.name.q_4 +
                a.rather.long.in.fact.very.long.name.q_2102 +
                a.rather.long.in.fact.very.long.name.q_2104 +
                a.rather.long.in.fact.very.long.name.q_3105 +
                a.rather.long.in.fact.very.long.name.q_3106 +
                a.rather.long.in.fact.very.long.name.q_4104 +
                a.rather.long.in.fact.very.long.name.q_6101 +
                a.rather.long.in.fact.very.long.name.q_6103 +
                a.rather.long.in.fact.very.long.name.q_7104 +
                a.rather.long.in.fact.very.long.name.q_3109 +
                a.rather.long.in.fact.very.long.name.q_4103 +
                a.rather.long.in.fact.very.long.name.q_2111 +
                a.rather.long.in.fact.very.long.name.q_3107 +
                a.rather.long.in.fact.very.long.name.q_3101 +
                a.rather.long.in.fact.very.long.name.q_3104 +
                a.rather.long.in.fact.very.long.name.q_7107,
                data = a.rather.long.in.fact.very.long.name.for.the.dataframe, minspan=-1)
print(summary(a.rather.long.in.fact.very.long.name.for.the.modelA, digits = 2))
plot(a.rather.long.in.fact.very.long.name.for.the.modelA)
plotmo(a.rather.long.in.fact.very.long.name.for.the.modelA)

a.rather.long.in.fact.very.long.name.for.the.modelC <-
        earth(x = a.rather.long.in.fact.very.long.name.for.the.dataframe[,-1],
          y = a.rather.long.in.fact.very.long.name.for.the.response,
                  degree = 3, minspan=-1)
print(summary(a.rather.long.in.fact.very.long.name.for.the.modelC, digits = 2))
plot(a.rather.long.in.fact.very.long.name.for.the.modelC)
plotmo(a.rather.long.in.fact.very.long.name.for.the.modelC)

data(etitanic)
a <- earth(survived ~ pclass+sex+age, data=etitanic, degree=2)
print(summary(a))
plotmo(a, trace=Trace, caption="plotmo with facs: pclass+sex+age")
plotmo(a, trace=Trace, clip=FALSE, degree2=FALSE, caption="plotmo (no degree2) with facs: pclass+sex+age")
plotmo(a, trace=Trace, clip=FALSE, grid.levels=list(pclass="2n", sex="ma"),
       caption="plotmo with grid.levels: pclass+sex+age")
# in above tests, all degree2 terms use facs
# now build a model with some degree2 term that use facs, some that don't
a <- earth(survived ~ pclass+age+sibsp, data=etitanic, degree=2)
print(summary(a))
plotmo(a, caption="plotmo with mixed fac and non-fac degree2 terms", border=NA)
plotmo(a, caption="plotmo with mixed fac and non-fac degree2 terms and grid.levels",
       grid.levels=list(pclass="2n", age=20), # test partial matching of grid levels, and numeric preds
       ticktype="d", nticks=2)

# check detection of illegal grid.levels argument
try(plotmo(a, grid.levels=list(pcla="1", pclass="2")))  # Expect error
try(plotmo(a, grid.levels=list(pclass="1", pcla="2")))  # Expect error
try(plotmo(a, grid.levels=list(pcla=1)))                # Expect error
try(plotmo(a, grid.levels=list(pcla=c("ab", "cd"))))    # Expect error
try(plotmo(a, grid.levels=list(pcla=NA)))               # Expect error
try(plotmo(a, grid.levels=list(pcla=Inf)))              # Expect error
try(plotmo(a, grid.levels=list(pcla=9)))                # Expect error
try(plotmo(a, grid.levels=list(age="ab")))              # Expect error
try(plotmo(a, grid.levels=list(age=NA)))                # Expect error
try(plotmo(a, grid.levels=list(age=Inf)))               # Expect error
try(plotmo(a, grid.lev=list(age=list(1,2))))            # Expect error

# more-or-less repeat above, but with glm models
a <- earth(survived ~ pclass+age+sibsp, data=etitanic, degree=2, glm=list(family=binomial))
print(summary(a))
plotmo(a, ylim=c(0, 1), caption="plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, ylim=c(0, 1), caption="plotmo glm with mixed fac and non-fac degree2 terms and grid.levels",
       grid.levels=list(pcl="2nd")) # test partial matching of variable name in grid levels
plotmo(a, type="earth", ylim=c(0, 1), caption="type=\"earth\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, type="link", ylim=c(0, 1), clip=FALSE, caption="type=\"link\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, type="class", ylim=c(0, 1), caption="type=\"class\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, ylim=c(0, 1), caption="default type (\"response\")\nplotmo glm with mixed fac and non-fac degree2 terms")
# now with different type2's
plotmo(a, do.par=FALSE, type2="persp",   theta=-20, degree1=FALSE, grid.levels=list(pclass="2nd"))
plotmo(a, do.par=FALSE, type2="contour", degree1=FALSE, grid.levels=list(pclass="2nd"))
plotmo(a, do.par=FALSE, type2="image",   degree1=FALSE, grid.levels=list(pclass="2nd"),
       col.response=as.numeric(etitanic$survived)+2, pch.response=20)
plotmo(a, do.par=FALSE, type="earth", type2="image", degree1=FALSE, grid.levels=list(pclass="2"))

# test vector main

a20 <- earth(O3 ~ humidity + temp + doy, data=ozone1, degree=2, glm=list(family=Gamma))

set.seed(1) # needed for nrug
plotmo(a20, nrug=-1)

set.seed(1) # needed for nrug
plotmo(a20, nrug=-1, caption="Test plotmo with a vector main",
       main=c("Humidity", "Temperature", "Day of year", "Humidity: Temperature", "Temperature: Day of Year"))

set.seed(1) # needed for nrug
cat("Expect warning below\n")
plotmo(a20, nrug=-1, caption="Test plotmo with a vector main, missing double titles",
       main=c("Humidity", "Temperature", "Day of year", "Humidity: Temperature"))

set.seed(1) # needed for nrug
cat("Expect warning below\n")
plotmo(a20, nrug=-1, caption="Test plotmo with a vector main, missing single titles",
       main=c("Humidity", "Temperature"))

aflip <- earth(O3~vh + wind + humidity + temp, data=ozone1, degree=2)

# test all1 and all2, with and without degree1 and degree2
plotmo(aflip, all2=T, caption="all2=T")
plotmo(aflip, all2=T, degree2=c(4, 2), caption="all2=T, degree2=c(4, 2)")
plotmo(aflip, all1=T, caption="all1=T")
plotmo(aflip, all1=T, degree1=c(3,1), degree2=NA, caption="all1=T, degree1=c(3,1), degree2=NA")

try(plotmo(aflip, no.such.arg=9)) # expect Error: plotmo: illegal argument "no.such.arg"
try(plotmo(aflip, degree1="all")) # Expect Error: degree1="all" is no longer legal, use all1=TRUE instead
try(plotmo(aflip, degree1="a"))   # Expect Error: degree1="all" is no longer legal, use all1=TRUE instead
try(plotmo(aflip, degree1="x"))   # Expect Error: degree1 must be an index vector (numeric or logical)
try(plotmo(aflip, degree2="all")) # Expect Error: degree2="all" is no longer legal, use all2=TRUE instead
try(plotmo(aflip, ycolumn=1))     # Expect Error: ycolumn is no longer legal, use nresponse instead
try(plotmo(aflip, title="abc"))   # Expect Error: "title" is illegal, use "caption" instead
try(plotmo(aflip, ticktype="d", ntick=3, tic=3, tick=9)) # expect Error : duplicated arguments "ticktype" "tic" "tick"
try(plotmo(aflip, ticktype="d", ntick=3, tic=3)) # expect Error : duplicated arguments "ticktype" "tic"
try(plotmo(aflip, ticktype="s", nt=3)) # expect Error : nticks is illegal with ticktype="simple"
try(plotmo(aflip, tic="s", nt=3)) # expect Error : nticks is illegal with ticktype="simple"
try(plotmo(aflip, tic="s", nt=3)) # expect Error : nticks is illegal with ticktype="simple"
try(plotmo(aflip, adj=8, adj=9)) # Error : duplicated arguments "adj" "adj"
try(plotmo(aflip, adj1=8, adj2=9)) # Error : plotmo: illegal argument "adj1"
try(plotmo(aflip, yc=8, x2=9)) # expect Error : "ycolumn" is no longer legal, use "nresponse" instead
try(plotmo(aflip, ticktype="d", ntick=3, ti=3)) # Error : "title" is illegal, use "caption" instead ("ti" taken to mean "title")
try(plotmo(aflip, ticktype="d", ntick=3, title=3)) # Error : "title" is illegal, use "caption" instead
try(plotmo(aflip, ticktype="d", ntick=3, tit=3, titl=7)) # Error : "title" is illegal, use "caption" instead ("tit" taken to mean "title")
try(plotmo(aflip, zlab="abc")) # expect Error : "zlab" is illegal, use "ylab" instead
try(plotmo(aflip, z="abc")) # expect Error : "zlab" is illegal, use "ylab" instead ("z" taken to mean "zlab")
try(plotmo(aflip, degree2="abc")) # expect Error : degree2 must be an index vector (numeric or logical)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

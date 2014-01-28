# test.earth.big: test earth with a biggish model

library(earth)
source("fast.postscript.R")
if(!interactive())
    fast.postscript(paper="letter")
options(warn=1) # print warnings as they occur
options(digits = 3)
set.seed(777)
N <- 1e4

ran <- function() runif(N, min=-1, max=1)

x <- cbind(
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(),
    ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran(), ran())

colnames(x) <- c(
    "x01", "x02", "x03", "x04", "x05", "x06", "x07", "x08", "x09", "x10",
    "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20",
    "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30",
    "x31", "x32", "x33", "x34", "x35", "x36", "x37", "x38", "x39", "x40",
    "x41", "x42", "x43", "x44", "x45", "x46", "x47", "x48", "x49", "x50",
    "x51", "x52", "x53", "x54", "x55", "x56", "x57", "x58", "x59", "x60",
    "x61", "x62", "x63", "x64", "x65", "x66", "x67", "x68", "x69", "x70",
    "x71", "x72", "x73", "x74", "x75", "x76", "x77", "x78", "x79", "x80",
    "x81", "x82", "x83", "x84", "x85", "x86", "x87", "x88", "x89", "x90",
    "x91", "x92", "x93", "x94", "x95", "x96", "x97", "x98", "x99", "x100")

func1 <- function(x)
{
     sin(6 * x[,1]) + x[,2] + 4 * x[,3] * x[,3]
}
y <- func1(x)
degree=1
nk=51
trace=0
cat("test 1: nrow(x)", nrow(x), "ncol(x)", ncol(x), "degree", degree, "nk", nk, "trace", trace, "\n")
a1 <- earth(x, y, nk=51, degree=degree, minspan=0, trace=trace)
print(summary(a1))
print(a1)
plot(a1)
plotmo(a1)

func2 <- function(x)
{
     x[,1] + x[,2] + x[,3] * x[,4]
}
cat("test 2: nrow(x)", nrow(x), "ncol(x)", ncol(x), "degree", degree, "nk", nk, "trace", trace, "\n")
y <- func2(x)
trace=4
a2 <- earth(x, y, nk=51, degree=degree, minspan=0, trace=trace)
print(summary(a2))

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}

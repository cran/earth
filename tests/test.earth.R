# test.earth.R
# The intention here is to check for porting problems by building
# a model with tracing enabled, and to call each earth routine once.
# For more comprehensive tests see earth\src\tests.
# Not tested: mars.to.earth, plot.earth, plot.earth.models, plotmo
library(earth)
options(digits=3)
data(trees)
fit1 <- earth(Volume ~ ., data = trees, degree = 2, trace = 4)
print(fit1, digits = 3)
print(summary(fit1, digits = 3))
print(evimp(fit1, trim=FALSE))
fit2 <- update(fit1, nprune = 3, trace = 0)
print(fit2, digits = 3)
print(earth:::get.nterms.per.degree(fit1))
print(earth:::get.nused.preds.per.subset(fit1$dirs, fit1$which.terms))
print(variable.names(fit1))
print(deviance(fit1))
print(earth:::reorder.earth(fit1))
print(predict(fit1, c(70, 12)))
fit3 <- earth(Volume ~ ., data = trees, degree = 1, minspan = 1)
print(fit3, digits = 4)
print(summary(fit3, digits = 3, style = "pmax"))

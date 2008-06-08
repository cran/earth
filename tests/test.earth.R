# test.earth.R
# The intention here is to check for porting problems by building
# a model with tracing enabled, and to call each earth routine once.
# For more comprehensive tests see earth\src\tests.
# Not tested: mars.to.earth, plot.earth, plot.earth.models, plotmo
library(earth)
options(digits=3)
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2, trace = 4)
print(a, digits = 3)
print(summary(a, digits = 3))
print(evimp(a, trim=FALSE))
a <- update(a, nprune = 5, trace = 0)
print(a, digits = 3)
print(earth:::get.nterms.per.degree(a))
print(earth:::get.nused.preds.per.subset(a$dirs, a$which.terms))
print(variable.names(a)[1:3])
print(deviance(a))
print(earth:::reorder.earth(a))
print(predict(a, c(5710, 4, 28, 40, 2693, -25, 87, 250, 33)))
a <- earth(O3 ~ ., data = ozone1, degree = 1, minspan = 1)
print(a, digits = 6)
print(summary(a, digits = 6, style = "pmax"))

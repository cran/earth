# test.earth.R
# Check for porting problems by building a few simple models.
# For much more comprehensive tests see earth\inst\slowtests.
library(earth)
options(digits=3) # prevent floating point implementation issues across machines
data(trees)
earth.mod <- earth(Volume~., data=trees)
print(summary(earth.mod))
allowed.func <- function(degree, pred, parents, namesx)
{
    namesx[pred] != "Height"  # disallow "Height"
}
earth.mod2 <- earth(Volume~., data=trees, allowed=allowed.func)
print(summary(earth.mod2))

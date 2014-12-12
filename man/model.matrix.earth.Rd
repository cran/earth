\name{model.matrix.earth}
\alias{model.matrix.earth}
\title{Get the earth basis matrix}
\description{
Get the basis matrix of an \code{\link{earth}} object.
}
\usage{
\method{model.matrix}{earth}(object = stop("no 'object' arg"),
    x = NULL, subset = NULL, which.terms = NULL,
    \dots,
    env = parent.frame(),
    trace = 0,
    Callers.name = "model.matrix.earth")
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{x}{
    An input matrix with the same number of columns as the \code{x} matrix
    used to construct the original \code{\link{earth}} object.
    Default is NULL, meaning use the original \code{x} matrix after
    taking the original \code{subset}, if any.
  }
  \item{subset}{
    Which rows to use in \code{x}.
    Default is NULL, meaning use all of \code{x}.
  }
  \item{which.terms}{
    Which terms to use.
    Default is NULL, meaning use \code{object$selected.terms}.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
  \item{env}{
     For internal use.
  }
  \item{trace}{
     Default 0. Set to non-zero to see which data \code{model.matrix.earth} is using.
  }
  \item{Callers.name}{
     For internal use (used by earth in trace messages).
  }
}
\value{
  A basis matrix \code{bx} of the same form returned by \code{\link{earth}}.\cr

  If \code{x}, \code{subset}, and \code{which.terms} are all NULL, this
  function returns the object's \code{bx}.  In this case, it is perhaps easier
  to simply use \code{object$bx}.\cr

The format of \code{bx} is described in \code{\link{earth}}.
The matrix \code{bx} can be used
as the input matrix to \code{\link{lm}} or \code{\link{glm}},
as shown below in the example.
In fact, that is what earth does internally after the pruning pass ---
it calls \code{\link{lm.fit}},
and additionally \code{\link{glm}} if earth's \code{glm} argument is used.
}
\seealso{
  \code{\link{earth}}
}
\examples{
data(trees)
earth.mod <- earth(Volume ~ ., data = trees)
summary(earth.mod, decomp = "none") # "none" to print terms in same seq as lm.mod below

# yields:
#  Call: earth(formula=Volume~., data=trees)
#
#               coefficients
#  (Intercept)         37.93
#  h(Girth-16)          7.40
#  h(16-Girth)         -3.92
#  h(Height-75)         0.48
#
#  Selected 4 of 6 terms, and 2 of 2 predictors
#  Termination condition: RSq changed by less than 0.001 at 6 terms
#  Importance: Girth, Height
#  Number of terms at each degree of interaction: 1 3 (additive model)
#  GCV 12.7    RSS 236    GRSq 0.954    RSq 0.971

bx <- model.matrix(earth.mod)         # equivalent to bx <- earth.mod$bx
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
summary(lm.mod)                       # yields same coeffs as above summary
                                      # displayed t values are not meaningful

# yields (note that the displayed stderrs and p values are meaningless):
#
#  Call:
#  lm(formula = trees$Volume ~ bx[, -1])
#
#  Residuals:
#     Min     1Q Median     3Q    Max
#  -5.272 -2.228  0.373  2.016  4.947
#
#  Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
#  (Intercept)            37.927      1.291   29.38  < 2e-16
#  bx[, -1]h(16-Girth)    -3.919      0.275  -14.25  4.4e-14
#  bx[, -1]h(Girth-16)     7.401      0.727   10.18  9.8e-11
#  bx[, -1]h(Height-75)    0.484      0.176    2.76     0.01
#
#  Residual standard error: 2.96 on 27 degrees of freedom
#  Multiple R-squared:  0.971,  Adjusted R-squared:  0.968
#  F-statistic:  300 on 3 and 27 DF,  p-value: <2e-16
}
\keyword{models}

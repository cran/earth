\name{model.matrix.earth}
\alias{model.matrix.earth}
\title{Get the "earth" basis matrix}
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
a <- earth(Volume ~ ., data = trees)
summary(a, decomp = "none") # "none" to print terms in same seq as a.lm below

# yields:
#  Call: earth(formula=Volume~., data=trees)
#
#               coefficients
#  (Intercept)        27.246
#  h(Girth-14)         6.177
#  h(14-Girth)        -3.266
#  h(Height-72)        0.491
#
#  Selected 4 of 6 terms, and 2 of 2 predictors
#  Importance: Girth, Height
#  Number of terms at each degree of interaction: 1 3 (additive model)
#  GCV 10.6    RSS 197    GRSq 0.962    RSq 0.976

bx <- model.matrix(a)               # equivalent to bx <- a$bx
a.lm <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
summary(a.lm)                       # yields same coeffs as above summary
                                    # displayed t values are not meaningful
# yields:
#  Call:
#  lm(formula = trees$Volume ~ bx[, -1])
#
#  Residuals:
#     Min     1Q Median     3Q    Max
#  -4.882 -1.770  0.281  1.646  4.983
#
#  Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
#  (Intercept)            27.246      1.123   24.26  < 2e-16
#  bx[, -1]h(Girth-14)     6.177      0.354   17.44  3.2e-16
#  bx[, -1]h(14-Girth)    -3.266      0.335   -9.76  2.4e-10
#  bx[, -1]h(Height-72)    0.491      0.123    3.99  0.00045
#
#  Residual standard error: 2.7 on 27 degrees of freedom
#  Multiple R-squared: 0.976,   Adjusted R-squared: 0.973
#  F-statistic:  361 on 3 and 27 DF,  p-value: <2e-16
}
\keyword{models}

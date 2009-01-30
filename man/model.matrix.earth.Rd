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
  A \code{bx} matrix of the same form returned by \code{\link{earth}}.\cr

  If \code{x}, \code{subset}, and \code{which.terms} are all NULL, this
  function returns the object's \code{bx}.  In this case, it is perhaps easier
  to simply use \code{object$bx}.\cr

The format of \code{bx} is described in \code{\link{earth}}.
The basis matrix \code{bx} can be used, for example,
as the input matrix to \code{\link{lm}} or \code{\link{glm}} as shown below.
In fact, that is what earth does internally --- 
it calls \code{\link{lm.fit}}, 
and additionally \code{\link{glm}} if earth's \code{glm} argument is used.
}
\seealso{
  \code{\link{earth}},
}
\examples{
data(trees)
a <- earth(Volume ~ ., data = trees)
summary(a, decomp = "none") # "none" to print terms in same seq as a.lm below

# yields:
#  Call: earth(formula = Volume ~ ., data = trees)
#
#                Volume
#  (Intercept)   23.208
#  h(Girth-12.9)  5.746
#  h(12.9-Girth) -2.866
#  h(Height-76)   0.718
#
#  Selected 4 of 5 terms, and 2 of 2 predictors
#  Estimated importance: Girth Height
#  Number of terms at each degree of interaction: 1 3 (additive model)
#  GCV 11.48697   RSS 213.4354   GRSq 0.958859   RSq 0.9736697

bx <- model.matrix(a)               # equivalent to bx <- a$bx
a.lm <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
summary(a.lm)                       # yields same coeffs as above summary
                                    # displayed p values are not meaningful
# yields:
#  Call:
#  lm(formula = trees$Volume ~ bx[, -1])
#
#  Residuals:
#      Min      1Q  Median      3Q     Max
#  -5.2900 -1.9782  0.0712  1.9209  4.1267
#
#  Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
#  (Intercept)             23.208      0.997   23.28  < 2e-16
#  bx[, -1]h(Girth-12.9)    5.746      0.297   19.35  < 2e-16
#  bx[, -1]h(12.9-Girth)   -2.866      0.444   -6.45  6.5e-07
#  bx[, -1]h(Height-76)     0.718      0.175    4.10  0.00034
#
#  Residual standard error: 2.81 on 27 degrees of freedom
#  Multiple R-squared: 0.974,   Adjusted R-squared: 0.971
#  F-statistic:  333 on 3 and 27 DF,  p-value: <2e-16
}
\keyword{models}

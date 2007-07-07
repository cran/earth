\name{model.matrix.earth}
\alias{model.matrix.earth}
\title{Get the 'earth' basis matrix}
\description{
Get the basis matrix of an \sQuote{earth} object.
}
\usage{
\method{model.matrix}{earth}(object, x, subset, which.terms, \dots)
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
    Default is NULL, meaning use \code{object$which.terms}.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
  A \code{bx} matrix of the same form returned by \code{\link{earth}}.\cr

  If \code{x}, \code{subset}, and \code{which.terms} are all NULL, this
  function returns the object's \code{bx}.  In this case, it is perhaps easier
  to simply use \code{object$bx}.\cr

The format of \code{bx} is described in \code{\link{earth}}.
The basis matrix \code{bx} can be used, for example,
as the input matrix to \code{\link{lm}} as shown below.
}
\seealso{
  \code{\link{earth}},
  \code{\link{get.nterms.per.degree}},
  \code{\link{get.nused.preds.per.subset}}
}
\examples{
data(trees)
a <- earth(Volume ~ ., data = trees)
summary(a, decomp = "none") # "none" to print terms in same seq as a.lm below

# yields:
#    Call:
#    earth(formula = Volume ~ ., data = trees)
#    
#    Expression:
#      26.64062 
#      +  6.049141 * pmax(0,  Girth -   13.7) 
#      -  3.202820 * pmax(0,   13.7 -  Girth) 
#      + 0.6079078 * pmax(0, Height -     75) 
#    
#    Number of cases: 31
#    Selected 4 of 5 terms, and 2 of 2 predictors
#    Number of terms at each degree of interaction: 1 3 (additive model)
#    GCV: 10.56081     RSS: 196.2266     GRSq: 0.962176    RSq: 0.9757927

bx <- model.matrix(a)               # equivalent to bx <- a$bx
a.lm <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
summary(a.lm)                       # yields same coeffs as above summary
                                    # displayed p values are not meaningful
# yields:
#    Call:
#    lm(formula = trees$Volume ~ bx[, -1])
#    
#    Residuals:
#          Min        1Q    Median        3Q       Max 
#    -4.691469 -1.987719  0.008402  1.421682  4.820064 
#    
#    Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)
#    (Intercept)            26.6406     1.0019  26.590  < 2e-16
#    bx[, -1]h(Girth-13.7)   6.0491     0.3344  18.087  < 2e-16
#    bx[, -1]h(13.7-Girth)  -3.2028     0.3502  -9.146  9.3e-10
#    bx[, -1]h(Height-75)    0.6079     0.1531   3.972 0.000477
#    
#    Residual standard error: 2.696 on 27 degrees of freedom
#    Multiple R-Squared: 0.9758,   Adjusted R-squared: 0.9731 
#    F-statistic: 362.8 on 3 and 27 DF,  p-value: < 2.2e-16 

}
\keyword{models}

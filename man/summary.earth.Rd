\name{summary.earth}
\alias{summary.earth}
\alias{print.summary.earth}
\title{Summary method for 'earth' objects}
\description{
Summary method for 'earth' objects.
}
\usage{
\method{summary}{earth}(object = stop("no 'object' arg"),
        digits = getOption("digits"), \dots)

\method{print}{summary.earth}(x = stop("no 'x' arg"), digits, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument for \code{summary.earth}.
  }
  \item{x}{
    A \code{\link{summary.earth}} object.
    This is the only required argument for \code{print.summary.earth}.
  }
  \item{digits}{
    The number of significant digits.
    For \code{summary.earth} the default is \code{getOption("digits")}.
    For \code{print.summary.earth} the default is the \code{$digits} component of \code{object}.
  }
  \item{\dots}{
    Extra arguments are passed to \code{\link{format.earth}}.
    Use \code{decomp="none"} to order the terms as created by the forward.pass.
  }
}
\value{
The value is the same as that returned by \code{\link{earth}}
but with the following extra components.

  \item{string}{
    String created by \code{\link{format.earth}}.
    For multiple response models, a vector of strings.
  }
  \item{digits}{
    The \code{digits} argument, passed on to \code{print.summary.earth}.
  }
}
\seealso{
  \code{\link{earth}},
  \code{\link{format.earth}}
}
\examples{
a <- earth(Volume~ ., data = trees)
summary(a, digits = 2)

# yields:
#    Call:
#    earth(formula = Volume ~ ., data = trees)
#    
#    Expression:
#      23 
#      +  5.7 * pmax(0,  Girth -     13) 
#      -  2.9 * pmax(0,     13 -  Girth) 
#      + 0.72 * pmax(0, Height -     76) 
#    
#    Number of cases: 31
#    Selected 4 of 5 terms, and 2 of 2 predictors
#    Number of terms at each degree of interaction: 1 3 (additive model)
#    GCV: 11          RSS: 213         GRSq: 0.96      RSq: 0.97
}
\keyword{models}

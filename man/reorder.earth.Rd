\name{reorder.earth}
\alias{reorder.earth}
\title{Order the terms in an 'earth' model}
\description{
Given an \sQuote{earth} object,
return an index vector ordered using the specified \code{decomp}.
The index vector is suitable for indexing into \code{object$selected.terms}
and \code{object$coefficients}.
}
\usage{
\method{reorder}{earth}(x = stop("no 'x' arg"),
             which.terms = x$selected.terms, decomp = c("anova", "none"),
             degree = 99, min.degree  = 0, \dots)
}
\arguments{
  \item{x}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{which.terms}{
    Which terms to use.
    Default is \code{x$which.terms}.
  }
  \item{decomp}{
        One of\cr
        \code{"none"}   order the terms as created during the earth forward pass
                  i.e. no reordering.\cr
        \code{"anova"}  (default) order the terms using the \sQuote{anova decomposition}
                  i.e. in increasing order of interaction.\cr
  }
  \item{degree}{
    Maximum order of interaction.
    Terms of order greater than \code{degree} will be dropped.
    Default is 99 which is effectively all.
    Set \code{degree=0} to return just the intercept.
  }
  \item{min.degree}{
    Minimum order of interaction.
    Terms of order less than \code{min.degree} will be dropped.
    Default is 0, meaning all.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
An index vector suitable for indexing into \code{object$selected.terms}
and \code{object$coefficients}, and ordered using the specified \code{decomp}.
}
\seealso{
  \code{\link{earth}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
reorder(a, decomp = "none")

# yields:
# [1]  1  2  3  4  5  6  7  8  9 10 11

reorder(a)   # defaults to decomp = "anova"

# yields:
# [1]  1  2 11  8  3  4  9  5  6  7 10

a$selected.terms[reorder(a)]

# yields:
# [1]  1  2 19 13  6  7 14  9 10 11 17
}
\keyword{models}
\keyword{regression}

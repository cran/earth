\name{get.nterms.per.degree}
\alias{get.nterms.per.degree}
\title{Number of terms at each degree of interaction in an 'earth' object}
\description{
Get the number of terms at each degree of interaction in an \sQuote{earth} object.
}
\usage{
get.nterms.per.degree(object, which.terms = object$selected.terms)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
  }
  \item{which.terms}{
    Typically the \code{selected.terms} component of the \code{\link{earth}} object.
  }
}
\value{
  A vector showing the number of terms
  at each degree of interaction in the model
  specified by \code{which.terms}.
  The first entry is the intercept entry, and is always 1 because
  there is always one intercept.
  The second entry is the number of additive terms,
  and so on.
}
\seealso{
  \code{\link{earth}},
  \code{\link{get.nused.preds.per.subset}},
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
get.nterms.per.degree(a)  # yields 1 5 5
}
\keyword{models}

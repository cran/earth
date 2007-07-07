\name{get.nused.preds.per.subset}
\alias{get.nused.preds.per.subset}
\title{Number of predictors in an 'earth' object}
\description{
Get the number of predictors in a model which is a subset of an \sQuote{earth} model.
}
\usage{
get.nused.preds.per.subset(dirs, which.terms)
}
\arguments{
  \item{dirs}{
    \code{dirs} component of an \code{\link{earth}} object.
  }
  \item{which.terms}{
   Typically the \code{selected.terms} or \code{prune.terms}
   component of an \code{\link{earth}} object.
  }
}
\value{
  The number of predictors used in the model specified by \code{which.terms}.
  If \code{which.terms} is a matrix such as \code{prune.terms}, the value
  is a vector.
}
\seealso{
  \code{\link{earth}},
  \code{\link{get.nterms.per.degree}},
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
get.nused.preds.per.subset(a$dirs, a$selected.terms)

# yields:
#    [1] 8     # there are 8 predictors in selected.terms

get.nused.preds.per.subset(a$dirs, a$prune.terms)

# yields:
#    [1] 0 1 2 3 4 5 5 6 6 8 8 8 8 8 8 8 8 8 8 8 8

}
\keyword{models}

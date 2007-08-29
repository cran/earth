\name{predict.earth}
\alias{predict.earth}
\title{Predict with an 'earth' model}
\description{
Predict with an 'earth' model.
}
\usage{
\method{predict}{earth}(object, newdata, \dots, type=c("response", "terms"))
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{newdata}{
    Make predictions using the \code{newdata}.
    Can be a dataframe, a matrix, or a vector with the same length as a row of the input matrix x.
    Default is NULL, meaning return \code{object}'s fitted values.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
  \item{type}{
     Type of prediction.  One of:\cr
\code{"response"} (default)\cr
\code{"terms"} return a response suitable for \code{\link{termplot}}.  Only
the additive terms and the first response (for multi-response models) are returned.\cr
  }
}
\value{
If \code{type="response"},  the predicted values.\cr
If \code{type="terms"},  a matrix with each column showing the contribution of a predictor.
}
\seealso{
  \code{\link{earth}},
  \code{\link{predict}}
}
\examples{
data(trees)
a <- earth(Volume ~ ., data = trees)
predict(a)           # same as a$fitted.values
predict(a, c(10,80)) # yields 17.82973
}
\keyword{models}

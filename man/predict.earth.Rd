\name{predict.earth}
\alias{predict.earth}
\title{Predict with an "earth" model}
\description{
Predict with an \code{\link{earth}} model.
}
\usage{
\method{predict}{earth}(object = stop("no 'object' arg"), newdata = NULL,
                        type = c("link", "response", "earth", "terms"),
                        trace = FALSE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{newdata}{
    Make predictions using \code{newdata}, which
    can be a dataframe, a matrix, or a vector with length equal to a multiple of the row length
    of the original input matrix \code{x}.
    Default is NULL, meaning return the fitted values from the training set.
  }
  \item{type}{
     Type of prediction.
     One of \code{"link"} (default), \code{"response"}, \code{"earth"}, or \code{"term"}.\cr\cr
        The \code{"link"}, \code{"response"}, and \code{"earth"} choices all return the same value unless
        the earth model has a GLM component (i.e., unless the \code{glm} argument
        was used in the original call to \code{\link{earth}}).
        See \code{\link{predict.glm}} for a description of the \code{"link"} and \code{"response"} options.
        (In brief: for logistic models use \code{type="response"} to return probabilities.) \cr\cr
        \code{"terms"} returns an \code{"earth"} response suitable for \code{\link{termplot}}.
        Only the additive terms and the first response (for multi-response models) are returned.
        Also, \code{"terms"} always returns the earth (not GLM) terms.\cr\cr
  }
  \item{trace}{
     Default FALSE. Set to TRUE to see which data, subset, etc.  \code{predict.earth} is using.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
If \code{type="terms"},  a matrix with each column showing the contribution of a predictor.
Else the predicted values (a matrix for multiple response models).\cr
}
\seealso{
  \code{\link{earth}},
  \code{\link{predict}}
}
\examples{
data(trees)
a <- earth(Volume ~ ., data = trees)
predict(a)           # same as a$fitted.values
predict(a, c(10,80)) # yields 17.76888
}
\keyword{models}

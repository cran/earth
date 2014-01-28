\name{predict.earth}
\alias{predict.earth}
\title{Predict with an "earth" model}
\description{
Predict with an \code{\link{earth}} model.
}
\usage{
\method{predict}{earth}(object = stop("no 'object' arg"), newdata = NULL,
        type = c("link", "response", "earth", "class", "terms"),
        thresh = .5, trace = FALSE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{newdata}{
    Make predictions using \code{newdata}, which
    can be a data frame, a matrix, or a vector with length equal to a multiple of the number of columns
    of the original input matrix \code{x}.
    Note that this is more flexible than the predict methods for most R models.
    NAs are allowed (and the predicted value will be NA unless the NAs are in variables
    that are unused in the earth model).
    Default is NULL, meaning return values predicted from the training set.
  }
  \item{type}{
     Type of prediction.
     One of \code{"link"} (default), \code{"response"}, \code{"earth"}, \code{"class"}, or \code{"terms"}.
     See the \bold{Note} below.
  }
  \item{thresh}{
    Threshold, a value between 0 and 1 when predicting a probability.
    Only applies when \code{type="class"}.
    Default is .5.
    See the \bold{Note} below.
  }
  \item{trace}{
     Default FALSE. Set to TRUE to see which data, subset, etc.  \code{predict.earth} is using.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
The predicted values (a matrix for multiple response models).
Except if \code{type="terms"}, then a matrix with each column showing the contribution of a predictor.
}
\note{
\bold{Predicting with standard earth models}

    Use the default \code{type="link"}, or possibly \code{type="class"}.

    Actually, the \code{"link"}, \code{"response"}, and \code{"earth"}
    choices all return the same value unless the \code{glm} argument
    was used in the original call to \code{\link{earth}}.

\bold{Predicting with earth-GLM models}

    This section applies to earth models with a GLM component, i.e.,
    when the \code{glm} argument was used
    in the original call to \code{\link{earth}}.

    The \code{"link"} and \code{"response"} options:
    see \code{\link{predict.glm}} for a description of these.
    In brief: for logistic models
    use \code{type="link"} to get log-odds and
    \code{type="response"} to get probabilities.

    Use option \code{"earth"} to get the linear fit (this gives the prediction you would get
    if your original call to earth had no \code{glm} argument).

\bold{Predicting with "class"}

    Use option \code{"class"} to get the predicted class.
    With option \code{"class"}, this function first makes predictions with
    \code{type="response"} and then assigns the predicted values to classes as follows:

    (i) When \code{y} is a \emph{logical}, predict TRUE if
    the predicted probability is greater than \code{thresh}.

    (ii) When \code{y} is a \emph{numeric}, predict TRUE if
    the predicted value is greater than \code{thresh}.
    Actually, this is identical to the above case,
    although \code{thresh} here may legitimately be a value
    outside the 0...1 range.

    (iii) When \code{y} is a \emph{two level factor},
    predict the second level if its probability is more than \code{thresh}.
    In other words, with the default \code{thresh=.5} predict the most probable level.

    (iv) When \code{y} is a \emph{three or more level factor},
    predict the most probable level (and \code{thresh} is ignored).

\bold{Predicting with "terms"}

    The \code{"terms"} option returns a \code{"link"} response suitable for \code{\link{termplot}}.
    Only the additive terms and the first response (for multi-response models) are returned.
    Also, \code{"terms"} always returns the earth terms, and ignores the GLM component
    of the model, if any.
}
\seealso{
  \code{\link{earth}},
  \code{\link{predict}}
}
\examples{
data(trees)
a <- earth(Volume ~ ., data = trees)
predict(a)           # same as a$fitted.values
predict(a, c(10,80)) # yields 18.11
}
\keyword{models}

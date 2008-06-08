\name{residuals.earth}
\alias{residuals.earth}
\alias{resid.earth}
\title{Residuals for an "earth" model}
\description{
Residuals for an \code{\link{earth}} model.
}
\usage{
\method{residuals}{earth}(object = stop("no 'object' arg"),
          type = NULL, warn = TRUE, \dots)

\method{resid}{earth}(object = stop("no 'object' arg"),
      type = NULL, warn = TRUE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{type}{
     One of\cr
        \code{"earth"} (default) return earth residuals (from the \code{\link{lm}} fit on \code{bx})\cr
        \code{"deviance"}
            Return the earth \code{\link{lm}} residuals unless the object has a \code{\link{glm}} component,
            in which case return the \code{\link{glm}} deviance residuals.\cr
        \code{"pearson"}\cr
        \code{"working"}\cr
        \code{"response"}\cr
        \code{"partial"}
            Return the corresponding \code{\link{glm}} residuals
           (from the \code{\link{glm}} fit on \code{bx}).
            Can only be used if the earth model has a \code{\link{glm}} component.
  }
  \item{warn}{
        This function gives warnings when the results are not what you may expect.
        Use \code{warn=FALSE} to turn of just these warnings.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
The residual values (will be a matrix for multiple response models).
}
\seealso{
  \code{\link{earth}}\cr
  \code{\link{residuals}}\cr
  \code{\link{resid}} identical to \code{\link{residuals}}
}
\examples{
data(etitanic)
a <- earth(pclass ~ ., data=etitanic, glm=list(family=binomial))
head(resid(a, warn=FALSE))      # earth residuals, a column for each response
head(resid(a, type="earth"))    # same
head(resid(a, type="deviance")) # GLM deviance residuals, a column for each response
}
\keyword{models}

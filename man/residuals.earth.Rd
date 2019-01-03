\name{residuals.earth}
\alias{residuals.earth}
\title{Residuals for an earth model}
\description{
Residuals of an \code{\link{earth}} model.
}
\usage{
\method{residuals}{earth}(object = stop("no 'object' argument"),
          type = NULL, warn = TRUE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{type}{
     One of:\cr
\cr
        \bold{\code{"earth"}} (default) Residuals from the \code{\link{lm}} fit on \code{bx}.\cr
        \bold{\code{"response"}}
            Residuals as above, but for earth-glm models
            return the \code{\link{glm}} response residuals.\cr
        \bold{\code{"standardize"}} Residuals divided by \code{se * sqrt(1 - h_ii)}.
See the \code{standardize} argument of \code{\link{plot.earth}}.\cr
        \bold{\code{"delever"}} Residuals divided by \code{sqrt(1 - h_ii)}.
See the \code{delever} argument of \code{\link{plot.earth}}.\cr
\cr
The following options are for earth-glm models.
They return the GLM residuals (from the \code{\link{glm}} fit on \code{bx}).
See \code{\link{residuals.glm}} for details:
\cr
\cr
        \bold{\code{"deviance"}}\cr
        \bold{\code{"pearson"}}\cr
        \bold{\code{"working"}}\cr
        \bold{\code{"partial"}}\cr
\cr
The following options for earth-glm models are redundant.
They are provided for compatibility with older versions of earth or other functions:
\cr
\cr
        \bold{\code{"glm.response"}} same as \code{"response"}\cr
        \bold{\code{"glm.deviance"}} same as \code{"deviance"}\cr
        \bold{\code{"glm.pearson"}} same as \code{"pearson"}\cr
        \bold{\code{"glm.working"}} same as \code{"working"}\cr
        \bold{\code{"glm.partial"}} same as \code{"partial"}
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
earth.mod <- earth(pclass ~ ., data=etitanic, glm=list(family=binomial))
head(resid(earth.mod, warn=FALSE))      # earth residuals, a column for each response
head(resid(earth.mod, type="response")) # GLM response resids, a column for each response
}
\keyword{models}

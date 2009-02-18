\name{print.evimp}
\alias{print.evimp}
\title{Print an evimp object}
\description{
    Print an \code{\link{evimp}} object.
    Just a simple utility that drops the class field when printing.
}
\usage{
\method{print}{evimp}(x = stop("no 'x' arg"), \dots)
}
\arguments{
  \item{x}{
     An \code{\link{evimp}} object.
  }
  \item{\dots}{
       Extra arguments passed to print functions.
  }
}
\seealso{
  \code{\link{earth}},
  \code{\link{evimp}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data=ozone1, degree=2)
evimp(a) # implicitly calls print.evimp
}
\keyword{models}

\name{update.earth}
\alias{update.earth}
\title{Update an "earth" model}
\description{
Update an \code{\link{earth}} model.
}
\usage{
\method{update}{earth}(object = stop("no 'object' arg"),
       formula. = NULL, ponly = FALSE, \dots, evaluate = TRUE)
}
\arguments{
  \item{object}{ The earth object}
  \item{formula.}{
     The \code{formula.} argument is treated like earth's \code{formula} argument.\cr
  }
  \item{ponly}{
     Force pruning only, no forward pass.
     Default is FALSE, meaning \code{update.earth} decides automatically if a forward pass is needed.
     See note below.
  }
  \item{\dots}{
     Arguments passed on to \code{\link{earth}}.
  }
  \item{evaluate}{
     If TRUE (default) evaluate the new call, else return the call.
     Mostly for compatibility with the generic \code{\link{update}}.
  }
}
\note{
If only the following arguments are used, a forward pass
is unnecessary, and \code{update.earth} will perform only the pruning pass.
This is usually much faster for large models.
\preformatted{
     object
     glm
     trace
     nprune
     pmethod
     Get.crit
     Eval.model.subsets
     Print.pruning.pass
     Force.xtx.prune
     Use.beta.cache
}
This automatic determination to do a forward pass can be overridden
with the \code{ponly} argument.
If \code{ponly=TRUE} the forward pass will be skipped and only the pruning pass will be executed.
This is useful for doing a pruning pass with new data, typically with \code{penalty=0}.
A similar use of \code{ponly=TRUE} is to do subset
selection with a different \code{penalty} from that used to build the original model.
With \code{trace=1}, \code{update.earth} will tell you if earth's forward pass was skipped.

If you used \code{keepxy=TRUE} in your original call to \code{earth}, then
\code{update.earth} will use the saved values of \code{x}, \code{y}, etc.,
unless you specify otherwise by arguments to \code{update.earth}.
It can be helpful to set \code{trace=1} to see which \code{x}, etc.
\code{update.earth} uses.
}
\value{
   The value is the same as that returned by \code{\link{earth}}.
   If \code{object} is the only parameter then no changes are made
   --- the returned value will be the same as the original \code{object}.
}
\seealso{
  \code{\link{earth}}
}
\examples{
data(ozone1)
(a <- earth(O3 ~ ., data = ozone1, degree = 2))

# yields:
#    Selected 11 of 21 terms, and 8 of 9 predictors
#    Importance: temp, humidity, ibt, doy, dpg, ibh, vis, wind, vh-unused
#    Number of terms at each degree of interaction: 1 5 5
#    GCV 13.5    RSS 3768    GRSq 0.791    RSq 0.822

update(a, formula = O3 ~ . - temp) # requires forward pass and pruning

# yields:
#    Selected 14 of 21 terms, and 8 of 8 predictors
#    Importance: ibt, humidity, doy, vh, ibh, dpg, wind, vis
#    Number of terms at each degree of interaction: 1 5 8
#    GCV 13.1    RSS 3499    GRSq 0.796    RSq 0.834

update(a, nprune = 8)              # requires only pruning

# yields:
#    Selected 8 of 21 terms, and 6 of 9 predictors
#    Importance: temp, humidity, ibt, doy, dpg, ibh-unused, vis, vh-unused, wind-unused
#    Number of terms at each degree of interaction: 1 5 2
#    GCV 15.1    RSS 4435    GRSq 0.766    RSq 0.79

update(a, penalty=1, ponly=TRUE)   # pruning pass only with a new penalty

# yields:
#    Selected 13 of 21 terms, and 8 of 9 predictors
#    Importance: temp, humidity, ibt, doy, dpg, ibh, vis, wind, vh-unused
#    Number of terms at each degree of interaction: 1 5 7
#    GCV 12.5    RSS 3671    GRSq 0.805    RSq 0.826
}
\keyword{models}

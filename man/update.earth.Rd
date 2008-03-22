\name{update.earth}
\alias{update.earth}
\title{Update an 'earth' model}
\description{
Update an \sQuote{earth} model.
}
\usage{
\method{update}{earth}(object = stop("no 'object' arg"),
       formula. = NULL, ponly = FALSE, \dots, evaluate = TRUE)
}
\arguments{
  \item{object}{}
  \item{formula.}{
     These arguments are the same as \code{\link{earth}}.
     The \code{formula.} argument is treated like earth's \code{formula} argument.\cr
  }
  \item{ponly}{
     Force only pruning, no forward pass.
     Default is FALSE, meaning decide automatically if a forward pass is needed.
     See note below.
  }
  \item{\dots}{
     Arguments passed on to \code{\link{earth}}.
  }
  \item{evaluate}{
     If TRUE (default) evaluate the new call else return the call.
     Mostly for compatibility with the generic \code{\link{update}}.
  }
}
\note{
If only the following arguments are used, a forward pass
is unnecessary, and \code{update.earth} will perform only the pruning pass.
This is usually much faster for large models.
\preformatted{
     trace
     pmethod
     nprune
     Get.crit
     Eval.model.subsets
     Print.pruning.pass
     Force.xtx.prune
     Use.beta.cache
}
This automatic determination to do a forward pass can be overriden
with the \code{ponly} argument.
If \code{ponly=TRUE} the forward pass will be skipped, and only the pruning pass will be executed.
This is useful for doing a pruning pass with new data, typically with \code{penalty=1}.
A similar use of \code{ponly=TRUE} is to do subset 
selection with a different \code{penalty} from that used to build the original model.
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
#    Number of terms at each degree of interaction: 1 5 5
#    GCV: 13.43422   RSS: 3762.234   GRSq: 0.7913151   RSq: 0.8218252

update(a, formula = O3 ~ . - temp) # requires forward pass and pruning

# yields:
#    Selected 15 of 21 terms, and 8 of 8 predictors
#    Number of terms at each degree of interaction: 1 5 9
#    GCV: 13.14588   RSS: 3443.265   GRSq: 0.7957941   RSq: 0.8369311

update(a, nprune = 8)              # requires only pruning

# yields:
#    Selected 8 of 21 terms, and 6 of 9 predictors
#    Number of terms at each degree of interaction: 1 5 2
#    GCV: 15.07944   RSS: 4433.915   GRSq: 0.7657586   RSq: 0.7900152

update(a, penalty=1, ponly=TRUE)   # pruning pass only with a new penalty

# yields:
#    Selected 13 of 21 terms, and 8 of 9 predictors
#    Number of terms at each degree of interaction: 1 5 7
#    GCV: 12.5         RSS: 3675         GRSq: 0.805      RSq: 0.826
}
\keyword{models}

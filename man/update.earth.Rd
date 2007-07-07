\name{update.earth}
\alias{update.earth}
\title{Update an 'earth' model}
\description{
Update an \sQuote{earth} model.
}
\usage{
\method{update}{earth}(object = stop("no 'object' arg"),
       formula. = NULL, \dots, evaluate = TRUE)
}
\arguments{
  \item{object}{}
  \item{formula.}{}
  \item{\dots}{
     These arguments are the same as \code{\link{earth}}.
     The \code{formula.} argument is treated like earth's \code{formula} argument.\cr

  }
  \item{evaluate}{
     If TRUE (default) evaluate the new call else return the call.
  }
}
\note{
If only the following arguments are used, a forward pass
is unnecessary, and \code{update.earth} will perform only the pruning pass.
This is usually much faster for large models.
\preformatted{
     trace
     pmethod
     ppenalty
     nprune
     Get.crit
     Eval.model.subsets
     Print.pruning.pass
     Force.xtx.prune
}
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
#    Selected 14 of 21 terms, and 8 of 9 predictors
#    Number of terms at each degree of interaction: 1 6 7
#    GCV: 13.74837   RSS: 3662.576   GRSq: 0.7864352   RSq: 0.8265448

update(a, formula = O3 ~ . - temp) # requires forward pass and pruning

# yields:
#    Selected 15 of 21 terms, and 8 of 8 predictors
#    Number of terms at each degree of interaction: 1 5 9
#    GCV: 13.11272   RSS: 3434.579   GRSq: 0.7963093   RSq: 0.8373425

update(a, nprune = 8)              # requires only pruning

# yields:
#    Selected 8 of 21 terms, and 6 of 9 predictors
#    Number of terms at each degree of interaction: 1 5 2
#    GCV: 15.8959    RSS: 4673.984   GRSq: 0.7530759   RSq: 0.7786458
}
\keyword{models}

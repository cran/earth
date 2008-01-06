\name{mars.to.earth}
\alias{mars.to.earth}
\title{Convert a 'mars' to an 'earth' object}
\description{
Convert a \sQuote{mars} object to an \sQuote{earth} object
}
\usage{
mars.to.earth(object)
}
\arguments{
  \item{object}{
    A \sQuote{mars} object, created using \code{\link[mda]{mars}} in the \code{mda} package.
  }
}
\value{
  The value is the same format as that returned by \code{\link{earth}} but
  with skeletal versions of \code{rss.per.subset}, \code{gcv.per.subset}, and \code{prune.terms}.\cr

  You can fully initialise these components by calling \code{\link{update.earth}}
  after \code{mars.to.earth}, but if you do this \code{selected.terms} may change.
  However with \code{pmethod="backward"} a change is unlikely ---
  \code{selected.terms} would change only if GCVs are so close that numerical errors
  have an effect.
}
\note{
  Perhaps the most notable difference between 
  earth and mars objects is that \code{mars} returns the
  MARS basis matrix in a field called "\code{x}"
  whereas \code{earth} returns "\code{bx}" with only the selected terms.
  Also, \code{earth} returns "\code{dirs}" rather than "\code{factors}",
  and in \code{earth} this matrix can have entries of value 2 for linear predictors.

  For details of other differences between earth and mars objects,
  see the comments in the source code of \code{mars.to.earth()}.
}
\seealso{
  \code{\link{earth}},
  \code{\link[mda]{mars}}
}
\examples{
\dontrun{
# The example is not run with example("mars.to.earth") because that would
# require making the earth package dependent on the mda package.

library(mda)
a <- mars(trees[,-3], trees[,3])
a <- mars.to.earth(a)
summary(a, digits = 2) # the standard earth functions can now be used

# yields (note the reconstructed call):
#    Call:
#    earth(x = trees[, -3], y = trees[, 3])
#
#    Expression:
#      20
#      +  5.4 * pmax(0,  Girth -     12)
#      -  2.6 * pmax(0,     12 -  Girth)
#      + 0.72 * pmax(0, Height -     76)
#
#    Number of cases: 31
#    Selected 4 of 5 terms, and 2 of 2 predictors
#    Number of terms at each degree of interaction: 1 3 (additive model)
#    GCV: 13     RSS: 251     GRSq: 0.95     RSq: 0.97

}
}
\keyword{models}

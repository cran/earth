\name{mars.to.earth}
\alias{mars.to.earth}
\title{Convert a "mars" object to an "earth" object}
\description{
Convert a \code{\link[mda]{mars}} object to an \code{\link{earth}} object
}
\usage{
mars.to.earth(object)
}
\arguments{
  \item{object}{
    A "mars" object, created using \code{\link[mda]{mars}} in the \code{mda} package.
  }
}
\value{
  The value is the same format as that returned by \code{\link{earth}} but
  with skeletal versions of \code{rss.per.subset},
  \code{gcv.per.subset}, and \code{prune.terms}.\cr

  You can fully initialize these components by calling \code{\link{update.earth}}
  after \code{mars.to.earth}, but if you do this \code{selected.terms} may change.
  However with \code{pmethod="backward"} a change is unlikely ---
  \code{selected.terms} would change only if GCVs are so close that numerical errors
  have an effect.
}
\note{
  Perhaps the most notable difference between
  \code{mars} and \code{earth} objects is that \code{mars} returns the
  MARS basis matrix in a field called "\code{x}"
  whereas \code{earth} returns "\code{bx}" with only the selected terms.
  Also, \code{earth} returns "\code{dirs}" rather than "\code{factors}",
  and in \code{earth} this matrix can have entries of value 2 for linear predictors.

  For details of other differences between \code{mars} and \code{earth} objects,
  see the comments in the source code of \code{mars.to.earth()}.

\bold{Weights}

  Note that the \code{w} argument is is actually silently ignored by \code{mars}.
  The equivalent \code{earth} argument \code{weights} is also not yet supported, 
  and you will get a warning.

  \code{mars} normalizes \code{wp} to (euclidean) length 1;
  \code{earth} normalizes
  \code{wp} to length equal to the number of responses, i.e., the number
  of columns in \code{y}.  This change was made so an all 1s \code{wp}
  (or in fact any all constant \code{wp}) is equivalent to using no \code{wp}.

  If the original call to \code{mars} used the \code{w} or \code{wp} arguments,
  \code{mars.to.earth} will run \code{\link{update.earth}} to force consistency.
  This could modify the model, so a warning is issued.
}
\seealso{
  \code{\link{earth}},
  \code{\link[mda]{mars}}
}
\examples{
if (require(mda)) {
    a <- mars(trees[,-3], trees[,3])
    a <- mars.to.earth(a)
    summary(a, digits = 2) # the standard earth functions can now be used

    # yields (note the reconstructed call):
    #    Call: earth(x=trees[, -3], y=trees[, 3])
    #
    #                    y
    #   (Intercept)  19.76
    #   h(Girth-12)   5.40
    #   h(12-Girth)  -2.56
    #   h(Height-76)  0.72
    #
    #   Selected 4 of 5 terms, and 2 of 2 predictors
    #   Importance: Girth, Height
    #   Number of terms at each degree of interaction: 1 3 (additive model)
    #   GCV 13    RSS 251    GRSq 0.95    RSq 0.97
}
}
\keyword{models}

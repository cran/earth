\name{plot.varmod}
\alias{plot.varmod}
\title{Plot a varmod object}
\description{
Plot a variance model (a \code{varmod} object).

Typically you call this function for a variance model
embedded in an \code{earth} model.
}
\usage{
\method{plot}{varmod}(x = stop("no 'x' arg"), which = 1:4,
  do.par = length(which) > 1, info=FALSE,
  cex = NULL, caption = if(do.par) NULL else "", main = NULL,
  col.line = "lightblue", col.min.sd = col.line, \dots)
}
\arguments{
  \item{x}{
    A \code{varmod} object.
    Typically this is embedded in a parent \code{earth} object,
    and so you invoke this function with \code{plot(earth.mod$varmod)}.
    The \code{varmod.method} argument must have been specified
    when building the \code{earth} model.
  }
  \item{which}{
    Which plots to plot.  Default is 1:4 meaning all.
The term \emph{parent} below refers to the \code{earth} model in which
the \code{varmod} is embedded.\cr
1) fitted vs parent fitted\cr
2) fitted vs parent first predictor\cr
3) residuals vs fitted\cr
4) model selection graph (only when \code{varmod.method="earth"} or \code{"x.earth"}).\cr
  }
  \item{do.par}{
    Call \code{par()} for global settings as appropriate.
    Default is \code{TRUE},
    Set to \code{FALSE} if you want to append figures to an existing plot.
  }
  \item{info}{
    Plot some additional information, including lowess fits in the first two plots.
  }
  \item{cex}{
    Character expansion.
  }
  \item{caption}{
    Default is NULL, meaning automatically generate an overall caption.
  }
  \item{main}{
    Default is NULL, meaning automatically generate a title for each graph.
  }
  \item{col.line}{
    Color of lines in the plots.
    Default is \code{lightblue}.
  }
  \item{col.min.sd}{
    Color of the \code{min.sd} dotted horizontal line.
    Default is \code{col.line}.  Use \code{0} to not plot this line.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
 }
}
\note{
The horizontal pale blue dotted line in the first two plots
shows the value of \code{min.sd}.
See \code{\link{earth}}'s \code{varmod.clamp} argument.
}
\examples{
data(ozone1)

set.seed(1) # optional, for cross validation reproducibility

# note: should really use ncross=30 below but for a quick demo we don't

earth.mod <- earth(O3~temp, data=ozone1, nfold=10, ncross=3, varmod.method="lm")

plot(earth.mod$varmod) # plot the embedded variance model (this calls plot.varmod)
}
\seealso{
  \code{\link{varmod}}
}

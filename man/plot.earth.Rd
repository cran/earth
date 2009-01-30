\name{plot.earth}
\alias{plot.earth}
\title{Plot an "earth" object}
\description{
Plot an \code{\link{earth}} object.
The plot shows model selection, cumulative distribution
of the residuals, residuals versus fitted values, and the residual QQ plot.
}
\usage{
\method{plot}{earth}(x = stop("no 'x' arg"),
     which = 1:4, ycolumn = 1,
     caption = if(do.par) NULL else "",
     col.rsq = "lightblue", col.loess = col.rsq,
     col.qq = col.rsq, col.grid = "grey",
     col.vline = "grey", lty.vline = 3,
     col.legend = 1, col.npreds = 1,
     nresiduals = 1000, cum.grid = "percentages", rlim = c(-1,-1),
     jitter = 0, id.n = 3, labels.id = rownames(residuals(x, warn=FALSE)),
     legend.pos = NULL, do.par = TRUE,
     main = NULL, pch = 1, \dots)
}
\arguments{
  \item{x}{
    An \code{\link{earth}} object.
    This is the only required argument.
    (This argument is called "x" for consistency with the generic \code{\link{plot}}.)
  }
  \item{which}{
    Which plots to plot. Default is \code{1:4}, meaning all.\cr
1)  model selection (GRSq plot)\cr
2)  cumulative distribution of absolute values of residuals\cr
3)  residuals versus fitted values\cr
4)  QQ.
  }
  \item{ycolumn}{
     Specify which column of the response to plot if the model has multiple responses.
     Default is 1.
     This argument does not affect the model selection plot which is always across all responses.
  }
  \item{caption}{
    Overall caption.  The default value is
    \code{if(do.par) NULL else ""}. Values are:\cr
       \code{"string"}  string\cr
       \code{""}  no caption\cr
       \code{NULL} generate a caption from the \code{$call} component of the \code{earth} object.
  }
  \item{}{
    \emph{For all the col arguments, 0 means don't plot the corresponding graph element.}
  }
  \item{col.rsq}{
    Color of RSq line in model selection plot.
    Default is \code{"lightblue"}.
  }
  \item{col.loess}{
    Color of \code{\link{loess}} line in residuals plot.
    Default is \code{col.rsq}.
    Generating the loess line occasionally causes warnings such as
    "Warning: pseudoinverse used".
    To get rid of these warnings, set \code{col.loess=0}
  }
  \item{col.qq}{
    Color of QQ line.
    Default is \code{col.rsq}.
  }
  \item{col.grid}{
    Color of grid lines in cumulative distribution plot.  Default is \code{"grey"}.
  }
  \item{col.vline}{
    Color of vertical line at best model in model plot.  Default is \code{"grey"}.
  }
  \item{lty.vline}{
    Line type of vertical line at best model in model plot.  Default is \code{3}.
  }
  \item{col.legend}{
    Color of legend (inside plot area) of model plot.
    Default is \code{1}, meaning draw a legend.
    Use \code{0} for no legend.
  }
  \item{col.npreds}{
     Color of the "number of predictors" plot within the model plot.
     Default is \code{1}.
     Use \code{0} for no "number of predictors" plot.
  }
  \item{nresiduals}{
    Maximum number of residuals to plot.
    Use \code{-1} for all.
    Default is \code{1000} (not all to reduce over-plotting).
    A systematic sample of size \code{nresiduals} is taken but
    the largest few residuals are always included.
  }
  \item{cum.grid}{
    Specify grid on cumulative distribution graph.
    Values are:\cr
    \code{"none"} no grid on cumulative distribution plot\cr
    \code{"grid"}  add grid\cr
    \code{"percentages"}  (default) add grid and percentage labels to quantile lines.
  }
  \item{rlim}{
    Two element vector \code{c(min,max)} specifying min and max
    values on the y axis of the RSq plot.
    Default is \code{c(-1,-1)}.

    Special value \code{min=-1} means the minimum y axis value
    is the smallest GRSq or RSq value excluding the intercept values.

    Special value \code{max=-1} means the maximum y axis value
    is the largest GRSq or RSq value.
  }
  \item{jitter}{
    Jitter applied to GRSq and RSq values to minimize over-plotting.
    Default is \code{0}, meaning no jitter.
    A typical useful value is 0.01.
  }
  \item{id.n}{
    Number of largest residuals to be labeled.  Default is \code{3}.
  }
  \item{labels.id}{
    Residual names.  Default is \code{rownames(residuals(x))}.
  }
  \item{legend.pos}{
    \code{NULL} (default) means position legend automatically.  Else specify \code{c(x,y)}
    in user coordinates.\cr

     \emph{The following settings are related to \code{par()} and are included so you can override the defaults.}
  }
  \item{do.par}{
    Call \code{par()} for global settings as appropriate.
    Default is \code{TRUE},
    which sets \code{mfrow, mar=c(4,4,2,1), mgp=c(1.6,0.6,0), cex=0.7}.
    Set to \code{FALSE} if you want to append figures to an existing plot.
  }
  \item{main}{
    Title of each plot.
    Default is \code{NULL}, meaning generate figure headings automatically.
  }
  \item{pch}{
    Plot character for QQ and residuals plot. Default is \code{1}.
  }
  \item{\dots}{
    Extra arguments passed to plotting functions.
  }
}
\note{
\bold{Why doesn't \code{plot.earth} print GLM information?}

It's just too much to display.
You can instead call \code{\link{plot}} on the \code{glm.list} in the earth model like this:
\preformatted{
data(etitanic)
a <- earth(survived ~ ., data=etitanic, glm=list(family=binomial))
par(mfrow=c(2,2))
plot(a$glm.list[[1]])
}
\bold{I want to add lines or points to the RSq plot,
and am having trouble getting my axis scaling right. Help?}

Use \code{do.par=FALSE}.
With \code{do.par=FALSE}, the axis scales match the axis labels.
With \code{do.par=TRUE}, \code{plot.earth} restores the
\code{\link{par}} parameters and axis scales to
what they were before calling \code{plot.earth}.
This usually means that the x- and y-axis scales are both 0 to 1.
}
\seealso{
  \code{\link{earth}},
  \code{\link{plot.earth.models}},
  \code{\link{plotmo}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(a, cum.grid = "none")
}
\keyword{models}

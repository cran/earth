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
    which = 1:4, nresponse = 1,
    caption = if(do.par) NULL else "",
    col.grsq = 1, lty.grsq = 1, col.rsq = "lightblue",
    lty.rsq = 5, col.vline = col.grsq, lty.vline = 3,
    col.npreds = if(is.null(x$cv.oof.rsq.tab)) 1 else 0, lty.npreds = 2,
    col.mean.oof.rsq = "palevioletred", col.oof.rsq = "mistyrose2",
    col.oof.vline = col.mean.oof.rsq, col.oof.labs = 0,
    col.pch.max.oof.rsq = 0, col.pch.cv.rsq = 0,
    col.mean.infold.rsq = 0, col.infold.rsq = 0, col.sel.grid = 0,
    ylim = c(-1,-1),
    col.legend = 1, cex.legend = NULL, legend.pos = NULL,
    col.cum.grid = "lightgray", cum.grid = "percentages",
    id.n = 3, labels.id = rownames(residuals(x, warn=FALSE)),
    col.residuals = 1, col.loess = col.rsq, nresiduals = 1000,
    col.qq = col.rsq,
    do.par = TRUE, main = NULL, pch = 1, rlim = NA, col.grid = NA, \dots)
}
\arguments{
  \item{x}{
    An \code{\link{earth}} object.
    This is the only required argument.
    (The argument is called "x" for consistency with the generic \code{\link{plot}}.)
  }
  \item{which}{
    Which plots to plot. Default is \code{1:4}, meaning all.\cr
1) model selection (GRSq versus number of terms)\cr
2) cumulative distribution of absolute values of residuals\cr
3) residuals versus fitted values\cr
4) QQ plot of residuals
  }
  \item{nresponse}{
     Specify which column of the response to plot if the model has multiple responses.
     Default is 1.
     This argument does not affect the Model Selection plot which is always across all responses.\cr
[TODO There is an issue in the handling of \code{nresponse} for multiple level factor responses.
Does \code{nresponse} refer to the column in the observed or predicted response?]
  }
  \item{caption}{
    Overall caption.  The default value is
    \code{if(do.par) NULL else ""}. Values are:\cr
       \code{"string"}  string\cr
       \code{""}  no caption\cr
       \code{NULL} generate a caption automatically.
  }
  \item{col.grsq}{
    Color of GRSq line (in the Model Selection plot).
    Default 1.
     Use \code{0} for no GRsq line.
  }
  \item{lty.grsq}{
    Line type of GRSq line (in the Model Selection plot).
    Default 1.
  }
  \item{col.rsq}{
    Color of RSq line (in the Model Selection plot).
    Default is \code{"lightblue"}.
    Use \code{0} for no Rsq line.
  }
  \item{lty.rsq}{
    Line type of RSq line (in the Model Selection plot).
    Default 5.
  }
  \item{col.vline}{
    Color of the vertical line at selected model (in the Model Selection plot).
    Default is \code{col.grsq}.
    This will be at the maximum GRSq unless \code{pmethod="none"}.
    Use \code{0} for no vertical line.
  }
  \item{lty.vline}{
    Line type of vertical line at selected model (in the Model Selection plot).
    Default is \code{3}.
  }
  \item{col.npreds}{
     Color of the "number of predictors" plot (in the Model Selection plot).
     The default displays the number of predictors unless the \code{oof.rsq}'s are displayed.
     Use \code{0} for no "number of predictors" plot.
  }
  \item{lty.npreds}{
     Line type of the "number of predictors" plot (in the Model Selection plot).
     Default is \code{2}.
  }
  \item{col.mean.oof.rsq}{
    Color of mean out-of-fold RSq for each number of terms (in the Model Selection plot).
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default \code{"palevioletred"}.
    Use \code{0} to not plot this line.
  }
  \item{col.oof.rsq}{
    Color of out-of-fold RSq lines for each fold (in the Model Selection plot).
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default is \code{"mistyrose2"}.
    Use \code{0} to not plot these lines.
    May be a vector of colors, which will be recycled if necessary.
  }
  \item{col.oof.vline}{
    Color of vertical line at the maximum \code{oof.rsq} (in the Model Selection plot).
    Default is \code{col.mean.oof.rsq}.
  }
  \item{col.oof.labs}{
    Color of fold number labels on the \code{oof.rsq} lines.
    Default is \code{0}, no labels.
  }
  \item{col.pch.max.oof.rsq}{
    Color of point plotted on the \code{oof.rsq} line to indicate the
    maximum \code{oof.rsq} for that fold.
    Default \code{0}, point not plotted.
  }
  \item{col.pch.cv.rsq}{
    Color of point plotted on the \code{oof.rsq} line to indicate the \code{cv.rsq}.
    for that fold
    (i.e., it is plotted at the number of terms selected by the in-fold GCV).
    Default \code{0}, point not plotted.
  }
  \item{col.mean.infold.rsq}{
    Color of mean in-fold RSq for each number of terms (in the Model Selection plot).
    Default \code{0}, line not plotted.
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
  }
  \item{col.infold.rsq}{
    Color of in-fold RSq lines for each fold (in the Model Selection plot).
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default is \code{0}, lines not plotted.
  }
  \item{col.sel.grid}{
    Color of grid lines in the Model Selection graph.
    Default is \code{0}, no grid.
    Try something like \code{"lightgray"}, \code{"linen"}, or \code{"seashell"}.
    See also \code{col.cum.grid}, for the grid in the Cumulative Distribution plot.
  }
  \item{ylim}{
    Two element vector \code{c(min,max)} specifying min and max
    values on the y axis in the Model Selection plot.
    This is useful if you want to focus on a specific region of the curve.
    Default is \code{c(-1,-1)}.\cr
    Special value \code{min=-1} means the minimum y axis value
    is the smallest GRSq or RSq value excluding the intercept values.\cr
    Special value \code{max=-1} means the maximum y axis value
    is the largest GRSq or RSq value.
  }
  \item{col.legend}{
    Legend color.
    Default is \code{1}.
    Use \code{0} for no legend.
  }
  \item{cex.legend}{
     Legend cex.
     Default is NULL, meaning choose automatically.
  }
  \item{legend.pos}{
    Legend position.
    Default \code{NULL}, meaning position the legend automatically.
    Else specify \code{c(x,y)} in user coordinates,
    or use \code{"topleft"} etc. as explained in \code{\link{legend}}.
  }
  \item{col.cum.grid}{
    Color of grid lines in the Cumulative Distribution plot.
    Default is \code{"lightgray"}.
    See also \code{col.sel.grid}, for the grid in the Model Selection plot.
  }
  \item{cum.grid}{
    Specify grid type in the Cumulative Distribution plot.
    Values are:\cr
    \code{"none"} no grid on Cumulative Distribution plot\cr
    \code{"grid"}  add grid\cr
    \code{"percentages"}  (default) add grid and percentage labels to quantile lines.
  }
  \item{id.n}{
    The largest \code{id.n} residuals will be labeled in the plot.
    Default is \code{3}.
  }
  \item{labels.id}{
    Residual names.  Default is \code{rownames(residuals(x))}.
    Only used if \code{id.n > 0}.
  }
  \item{col.loess}{
    Color of \code{loess} line in the Residuals plot.
    Default is \code{col.rsq}.
    (Actually \code{\link{lowess}} is used rather than \code{\link{loess}},
    to avoid ugly warnings issued by \code{loess}.
    The arguments is so named for backward compatibility.)
  }
  \item{col.residuals}{
    Color of the residual points in the Residuals and QQ plots.
    Default is \code{1}.
    Can be vectorized, for example \code{col.residuals=iris$Species}.
  }
  \item{nresiduals}{
    Maximum number of residuals to plot.
    Use \code{-1} for all.
    Default is \code{1000} (not all to reduce over-plotting).
    A systematic sample of size \code{nresiduals} is taken but
    the largest few residuals are always included.
    This parameter applies to the Residuals and Normal QQ plots;
    the Cumulative Distribution plot uses all residuals.
  }
  \item{col.qq}{
    Color of QQ line (in the QQ plot).
    Default is \code{col.rsq}.
    Use \code{0} for no QQ line.
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
    Plot character in the QQ and Residuals plot. Default is \code{1}.
  }
  \item{rlim}{
    Deprecated.  Please use \code{ylim} instead.
  }
  \item{col.grid}{
    Deprecated.  Please use \code{col.cum.grid} instead.
  }
  \item{\dots}{
    Extra arguments passed to plotting functions.
  }
}
\note{
For details on interpreting the graphs,
please see the \code{earth} package vignette
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the \code{earth} package}}.

Note that cross-validation data will not be displayed unless
both \code{nfold} and \code{keepxy} were used in the original call to
\code{earth}.
}
\seealso{
  \code{\link{earth}},
  \code{\link{plot.earth.models}},
  \code{\link{plotd}},
  \code{\link[plotmo]{plotmo}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(a)
}
\keyword{models}

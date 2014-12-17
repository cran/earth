\name{plot.earth}
\alias{plot.earth}
\title{Plot an earth object}
\description{
Plot an \code{\link{earth}} object.
The plot shows model selection, cumulative distribution
of the residuals, residuals versus fitted values, and the residual QQ plot.
}
\usage{
\method{plot}{earth}(x = stop("no 'x' arg"), which = 1:4,
    info = FALSE, delever = FALSE, pearson = FALSE, level = 0, versus = NULL,
    nresponse = 1, npoints = 1000, id.n = 3,
    labels.id = rownames(residuals(object, warn=FALSE)),
    center = TRUE, loess.f = .5,
    do.par = length(which) > 1, xlim = NULL, ylim = NULL,
    main = NULL, cex.main = 1.2, caption = if(do.par) NULL else "",
    xlab = NULL, ylab = NULL,
    pch = 20, col.line = "lightblue",
    col.loess = NULL, lwd.loess = NULL, col.cv = "red",
    col.qq = col.line, col.grid = "lightgray", col.points = 1, cex.points = NULL,
    shade.pints = "mistyrose2", shade.cints = "mistyrose4",
    cum.grid = "percentages",
    col.rsq = NA, col.residuals = NA, nresiduals = NA,
    # following are passed to plot.model.selection
    legend.pos = NULL, cex.legend = NULL, col.grsq = 1, col.infold.rsq = 0,
    col.mean.infold.rsq = 0,  col.mean.oof.rsq = "palevioletred",
    col.npreds = if(is.null(object$cv.oof.rsq.tab)) 1 else 0,
    col.oof.labs = 0, col.oof.rsq = "mistyrose2", col.oof.vline = col.mean.oof.rsq,
    col.pch.cv.rsq = 0, col.pch.max.oof.rsq = 0, col.sel.grid = 0,
    col.vline = col.grsq, col.vseg = 0, lty.grsq = 1, lty.npreds = 2,
    lty.rsq = 5, lty.vline = 3, col.legend = NA,
    \dots)
}
\arguments{
  There are numerous arguments, but most are merely for
  embellishing the plot.
  The important arguments are early in the argument list.

  \item{x}{
    An \code{\link{earth}} object.
    This is the only required argument.
    (The argument is called "x" for consistency with the generic \code{\link{plot}}.)
  }
  \item{which}{
    Which plots to plot. Default is \code{1:4}.\cr
1) Model selection (GRSq versus number of terms)\cr
2) Cumulative distribution of abs residuals\cr
3) Residuals versus fitted\cr
4) QQ plot of residuals\cr
5) Abs residuals versus fitted\cr
6) Abs residuals versus log fitted\cr
7) Cube root of the squared residuals versus log fitted\cr
8) Log abs residuals versus log fitted\cr
  }
  \item{info}{
    Print some additional information in the Residuals plot.
    Default is \code{FALSE}.\cr
    (1) Plot the distribution of the residuals along the bottom of the plot.\cr
    (2) Show the Spearman Rank Correlation of the absolute residuals
    with the fitted values.
    (The delevered and/or pearson residuals are used if the corresponding
    arguments are set.)
    Note that this correlation is invariant to monotone transforms of the residuals.\cr
    (3) Only for \code{which=5} or \code{8}.
    Regress the absolute residuals against the fitted values
    and show the regression slope.
    Robust linear regression is used via \code{\link[MASS]{rlm}} in the MASS package.
    % If we assume that the "power of the mean model" holds:\cr
    % \code{sd(residuals) = const * meanresponse ^ lambda}\cr
    % then the slope of this line estimates the exponent
    % (a large sample is necessary).
  }
  \item{delever}{
    Default is \code{FALSE}.
    Divide the residuals by \code{sqrt(1 - h_ii)},
    where \code{h_ii} are the diagonal entries of the hat matrix.
    The hat matrix here is from the linear fit on earth's basis matrix \code{bx}.
    (Note: This is sometimes called \dQuote{(internal) studentization}
    or \dQuote{standardization}, but those terms are
    also used in other senses.)
  }
  \item{pearson}{
    Default is \code{FALSE}.
    Divide each residual by its estimated standard deviation.
    Requires that the model was built with the \code{varmod.method} argument,
    because we need the residual model to get the standard deviations.
    Does not include de-levering unless that argument is set.
    (Note: Some authors use \dQuote{pearson residuals} in a different sense.)
  }
  \item{level}{
    Plot prediction intervals at the given confidence level.
    Default is \code{0} meaning do not plot the bands.
    A typical value is \code{level=.95}.\cr
    Requires that the model was built with the \code{varmod.method} argument.\cr
    The color of the bands is determined by the arguments \code{shade.pints} and
    \code{shade.cints}.
  }
  \item{versus}{
    What do we plot the residuals against?  One of:\cr\cr
\code{NULL} Default. Plot the residuals versus the fitted values
(the log fitted values when \code{which=6} to \code{8}).\cr
\cr
\code{""} Residuals versus the predictors (multiple plots if multiple predictors).\cr
\cr
\code{"*"} Residuals versus the \code{earth} basis functions.\cr
\cr
\code{string} Character vector specifying which predictors to plot against.
e.g. \code{versus=c("wind", "vis")}.
Note: \code{\link{grep}} is used for matching.  Thus \code{"wind"} will match
all variables that have \code{"wind"} in their names. Use \code{"^wind$"}
to match only the variable named \code{"wind"}.
Note that \code{versus=""} above is actually just a special case,
since the regex \code{""} matches anything.\cr\cr
  }
  \item{nresponse}{
     Specify which column of the response to plot if the model has multiple responses.
     Default is 1.
     This argument does not affect the Model Selection plot which is always across all responses.
% [TODO There is an issue in the handling of \code{nresponse} for multiple level factor responses.
% Does \code{nresponse} refer to the column in the observed or predicted response?]
  }
  \item{npoints}{
    Maximum number of points to display in residual and other plots.
    Use \code{-1} for all.
    Default is \code{1000} (not all to reduce over-plotting).
    A systematic sample of size \code{nresiduals} is taken but
    the largest few residuals are always included.
  }
  \item{id.n}{
    The largest \code{id.n} residuals will be labeled in the plot.
    Default is \code{3}.
  }
  \item{labels.id}{
    Residual names.  Default is \code{rownames(residuals(x))}.
    Only used if \code{id.n > 0}.
  }
  \item{center}{
    Default is TRUE, meaning center the horizontal axis in the residuals plot,
    so asymmetry in the residual distribution is more obvious.
  }
  \item{loess.f}{
    Smoothing parameter for the light blue smooth line in the residuals plot.
    Default \code{.5}.
    Lower values make the line bumpier.
    This argument is passed as \code{f} in the internal call to \code{\link{lowess}}.
  }
  \item{do.par}{
    If \code{TRUE}, start a new page and call \code{\link{par}} as appropriate.\cr
    Default is \code{length(which) > 1}, i.e., call \code{par} if
    drawing more than one plot.\cr
    Use \code{FALSE} to use the current graphics settings,
    for example if you want to add plots to an existing  multifigure page.\cr
    The value \code{2} means act like \code{TRUE} but do not restore
    the \code{par} settings to their original state when \code{plot.earth} exits.
    This is useful if you want to add a few more plots on the same page
    after running \code{plot.earth}.
  }
  \item{xlim}{
    Horizontal axis limits.
    The default is NULL, meaning automatic.
  }
  \item{ylim}{
    Vertical axis limits.
    The default is NULL, meaning automatic.\cr
    If \code{which} has more than one element,
    this argument applies only to the Model Selection plot.\cr
    In the model selection plot,
    the special value \code{min=-1} means the minimum y axis value
    is the smallest GRSq or RSq value excluding the intercept values.
    The special value \code{max=-1} means the maximum y axis value
    is the largest GRSq or RSq value.
  }
  \item{main}{
    Heading of each plot.
    Default is \code{NULL}, meaning generate the headings automatically.
  }
  \item{cex.main}{
    Size of plot titles.  Default is \code{1.2}.
    Used only if \code{do.par} is \code{TRUE} (default).
  }
  \item{caption}{
    Overall caption.  The default value is
    \code{if(do.par) NULL else ""}. One of:\cr
       \code{"string"}  string\cr
       \code{""}  no caption\cr
       \code{NULL} generate a caption automatically.
  }
  \item{xlab}{
    Default is \code{NULL}, meaning label the x axis automatically.
  }
  \item{ylab}{
    Default is \code{NULL}, meaning label the y axis automatically.
  }
  \item{pch}{
    Default is \code{20}.
  }
  \item{col.line}{
    Color of lines such as the RSq line in the Model Selection plot.
    Default is \code{"lightblue"}
    Use \code{0} for no RSq line.
  }
  \item{col.loess}{
    Color of \code{loess} line in the residuals plot.
    Default is \code{NULL}, meaning automatic.
    (Actually \code{\link{lowess}} is used rather than \code{\link{loess}},
    to avoid ugly warnings issued by \code{loess}.
    The arguments is so named for backward compatibility.)
  }
  \item{lwd.loess}{
    Line width of the \code{loess} line.
    Default is \code{NULL}, meaning automatic.
  }
  \item{col.cv}{
    Default is \code{"red"}.
    Color of cross validation line in the residuals plot.
    This is the residual of the mean out-fold-predicted value.
  }
  \item{col.qq}{
    Color of QQ line in the QQ plot.
    Default is \code{col.line}.
    Use \code{0} for no QQ line.
  }
  \item{col.grid}{
    Color of grid lines in the various plots.
    Default is \code{"lightgray"}.  Use \code{0} for no grid.
  }
  \item{col.points}{
    Color of plotted points.
    Default is \code{1}
  }
  \item{cex.points}{
    Size of the plotted points.
    Default is \code{NULL}, meaning automatic.
  }
  \item{shade.pints}{
    Only used if the \code{level} argument was used.
    Color of the prediction intervals.
    Default is \code{"mistyrose2"}
  }
  \item{shade.cints}{
    Only used if the \code{level} argument was used.
    Color of the confidence intervals.
    Default is \code{"mistyrose4"}
  }
  \item{cum.grid}{
    Specify grid type in the Cumulative Distribution plot.
    One of:\cr
    \code{"none"} no grid on Cumulative Distribution plot\cr
    \code{"grid"}  add grid\cr
    \code{"percentages"}  (default) add grid and percentage labels to quantile lines.
  }
  \item{col.rsq}{
    Deprecated, please use \code{col.line} instead.
    Default is \code{NA}.
  }
  \item{col.residuals}{
    Deprecated, please use \code{col.points} instead.
    Default is \code{NA}.
  }
  \item{nresiduals}{
    Deprecated, please use \code{npoints} instead.
    Default is \code{NA}.
\cr\cr
    \bold{The following are passed to \code{plot.model.selection}.}
  }
  \item{legend.pos}{
    Legend position.
    Default is \code{NULL}, meaning automatic.
    Else specify \code{c(x,y)} in user coordinates,
    or use \code{"topleft"} etc. as explained in \code{\link{legend}}.
  }
  \item{cex.legend}{
     Legend cex.
     Default is \code{NULL}, meaning automatic.
  }
  \item{col.grsq}{
    Color of GRSq line in the Model Selection plot.
    Default is \code{1}.
    Use \code{0} for no GRSq line.
  }
  \item{col.infold.rsq}{
    Color of in-fold RSq lines for each fold in the Model Selection plot.
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default is \code{0}, lines not plotted.
  }
  \item{col.mean.infold.rsq}{
    Color of mean in-fold RSq for each number of terms in the Model Selection plot.
    Default is \code{0}, line not plotted.
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
  }
  \item{col.mean.oof.rsq}{
    Color of mean out-of-fold RSq for each number of terms in the Model Selection plot.
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default is \code{"palevioletred"}.
    Use \code{0} to not plot this line.
  }
  \item{col.npreds}{
     Color of the "number of predictors" plot in the Model Selection plot.
     The default displays the number of predictors unless the \code{oof.rsq}'s are displayed.
     Use \code{0} for no "number of predictors" plot.
  }
  \item{col.oof.labs}{
    Color of fold number labels on the \code{oof.rsq} lines.
    Default is \code{0}, no labels.
  }
  \item{col.oof.rsq}{
    Color of out-of-fold RSq lines for each fold in the Model Selection plot.
    Applies only if \code{nfold} and \code{keepxy} were used
    in the original call to \code{earth}.
    Default is \code{"mistyrose2"}, a pale pink.
    Use \code{0} to not plot these lines.
    May be a vector of colors, which will be recycled if necessary.
  }
  \item{col.oof.vline}{
    Color of vertical line at the maximum \code{oof.rsq} in the Model Selection plot.
    Default is \code{col.mean.oof.rsq}.
  }
  \item{col.pch.cv.rsq}{
    Color of point plotted on the \code{oof.rsq} line to indicate the \code{cv.rsq}.
    for that fold
    (i.e., it is plotted at the number of terms selected by the in-fold GCV).
    Default is \code{0}, point not plotted.
  }
  \item{col.pch.max.oof.rsq}{
    Color of point plotted on the \code{oof.rsq} line to indicate the
    maximum \code{oof.rsq} for that fold.
    Default is \code{0}, point not plotted.
  }
  \item{col.sel.grid}{
    Color of grid lines in the Model Selection graph.
    Default is \code{0}, no grid.
    Try something like \code{"lightgray"}, \code{"linen"}, or \code{"seashell"}.
    See also \code{col.cum.grid}, for the grid in the Cumulative Distribution plot.
    (TODO Is this argument now subsumed by col.grid?)
  }
  \item{col.vline}{
    Color of the vertical line at selected model in the Model Selection plot.
    Default is \code{col.grsq}.
    This will be at the maximum GRSq unless \code{pmethod="none"}.
    Use \code{0} for no vertical line.
  }
  \item{col.vseg}{
    Default is \code{0}.
    Color of triangular marker at top of vertical line for best GRSq.
  }
  \item{lty.grsq}{
    Line type of GRSq line in the Model Selection plot.
    Default is \code{1}
  }
  \item{lty.npreds}{
    Line type of the "number of predictors" plot in the Model Selection plot.
    Default is \code{2}.
  }
  \item{lty.rsq}{
    Line type of RSq line in the Model Selection plot.
    Default is \code{5}.
  }
  \item{lty.vline}{
    Line type of vertical line at selected model in the Model Selection plot.
    Default is \code{3}.
  }
  \item{col.legend}{
    Legend color.
    Default is \code{1}.
    Use \code{0} for no legend.
  }
  \item{\dots}{
    Unused, but provided for generic/method consistency.
  }
}
\note{
For details on interpreting the graphs,
please see the \code{earth} package vignettes
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}
and
\dQuote{\href{../doc/earth-varmod.pdf}{Variance models in earth}}.

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
earth.mod <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(earth.mod)
}
\keyword{models}

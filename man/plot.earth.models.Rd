\name{plot.earth.models}
\alias{plot.earth.models}
\title{Compare earth models by plotting them.}
\description{
Compare \code{\link{earth}} models by plotting them.
}
\usage{
\method{plot}{earth.models}(x = stop("no 'x' arg"), which = c(1:2),
    caption = "", jitter = 0,
    col.grsq = discrete.plot.cols(length(x)), lty.grsq = 1,
    col.rsq = 0, lty.rsq = 5,
    col.vline = col.grsq, lty.vline = 3,
    col.npreds = 0,  lty.npreds  = 2, col.sel.grid = 0,
    ylim = c(0,1),
    col.legend = 1, cex.legend = NULL, legend.pos = NULL, legend.text = NULL,
    col.cum = NULL, do.par = TRUE,
    main = "Model Comparison", cex.main = 1.1, \dots)
}
\arguments{
  \item{x}{
    A list of one or more \code{\link{earth}} objects,
    or a single \code{\link{earth}} object.
    This is the only required argument.
    (This argument is called 'x' for consistency with the generic \code{\link{plot}}.)
  }
  \item{which}{
    Which plots to plot: 1 model, 2 cumulative distribution of residuals.
    Default is \code{1:2}, meaning both.
  }
  \item{caption}{
    Overall caption. Values:\cr
       \code{"string"} string\cr
       \code{""} (default) no caption\cr
       \code{NULL} generate a caption from the \code{$call} component of the \code{earth} objects.
  }
  \item{jitter}{
    Jitter applied to GRSq and RSq values to minimize over-plotting.
    Default is \code{0}, meaning no jitter.
    A typical useful value is 0.01.
\cr
\cr
    \emph{For the col arguments below, 0 means do not plot the corresponding graph element.
          You can use vectors of colors.}
  }
  \item{col.grsq}{
    Vector of colors for the GRSq plot.
    The default is \code{discrete.plot.cols(length(x))} which is vector
    of distinguishable colors, the first three of which are also distinguishable on
    a monochrome printer.
    You can examine the colors using\cr
    \code{earth:::discrete.plot.cols()}.
  }
  \item{lty.grsq}{
    Line type for the GRSq plot.
    Default is \code{1}.
  }
  \item{col.rsq}{
    Vector of colors for the RSq plot.
    Default is \code{0}, meaning no RSq plot.
  }
  \item{lty.rsq}{
    Line type for the RSq plot.
    Default is \code{5}.
  }
  \item{col.vline}{
    A vertical line is drawn for each object
    to show which model size was chosen for that object.
    The color of the line is \code{col.vline}.
    Default is \code{col.grsq}.
  }
  \item{lty.vline}{
    Line type of vertical lines (a vertical line is drawn to show the selected model for each object).
    Can be a vector.
    Default is \code{3}.
  }
  \item{col.npreds}{
    Vector of colors for the "number of predictors" plot within the model selection plot.
    Default is \code{0}, meaning no "number of predictors" plot.
    The special value \code{NULL} means borrow \code{col.grsq}
    (or \code{col.line} if \code{col.grsq} is \code{NULL}).
  }
  \item{lty.npreds}{
     Line type of the "number of predictors" plot (in the Model Selection plot).
     Default is \code{2}.
  }
  \item{col.sel.grid}{
    Color of grid lines in the Model Selection graph.
    Default is \code{0}, no grid.
    Try something like \code{"lightgray"}, \code{"linen"}, or \code{"seashell"}.
  }
  \item{ylim}{
    Two element vector \code{c(min,max)} specifying min and max
    values on the y axis of the RSq/GRSq plot.
    Default is \code{c(0,1)}.\cr
    The special value min=\code{-1} means the minimum y axis value
    is the smallest GRSq or RSq, excluding the intercept.\cr
    The special value max=\code{-1} means the maximum y axis value
    is the largest GRSq or RSq.
  }
  \item{col.legend}{
    Deprecated.   Use \code{legend.pos=NA} for no legend.
  }
  \item{cex.legend}{
     Legend cex.
     Default is NULL, meaning choose automatically.
  }
  \item{legend.pos}{
    Legend position.
    Default \code{NULL}, meaning position the legend automatically.
    Use NA for no legend.
    Specify \code{c(x,y)} in user coordinates,
    or use \code{"topleft"} etc. as explained in \code{\link{legend}}.
  }
  \item{legend.text}{
    Vector of strings to use as legend text.
    The special value \code{NULL} (default) means generate the legend text automatically
    from \code{call$formula}.
  }
  \item{col.cum}{
    Vector of colors for the cumulative distribution plot.
    The special value \code{NULL} (default) means borrow \code{col.grsq}
    (or \code{col.line} if \code{col.grsq} is \code{NULL}).
\cr
\cr
     \emph{The following settings are related to \code{par()} and are included so you can override the defaults.}
  }
  \item{do.par}{
    Call \code{par()} for global settings as appropriate.
    Default is \code{TRUE},
    which sets \code{mfrow, mar=c(4,4,2,3), mgp=c(1.6,0.6,0), cex=0.7}.
    Set to \code{FALSE} if you want to append figures to an existing plot.
  }
  \item{main}{
    Title of each plot.
    Default is \code{NULL}, meaning generate figure headings automatically.
  }
  \item{cex.main}{
    Size of plot titles.  Default is \code{1.1}.
    Used only if \code{do.par} is \code{TRUE} (default).
  }
  \item{\dots}{
    Unused, but provided for generic/method consistency.
  }
}
\note{
   This function ignores GLM and cross-validation components of the earth model, if any.
}
\seealso{
  \code{\link{earth}},
  \code{\link{plot.earth}},
  \code{\link{plot.earth.models}},
  \code{\link{plotd}},
  \code{\link[plotmo]{plotmo}}
}
\examples{
data(ozone1)
a1 <- earth(O3 ~ .,          data = ozone1, degree = 2)
a2 <- earth(O3 ~ .-wind,     data = ozone1, degree = 2)
a3 <- earth(O3 ~ .-humidity, data = ozone1, degree = 2)
plot.earth.models(list(a1,a2,a3), ylim=c(.65,.85))
}
\keyword{models}

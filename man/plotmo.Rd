\name{plotmo}
\alias{plotmo}
\title{Plot the predicted model response}
\description{
Plot the predicted model response when varying one
or two predictors while holding all the other predictors constant.
}
\usage{
plotmo(object = stop("no 'object' arg"),
       degree1 = TRUE, degree2 = TRUE, ycolumn = 1,
       caption = if(do.par) NULL else "",
       ylim = NULL, clip = TRUE, inverse.func = NULL,
       col.response = 0, pch.response = 1, trace = FALSE, 
       grid.func = median, grid.levels = NULL,
       ndegree1 = 500, lty.degree1 = 1, col.degree1 = 1,
       se = 0, col.shade = "lightblue", col.se = 0, lty.se = 2,
       func = NULL, col.func = "pink", pch.func = 20, nrug = 0,
       type2 = "persp", ngrid = 20,
       col.persp = "lightblue", col.image = grey(0:9/10),
       do.par = TRUE, main = NULL, theta = NA, phi = 30, shade = 0.5,
       ticktype = "simple", xlab = "", ylab = "", cex = NULL, \dots)
}
\arguments{
To start off, look at the arguments
  \code{object},
  \code{degree1},
  \code{degree2},
  \code{ylim},
  \code{clip}, and
  \code{do.par}.

  \item{object}{
        Model object.
        This is the only required argument.
  }
  \item{degree1}{
     Index vector specifying main effect plots to include.
     Default is \code{TRUE}, meaning all.
     Perhaps the easiest way to use this argument (and \code{degree2}) is to
     first plot all figures to see how the figures are numbered, then
     replot, using \code{degree1} to select the figures you want.
  }
  \item{degree2}{
     Index vector specifying interaction plots to include.
     Default is \code{TRUE}, meaning all.
  }
  \item{ycolumn}{
     Specify which column of the response to plot if the model has multiple responses.
     Default is 1.
  }
  \item{caption}{
    Overall caption.  The default value is
    \code{if(do.par) NULL else ""}. Values are:\cr
       \code{"string"}  string\cr
       \code{""}  no caption\cr
       \code{NULL}  generate a caption from \code{object$call} and the response name.
  }
  \item{trace}{
     Set to \code{TRUE} to trace operation.
     Useful for tracking down error messages.
     Default is \code{FALSE}.
  }
  \item{ylim}{Values are:\cr

     \code{NULL} (default) all y axes have same limits
     (where "y" is actually "z" on degree2 plots).
     The limits are the min and max values of y across all (degree1 and degree2) plots.
     If \code{col.response!=0} then the original response is
     included in the min and max calculation.\cr

     \code{NA} each graph has its own y limits.\cr

     \code{c(ymin,ymax)} all graphs have the specified y limits.\cr

  }
  \item{clip}{
     Plot only values that are in the range of the response of the original data.
     Default is \code{TRUE}.
  }
  \item{inverse.func}{
     Function applied to the predicted response before plotting.
     Default is \code{NULL}, meaning do not apply a function.
     For example, you could use \code{inverse.func=exp} if your
     model formula is \code{log(y)~x}.
     Note, however, that is usually unnecessary
     to use \code{inverse.func} for models such as \code{\link{glm}}
     with \code{family=binomial}, because the plotted response is already a probability.
  }
  \item{col.response}{
     Colour of response values (actually, the response sites in degree2 plots).
     Here "response" refers to the original data used to build the model.
     Default is 0, meaning don't plot the response.
     Ignored by \code{persp} plots.
  }
  \item{pch.response}{
     Plot character for \code{col.response}.
     Default is 1.
  }
  \item{grid.func}{
     Function applied to x columns to calculate the values for numeric predictors
     not on the axes.
     (For factor predictors, see the \code{grid.levels} argument below).
     Default is \code{\link{median}}.
     Examples:
\preformatted{plotmo(a, grid.func = mean) # mean instead of median
grid.func <- function(x) quantile(x)[2]   # 25\% quantile
plotmo(a, grid.func = grid.func)
}
  }
  \item{grid.levels}{
     This argument only applies if model x-matrix has factors.
     It has the same purpose as the \code{grid.func} argument above but applies to
     factor predictors, not numeric predictors.
     It is a list specifying which factor level to use for factor predictors
     not on the axes.  If specified, plotmo will print messages reminding you of how the levels are set.
     Default is \code{NULL} meaning use the first level of each factor.
\preformatted{a <- earth(survived ~ pclass+sex+age, data=etitanic, degree=2)
plotmo(a, degree2=0, grid.levels=list(pclass="2nd", sex="male"))
}
\bold{The following arguments are for degree1 (main effect) plots}
  }
  \item{ndegree1}{
     Maximum number of points in each degree1 plot.
     Default is \code{500}.
     Special value \code{-1} means use \code{nrow(x)}.
  }
  \item{lty.degree1}{
     Line type of degree 1 plots.  Default is \code{1}.
  }
  \item{col.degree1}{
     Colour of degree 1 plots.  Default is \code{1}.
  }
  \item{nrug}{
     Number of points in (jittered) rug.
     Default is \code{0} meaning no rug.
     Special value \code{-1} for all i.e. \code{nrow(x)}.
  }
  \item{se}{
     Draw standard error bands at plus and minus \code{se} times the pointwise standard errors.
     Default is \code{0}, meaning no standard error bands.
     A common value is \code{2}.
     The predict method for \code{object} must be callable with \code{se.fit=TRUE}
     (examples are \code{\link{predict.lm}} but not \code{\link{predict.earth}}).
  }
  \item{col.shade}{
     Colour of \code{se} shading.  Default is \code{"lightblue"}.
     Set to \code{0} for no shading.
  }
  \item{col.se}{
     Colour of \code{se} lines.  Default is \code{0} meaning no lines just shading.
  }
  \item{lty.se}{
     Line type of \code{se} lines.  Default is \code{2}.
  }
  \item{func}{
     Superimpose \code{func(x)} if \code{func} is not \code{NULL}.
     Default is \code{NULL}.
     This is useful if you are comparing the model to a known function.
     Note that \code{func} is called with a single argument which
     is a dataframe with columns in the same order as the predictors
     in the \code{formula} or \code{x} used to build the model.
     Set \code{trace=TRUE} to see the column names and first few rows of this dataframe.
  }
  \item{col.func}{
     Colour of \code{func} points.
     Default is \code{"pink"}.
  }
  \item{pch.func}{
     Plot character for \code{func} points.
     Default is 20.\cr

  \bold{The following arguments are for degree2 plots}
  }
  \item{type2}{
     Degree2 plot type.
     One of \code{"persp"} (default), \code{"contour"}, or \code{"image"}.
  }
  \item{col.persp}{
     Colour of \code{\link{persp}} surface. Default is \code{"lightblue"}.
     Set to 0 for no colour.
  }
  \item{col.image}{
     Colours of \code{\link{image}} plot. Default is \code{grey(0:9/10)}.
     The default excludes \code{grey(1)} because that is the "colour" of clipped
     values, see \code{clip}.
  }
  \item{ngrid}{
     Grid side length for degree2 plots. Default is \code{20}.\cr

  \bold{The following settings are related to \code{par()} and are included so you can override the defaults.}
  }
  \item{do.par}{
    Call \code{par()} for global settings as appropriate.
    Default is \code{TRUE}.
    Set to \code{FALSE} if you want to append figures to an existing plot.
  }
  \item{main}{
    Title of each plot.
    Default is \code{NULL}, meaning generate figure headings automatically.
  }
  \item{theta}{
     Rotation argument for \code{\link{persp}}.
     Default is \code{NA}, meaning automatically rotate each graph so the highest corner is furthest away.
     Use \code{theta=-45} for x and y increasing as you move into the paper.
  }
  \item{phi}{
     Phi argument for \code{\link{persp}}. Default is \code{30}.
  }
  \item{shade}{
     Shade arg for \code{\link{persp}}. Default is \code{0.5}.
  }
  \item{ticktype}{
     Ticks on \code{\link{persp}} plot.  One of \code{simple} or \code{detailed}.
     Default is \code{"simple"}.
  }
  \item{xlab}{
     Horizontal axis label on degree1 plots
     (for degree2 plots the abscissa labels are always the predictor names).
     Default is \code{""}, meaning none, which gives more plottable area.
     Set to \code{NULL} to use the predictor names as labels.
     If you use \code{NULL}, you may want to set \code{main=""} to avoid
     redundant labeling.
  }
  \item{ylab}{
     Vertical axis label.
     Default is \code{""}, meaning none, which gives more plottable area.
  }
  \item{cex}{
    Character expansion.
  }
  \item{\dots}{
     Extra arguments passed to plotting functions called by \code{plotmo}.
     Using arguments here may cause warnings which can often be safely ignored.
  }
}
\details{
\code{Plotmo} is a general purpose model plotting function
(but comes with the \code{\link{earth}} package).
It is intended for models with quantitative responses.

\bold{The general idea}

Plotmo plots a degree1 plot by changing one predictor while holding 
all other predictors at a constant value.
For degree2 plots, plotmo changes two variables
while holding all other predictors at a constant value.

The question arises: what should that constant value be?
For numerical predictors it is the median of the predictor values in the input matrix x;
for factor predictors it is the first level of the factor.
You can change those defaults with the \code{grid.func} and \code{grid.levels} arguments.

Each graph shows only a thin slice of the data because most variables are constant.
You should be aware of that when interpreting the graph.

\bold{Details of operation}

Let's say the model \code{object} has three predictors,
\code{x1}, \code{x2}, and \code{x3} (all numeric)
and \code{plotmo} is about to plot the degree1 plot for \code{x2}.
\code{Plotmo} first builds an input matrix
with \code{ndegree1} rows and with column names \code{x1}, \code{x2}, and \code{x3}.
It sets all entries for the \code{x1} column to \code{x1}'s median value (actually,
the value returned by \code{grid.func} applied to \code{x1}).
Likewise for the \code{x3} column.
It sets the \code{x2} column to an equally spaced sequence of values
from \code{min(x2)} to \code{max(x2)}.
Finally, it calls \code{\link{predict}}\code{(type="response")} with the newly created input matrix,
and plots the predicted values against the sequence of \code{x2} values.

Operation is similar for degree2 plots: all columns of the input
matrix for \code{\link{predict}} are set to their medians except for the columns of the two
predictors being plotted.

Note that \code{plotmo} calls \code{\link{predict}} with new data and
\code{type="response"}, whereas \code{\link{termplot}}
calls \code{\link{predict}} with \code{type="terms"}.

\bold{Limitations}

NAs are not allowed.
To prevent confusing error messages from \code{plotmo},
remove NAs before building your model.
(To be safe, you can use \code{na.action=na.fail} when building your model
so you get an error message if you inadvertently have a NA).

Weights are currently ignored, with a warning.

Factor predictors are not supported on degree2 graphs (you will get a reminder message).

To avoid confusing error messages from plotmo,
keep the original formula you used to build the model simple.
By default (i.e. when using \code{get.plotmo.x.default} and \code{get.pairs.default}),
\code{plotmo} parses the input \code{\link{formula}} using \code{\link{gsub}}.
This crude approach is not infallible but works for the common formulae.
It determines which predictors are paired by looking for
forms such as \code{"x1:x2"} or \code{"x1*x2"} in the model formula.

Variable names containing \$ are not supported.
The work around is to build the model using temporary variables
or to use \code{\link{attach}}.

\code{Plotmo} can get confused by predictors in formulae which use indexing, such as \code{x[,3]}.
The symptom is usually a message along the lines\cr
\code{Error in model.frame: invalid type (list) for variable 'x[,3]'}.

A mesage like\cr
\code{Warning in model.frame.default: 'newdata' had 50 rows but variable(s) found have 31 rows}\cr
means that \code{\link{model.frame.default}} cannot find all the variables in the data frame
created by \code{plotmo}.

\bold{Minimum Requirements}

\code{Plotmo} requires the following of the model object.
These requirements are for default operation, which can be changed
as described in the next section.

1) \code{object} must have a \code{\link{predict}} method that supports \code{type=response}.

2) for standard error bands (see the \code{se} argument),
\code{object} must have a \code{\link{predict}}
method that can be called with \code{se.fit=TRUE}.

3) \code{object} must have the following two components
(which are searched for in the order given for each):

\itemize{
  \item \code{$x}, or \code{$call$formula}
  (\code{formula} is required for degree2 plots),
  or \code{$call$x}.

  \item \code{$y}, or \code{$call$formula}, or \code{$call$y}.
}

\bold{Extending plotmo}

\code{Plotmo} calls the following generic functions, all defined in
the file \code{plotmo.R}:

\code{plotmo.prolog}\cr
\code{get.plotmo.x}\cr
\code{get.plotmo.y}\cr
\code{plotmo.predict}\cr
\code{get.singles}\cr
\code{get.pairs}

Thus \code{plotmo} can be extended
by writing new method functions, although the default
functions may suffice for your object's class.
See the source comments for details.

\bold{FAQ}

I want to add lines or points to a plot created by \code{plotmo}.
and am having trouble getting my axis scaling right. Help?

Use \code{do.par=FALSE}.
With \code{do.par=FALSE}, the axis scales match the axis labels.
With \code{do.par=TRUE}, plot.earth restores the
\code{\link{par}} parameters and axis scales to
what they were before calling plot.earth.
This usually means that the x- and y-axis scales are both 0 to 1.
}
\author{
Stephen Milborrow.
Users are encouraged to send feedback --- use
milboATsonicPERIODnet \url{http://www.milbo.users.sonic.net}.
}
\seealso{
  \code{\link{termplot}},
  \code{\link{plot.earth}},
  \code{\link{plot.earth.models}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plotmo(a)

# example with some arguments:
# plotmo(a, caption = "example", ylim = NULL, degree1 = c(1,2,4),
#   degree2 = 4, col.response = 3, clip = FALSE, ticktype = "d", theta = -30)

# examples using functions other than earth:
#
# plotmo(lm(O3 ~ log(temp) + humidity*temp, data=ozone1), se=2)
#
# library(gam)
# data(airquality)
# airquality <- na.omit(airquality)  # plotmo doesn't know how to deal with NAs
# plotmo(gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality))
#
# library(mgcv)
# plotmo(gam(O3 ~ s(doy) + s(humidity,temp), data=ozone1), se=2, ylim=NA)
}
\keyword{models}
\keyword{regression}

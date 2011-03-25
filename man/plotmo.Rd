\name{plotmo}
\alias{plotmo}
\concept{partial dependence plot}
\title{Plot a model's response with varying predictor values}
\description{
Plot a model's response when varying one or two predictors while
holding all the other predictors constant.  A poor man's partial
dependence plot.
}
\usage{
plotmo(object = stop("no 'object' arg"),
       type=NULL, nresponse = NA, clip = TRUE, ylim = NULL,
       degree1 = TRUE, all1=FALSE, degree2 = TRUE, all2=FALSE,
       grid.func = median, grid.levels = NULL,
       col.response = 0, cex.response = 1, pch.response = 1,
       jitter.response=0, inverse.func = NULL,
       xflip=FALSE, yflip=FALSE, swapxy=FALSE,
       trace = FALSE, nrug = 0,
       se = 0, col.shade = "lightblue", col.se = 0, lty.se = 2,
       func = NULL, col.func = "pink", pch.func = 20,
       ngrid1 = 500, lty.degree1 = 1, col.degree1 = 1,
       type2 = "persp", ngrid2 = 20,
       col.image = grey(0:9/10), col.persp = "lightblue",
       theta = NA, phi = 30, shade = 0.5,
       do.par = TRUE, caption = NULL, main = NULL,
       xlab = "", ylab = "", cex = NULL, cex.lab = 1, ...)
}
\arguments{
To start off, look at the arguments
  \code{object},
  \code{type},
  and \code{clip}.

  \item{object}{
     Model object.
  }
  \item{type}{
     Type parameter passed to \code{\link{predict}}.
     For legal values see the \code{predict} method for
     your \code{object}; for example, see
     \code{\link[earth]{predict.earth}} or \code{\link{predict.rpart}}.
     By default,
     \code{plotmo} tries to automatically select a suitable value
     (usually \code{"response"}, set \code{trace=TRUE} to see the value).
  }
  \item{nresponse}{
     Which column to use when \code{predict} returns multiple columns.
     This can be a column index or column name
     (which may be abbreviated, partial matching is used).
     Ignored when \code{predict} returns a single column.
  }
  \item{clip}{
     Default is \code{TRUE}, meaning plot only predicted values that
     are in the range of the response of the original data.
     Use \code{FALSE} to plot all predicted values.
  }
  \item{ylim}{Three possibilities:\cr

     (i) \code{NULL} (default) all y axes have same limits
     (where ``y'' is actually ``z'' on \code{degree2} plots).
     The limits are the min and max values of the predicted response
     across all plots (after applying \code{clip}).
     If \code{col.response != 0} then the original response is
     also included in the min and max calculation.\cr

     (ii) \code{NA} each graph has its own y limits.\cr

     (iii) \code{c(ymin,ymax)} graphs have the specified y limits.\cr
\cr
  }
  \item{degree1}{
     Index vector specifying main effect plots to include.
     Default is \code{TRUE}, meaning all \code{degree1} plots (the \code{TRUE} gets recycled).
     Use \code{FALSE} (or \code{0} or \code{NA}) for no \code{degree1} plots.\cr
     The easiest way to use this argument (and \code{degree2}) is to
     first use the default (and possibly \code{all1=TRUE})
     to plot all figures.  This shows how the figures are numbered.
     Then replot using \code{degree1} to select the figures you want,
     for example, \code{degree1=c(1,3)}.
  }
  \item{all1}{
     Default is \code{FALSE}.
     Use \code{TRUE} to plot all predictors,
     not just those usually selected by \code{plotmo}.
     See ``What variables get plotted?'' below.
     The \code{all1} argument \emph{increases} the number of plots;
     the \code{degree1} argument \emph{reduces} the number of plots.
  }
  \item{degree2}{
     Index vector specifying interaction plots to include.
     Default is \code{TRUE}, meaning all \code{degree2} plots.
     Use \code{FALSE} for no \code{degree2} plots.
  }
  \item{all2}{
     Default is \code{FALSE}.
     Use \code{TRUE} to plot all pairs of predictors,
     not just those usually selected by \code{plotmo}.
  }
  \item{grid.func}{
     Function applied to columns of the x matrix to fix the values of
     variables not on the axes.
     Default is \code{\link{median}}.
     (This argument is ignored for factors.  The first level of
     factors is used. That can be changed with \code{grid.levels}.)
     Examples:
\preformatted{  plotmo(a, grid.func = mean) # mean instead of median
  grid.func <- function(x) quantile(x)[2] # 25\% quantile
  plotmo(a, grid.func = grid.func)
}
  }
  \item{grid.levels}{
     Default is \code{NULL}.
     Else a list of variables and their fixed value to be used
     when the variable is not on the axis.
     Supersedes \code{grid.func} for variables in the list.
     Names can be abbreviated, partial matching is used.
     Example:
\preformatted{  a <- earth(survived ~ pclass+sex+age, data=etitanic, degree=2)
  plotmo(a, grid.levels=list(sex="m", age=21))
}
     }
  \item{col.response}{
     Color of response points (response sites in \code{degree2} plots).
     Here ``response'' refers to the response in the original data
     used to build the model.
     Default is 0, don't plot the response.
     Can be a vector, for example, \code{col.response=as.numeric(survived)+1}.
  }
  \item{cex.response}{
     Relative size of response points.
     Default is 1.
     Applies only if \code{col.response!=0}.
  }
  \item{pch.response}{
     Plot character for response points.
     Default is 1.
     Applies only if \code{col.response!=0}.
  }
  \item{jitter.response}{
     Amount to jitter the response points
     (passed to \code{\link{jitter}} as the \code{factor} argument).
     Default \code{0}, no jitter.
     A typical useful value is \code{1}.
     Applies only if \code{col.response!=0}.
  }
  \item{inverse.func}{
     Function applied to the predicted response before plotting.
     Default is \code{NULL}, meaning do not apply a function.
     For example, you could use \code{inverse.func=exp} if your
     model formula is \code{log(y)~x}.
  }
  \item{xflip}{
     Default \code{FALSE}.
     Use \code{TRUE} to flip the direction of the x axis.
     This argument (and \code{yflip} and \code{swapxy}) is useful when comparing
     to a plot from another source and you want the axes to be the same.
     (Note that \code{xflip} and \code{yflip} cannot be used on the \code{persp} plots,
     a limitation of the \code{persp} function.)
  }
  \item{yflip}{
     Default \code{FALSE}.
     Use \code{TRUE} to flip the direction of the y axis of the \code{degree2} graphs.
  }
  \item{swapxy}{
     Default \code{FALSE}.
     Use \code{TRUE} to swap the x and y axes on the \code{degree2} graphs.
  }
  \item{trace}{
     Default is \code{FALSE}.
     Use \code{TRUE} to trace operation.
     Use values greater than \code{1} for more detailed tracing.\cr
\cr
\bold{The following arguments are for \code{degree1} (main effect) plots}
  }
  \item{nrug}{
     Number of points in (jittered) rug.
     Default is \code{0}, no rug.
     Special value \code{-1} for all, i.e., \code{nrow(x)}.
     Otherwise a random subset of \code{nrug} points is taken.
  }
  \item{se}{
     Draw standard error bands at plus and minus \code{se} times the pointwise standard errors.
     Default is \code{0}, no standard error bands.
     A typical value would be \code{2}.
     The predict method for \code{object} must support standard errors,
     i.e. be callable with \code{se.fit=TRUE}
     (for example \code{\link{predict.lm}} but not \code{\link{predict.earth}}).
     Example:
\preformatted{  fit <- lm(stack.loss~., stackloss)
  plotmo(fit, se=2, col.response=2, nrug=-1)
}
  }
  \item{col.shade}{
     Color of \code{se} shading.  Default is \code{"lightblue"}.
     Use \code{0} for no shading.
  }
  \item{col.se}{
     Color of \code{se} lines.  Default is \code{0}, no lines just shading.
  }
  \item{lty.se}{
     Line type of \code{se} lines.  Default is \code{2}.
  }
  \item{func}{
     Superimpose \code{func(x)} if \code{func} is not \code{NULL}.
     Default is \code{NULL}.
     This is useful if you are comparing the model to a known function.
     The \code{func} is called with a single argument which
     is a dataframe with columns in the same order as the variables
     in the \code{formula} or \code{x} used to build the model.
     Use \code{trace=TRUE} to see the column names and first few rows of this dataframe.
  }
  \item{col.func}{
     Color of \code{func} points.
     Default is \code{"pink"}.
  }
  \item{pch.func}{
     Plot character for \code{func} points.
     Default is 20.
  }
  \item{ngrid1}{
     Number of points in \code{degree1} plots.
     Default is \code{500}.
  }
  \item{lty.degree1}{
     Line type of \code{degree1} plots.  Default is \code{1}.
  }
  \item{col.degree1}{
     Color of \code{degree1} plots.  Default is \code{1}.\cr
\cr
  \bold{The following arguments are for \code{degree2} plots}
  }
  \item{type2}{
     \code{Degree2} plot type.
     One of \code{"\link{persp}"} (default), \code{"\link{contour}"}, or \code{"\link{image}"}.
  }
  \item{ngrid2}{
     Grid size for \code{degree2} plots (\code{ngrid2 x ngrid2} points are plotted).
     Default is \code{20}.\cr
     Note 1: the default will often be too small for \code{contour} and \code{image} plots.\cr
     Note 2: with large \code{ngrid2} values, \code{persp} plots look better with \code{border=NA}.
  }
  \item{col.image}{
     Colors of \code{\link{image}} plot. Default is \code{grey(0:9/10)}.
     The default excludes \code{grey(1)} (white) because that is
     the color of clipped values, see \code{clip}.
  }
  \item{col.persp}{
     Color of \code{\link{persp}} surface. Default is \code{"lightblue"}.
     Use 0 for no color.
  }
  \item{theta}{
     Rotation parameter for \code{\link{persp}}.
     Default is \code{NA}, meaning automatically rotate each graph
     so the highest corner is furthest away.
     (Use \code{trace=TRUE} to see the calculated value for \code{theta}.)
  }
  \item{phi}{
     Passed to \code{\link{persp}}. Default is \code{30}.
  }
  \item{shade}{
     Passed to \code{\link{persp}}. Default is \code{0.5}.\cr
\cr
  \bold{The following settings are related to \code{par()}.}
  }
  \item{do.par}{
    Default is \code{TRUE}, meaning
    call \code{\link{par}()} as appropriate for settings such as \code{mfrow}.
    Use \code{FALSE} if you don't want \code{plotmo} to start a new page.
  }
  \item{caption}{
    Overall caption.  The default is to automatically create a caption
    from the \code{call} and response name.\cr
  }
  \item{main}{
    A vector of titles, one for each plot.  Will be recycled if necessary.
    The default generates titles automatically.
    See also \code{caption}, for the overall title.
  }
  \item{xlab}{
     Horizontal axis label on \code{degree1} plots
     (for \code{degree2} plots the abscissa labels are always the variable names).
     Default is \code{""}, no label, which gives more plottable area.
     Use the special value \code{NULL} for the variable names as labels.
     If you use \code{NULL}, you may want to use \code{main=""} to avoid
     redundant labeling.
  }
  \item{ylab}{
     Vertical axis label.  Values as for \code{xlab}.
  }
  \item{cex}{
    Character expansion.
  }
  \item{cex.lab}{
    Relative size of axis labels and text. Default \code{1}.
  }
  \item{\dots}{
     Extra arguments are passed on to the \code{degree2} plotting functions.
  }
}
\details{
\code{Plotmo} can be used on a wide variety of regression models even
though it comes with the \code{\link{earth}} package.
It plots a \code{degree1} plot by calling \code{\link{predict}} to
predict the response when changing one variable while holding all
other variables at their median values.
For \code{degree2} plots, two variables are changed while holding others
at their medians.
For factors, the first level is used instead of the median.
You can change this value with the \code{grid.func} and
\code{grid.levels} arguments.

Each graph shows only a thin slice of the data because most variables
are fixed.  Please be aware of that when interpreting the graph ---
over-interpretation is a temptation.

\bold{Limitations}

Observation weights are currently ignored, with a warning.

NAs are not supported.
To prevent confusing error messages from functions called by \code{plotmo},
it is safest to remove NAs before building your model.
(However, \code{\link[rpart]{rpart}} models are treated specially by
\code{plotmo}, actually allowing NAs so you can use \code{plotmo} with
the default arguments for \code{rpart}.)

Keep the variable names in the original model formula simple.
Use temporary variables or \code{\link{attach}} rather than using
\code{$} and similar in names in formulas.

% By default (i.e., when using \code{get.plotmo.x.default} and \code{get.pairs.default},
% but not when \code{plotmo} is used on an \code{earth} model),
% \code{plotmo} parses the input \code{\link{formula}} using \code{\link{gsub}}.
% This crude approach is not infallible but works for the common formulas.
% It determines which predictors are paired by looking for
% forms such as \code{"x1:x2"} or \code{"x1*x2"} in the model formula.

\bold{What variables get plotted?}

\code{Plotmo} invokes object-specific methods to select which variables to plot.
This selection may leave out variables you want to see ---
in that case use \code{all1=TRUE} or \code{all2=TRUE}.
The set of variables plotted for some common classes is as follows:
\describe{
\item{\code{o earth}}{\code{degree1}: variables in additive (non interaction) terms\cr
\code{degree2}: variables appearing together in interaction terms.}

\item{\code{o rpart}}{\code{degree1}: variables used in the tree\cr
\code{degree2}: parent-child pairs.}

\item{\code{o randomForest}}{\code{degree1}: all variables\cr
\code{degree2}: pairs of the four most important variables (ranked
on the first column of \code{object$importance)}.}

\item{\code{o lm, glm, lda, etc.}}{\code{degree1}: all variables\cr
\code{degree2}: variables in formula terms
like \code{x1*x2} and \code{x1:x2}.}
}

\bold{Using \code{plotmo} on other objects}

Here are some examples which illustrate \code{plotmo} on
non-\code{earth} objects.
\preformatted{
    # use a small set of variables for illustration
    data(ozone1)
    oz <- ozone1[,c("O3", "humidity", "temp", "ibt")]

    lm.fit <- lm(O3 ~ humidity + temp*ibt, data=oz)     # linear model
    plotmo(lm.fit, se=2, col.response=2, nrug=-1)

    library(mgcv)                                       # GAM
    gam.fit <- gam(O3 ~ s(humidity) + s(temp) + s(ibt) + s(temp, ibt), data=oz)
    plotmo(gam.fit, se=2, all2=TRUE)

    library(rpart)
    plotmo(rpart(O3 ~ ., data=oz), all2=TRUE)

    library(randomForest)                               # randomForest
    rf.fit <- randomForest(O3~., data=oz)
    plotmo(rf.fit) # surfaces are smoother than the rpart model above
    partialPlot(rf.fit, oz, temp) # compare partial dependence plot

    library(MASS)                                       # QDA
    lcush <- data.frame(Type=as.numeric(Cushings$Type), log(Cushings[,1:2]))
    lcush <- lcush[1:21,]
    qda.fit <- qda(Type~., data=lcush)
    plotmo(qda.fit, type="class", all2=TRUE,
       type2="contour", ngrid2=100, nlevels=2, drawlabels=FALSE,
       col.response=as.numeric(lcush$Type)+1,
       pch.response=as.character(lcush$Type))

}
\code{Plotmo} has to make some assumptions about the model object.  If
the model function did not save the call or data with the object in a
standard fashion, \code{plotmo} cannot proceed and will issue an error.
Object-specific methods can be usually written to deal with such
issues, see the next section.

\bold{Extending plotmo}

\code{Plotmo} calls the following class methods, defined in \code{plotmo.methods.R}:
\describe{
\item{\code{plotmo.prolog}}{ called before plotting begins}
\item{\code{get.singles}}{ get the vector of variables for \code{degree2} plots}
\item{\code{get.pairs}}{ get the list of pairs for \code{degree2} plots}
\item{\code{get.plotmo.x}}{ get the model matrix \code{x}}
\item{\code{get.plotmo.y}}{ get the model response \code{y}}
\item{\code{plotmo.predict}}{ called for each sub-plot}
}
Thus \code{plotmo} can be extended by writing new methods, although
the defaults suffice for many objects.  See the source code
for details.

\bold{Alternatives}

An alternative approach is to use partial dependence plots.
(See the book by Hastie, Tibshirani, and
Friedman 2nd ed. section 10.13.2.)
\code{Plotmo} sets the ``other'' variables to their median
value, whereas in a partial dependence plot at each plotted point the
effect of the other variables is averaged.

There appears to be no general-purpose R function for partial
dependence plots similar to \code{plotmo}.
Averaging over the sample at every point is a slow process
unless the effect of the averaging can be
determined without actually doing the calculation.
That is not the case for most models, but it is for trees.
See \code{\link[randomForest]{partialPlot}} in the
\code{\link[randomForest]{randomForest}} package.

\code{\link[=termplot]{Termplot}} is effective where applicable, but
it can be used only on models with a \code{predict} method that
supports \code{type="terms"}, and it does not generate
\code{degree2} plots.

\bold{Common error messages}

\code{Warning in model.frame.default: 'newdata' had 50 rows but vars have 31 rows}

This message usually means that \code{\link{model.frame}} cannot
find all the variables in the data frame created by \code{plotmo}.
Make sure the variables you used to build the model are still
available when you call \code{plotmo}.  Try also simplifying
the formula used to create the model.

\code{Error in model.frame: invalid type (list) for variable 'x[,3]'}

\code{Plotmo} can get confused by variables in formulas
which use indexing, such as \code{x[,3]}.
The symptom is usually a message similar to the above.

\code{Error in match.arg(type): 'arg' should be one of ...}

The message is probably issued by the \code{predict} method for your \code{object}.
Set \code{type} to an appropriate value.

\code{Error: get.plotmo.x.default cannot get the x matrix}

This and similar messages mean that \code{plotmo} cannot get the data
it needs from the model \code{object}.
Typically this means that class methods need to be written for the
object, see ``Extending plotmo'' above (and contact the author).

\bold{FAQ}

\emph{I want to add lines or points to a plot created by \code{plotmo}.
and am having trouble getting my axis scaling right. Help?}

Use \code{do.par=FALSE}.
With \code{do.par=FALSE}, the axis scales match the axis labels.
With \code{do.par=TRUE}, \code{plotmo} restores the
\code{\link{par}} parameters and axis scales to
their values before \code{plotmo} was called.

\emph{The \code{persp} display is very jagged. How can I change that?}

Try using \code{clip=FALSE}.
The jaggedness is probably an artifact of the way  \code{\link{persp}} works.
You can also try increasing \code{ngrid2}.
}
\seealso{
  \code{www.milbo.org/rpart-plot/prp.pdf} has a section on \code{plotmo}\cr\cr
  \code{\link{plot.earth}},
  \code{\link{plot.earth.models}},
  \code{\link{plotd}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plotmo(a)
}
\keyword{models}
\keyword{regression}

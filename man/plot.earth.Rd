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
     col.rsq = "lightblue", col.loess = col.rsq,
     col.qq = col.rsq, col.grid = "grey",
     col.vline = "grey", lty.vline = 3,
     col.legend = 1, col.npreds = 1,
     nresiduals = 1000, cum.grid = "percentages", rlim = c(-1,-1),
     id.n = 3, labels.id = rownames(residuals(x, warn=FALSE)),
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
4)  QQ plot of residuals
  }
  \item{nresponse}{
     Specify which column of the response to plot if the model has multiple responses.
     Default is 1.
     This argument does not affect the model selection plot which is always across all responses.\cr
TODO there is an issue in the handling of \code{nresponse} for multiple
level factor responses.
Does \code{nresponse} refer to the column in the observed or predicted response?
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
  \item{id.n}{
    Number of largest residuals to be labeled.  Default is \code{3}.
  }
  \item{labels.id}{
    Residual names.  Default is \code{rownames(residuals(x))}.
  }
  \item{legend.pos}{
    \code{NULL} (default) means position legend automatically.
    Else specify \code{c(x,y)} in user coordinates,
    or use \code{"topleft"} etc. as explained in \code{\link{legend}}.\cr

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
\bold{Interpreting the \code{plot.earth} graphs}

For concreteness, the description below is based
on the graphs plotted by \code{example(plot.earth)}.
The graphs plotted by \code{plot.earth}, apart from the Model Selection graph, are
standard tools used in residual analysis and more information
can be found in most linear regression textbooks.

One should be wary of over-interpretation of the graphs,
since the residuals are measured on the training data rather than on new data.
In linear models that is usually not an issue, but for flexible models
like MARS the residuals measured on the training data give an
optimistic view of the model's predictive ability.

\bold{Nomenclature}.
% The \emph{response value} is the observed response value.
% The \emph{fitted value} is a synonym for the \emph{predicted value} on the training data.
The \emph{residuals} are the differences between the values predicted by the model
and the corresponding response values.
In this help page the residuals are all measured on the training data.
The \emph{residual sum of squares} (RSS) is the sum of the squared values of the residuals.
\emph{R-Squared} (RSq, also called the \emph{coefficient of determination}
or \emph{fraction of variance explained})
is a normalized form of the RSS,
and, depending on the model, varies from 0
(a model that always predicts the same value i.e. the mean observed response value)
to 1 (a model that perfectly predicts the responses in the training data).
The \emph{Generalized Cross Validation} (GCV)
is a form of the RSS penalized by the effective number of model parameters
(and divided by the number of observations).
More details can be found in the FAQ section of the \code{\link{earth}} help page.
The \emph{GRSq} normalizes the GCV in the same way that the RSq normalizes
the RSS.
The GCV and GRSq are measures of the generalization ability of the model,
i.e., how well the model would predict using data not in the training set.
There is some arbitrariness in their values since the effective
number of model parameters is a just an estimate in MARS models.

In the example \bold{Model Selection} graph,
the RSq and GRSq run together at first, but diverge as the number of terms increases.
This is typical behavior, and what we are seeing is an increased penalty
being applied to the GCV as the number of model parameters increases.
The vertical gray dotted line is positioned at the maximum GRSq
and indicates that the best model has 11 terms and uses all 8 predictors
(the number of predictors is shown by the black dotted line).
The graph also shows the number of predictors and terms we would need
if we were prepared to accept a lower GRSq (use the earth parameter \code{nprune}
to trim the model).

The \bold{Cumulative Distribution} graph
shows the cumulative distribution of the absolute values of residuals.
What we would ideally like to see is a graph that starts at 0 and shoots up quickly to 1.
In the example graph, the median absolute residual is about 2.2
(look at the vertical gray line for 50\%).
We see that 95\% of the absolute values of residuals are less than about 7.1
(look at the vertical gray line for 95\%).
So in the training data,
95\% of the time the predicted value is within 7.1 units of the observed value.

The \bold{Residuals vs Fitted} graph
shows the residual for each value of the predicted response.
By comparing the scales of the axes one can get an immediate
idea of the size of the residuals relative to the predicted values.

Ideally the residuals should show constant variance
i.e. the residuals should remain evenly spread out as the fitted values increase.
However, in the example graph we see heteroscedasity --- the residuals spread out in a "<" shape.
There is a decrease in the accuracy of the predictions as the predicted value increases.
% This indicates that variable selection in the MARS pruning pass
% will be biased towards large fitted values, i.e.,
% rather than treating the full range of the response fairly, the pruning pass
% in our example model will favor MARS terms that better fit large response values.
% (Variable selection is done using linear fits on the MARS basis functions.
% Variable selection is based on comparing GCVs and GCVs are derived from RSSs.
% In our example model, the bigger residuals at larger fitted values
% imply that larger fitted values have more ability to affect the RSS.)
To reduce the heteroscedasity,
we could refit the model after performing a transform on the response.
A \code{log} transform, for instance, would even out the residuals:\preformatted{
    a1 <- earth(log(O3) ~ ., data = ozone1, degree = 2)
    plot(a1)
}
Transforming the data may cause other problems,
such as mismatches to a known underlying physical model or
difficulties in interpretation,
so it's best to consult (or become) an expert on the type of data being modeled
(in this case, ozone pollution data --- an expert may say that taking the log is meaningless).

The pale blue line is a \code{\link{loess}} fit.
(Readers not familar with \code{loess} fits can think of them as fancy moving averages.)
In this instance it shows that the mean residual is more or less constant
except at low fitted values.
The end effect is possibly due to failure of the model in that region because
of smaller residuals, but cause and effect get tangled here.
Compare the residuals of the earth model to the linear model, and notice
how the pattern of residuals show that the earth model is more successful
at modeling non-linearities in the data:
\preformatted{
    a2 <- lm(O3 ~ ., data = ozone1)
    plot(a2, which=1)
}
One should always eyeball the residuals themselves rather than blindly
trusting the \code{loess} fit, which is itself an approximation.
However, in our example earth model the \code{loess} line appears reliable.
Also, heteroscedasity is not as big an issue with earth as it is with linear models,
where homoscedasity of the residuals is assumed for inference (such
as calculation of p-values).

Cases 192, 193, and 226 have the largest residuals and
fall suspiciously into a separate cluster.
(If overplotting makes the labels hard to read,
reduce the number of labels with the \code{id.n} argument of \code{plot.earth}.)
As a general rule, it is worthwhile investigating cases with large residuals.
Perhaps they should be excluded when building the model.
On the other hand, it is possible that they reveal something important about the data
that could warrant changes to the model.
In our example it is also worthwhile looking at cases
with \emph{small} residuals because of non-linearity in that region.
To see the example input matrix ordered on the magnitude of the residuals,
use \code{ozone1[order(abs(a$residuals)),]}.

Sometimes groups of residuals appear in a series of straight lines with slopes of -1.
This effect is slightly visible in the example graph.
These lines usually do not indicate a problem.
They are are formed when a set of plotted points has the same observed
value, commonly due to discretization in the measurement of the
observed response (e.g. caused by rounding to the nearest inch).

The \bold{Normal Q-Q} graph
compares the distribution of the residuals to a normal distribution.
If the residuals are distributed normally they will lie on the line.
Following R convention, the abscissa is the normal axis
and the ordinate is the residual axis;
some popular books have it the other way round.
In our example,
we see divergence from normality in the left tail ---
the left tail of the distribution is fatter than that of a normal distribution.
Once again, we see that cases 192, 193, and 226 have the largest residuals.

\bold{About \code{plot.earth} and earth-glm models}

"Earth-glm" models are models created with a \code{glm} argument to \code{earth}.
In earth-glm models, much of the analysis in the above section does not apply because
in these models the residuals are not assumed to have a normal distribution.

Note that the residuals plotted by \code{plot.earth}
are residuals from earth's call to \code{\link{lm}} after the pruning pass,
not \code{\link{glm}} residuals
That is, \code{plot.earth} ignores the \code{glm} part of the model, if any.

For earth-glm models, \code{\link{plotd}} can be useful.

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
  \code{\link{plotd}},
  \code{\link[plotmo]{plotmo}}
}
\examples{
data(ozone1)
a <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(a)
}
\keyword{models}

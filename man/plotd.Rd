\name{plotd}
\alias{plotd}
\title{Plot the distribution of predictions for each class}
\description{
Draw a plot of the distribution of the predicted values for each class.
}
\usage{
plotd(object, hist = FALSE, type="response", thresh = .5, xlim = NULL, 
    main=NULL, xlab="Predicted Value", ylab=if(hist) "Count" else "Density",
    lty = 1, col = c("grey60", 1, "lightblue", "brown", "pink", 2, 3, 4),
    borders = NULL, labels=FALSE, zero.line=FALSE,
    legend = TRUE, legend.names = NULL, legend.pos = NULL,
    legend.cex = .8, legend.bg = "white", legend.extra = FALSE,
    vline.col = 0, vline.thresh = thresh, vline.lty = 1,
    ...)
}
\arguments{
  To start off, look at the arguments \code{object}, \code{hist}, and \code{type}.
  \item{object}{
    model object.
    Typically a model which predicts a class or a class probability.
  }
  \item{hist}{
    \code{FALSE} (default) to call \code{\link{density}} internally\cr
    \code{TRUE} to call \code{\link{hist}} internally
  }
  \item{type}{
    type parameter passed on to \code{\link{predict}}.
    Default is \code{"response"}.
    See the predict method for your object to see other values;
    for example see \code{\link{predict.earth}}.
    Typically you would set \code{hist=TRUE} if \code{type="class"}.
  }
  \item{thresh}{
    Passed to \code{\link{predict.earth}} (only used if \code{type="class"}).
    Default is \code{.5}.
    Also used as a default for the \code{vline.thresh} argument (see below).
    Also used for multiple column numeric responses.
  }
  \item{xlim}{
    limits of the x axis.
    The default \code{NULL} means determine these limits automatically,
    else specify \code{c(xmin,xmax)}.
  }
  \item{main}{
    main graph title. Values:\cr
       \code{"string"} string\cr
       \code{""} no title\cr
       \code{NULL} (default) generate a title from the \code{$call} component of the model object.
  }
  \item{xlab}{
    x axis label.
    Default is \code{"Predicted Value"}.
  }
  \item{ylab}{
    y axis label.
    Default is \code{if(hist) "Count" else "Density"}.
  }
  \item{lty}{
    per class line types for the plotted lines.
    Default is 1 for all lines.
  }
  \item{col}{
    per class graph colors.
    The first few colors of the default are intended to be
    easily distinguishable on both color displays and monochrome printers.
  }
  \item{borders}{
    per class border colors for \code{\link{hist}} plot.
    Only used if \code{hist=TRUE}.
    The default \code{NULL} means use the values in the \code{col} argument.
  }
  \item{labels}{
    TRUE to draw counts or rounded densities on \code{\link{hist}} plot.
    Default is \code{FALSE}.
  }
  \item{zero.line}{
    parameter for \code{\link{plot.density}}.
    Only used if \code{hist=FALSE}.
    Default is \code{FALSE}.
  }
  \item{legend}{
    \code{TRUE} (default) to draw a legend
  }
  \item{legend.names}{
    class names in legend.
    The default \code{NULL} means determine these automatically.
  }
  \item{legend.pos}{
    position of the legend.
    The default \code{NULL} means try to position the legend automatically,
    else specify \code{c(x,y)}.
  }
  \item{legend.cex}{
    \code{cex} for \code{\link{legend}}.
    Default is \code{.8}.
    Make this smaller for a smaller legend.
  }
  \item{legend.bg}{
   \code{bg} color for \code{\link{legend}}.
   Default is \code{"white"}.
  }
  \item{legend.extra}{
    show number of occurrences of each class.
    Default is \code{FALSE}.
  }
  \item{vline.col}{
    color of vertical line.
    Default is 0, meaning no vertical line.
    The vertical line is intended to indicate class separation.
  }
  \item{vline.thresh}{
    x position of vertical line.  Default is \code{thresh}.
  }
  \item{vline.lty}{
    line type of vertical line.
    Default is 1.
  }
  \item{\dots}{
    extra arguments passed to \code{\link{density}} or \code{\link{hist}}.
    Arguments for \code{\link{hist}} which affect only plotting 
    are not allowed (you will get an error message).
  }
}
\note{
  This function estimates distributions with the
  \code{\link{density}} and \code{\link{hist}} functions,
  and also calls \code{\link{plot.density}} and \code{\link{plot.histogram}}.
  For an overview of the issues involved, see Venables and Ripley MASS section 5.6.

  This function is included in the \code{\link{earth}} package
  but can also be used with other models.  Example: \preformatted{
      data(etitanic)
      glm.model <- glm(survived ~ ., data=etitanic, family=binomial)
      plotd(glm.model)
  } or even: \preformatted{
      lm.model <- lm(survived ~ ., data=etitanic)
      plotd(lm.model)
  }
  \bold{How plotd works}

  This function calls \code{\link{predict}}
  with the data originally used to build
  the model, and with the \code{type} specified above.
  It then separates the predicted values into classes,
  where the class for each predicted value
  is determined by the class of the observed response.
  Finally, it calls \code{\link{density}} (or \code{\link{hist}}
  if \code{hist=TRUE}) for each class-specific set of values,
  and plots the results.

  \bold{Partitioning the response into classes}

  Considerable effort is made to partition the predicted response 
  into classes in a sensible way.
  The general idea is that the observed response is the reference
  used to partition the predicted response into classes.
  This section gives details.

  If the \emph{observed response is numeric or logical}, 
  \code{plotd} separates the predicted response values into two groups corresponding to
  \emph{zero-or-less} and \emph{greater-than-zero} observed responses.
  Thus \code{plotd} does what you would expect when the observed response is a
  \code{\link{logical}}, or if your model predicts \code{-1} and \code{1}.

  If the \emph{observed response is a factor and the predicted response 
  is a single numeric column}, \code{plotd} partitions the predicted 
  response into two classes --- the first level
  of the factor (in the observed response) versus all others.

  If the \emph{observed response is a factor and the predicted response 
  is a single \emph{factor} column}, \code{plotd} partitions the predicted
  response into the levels of the observed response.
  This occurs for earth models using predictions with \code{type="class"}.

  If the \emph{observed response is a factor and the predicted response 
  has multiple numeric columns}, 
  \code{plotd} assumes each column of the predicted response corresponds
  to a dummy variable for a factor level.
  This type of response occurs for earth models where the observed response 
  is a factor with two or more levels 
  (see the "Factors" section on the \code{\link{earth}} help page).
  The number of classes is taken to be the number of factor levels,
  or, equivalently, the number of columns in the predicted response.
  The predicted values for each factor are taken to be 
  those elements of each column of the predicted response
  which match the column's factor level in the observed response.
  The remaining elements in the column are not plotted.

  If the \emph{observed response has multiple numeric (or logical) columns 
  and the predicted response has the same number of numeric (or logical) columns}, 
  \code{plotd} assumes each column of the predicted response corresponds
  to a column in the observed response.
  This type of response occurs for \code{\link{lm}} or earth models where the observed response 
  has multiple numeric or logical columns.
  The number of classes is taken to be the number of columns in the response.
  In this situation, \code{plotd} selects the elements in a column of a predicted response
  matching elements in the corresponding column in the observed response
  which are greater than \code{thresh}.
  The remaining elements in the column are not plotted.
}
\seealso{
  \code{\link{earth}}, \code{\link{plot.earth}}, \code{\link{plotmo}}\cr
  \code{\link{density}}, \code{\link{plot.density}}\cr
  \code{\link{hist}}, \code{\link{plot.histogram}}

}
\examples{
old.par <- par(no.readonly=TRUE); par(mfrow=c(2,2)); 
par(mar=c(4, 3, 1.7, 0.5)); par(mgp=c(1.6, 0.6, 0))
data(etitanic)
a <- earth(survived ~ ., data=etitanic, degree=2, glm=list(family=binomial))

plotd(a, main="plotd default")

plotd(a, main="histogram", hist=TRUE, legend.pos=c(.25,220))

plotd(a, main="histogram, type=class", hist=TRUE, type="class",
      legend.pos=c(.14,500), legend.extra=TRUE)

par(old.par)
}
\keyword{models}

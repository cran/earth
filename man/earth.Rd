\name{earth}
\alias{earth}
\alias{earth.default}
\alias{earth.formula}
\concept{regression}
\concept{mars}
\concept{Friedman}
\title{Earth: Multivariate Adaptive Regression Splines}
\description{
Build a regression model
using the techniques in Friedman's papers \sQuote{Multivariate Adaptive Regression Splines}
and \sQuote{Fast MARS}.
}
\usage{
\method{earth}{formula}(formula, data, \dots)

\method{earth}{default}(x = stop("no 'x' arg"), y = stop("no 'y' arg"),
      weights = NULL, subset = NULL, na.action = na.fail,
      penalty = if(degree > 1) 3 else 2, trace = 0, keepxy = FALSE,
      degree = 1, nk = max(21, 2 * NCOL(x) + 1),
      thresh = 0.001, minspan = 0, newvar.penalty = 0,
      fast.k = 20, fast.beta = 1, fast.h = NULL,
      pmethod = "backward", ppenalty = penalty, nprune = NULL,
      Object  = NULL, Get.crit = get.gcv,
      Eval.model.subsets = eval.model.subsets,
      Print.pruning.pass = print.pruning.pass,
      Force.xtx.prune = FALSE, \dots)
}
\arguments{
All arguments are optional except \code{formula}, or \code{x} and \code{y}.
NAs are not allowed.\cr

To start off, look at the arguments
\code{formula},
\code{x},
\code{y},
\code{degree},
\code{nk}, and
\code{nprune}.

  \item{formula}{
     Model formula.
  }
  \item{data}{
    Data frame for \code{formula}.
  }
  \item{x}{
     Matrix containing the independent variables.
  }
  \item{y}{
     Vector containing the response variable,
     or in the case of multiple responses,
     a matrix whose columns are the response values for each variable.
     If the \code{y} values are very big or very small then you may get
     better results if you \code{\link{scale}} \code{y} first for
     better numerical stability in the forward pass.
  }
  \item{weights}{
     Weight vector (not yet supported).
  }
  \item{subset}{
     Index vector specifying which cases to use i.e. which rows in \code{x} to use.
     Default is NULL, meaning all.
  }
  \item{na.action}{
     NA action. Default is \code{na.fail}, and only \code{na.fail} is supported.
  }
  \item{penalty}{
    Generalized Cross Validation (GCV) penalty per knot.
    Default is \code{if(degree>1) 3 else 2}.
    A value of 0 penalises only terms, not knots.
    The value -1 is a special case, meaning no penalty, so GCV=RSS/n.
    Theory suggests values in the range of about 2 to 3.
    In practice for big data sets, larger values can be useful to force a smaller rmodel.
    See also \code{ppenalty}.
  }
  \item{trace}{
     Trace earth's execution.  Default is 0.  Values:

     0 none 1 overview 2 forward 3 pruning 4 more pruning 5 \dots
  }
  \item{keepxy}{
     Set to TRUE to retain \code{x}, \code{y}, and  \code{subset} in the returned value.
     Default is FALSE\cr

     \emph{The following arguments are for the forward pass}
  }
  \item{degree}{
     Maximum degree of interaction (Friedman's \eqn{mi}).
     Default is 1, meaning build an additive model.
  }
  \item{nk}{
     Maximum number of model terms before pruning.
     Includes the intercept.
     Default is \code{max(21,2*NCOL(x)+1)}.
     The number of terms created by the forward pass will be
     less than \code{nk} if there are linearly dependent terms
     which must be discarded,
     or if a forward stopping condition is reached before \code{nk} terms.
     See the section below on the forward pass.
  }
  \item{thresh}{
     Forward stepping threshold.
     Default is 0.001.
     This is one of the arguments used to decide when forward stepping
     should terminate.
     See the section below on the forward pass.
  }
  \item{minspan}{
    Minimum distance between knots.
    Set \code{trace>=2} to see the calculated value.  Values:\cr
    \code{<0} add to the internally calculated min span (i.e. decrease span).\cr
    \code{0} (default) use internally calculated min span as per
    Friedman's MARS paper section 3.8 with \eqn{alpha} = 0.05.
    Intended to increase resistance to runs of noise in the input data.\cr
    \code{>0} use instead of the internally calculated min span.
    Thus a value of 1 means consider all knots.
  }
  \item{newvar.penalty}{
     Penalty for adding a new variable in the forward pass
     (Friedman's \eqn{gamma}, equation 74 in the MARS paper).
     This argument can mitigate the effects of collinearity or concurvity
     in the input data.
     Default is 0.
     Useful non-zero values range from about 0.01 to 0.2 --- you will
     need to experiment.
  }
  \item{fast.k}{
     Maximum number of considered parent terms, as
     as described in Friedman's Fast MARS paper section 3.0.
     Default is 20.
     The special value -1 is equivalent to infinity, meaning no Fast MARS.
     Typical values, apart from -1, range
     from about 20 to 5, in steps of 5 or 10.
     In general, with a lower \code{fast.k} (say 5), \code{earth} is faster;
     with a higher \code{fast.k}, or with \code{fast.k} disabled (set to -1), \code{earth} builds
     a better model.
     However it is not unusual to get a slightly better model with a lower \code{fast.k},
     and you may need to experiment.
  }
  \item{fast.beta}{
     Fast MARS ageing coefficient, as described in the
     Fast MARS paper section 3.1.
     Default is 1.
     A value of 0 sometimes gives better results.
  }
  \item{fast.h}{
     Fast MARS \eqn{h}, as described in the Fast MARS paper section 4.0.
     (not yet implemented).\cr

     \emph{The following arguments are for the pruning pass.
     One strategy is to build a large model and then trim it
     back using \code{\link{update.earth}} with the arguments below.
     }
  }
  \item{pmethod}{
     Pruning method.
     Default is \code{"backward"}.
     One of: \code{backward none exhaustive forward seqrep}.
     If \code{y} has multiple columns, then only \code{backward} or \code{none}
     is allowed.
     Pruning can take a while if \code{exhaustive} is chosen and
     the model is big (more than about 30 terms).
     The current version of the \code{\link[leaps]{leaps}} package
     used during pruning does not allow user interrupts
     (i.e. you have to kill your R session to interrupt).
  }
  \item{ppenalty}{
     Like \code{penalty} but for the pruning pass.
     Default is \code{penalty}.
  }
  \item{nprune}{
     Maximum number of terms (including intercept) in the pruned model.
     Default is NULL, meaning all terms.
     Use this to reduce exhaustive search time, or to enforce a maximum model size.\cr

    \emph{The following arguments are for internal or advanced use}
  }
  \item{Object}{
     Earth object to be updated, for use by \code{\link{update.earth}}.
  }
  \item{Get.crit}{
     Criterion function for model selection during pruning.
     By default a function that returns the GCV.
     See the section below on the pruning pass.
  }
  \item{Eval.model.subsets}{
     Function used to evaluate model subsets
     --- see notes in source code.
  }
  \item{Print.pruning.pass}{
     Function used to print pruning pass results
     --- see notes in source code.
  }
  \item{Force.xtx.prune}{
Default is FALSE.
This argument pertains to subset evaluation in the pruning pass.
By default,
if \code{y} has a single column then \code{earth} calls the \code{\link[leaps]{leaps}} routines;
if \code{y} has multiple columns then \code{earth} calls \code{EvalSubsetsUsingXtx}.
The \code{\link[leaps]{leaps}} routines are based on the QR decomposition;
the \code{EvalSubsetsUsingXtx} routine is based on the inverse of X'X.
The \code{\link[leaps]{leaps}} routines are more accurate but do not support multiple responses.
Setting \code{Force.xtx.prune=TRUE} forces use of \code{EvalSubsetsUsingXtx}, even
if \code{y} has a single column.
}
  \item{\dots}{
     \code{earth.formula}: arguments passed to \code{earth.default}.

     \code{earth.default}: unused, but provided for generic/method consistency.
  }
}
\value{
  An object of class \sQuote{earth} which is a list with the components
  listed below.
  \emph{Term} refers to a term created during the forward pass
  (each line of the output from \code{\link{format.earth}} is a term).
  Term number 1 is always the intercept.
  \item{bx}{
     Matrix of basis functions applied to \code{x}.
     Each column corresponds to a selected term.
     Each row corresponds to a row in in the input matrix \code{x},
     after taking \code{subset}.
     See \code{\link{model.matrix.earth}} for an example of \code{bx} handling.
     Example \code{bx}:\preformatted{    (Intercept) h(Girth-12.9) h(12.9-Girth) h(Girth-12.9)*h(...
[1,]          1           0.0           4.6                   0
[2,]          1           0.0           4.3                   0
[3,]          1           0.0           4.1                   0
...}
}
  \item{dirs}{
     Matrix with \eqn{ij}-th element equal to 1 if term
     \eqn{i} has a factor of the form \eqn{x_j > c}, equal to \eqn{-1} if
     term \eqn{i} has a factor of the form \eqn{x_j \le c}, and to 0 if
     \eqn{x_j} is not in term \eqn{i}.
     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Example \code{dirs}:\preformatted{                       Girth Height
(Intercept)                0      0  #no factors in intercept
h(Girth-12.9)              1      0  #2nd term uses Girth
h(12.9-Girth)             -1      0  #3rd term uses Girth
h(Girth-12.9)*h(Height-76) 1      1  #4th term uses Girth and Height
...}
  }
  \item{cuts}{
     Matrix with \eqn{ij}-th element equal to the cut point
     for variable \eqn{j} in term \eqn{i}.
     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Example \code{cuts}:\preformatted{                           Girth Height
(Intercept)                  0.0      0  #intercept, no cuts
h(Girth-12.9)               12.9      0  #2nd term has cut at 12.9
h(12.9-Girth)               12.9      0  #3rd term has cut at 12.9
h(Girth-12.9)*h(Height-76)  12.9     76  #4th term has two cuts
...}
  }
  \item{selected.terms}{
     Vector of term numbers in the best model.
     Can be used as a row index vector into \code{cuts} and \code{dirs}.
     The first element \code{selected.terms[1]} is always 1, the intercept.
  }
  \item{prune.terms}{
     A matrix specifying which terms appear in which subsets.
     The row index of \code{prune.terms} is the model size
     (the model size is the number of terms in the model).
     Each row is a vector of term numbers for the best model of that size.
     An element is 0 if the term is not in the model, thus \code{prune.terms} is a
     lower triangular matrix, with dimensions \code{nprune x nprune}.
     The model selected by the pruning pass
     is at row \code{length(selected.terms)}.
     Example \code{prune.terms}:\preformatted{[1,]    1  0  0  0  0  0  0  #intercept-only model
[2,]    1  2  0  0  0  0  0  #best 2 term model uses terms 1,2.
[3,]    1  2  4  0  0  0  0  #best 3 term model uses terms 1,2,4
[4,]    1  2  9  8  0  0  0
...}
  }
  \item{rss}{
     Residual sum-of-squares (RSS) of the model
     (summed over all responses if \code{y} has multiple columns).
  }
  \item{rsq}{
    \code{1-rss/rss.null}.
     R-Squared of the model
     (calculated over all responses if \code{y} has multiple columns).
     A measure of how well the model fits the training data.
  }
  \item{gcv}{
     Generalised Cross Validation (GCV) of the model
     (summed over all responses if \code{y} has multiple columns).
     The GCV is calculated using \code{ppenalty} as are all GCVs in the returned value.
     For details of the GCV calculation, see
     equation 30 in Friedman's MARS paper and \code{earth:::get.gcv}.
  }
  \item{grsq}{
     \code{1-gcv/gcv.null}.
     An estimate of the predictive power of the model
     (calculated over all responses if \code{y} has multiple columns).
     Unlike \code{rsq}, \code{grsq} can be negative.
     A negative \code{grsq} would indicate
     a severely over parameterised model --- a model that
     would not generalise well
     even though it may be a good fit to the training data.
     Example of a negative \code{grsq}:\cr
     \code{earth(mpg~., data=mtcars, pmethod="none", trace=4)}
  }
  \item{rss.per.response}{
    A vector of the RSS for each response.
    Length is \code{ncol(y)}.
    The \code{rss} component above is  equal to \code{sum(rss.per.response)}.
  }
  \item{rsq.per.response}{
    A vector of the R-Squared for each response.
    Length is \code{ncol(y)}.
  }
  \item{gcv.per.response}{
    A vector of the GCV for each response.
    Length is \code{ncol(y)}.
    The \code{gcv} component above is equal to \code{sum(gcv.per.response)}.
  }
  \item{grsq.per.response}{
    A vector of the GRSq for each response.  A vector of length \code{ncol(y)}.
  }
  \item{rss.per.subset}{
     A vector of the RSS for each model subset generated by the pruning pass.
     Length is \code{nprune}.
     If \code{y} has multiple columns, the RSS is summed over all responses for each subset.
     The null RSS (i.e. the RSS of an intercept only-model) is \code{rss.per.subset[1]}.
     The \code{rss} above is \code{rss.per.subset[length(selected.terms)]}.
  }
  \item{gcv.per.subset}{
     A vector of the GCV for each model in \code{prune.terms}.
     Length is is \code{nprune}.
     If \code{y} has multiple columns, the GCV is summed over all responses for each subset.
     The null GCV (i.e. the GCV of an intercept-only model) is \code{gcv.per.subset[1]}.
     The \code{gcv} above is \code{gcv.per.subset[length(selected.terms)]}.
  }
  \item{fitted.values}{
    Fitted values.
    A matrix with dimensions \code{nrow(y) x ncol(y)}.
  }
  \item{residuals}{
    Residuals.
    A matrix with dimensions \code{nrow(y) x ncol(y)}.
  }
  \item{coefficients}{
    Regression coefficients.
    A matrix with dimensions \code{length(selected.terms) x ncol(y)}.
    Each column holds the least squares coefficients from regressing that
    column of \code{y} on \code{bx}.
    The first row holds the intercept coefficients.
  }
  \item{ppenalty}{
     The GCV penalty used during pruning.
     A copy of \code{earth's} \code{ppenalty} argument.
  }
  \item{call}{
     The call used to invoke \code{earth}.
  }
  \item{terms}{
     Model frame terms.
     This component exists only if the model was built using \code{earth.formula}.
  }
  \item{x}
  \item{y}
  \item{subset}{
     Copy of input arguments \code{x}, \code{y}, and  \code{subset}.
     These components exist only if \code{keepxy=TRUE}.
  }
}
\note{
\strong{Standard Model Functions}

Standard model functions such as \code{\link{case.names}}
are provided for \code{earth} objects and are not explicitly documented.

\strong{Other Implementations}

The results are similar to but not identical to other
Multivariate Adaptive Regression Splines implementations.
The differences stem from the forward pass where very small
implementation differences (or perturbations of the input data) can cause
rather different selection of terms and knots.
The backward passes give identical or near identical results,
given the same forward pass results.

The source code of \code{earth} is derived from \code{\link[mda]{mars}}
in the \code{mda} package written by
by Trevor Hastie and Robert Tibshirani.
See also \code{\link{mars.to.earth}}.
Earth's principal advantages over mda::mars are that it is
much faster and has plotting and printing methods.

The term \sQuote{MARS} is trademarked and licensed exclusively to
Salford Systems \url{http://www.salfordsystems.com}.
Their implementation uses an engine written by Friedman and
offers more features than \code{earth}.

\strong{Limitations}

The following aspects of MARS
are mentioned in Friedman's papers but not implemented in \code{earth}:

i) Piecewise cubic models\cr
ii) Specifying which predictors must enter linearly\cr
iii) Specifying which predictors can interact\cr
iv) Model slicing (\code{\link{plotmo}} goes part way)\cr
v) Handling missing values\cr
vi) Logistic regression \cr
vii) Special handling of categorical predictors\cr
viii) Fast MARS h parameter.

\strong{Large Models and Execution Time}

For a given set of input data,
the following can increase the speed of the forward pass:

i) decreasing \code{fast.k}\cr
ii) decreasing \code{nk}\cr
iii) decreasing \code{degree}\cr
iv) increasing \code{threshold}\cr
v) increasing \code{min.span}.

The backward pass is normally much faster than the forward pass,
unless \code{pmethod="exhaustive"}.
Reducing \code{npune} reduces exhaustive search time.
One strategy is to first build a large model
and then adjust pruning parameters using \code{\link{update.earth}}.

\strong{The Forward Pass}

The forward pass adds terms in pairs until the first of the
following conditions is met:

i) reach maximum number of terms \code{(nterms>=nk)}\cr
ii) reach DeltaRSq threshold \code{(DeltaRSq<thresh)} where
DeltaRSq is the difference in R-Squared caused by adding the current term pair\cr
iii) reach max RSq \code{(RSq>1-thresh)}\cr
iv) reach min GRSq \code{(GRSq< -10)}.

Set \code{trace>=2} to see the stopping condition.

The result of the forward pass is the set of terms defined by \code{$dirs}
and \code{$cuts}.
As a final step, the forward pass deletes linearly dependent terms, if any,
so all terms in \code{$dirs} and \code{$cuts} are independent.

Note that GCVs (via GRSq) are used during the forward pass only as one of the
stopping conditions and in trace prints.
Changing the GCV penalty does not change the knot positions.

\strong{The Pruning Pass}

The pruning pass is handed the sets of terms created by the forward pass
and works like this:
it determines the subset of terms (using \code{pmethod})
with the lowest RSS for each model size in \code{1:nprune}
(see \code{Force.xtx.prune} for some details).
It saves the RSS and term numbers for each such subset in \code{rss.per.subset}
and \code{prune.terms}.
It then applies the \code{Get.crit} function with \code{ppenalty}
to \code{rss.per.subset} to yield \code{gcv.per.subset}.
It chooses the model with lowest value in \code{gcv.per.subset},
and puts its term numbers into \code{selected.terms}.
Finally, it runs \code{\link{lm}} to determine the
\code{fitted.values}, \code{residuals}, and \code{coefficients}, by regressing
the response \code{y} on the \code{selected.terms} of \code{bx}.

Set \code{trace>=3} to trace the pruning pass.

By default \code{Get.crit} is \code{earth:::get.gcv}.
Alternative \code{Get.crit} functions can be defined.
See the source code of \code{get.gcv} for an example.

\strong{Testing on New Data}

This example demonstrates one way to train on 80\% of the data and
test on the remaining 20\%.
In practice a larger dataset should be used for splitting.
Also, remember that the test set should not be used for parameter tuning ---
use separate validation set(s) for that.
\preformatted{
train.subset <- sample(1:nrow(trees), .8 * nrow(trees))
test.subset <- (1:nrow(trees))[-train.subset]
a <- earth(Volume ~ ., data = trees[train.subset, ])
yhat <- predict(a, newdata = trees[test.subset, ])
y <- trees$Volume[test.subset]
print(1 - sum((y - yhat)^2) / sum((y - mean(y))^2)) # print R-Squared
}
\strong{Using earth with fda and mda}

Earth can be used with \code{\link[mda]{fda}} and \code{\link[mda]{mda}}
in the \code{\link[mda]{mda}} package.
Use \code{keepxy=TRUE} if you want to call \code{\link{plotmo}} later.
Use \code{trace>=3} to see the call to \code{earth} generated by \code{fda} or \code{mda}.
Example:
\preformatted{
library(mda)
(a <- fda(Species~., data=iris, method=earth, keepxy=TRUE))
plot(a)
plotmo(a$fit, ycolumn=1, ylim=c(-1.5,1.5), clip=FALSE)
plotmo(a$fit, ycolumn=2, ylim=c(-1.5,1.5), clip=FALSE)
}

\strong{Warning and Error Messages}

\code{Earth} prints most error and warning messages without
printing the \sQuote{call}.
If you are mystified by a warning message, try
setting \code{\link{options}(warn=2)}
and using \code{\link{traceback}}.
}
\author{
  Stephen Milborrow, derived from \code{mda::\link[mda]{mars}}
  by Trevor Hastie and Robert Tibshirani.

This is an early release and users are encouraged to send feedback --- use
milboATsonicPERIODnet.
}
\references{
  The primary references are
  the Friedman papers.  Readers may find the MARS section in Hastie, Tibshirani,
  and Friedman a more accessible introduction.  Faraway takes a hands-on approach,
  using the \code{\link[=ozone1]{ozone}} data to compare \code{mda::mars} with other techniques.
  (If you use Faraway's examples with \code{earth} instead of \code{mars}, use \code{$bx}
  instead of \code{$x}).
  Earth's pruning pass uses \code{\link[leaps]{leaps}} which is based on
  techniques in Miller.

  Faraway \emph{Extending the Linear Model with R}
  \url{http://www.maths.bath.ac.uk/~jjf23}

  Friedman  (1991) \emph{Multivariate Adaptive Regression Splines (with discussion)}
  Annals of Statistics 19/1, 1--141

  Friedman  (1993) \emph{Fast MARS}
  Stanford University Department of Statistics, Technical Report 110
  \url{http://www-stat.stanford.edu/research/index.html}

  Hastie, Tibshirani, and Friedman (2001) \emph{The Elements of Statistical Learning}
  \url{http://www-stat.stanford.edu/~hastie/pub.htm}

  Miller, Alan (1990, 2nd ed. 2002) \emph{Subset Selection in Regression}
}
\seealso{
  \code{\link{format.earth}},
  \code{\link{get.nterms.per.degree}},
  \code{\link{get.nused.preds.per.subset}},
  \code{\link{mars.to.earth}},
  \code{\link{model.matrix.earth}},
  \code{\link{ozone1}},
  \code{\link{plot.earth.models}},
  \code{\link{plot.earth}},
  \code{\link{plotmo}},
  \code{\link{predict.earth}},
  \code{\link{reorder.earth}},
  \code{\link{summary.earth}},
  \code{\link{update.earth}}
}
\examples{
a <- earth(Volume ~ ., data = trees)
summary(a, digits = 2)

# yields:
#    Call:
#    earth(formula = Volume ~ ., data = trees)
#
#    Expression:
#      23
#      +  5.7 * pmax(0,  Girth -     13)
#      -  2.9 * pmax(0,     13 -  Girth)
#      + 0.72 * pmax(0, Height -     76)
#
#    Number of cases: 31
#    Selected 4 of 5 terms, and 2 of 2 predictors
#    Number of terms at each degree of interaction: 1 3 (additive model)
#    GCV: 11     RSS: 213     GRSq: 0.96     RSq: 0.97
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

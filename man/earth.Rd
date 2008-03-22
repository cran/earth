\name{earth}
\alias{earth}
\alias{earth.default}
\alias{earth.formula}
\concept{regression}
\concept{mars}
\concept{Friedman}
\title{Multivariate Adaptive Regression Splines}
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
      nk = max(21, 2 * NCOL(x) + 1), degree = 1,
      linpreds = FALSE, allowed  = NULL,
      thresh = 0.001, minspan = 0, newvar.penalty = 0,
      fast.k = 20, fast.beta = 1,
      pmethod = "backward", nprune = NULL,
      Object  = NULL, Get.crit = get.gcv,
      Eval.model.subsets = eval.model.subsets,
      Print.pruning.pass = print.pruning.pass,
      Force.xtx.prune = FALSE, Use.beta.cache = TRUE, \dots)
}
\arguments{
All arguments are optional except \code{formula}, or \code{x} and \code{y}.
NAs are not allowed.\cr
\cr
To start off, look at the arguments
\code{formula},
\code{x},
\code{y},
\code{nk}, and
\code{degree}.
Many users will find that those arguments are all they need,
plus in some cases \code{nprune}, \code{penalty}, and \code{minspan}.

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
     If the \code{y} values are very big or very small you may get
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
    Generalised Cross Validation (GCV) penalty per knot.
    Default is \code{if(degree>1) 3 else 2}.
    A value of 0 penalises only terms, not knots.
    The value -1 is a special case, meaning no penalty, so GCV=RSS/n.
    Theory suggests values in the range of about 2 to 4.
    In practice, for big data sets larger values can be useful to force a smaller model.
  }
  \item{trace}{
     Trace earth's execution.  Default is 0.  Values:

     0 none 1 overview 2 forward pass 3 pruning 4 more pruning 5 \dots
  }
  \item{keepxy}{
     Set to TRUE to retain \code{x}, \code{y}, and  \code{subset} in the returned value.
     Default is FALSE\cr

     \strong{The following arguments are for the forward pass}
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
  \item{degree}{
     Maximum degree of interaction (Friedman's \eqn{mi}).
     Default is 1, meaning build an additive model (i.e. no interaction terms).
  }
  \item{linpreds}{
      \strong{[new Jan 2008]}\cr
     Index vector specifying which predictors should enter linearly,
     as in \code{\link{lm}}.\cr
     The default is FALSE, meaning all predictors enter
     in the standard MARS fashion i.e. in hinge functions.\cr
     A predictor's index in \code{linpreds} is the column number in the input matrix \code{x}.
     If you use the formula interface to \code{earth}, 
     earth's internal call to \code{\link{model.matrix}} may reorder or add columns.
     For example, you may have to exclude the response when counting columns.
     As a sanity check, use \code{trace >= 2} to see a list of numbered predictors.\cr
     An example is given in the FAQ section below.\cr
     Note: in the current implementation, the GCV penalty for predictors
     that enter linearly is the same as that for predictors with knots.
     That is not quite correct; linear terms should be penalised less.
  }
  \item{allowed}{
      \strong{[new Jan 2008]}\cr
     Function specifying which predictors can interact and how.
     Default is NULL, meaning all standard MARS terms are allowed.\cr
     Earth calls the \code{allowed} function just before adding a term.
     If \code{allowed} returns TRUE the term goes into the model as usual;
     if \code{allowed} returns FALSE the term is discarded.
     Your \code{allowed} function should have the following prototype\cr
     \code{allowed <- function(degree, pred, parents) \{...\}}\cr
     where\cr
        \code{degree} is the interaction degree of the candidate term.
              Will be 1 for additive terms.\cr
        \code{pred} is the index of the candidate predictor.
              For details of predictor indexing, see the description of
              the \code{linpreds} argument above.\cr
        \code{parents} is the candidate parent term's row in \code{dirs}.\cr
     An example is given in the FAQ section below.\cr
  }
  \item{thresh}{
     Forward stepping threshold.
     Default is 0.001.
     This is one of the arguments used to decide when forward stepping
     should terminate.
     See the section below on the forward pass.
  }
  \item{minspan}{
    Minimum distance between knots.\cr
    Use a value of 1 to consider all x values (which is good if the
    data are not noisy).\cr
    The default value is 0
    The value 0 is treated specially and
    means calculate the \code{minspan} internally as per
    Friedman's MARS paper section 3.8 with \eqn{alpha} = 0.05.
    Set \code{trace>=2} to see the calculated value.\cr
    (The \strong{default was 1 prior to Mar 2008} but was changed for
    compatibility with \code{mda::\link[mda]{mars}} and Friedman's papers).\cr
    This argument is intended to increase resistance to runs of noise in the input data.
    Higher values will increase smoothness in your model.\cr
    Note: predictor value extremes are ineligible for knots
    regardless of the \code{minspan} setting, as per the MARS paper equation 45.
  }
  \item{newvar.penalty}{
     Penalty for adding a new variable in the forward pass
     (Friedman's \eqn{gamma}, equation 74 in the MARS paper).
     Default is 0, meaning no penalty for adding a new variable.
     Useful non-zero values range from about 0.01 to 0.2 --- you will
     need to experiment.
     This argument can mitigate the effects of colinearity or concurvity
     in the input data, but anecdotal evidence is that it does not work very well.
     If you know two variables are strongly correlated then you
     would do better to delete one of them before building your model.
  }
  \item{fast.k}{
     Maximum number of parent terms considered at each step of the forward pass.
     Friedman invented this parameter to speed up the forward pass
     (see the Fast MARS paper section 3.0).
     Default is 20.
     Values of 0 or less are equivalent to infinity, meaning no Fast MARS.
     Typical values, apart from 0, are 20, 10, or 5.
     In general, with a lower \code{fast.k} (say 5), \code{earth} is faster;
     with a higher \code{fast.k}, or with \code{fast.k} disabled (set to 0),
     \code{earth} builds a better model.
     However it is not unusual to get a slightly better model with a lower \code{fast.k},
     and you may need to experiment.
  }
  \item{fast.beta}{
     Fast MARS ageing coefficient, as described in the
     Fast MARS paper section 3.1.
     Default is 1.
     A value of 0 sometimes gives better results.\cr

        \strong{The following arguments are for the pruning pass.}\cr
        One strategy is to build a large model and then trim it
        back using \code{\link{update.earth}} with the arguments below.
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
     (i.e. you have to kill your R session to interrupt;
      in Windows you can use \code{tskill}).
  }
  \item{nprune}{
     Maximum number of terms (including intercept) in the pruned model.
     Default is NULL, meaning all terms created by the forward pass
     (but typically not all terms will remain after pruning).
     Use this to reduce exhaustive search time, or
     to enforce an upper bound on the model size.\cr

    \strong{The following arguments are for internal or advanced use}
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
     Function to evaluate model subsets
     --- see notes in source code.
  }
  \item{Print.pruning.pass}{
     Function to print pruning pass results
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
  \item{Use.beta.cache}{
    Default is TRUE.
    Using the "beta cache" takes more memory but is faster
    (by 20\% and often much more for large models).
    The beta cache uses \code{nk * nk * ncol(x) * sizeof(double)} bytes.
    Set \code{Use.beta.cache=FALSE} to save memory.
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
     The GCV is calculated using \code{penalty} (as are all returned GCVs).
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
  \item{bx}{
     Matrix of basis functions applied to \code{x}.
     Each column corresponds to a selected term.
     Each row corresponds to a row in in the input matrix \code{x},
     after taking \code{subset}.
     See \code{\link{model.matrix.earth}} for an example of \code{bx} handling.
     For brevity, "h" is used instead of "pmax" in column names.
     Example \code{bx}:\preformatted{    (Intercept) h(Girth-12.9) h(12.9-Girth) h(Girth-12.9)*h(...
[1,]          1           0.0           4.6                   0
[2,]          1           0.0           4.3                   0
[3,]          1           0.0           4.1                   0
...}
}
\cr
  \item{dirs}{
     Matrix with one row per MARS term, and with with ij-th element equal to\cr

\code{0}  if predictor j is not in term i\cr
\code{-1} if a factor of the form \code{pmax(c - xj)} is in term i\cr
\code{1}  if a factor of the form \code{pmax(xj - c)} is in term i\cr
\code{2}  if predictor j enters term i linearly \strong{[new Jan 2008]}.\cr

     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Example \code{dirs}:\preformatted{                       Girth Height
(Intercept)                0  0 #no factors in intercept
h(Girth-12.9)              1  0 #2nd term uses Girth
h(12.9-Girth)             -1  0 #3rd term uses Girth
h(Girth-12.9)*h(Height-76) 1  1 #4th term uses Girth and Height
...}
  }
\cr
  \item{cuts}{
     Matrix with ij-th element equal to the cut point
     for predictor j in term i.
     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Example \code{cuts}:\preformatted{                           Girth Height
(Intercept)                    0   0  #intercept, no cuts
h(Girth-12.9)               12.9   0  #2nd term has cut at 12.9
h(12.9-Girth)               12.9   0  #3rd term has cut at 12.9
h(Girth-12.9)*h(Height-76)  12.9  76  #4th term has two cuts
...}
  }
\cr
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
     Example \code{prune.terms}:\preformatted{[1,]    1  0  0  0  0  0  0 #intercept-only model
[2,]    1  2  0  0  0  0  0 #best 2 term model uses terms 1,2
[3,]    1  2  4  0  0  0  0 #best 3 term model uses terms 1,2,4
[4,]    1  2  9  8  0  0  0 #and so on
...}
  }
\cr
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
    A vector of the GRSq for each response.  Length is \code{ncol(y)}.
  }
  \item{rss.per.subset}{
     A vector of the RSS for each model subset generated by the pruning pass.
     Length is \code{nprune}.
     If \code{y} has multiple columns, the RSS is summed over all responses for each subset.
     The null RSS (i.e. the RSS of an intercept only-model) is \code{rss.per.subset[1]}.
     The \code{rss} above is\cr
\code{rss.per.subset[length(selected.terms)]}.
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
  \item{penalty}{
     The GCV penalty used during pruning.
     A copy of \code{earth's} \code{penalty} argument.
  }
  \item{ppenalty}{
     Deprecated.
     A copy of \code{earth's} \code{penalty} argument.
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
\strong{Other implementations}

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
Earth's principal advantages over \code{mda::mars} are that it is
much faster and provides plotting and printing methods.

The term \sQuote{MARS} is trademarked and licensed exclusively to
Salford Systems \url{http://www.salfordsystems.com}.
Their implementation uses an engine written by Friedman and
offers more features than \code{earth}.

StatSoft also have an implementation which they call MARSplines
\url{http://www.statsoft.com/textbook/stmars.html}.

\strong{Limitations}

The following aspects of MARS
are mentioned in Friedman's papers but not implemented in \code{earth}:

i) Piecewise cubic models\cr
ii) Model slicing (\code{\link{plotmo}} goes part way)\cr
iii) Handling missing values\cr
iv) Logistic regression \cr
v) Special handling of categorical predictors\cr
vi) Fast MARS h parameter\cr
vii) Cross validation to determine \code{penalty}\cr
viii) Anova tables with sigma and other information.

\strong{The forward pass}

Understanding the details of the forward and pruning passes
will help you understand earth's return value and
the admittedly large number of arguments.
The result of the forward pass is \code{bx} and
the set of terms defined by \code{dirs}
and \code{cuts} (these are all fields in earth's return value,
but the \code{bx} here includes all terms before trimming
back to \code{selected.terms}).

The forward pass adds terms in pairs until the first of the
following conditions is met:

i) reach maximum number of terms \code{(nterms >= nk)}\cr
ii) reach DeltaRSq threshold \code{(DeltaRSq < thresh)} where
DeltaRSq is the difference in R-Squared caused by adding the current term pair\cr
iii) reach max RSq \code{(RSq > 1-thresh)}\cr
iv) reach min GRSq \code{(GRSq < -10)}.

Set \code{trace>=1} to see the stopping condition.

Note that GCVs (via GRSq) are used during the forward pass only as one of the
stopping conditions and in \code{trace} prints.
Changing the \code{penalty} argument does not change the knot positions.

The various stopping conditions mean that the actual number of terms
created by the forward pass may be less than \code{nk}.
There are other
reasons why the actual number of terms may be less than \code{nk}:
(i) the forward pass discards one side of each term pair
if it adds nothing to the model ---
but the forward pass counts terms as if they were actually created in pairs,
and, (ii) as a final step, the forward pass deletes linearly dependent terms, if any,
so all terms in \code{dirs} and \code{cuts} are independent.
And remember that the pruning pass will further discard terms.

\strong{The pruning pass}

The pruning pass is handed the sets of terms created by the forward pass.
Its job is to find the subset of those terms that gives the lowest GCV.
The following description of the pruning pass
explains how various fields in earth's returned value are generated.

The pruning pass works like this:
it determines the subset of terms (using \code{pmethod})
with the lowest RSS for each model size in \code{1:nprune}
(see the \code{Force.xtx.prune} argument above for some details).
It saves the RSS and term numbers for each such subset in \code{rss.per.subset}
and \code{prune.terms}.
It then applies the \code{Get.crit} function with \code{penalty}
to \code{rss.per.subset} to yield \code{gcv.per.subset}.
It chooses the model with the lowest value in \code{gcv.per.subset},
and puts its term numbers into \code{selected.terms}.
Finally, it runs \code{\link{lm}} to determine the
\code{fitted.values}, \code{residuals}, and \code{coefficients}, by regressing
the response \code{y} on the \code{selected.terms} of \code{bx}.

Set \code{trace>=3} to trace the pruning pass.

By default \code{Get.crit} is \code{earth:::get.gcv}.
Alternative \code{Get.crit} functions can be defined.
See the source code of \code{get.gcv} for an example.

\strong{Large models: execution time and memory}

For a given set of input data,
the following can increase the speed of the forward pass:

i) decreasing \code{fast.k}\cr
ii) decreasing \code{nk}\cr
iii) decreasing \code{degree}\cr
iv) increasing \code{threshold}\cr
v) increasing \code{min.span}.

The backward pass is normally much faster than the forward pass,
unless \code{pmethod="exhaustive"}.
Reducing \code{nprune} reduces exhaustive search time.
One strategy is to first build a large model
and then adjust pruning parameters such as \code{nprune} using \code{\link{update.earth}}.

The following very rough rules of thumb apply for large models.
Using \code{minspan=1} instead of the default \code{0} will increase times by 20 to 50\%.
Using \code{fast.k=5} instead of the default \code{20} can give you substantial speed gains
but will sometimes give you a much smaller \code{GRSq}.
Using an \code{allowed} function slows down model building by about 10\%.

Earth does not impose specific limits on the model size.
Model size is limited only by the amount of memory on your system and your patience.
An example: the earth test suite builds a model using \code{earth.default}
with a 1e4 by 100 input matrix with \code{nk=21} and \code{degree=1}.
Windows XP reports that the peak memory use when building this model is 47 MBytes.
Using the formula interface to earth pushes memory to 62 MBytes.
Increasing the number of rows in the input matrix to 1e5 pushes memory to 240 MBytes.
These numbers are for earth version 1.3 (Mar 2008),
which requires less memory than previous versions.

Increasing the \code{degree} does not change the memory requirement but greatly increases the running time.
Earth uses more memory if \code{trace>=2}.
Earth uses more memory if any elements of the arguments \code{x} and \code{y} are not \link{double},
because it must convert them to double internally.

\strong{Migrating from mda::mars}

Using earth is usually just a matter
of changing the call from "\code{mars}" to "\code{earth}".
But there are a few argument differences and
earth will issue a warning if you give it a mars-only argument.
See also \link{mars.to.earth}.

The resulting model will be similar but not identical because
of small implementation differences which are magnified
by the inherent instability of the MARS forward pass.

If you are further processing the output of earth you will need to
consider differences in the returned value.  The header of the
source file \code{mars.to.earth.R} describes these.
Perhaps the most important is that \code{mars} returns the
MARS basis matrix in a field called "\code{x}"
whereas \code{earth} returns "\code{bx}" with only the selected terms.
Also, \code{earth} returns "\code{dirs}" rather than "\code{factors}",
and in \code{earth} this matrix can have entries of value 2 for linear predictors.

Earth does not yet support the \code{mda::mars} \code{w} argument 
However, this argument is is actually ignored internally by \code{mda::mars}.
Support for the \code{mda::mars} argument \code{wp} will be added
to earth in an upcoming release.

\strong{Multiple response models}

If \code{y} has K columns then earth builds K simultaneous models.
Each model has the same set of basis functions
(the same \code{bx}, \code{selected.terms}, \code{dirs} and \code{cuts})
but different coefficients (the returned \code{coefficients} will have K columns).
The models are built and pruned as usual but with the GCVs
and RSSs averaged across all K responses.

Since earth attempts to optimise for all models simultaneously,
the results will not be as "good" as building the models independently
i.e. the GCV of the combined model will usually not be as good as
the GCVs for independently built models.
However, the combined model may be a better model in other senses,
depending on what you are trying to achieve.
For example, it could be useful for earth to select
the set of predictors that is best across \emph{all} responses.

For more details on using residual errors averaged over multiple responses see section 4.1 of
Hastie, Tibshirani, and Buja \emph{Flexible Discriminant Analysis by Optimal Scoring},
JASA, December 1994 \url{http://www-stat.stanford.edu/~hastie/Papers/fda.pdf}.

\strong{Using earth with fda and mda}

\code{Earth} can be used with \code{\link[mda]{fda}} and \code{\link[mda]{mda}}
in the \code{mda} package.  Earth will generate a multiple response model,
as described above.
Use the \code{fda/mda} argument \code{keep.fitted=TRUE} if
you want to call \code{\link{plot.earth}} later
(actually only necessary for large datasets, see the description of
\code{keep.fitted} in \code{\link[mda]{fda}}).
Use the \code{earth} argument \code{keepxy=TRUE} if you want to call
\code{\link{update.earth}} or \code{\link{plotmo}} later.
Example:
\preformatted{
library(mda)
(a <- fda(Species~., data=iris, keep.fitted=TRUE, method=earth, keepxy=TRUE))
plot(a)
summary(a$fit)  # examine earth model embedded in fda model
plot(a$fit)
plotmo(a$fit, ycolumn=1, ylim=c(-1.5,1.5), clip=FALSE)
plotmo(a$fit, ycolumn=2, ylim=c(-1.5,1.5), clip=FALSE)
}
\strong{Standard model functions}

Standard model functions such as \code{\link{case.names}}
are provided for \code{earth} objects and are not explicitly documented.\cr\cr

\strong{FREQUENTLY ASKED QUESTIONS}
\cr\cr
\strong{What are your plans for earth?}

The next set of releases will support logistic regression
and functions for model validation.

\strong{How can I train on one set of data and test on another?}

The example below demonstrates one way to train on 80\% of the data and
test on the remaining 20\%.
\preformatted{
train.subset <- sample(1:nrow(trees), .8 * nrow(trees))
test.subset <- (1:nrow(trees))[-train.subset]
a <- earth(Volume ~ ., data = trees[train.subset, ])
yhat <- predict(a, newdata = trees[test.subset, ])
y <- trees$Volume[test.subset]
print(1 - sum((y - yhat)^2) / sum((y - mean(y))^2)) # print R-Squared
}
In practice a dataset larger than the one in the example should be used for splitting.
The model variance is too high with this small set --- run the example a few
times to see how the model changes as \code{\link[basic]{sample}}
splits the dataset differently on each run.
Also, remember that the test set should not be used for parameter tuning
because you will be optimising for the test set ---
use GCVs or separate validation sets for that.

\strong{How can I establish variable importance?}

Establishing predictor importance is in general a tricky and even controversial problem.

Running \code{\link{plotmo}} with \code{ylim=NULL} (the default)
gives an idea of which predictors make the largest changes to the predicted value.

You can use \code{\link{drop1}}
(assuming you are using the formula interface to earth).
Calling \code{drop1(my.earth.model)} will delete each predictor in turn from your model,
rebuild the model from scratch each time, and calculate the GCV each time.
You will get warnings that the earth library function \code{extractAIC.earth} is
returning GCVs instead of AICs --- but that is what you want so you can
ignore the warnings.
The column labeled \code{AIC} in the printed response
from \code{\link{drop1}} will actually be a column of GCVs not AICs.
The \code{Df} column is not much use in this context.
Note that \code{\link{drop1}} drops \emph{predictors} from the model
while earth's pruning pass drops \emph{terms}.

You will get lots of output from \code{\link{drop1}} if you built your original earth
model with \code{trace>0}.
You can set \code{trace=0} by updating your model before calling \code{\link{drop1}}.
Do it like this:\cr
\code{my.model <- \link{update.earth}(my.model, trace=0)}.

You can also use \code{\link[caret]{varImp}} in the \code{caret} package,
which calculates for each variable the effect
of all MARS terms containing that variable.

Remember that these techniques only tell you how important
a variable is with the other variables already in the model.

MARS models have a high variance --- if the data changes a little,
the set of basis terms created by the forward pass can change a lot.
So the above estimates of predictor importance are unreliable
because they can vary with even slightly different training data.
Furthermore, if two predictors are closely related, the forward pass will
somewhat arbitrarily choose one over the other.
The chosen predictor will incorrectly appear more important.
There are better ways of measuring variable importance (using
resampling) but they are not yet implemented in the \code{earth} package.

\strong{Which predictors were added to the model first?}

You can see the forward pass adding terms with \code{trace=2} or higher.
But remember, pruning will usually remove some of the terms.
You can also use
\preformatted{
summary(my.model, decomp="none")
}
which will list the basis functions remaining after pruning,
in the order they were added by the forward pass.

\strong{Which predictors are actually in the model?}

Well, you can always use \code{\link{summary.earth}}.
Or use the following function which will give you a
list of the predictors in the model:
\preformatted{
get.used.pred.names <- function(obj)   # obj is an earth object
  names(which(apply(obj$dirs[obj$selected.terms,,drop=FALSE],2,any)))
}
\strong{Why do I get fewer terms than \code{nk}, even with \code{prune="none"}?}

There are several conditions that can terminate the forward pass,
and reaching \code{nk} is just one of them.
See the section above on the forward pass.

\strong{Why do I get fewer terms than \code{nprune}?}

The pruning pass selects a model with the lowest GCV
that has \code{nprune} or fewer terms.
Thus the \code{nprune} argument specifies the \emph{maximum}
number of permissible terms in the final pruned model.

You can work around this because you will get exactly \code{nprune}
terms if you specify \code{pmethod="none"}.
Compare the output of these two examples:
\preformatted{
earth(Volume ~ ., data = trees, trace=3)
earth(Volume ~ ., data = trees, trace=3, pmethod="none")
}
Another way to get exactly \code{nprune} terms is to specify \code{penalty = -1}.
This special value of \code{penalty} causes earth to set the GCV to \code{RSS/nrow(x)}.
Since the training RSS always decreases with more terms,
the pruning pass will choose the maximum allowable number of terms.
An example:
\preformatted{
earth(Volume ~ ., data = trees, trace=3, penalty=-1)
}
\strong{Can you give an example of the \code{linpreds} argument?}

With \code{linpreds} you can specify which predictors should enter linearly,
instead of in hinge functions.
Starting with
\preformatted{
a1 <- earth(Volume ~ ., data = trees)
plotmo(a1)
}
we see in the \code{\link{plotmo}} graphs that \code{Height}
is not as important as \code{Girth}.
For collaborative evidence that \code{Girth} is a more reliable
indicator of \code{Volume} you can use \code{\link{pairs}}:
\preformatted{
pairs(trees, panel = panel.smooth)
}
Since we want the simplest model that describes the data,
we can specify that \code{Height} should enter linearly:
\preformatted{
a2 <- earth(Volume ~ ., data = trees, linpreds = 2)  # 2 is Height column
summary(a2)
}
which yields
\preformatted{
Expression:
  -7.41
  + 0.418 * Height
  +  5.86 * pmax(0,  Girth -   12.9)
  -  2.41 * pmax(0,   12.9 -  Girth)
}
In this example, the second simpler model has almost the same RSS as the first model.
We can make both \code{Girth} and \code{Height} enter linearly with
\preformatted{
a2 <- earth(Volume ~ ., data = trees, linpreds = c(1,2))
}
or with
\preformatted{
a3 <- earth(Volume ~ ., data = trees, linpreds = TRUE)
}
But specifying that all predictors should enter linearly is not really a useful thing to do.
In our simple example, the all-linear MARS model is the
same as a standard linear model (compare the \code{\link{summary}} for each)
\preformatted{
a4 <- lm(Volume ~ ., data = trees)
}
but in general that will not be true.
The \code{lindeps} argument does not stipulate that a predictor \emph{must} enter the model,
only that if it enters it should enter linearly.
The forward pass will stop if it realizes
that adding more predictors does not improve the model.

\strong{Can you give an example of the \code{allowed} argument?}

You can specify how variables are allowed to enter MARS terms
with the \code{allowed} argument.
The interface is flexible but requires a bit of programming.
We start with a simple example, which completely excludes one predictor
from the model:
\preformatted{
example1  <- function(degree, pred, parents)   # returns TRUE if allowed
{
    pred != 2  # disallow predictor 2, which is "Height"
}
a1 <- earth(Volume ~ ., data = trees, allowed = example1)
print(summary(a1))
}
But that's not much use, because it's simpler to exclude the predictor
from the input matrix when invoking earth:
\preformatted{
a1 <- earth(Volume ~ . - Height, data = trees)
}
Here is a more useful example. It prevents the specified predictor
from being used in interaction terms.
The \code{parents} argument is the candidate parent's row in the \code{dirs} matrix
(\code{dirs} is described in the Value section above).
Each entry of \code{parents} is 0, 1, -1, or 2, and you index
\code{parents} on the predictor index.
Thus \code{parents[pred]} is 0 if \code{pred} is not in the parent term.

The example is artificial because it's unlikely you
would want to single out humidity from interactions in the ozone data.
\preformatted{
example2 <- function(degree, pred, parents)
{
    # disallow humidity in terms of degree > 1
    # 3 is the "humidity" column in the input matrix
    if (degree > 1 && (pred == 3 || parents[3]))
        return(FALSE)
    TRUE
}
a <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example2)
print(summary(a))
}
The following example allows only the specified predictors
in interaction terms:
\preformatted{
example3 <- function(degree, pred, parents)
{
    # allow only humidity and temp in terms of degree > 1
    # 3 and 4 are the "humidity" and "temp" columns
    allowed.set = c(3,4)
    if (degree > 1 &&
           (all(pred != allowed.set) || any(parents[-allowed.set])))
        return(FALSE)
    TRUE
}
a <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example3)
print(summary(a))
}
The basic MARS model building strategy is always 
applied even when there is an allowed function.
For example, earth considers a term for addition only
if all factors of that term except the new one are already in a term.
This means that an allowed function that inhibits, say, all degree 2 
terms will also effectively inhibit higher degrees too because
there will be no degree 2 terms for earth to extend to degree 3.

You can expect model building to be about 10\% slower with an \code{allowed} function.

\strong{How does \code{summary.earth} order terms?}

With \code{decomp="none"},
the terms are ordered as created by the forward pass.

With the default \code{decomp="anova"},
the terms are ordered in increasing order of interaction.
In detail:\cr
(i) terms are sorted first on degree of interaction\cr
(ii) then terms with a \code{linpreds} linear factor before standard terms\cr
(iii) then on the predictors in the factors (in
the order of the columns in the input matrix)\cr
(ii) and finally on increasing knot values.

It's actually \code{earth:::reorder.earth} that does the ordering.

A slightly different and less predictable ordering was used prior to earth version 1.2 (Jan 2008).
The change also affects the ordering of \code{degree2} \code{\link{plotmo}} plots.

\strong{Why \code{pmax} and not \code{max} in the output from \code{summary.earth}?}

With \code{\link{pmax}} the earth equation is an R expression
that can handle multiple cases.
Thus the expression is consistent with the
way \code{\link{predict.earth}} works --- you can give \code{\link{predict}}
multiple cases (i.e. multiple rows in the input matrix)
and it will return a vector of predicted values.

\strong{What about boosting MARS?}

If you want to boost, use boosted trees rather than boosted MARS --- you
will get better results.

More precisely, although boosted MARS gives better results than plain MARS,
if you would like to improve prediction performance (at the cost
of a more complicated and less interpretable model)
you will usually get better results with
boosted \code{\link[rpart]{rpart}} trees than with boosted MARS.
See Gillian Ward (2007) \emph{Statistics in Ecological Modeling:
Presence-Only Data and Boosted Mars (Doctoral Thesis)}
\url{http://www-stat.stanford.edu/~hastie/THESES/Gill\_Ward.pdf}.

This could change as the state of the art advances.

\strong{What is a GCV, in simple terms?}

GCVs are important for MARS because the pruning pass
uses GCVs to evaluate model subsets.

In general terms,
when testing a model
(not necessarily a MARS model)
we want to test \emph{generalisation} performance and so
want to measure error on independent data i.e. not on the training data.
Often a decent set of independent data is unavailable and so we
resort to cross validation or leave-one-out methods.
But that can be painfully slow.
For certain forms of model we can use a formula to
approximate the error that would determined by leave-one-out
validation --- that approximation is the GCV.
The formula adjusts (i.e. increases) the mean training RSS to take into
account the flexibility of the model.
Summarising, the GCV approximates the mean RSS that would
be measured on independent data.
Even when the approximation is not that good it is usually
good enough for comparing models during pruning.

GCVs were introduced by Craven and Wahba, and extended by Friedman for MARS.
See Hastie et al. p216 and the Friedman MARS paper.
GCV stands for "Generalised Cross Validation", a perhaps misleading term.

The \code{GRSq} measure used in the earth package standardises the raw GCV,
in the same way that R-Squared standardises the RSS.

\strong{If GCVs are so important, why don't linear models use them?}

The more flexible a model, the more its propensity to overfit the training data.
(An overfitted model fits the training data well but will not
give good prediction results on new data.)
Linear models are constrained, with usually only a few parameters, 
and don't have the tendency to overfit the data like more flexible models such as MARS.
This means that for linear models the RSS on the data used to build the model is usually an
adequate measure of generalisation ability, assuming enough data.

This is no longer true if you do automatic variable selection on linear models,
because the process of selecting variables increases the flexibility
of the model. Hence the AIC --- as used in, say, \code{\link{drop1}}.
The GCV, AIC, and friends are means to the same end.
Depending on what information is available during model building.
we use one of these statistics to estimate model generalisation performance 
for the purpose of selecting a model.

\strong{What happened to \code{get.nterms.per.degree}, 
\code{get.nused.preds.per.subset}, and \code{reorder.earth}?}

From release 1.3.0, some earth functions are no longer public,
to help simplify the user interface.
The functions are still available (and stable) if you need them --- 
use for example \code{earth:::reorder.earth()}.

\strong{What happened to the \code{ppenalty} argument?}

This was removed (release 1.3.1) because it is no longer needed.
The \code{ponly} argument of \code{\link{update.earth}}
is a more flexible way of achieving the same end.
}
\author{
  Stephen Milborrow, derived from \code{mda::\link[mda]{mars}}
  by Trevor Hastie and Robert Tibshirani.

Users are encouraged to send feedback --- use
milboATsonicPERIODnet \url{http://www.milbo.users.sonic.net}.
}
\references{
  The primary references are
  the Friedman papers.  Readers may find the MARS section in Hastie, Tibshirani,
  and Friedman a more accessible introduction.  Faraway takes a hands-on approach,
  using the \code{\link[=ozone1]{ozone}} data to compare \code{mda::mars} with other techniques.
  (If you use Faraway's examples with \code{earth} instead of \code{mars}, use \code{$bx}
  instead of \code{$x}.)
  Earth's pruning pass uses the \code{\link[leaps]{leaps}} package which is based on
  techniques in Miller.

  Faraway (2005) \emph{Extending the Linear Model with R}
  \url{http://www.maths.bath.ac.uk/~jjf23}

  Friedman (1991) \emph{Multivariate Adaptive Regression Splines (with discussion)}
  Annals of Statistics 19/1, 1--141

  Friedman (1993) \emph{Fast MARS}
  Stanford University Department of Statistics, Technical Report 110
  \url{http://www-stat.stanford.edu/research/index.html}

  Hastie, Tibshirani, and Friedman (2001) \emph{The Elements of Statistical Learning}
  \url{http://www-stat.stanford.edu/~hastie/pub.htm}

  Miller, Alan (1990, 2nd ed. 2002) \emph{Subset Selection in Regression}
}
\seealso{
  Start with \code{\link{summary.earth}} and \code{\link{plotmo}}.

  \code{\link{format.earth}},
  \code{\link{mars.to.earth}},
  \code{\link{model.matrix.earth}},
  \code{\link{ozone1}},
  \code{\link{plot.earth.models}},
  \code{\link{plot.earth}},
  \code{\link{plotmo}},
  \code{\link{predict.earth}},
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
#    GCV: 11          RSS: 213         GRSq: 0.96      RSq: 0.97
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

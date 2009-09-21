\name{earth}
\alias{earth}
\alias{earth.default}
\alias{earth.formula}
\alias{earth.fit}
\concept{regression}
\concept{mars}
\concept{Friedman}
\title{Multivariate Adaptive Regression Splines}
\description{
Build a regression model
using the techniques in Friedman's papers "Multivariate Adaptive Regression Splines"
and "Fast MARS".
}
\usage{
\method{earth}{formula}(formula = stop("no 'formula' arg"),
   data, weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
   na.action = na.fail, glm = NULL, trace = 0,
   keepxy = FALSE, nfold=0, stratify=TRUE, \dots)

\method{earth}{default}(x = stop("no 'x' arg"), y = stop("no 'y' arg"),
    weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
    na.action = na.fail, glm = NULL, trace = 0,
    keepxy = FALSE, nfold=0, stratify=TRUE, \dots)

\method{earth}{fit}(x = stop("no 'x' arg"), y = stop("no 'y' arg"),
    weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
    na.action = na.fail, glm = NULL, trace = 0,
    nk = max(21, 2 * ncol(x) + 1), degree = 1,
    penalty = if(degree > 1) 3 else 2, thresh = 0.001,
    minspan = 0, newvar.penalty = 0, fast.k = 20, fast.beta = 1,
    linpreds = FALSE, allowed = NULL,
    pmethod = "backward", nprune = NULL, Object = NULL,
    Get.crit = get.gcv, Eval.model.subsets = eval.model.subsets,
    Print.pruning.pass = print.pruning.pass, Force.xtx.prune = FALSE,
    Use.beta.cache = TRUE, \dots)
}
\arguments{
% % All arguments are optional except \code{formula}, or \code{x} and \code{y}.
% % NAs are not allowed.\cr
% \cr
To start off, look at the arguments
\code{formula},
\code{data},
\code{x},
\code{y},
\code{nk}, and
\code{degree}.
Many users will find that those arguments are all they need,
plus in some cases
\code{keepxy},
\code{nprune},
\code{penalty},
\code{minspan}, and
\code{trace}.
For GLM models, use the \code{glm} argument.
For cross validation, use the \code{nfold} argument.

  \item{formula}{
     Model formula.
  }
  \item{data}{
    Data frame for \code{formula}.
  }
  \item{x}{
     Matrix or dataframe containing the independent variables.
  }
  \item{y}{
     Vector containing the response variable, or, in the case of multiple responses,
     a matrix or dataframe whose columns are the values for each response.
  }
  \item{subset}{
     Index vector specifying which cases to use, i.e., which rows in \code{x} to use.
     Default is NULL, meaning all.
  }
  \item{weights}{
     Weight vector (not yet supported).
%      Vector of case weights.
%      Default is NULL, meaning no case weights.
%      If specified, \code{weights} must have length equal to \code{nrow(x)}
%      before applying \code{subset}.\cr
% %      You do not need to normalize \code{weights} --- earth does that.
% %      \code{wp} is normalized so that its elements sum to
% %      \code{nrow(x)} and all residuals calculated in the forward and pruning
% %      pass are weighted by the the corresponding element of \code{weights}.
% %      Zero values are converted to a small non zero value.
  }
  \item{wp}{
     Vector of response weights.
     Default is NULL, meaning no response weights.
     If specified, \code{wp} must have an element for each column of
     \code{y} (after \code{\link[=factor]{factors}}, if any, have been expanded).\cr
     Note for \code{mda::mars} users:
     earth's internal normalization of \code{wp} is different from \code{mars}.
     Earth uses \code{wp<-sqrt(wp/mean(wp))}
     and \code{mars} uses \code{wp<-sqrt(wp/sum(wp))}.
     Thus in earth, a \code{wp} with all elements equal is equivalent to no \code{wp}.
     For models built with \code{wp}, multiply the GCV
     calculated by \code{mars} by \code{length(wp)} to compare it to earth's GCV.
%      Earth uses \code{wp} as follows:
%      It normalizes \code{wp} as follows \code{wp <- sqrt(wp / mean(wp))}
%      It multiplies each column of \code{y} by the corresponding
%      element of the normalized \code{wp}.
%      It runs its internal MARS machine with the modified \code{y}.
%      Finally, after running the machine it divides each column of
%      \code{fitted.values}, \code{residuals}, and \code{coefficients} by the
%      corresponding element of the normalized \code{wp}.
  }
  \item{scale.y}{
    \code{\link[=scale]{Scale}} \code{y} in the forward pass for better numeric stability
    \bold{[added Jan 2009]}.
    Scaling here means subtract the mean and divide by the standard deviation.
    Default is \code{NCOL(y)==1},
    i.e., scale \code{y} unless \code{y} has multiple columns.
  }
  \item{na.action}{
     NA action. Default is \code{na.fail}, and only \code{na.fail} is supported.
     (Why? Because adding support to \code{earth} for other NA actions is easy, but making
     sure that is handled correctly internally in \code{predict}, \code{plotmo}, etc. is tricky.
     It is more reliable to make the user remove NAs before calling earth.)
  }
  \item{glm}{
     NULL (default) or a list of arguments to pass on to \code{\link{glm}}.
     See the documentation of \code{\link{glm}} for a description of these arguments
     (but not all of glm's arguments are supported, earth will give an error message
     if you use a glm argument it does not support).\cr
     See also the "Generalized linear models" section below. Example:\cr
     \code{earth(y~x, glm=list(family=binomial))}
  }
  \item{trace}{
     Trace earth's execution.  Default is 0.  Values:

     0 none .5 cross validation 1 overview 2 forward pass 3 pruning 4 model mats, memory use, more pruning, etc. 5 \dots
  }
  \item{keepxy}{
     Set to TRUE to retain the following in the returned value: \code{x} and \code{y} (or \code{data}),
     \code{subset}, and \code{weights}.
     Default is FALSE.
     The function \code{\link{update.earth}} and friends will use these instead of searching for them
     in the environment at the time \code{\link{update.earth}} is invoked.
     This argument also affects the amount of data kept when
     the \code{nfold} argument is used
     (see \code{cv.list} in the "Value" section below).\cr

     \bold{The following arguments are for the forward pass.}
  }
  \item{nk}{
     Maximum number of model terms before pruning, i.e., the
     maximum number of terms created by the forward pass.
     Includes the intercept.
     Default is \code{max(21,2*NCOL(x)+1)}.
     The number of terms created by the forward pass will be
     less than \code{nk} if a forward stopping condition is reached before \code{nk} terms,
     or if the forward pass drops one side of a hinge pair to prevent linear dependencies.
     See the "Forward pass" section below.
  }
  \item{degree}{
     Maximum degree of interaction (Friedman's \eqn{mi}).
     Default is 1, meaning build an additive model (i.e., no interaction terms).
  }
  \item{penalty}{
    Generalized Cross Validation (GCV) penalty per knot.
    Default is \code{if(degree>1) 3 else 2}.
    A value of 0 penalizes only terms, not knots.
    The value -1 is treated specially to mean no penalty, so GCV=RSS/n.
    Theory suggests values in the range of about 2 to 4.
    In practice, for big data sets larger values can be useful to force a smaller model.
    The FAQ section below has some information on GCVs.
  }
  \item{thresh}{
     Forward stepping threshold.
     Default is 0.001.
     This is one of the arguments used to decide when forward stepping
     should terminate.
     See the "Forward pass" section below.
  }
  \item{minspan}{
    Minimum distance between knots.\cr
    Use a value of 1 to consider all x values (which is good if the
    data are not noisy).\cr
    The default value is 0.
    The value 0 is treated specially and
    means calculate the \code{minspan} internally as per
    Friedman's MARS paper section 3.8 with \eqn{alpha} = 0.05.
    Set \code{trace>=2} to see the calculated value.\cr
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
     This argument can mitigate the effects of collinearity or concurvity
     in the input data, but anecdotal evidence is that it does not work very well.
     If you know two variables are strongly correlated then you
     would do better to delete one of them before calling earth.
  }
  \item{fast.k}{
     Maximum number of parent terms considered at each step of the forward pass.
     Friedman invented this parameter to speed up the forward pass
     (see the Fast MARS paper section 3.0).
     Default is 20.
     Values of 0 or less are treated specially
     (as being equivalent to infinity), meaning no Fast MARS.
     Typical values, apart from 0, are 20, 10, or 5.
     In general, with a lower \code{fast.k} (say 5), \code{earth} is faster;
     with a higher \code{fast.k}, or with \code{fast.k} disabled (set to 0),
     \code{earth} builds a better model.
     However it is not unusual to get a slightly better model with a lower \code{fast.k},
     and you may need to experiment.
  }
  \item{fast.beta}{
     Fast MARS aging coefficient, as described in the
     Fast MARS paper section 3.1.
     Default is 1.
     A value of 0 sometimes gives better results.
  }
  \item{linpreds}{
     Index vector specifying which predictors should enter linearly, as in \code{\link{lm}}.\cr
     The default is FALSE, meaning all predictors enter
     in the standard MARS fashion, i.e., in hinge functions.\cr
     A predictor's index in \code{linpreds} is the column number in the input matrix \code{x}
     after factors have been expanded.
     Examples are given in the FAQ section below.\cr
     Note: in the current implementation, the GCV penalty for predictors
     that enter linearly is the same as that for predictors with knots.
     That is not quite correct; linear terms should be penalized less.
  }
  \item{allowed}{
     Function specifying which predictors can interact and how.
     Default is NULL, meaning all standard MARS terms are allowed.\cr
     Earth calls the \code{allowed} function just before adding a term
     in the forward pass.
     If \code{allowed} returns TRUE the term goes into the model as usual;
     if \code{allowed} returns FALSE the term is discarded.
     Examples are given in the FAQ section below.\cr
     Your \code{allowed} function should have the following prototype\cr
     \code{allowed <- function(degree, pred, parents, namesx, first) \{.....\}}\cr
     where\cr
        \code{degree} is the interaction degree of the candidate term.
              Will be 1 for additive terms.\cr
        \code{pred} is the index of the candidate predictor.
             A predictor's index in \code{linpreds} is the column number in the input matrix \code{x}
             after factors have been expanded.\cr
        \code{parents} is the candidate parent term's row in \code{dirs}.\cr
        \code{namesx} is optional and if present is the column names of \code{x} after factors have been expanded.\cr
        \code{first} is optional and if present is TRUE the first time your \code{allowed} function is invoked for
        the current model, and thereafter FALSE.
  }
  \item{nfold}{
     Number of cross validation folds.
     Default is 0, i.e., no cross validation.
     If greater than 1, earth first builds a standard model as usual, with all the data.
     It then builds \code{nfold} cross validated models,
     measuring R-Squared on the left-out data each time.
     The final cross validation R-Squared is the mean of these R-Squareds.
     If a binomial or poisson model (see the \code{glm} argument),
     then further statistics are calculated.
     See the "Cross validation" section below for details.
  }
  \item{stratify}{
     Only applies if \code{nfold>1}.
     Default is TRUE.
     Stratify the cross validation samples so that, for each column of the response
     \code{y} (after factors have been expanded),
     an approximately equal number of cases with a non-zero response
     occur in each cross validation subset.
     That means that if \code{y} is a factor, there will be
     approximately equal numbers of each factor level in each fold
     (see the "Factors" section below).
     We say "approximately equal" because the number of occurrences of a factor
     level may not be exactly divisible by the number of folds.
\cr\cr
     \bold{The following arguments are for the pruning pass.}
  }
  \item{pmethod}{
     Pruning method.
     Default is \code{"backward"}.
     One of: \code{backward none exhaustive forward seqrep}.
     Use \code{none} to retain all the terms created by the forward pass.
     If \code{y} has multiple columns, then only \code{backward} or \code{none}
     is allowed.
     Pruning can take a while if \code{exhaustive} is chosen and
     the model is big (more than about 30 terms).
     The current version of the \code{\link[leaps]{leaps}} package
     used during pruning does not allow user interrupts
     (i.e., you have to kill your R session to interrupt;
     in Windows hit Ctl-Alt-Delete or from the command line use \code{tskill}).
  }
  \item{nprune}{
     Maximum number of terms (including intercept) in the pruned model.
     Default is NULL, meaning all terms created by the forward pass
     (but typically not all terms will remain after pruning).
     Use this to reduce exhaustive search time, or
     to enforce an upper bound on the model size.\cr

    \bold{The following arguments are for internal or advanced use.}
  }
  \item{Object}{
     Earth object to be updated, for use by \code{\link{update.earth}}.
  }
  \item{Get.crit}{
     Criterion function for model selection during pruning.
     By default a function that returns the GCV.
     See the "Pruning pass" section below.
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
    The \code{leaps} routines are more accurate but do not support multiple responses
    (\code{leaps} is based on the QR decomposition and
    \code{EvalSubsetsUsingXtx} is based on the inverse of X'X).
    Setting \code{Force.xtx.prune=TRUE} forces use of \code{EvalSubsetsUsingXtx}, even
    if \code{y} has a single column.
  }
  \item{Use.beta.cache}{
    Default is TRUE.
    Using the "beta cache" takes more memory but is faster
    (by 20\% and often much more for large models).
    The beta cache uses \code{nk * nk * ncol(x) * sizeof(double)} bytes.
    Set \code{Use.beta.cache=FALSE} to save memory.
    (The beta cache is an innovation in this implementation of MARS
    and does not appear in Friedman's papers. It is not related to
    the \code{fast.beta} argument.)
  }
  \item{\dots}{
     Dots are passed on to \code{earth.fit}.
 }
}
\value{
  An object of class "earth" which is a list with the components
  listed below.
  \emph{Term} refers to a term created during the forward pass
  (each line of the output from \code{\link{format.earth}} is a term, and
  \code{dirs} defines which predictors are in which terms).
  Term number 1 is always the intercept.
  \item{rss}{
     Residual sum-of-squares (RSS) of the model (summed over all responses
     if \code{y} has multiple columns).
  }
  \item{rsq}{
    \code{1-rss/rss.null}.
     R-Squared of the model (calculated over all responses).
     A measure of how well the model fits the training data.
     Note that \code{rss.null} is \code{sum((y - mean(y))^2)}.
  }
  \item{gcv}{
     Generalized Cross Validation (GCV) of the model (summed over all responses)
     The GCV is calculated using the \code{penalty} argument.
     For details of the GCV calculation, see
     equation 30 in Friedman's MARS paper and \code{earth:::get.gcv}.
  }
  \item{grsq}{
     \code{1-gcv/gcv.null}.
     An estimate of the predictive power of the model (calculated over all responses).
     Unlike \code{rsq}, in MARS models \code{grsq} can be negative.
     A negative \code{grsq} would indicate
     a severely over parameterized model --- a model that
     would not generalize well
     even though it may be a good fit to the training data.
     Watch the GRSq take a nose dive in this example:\cr
     \code{earth(mpg~., data=mtcars, pmethod="none", trace=3)}
  }
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
% \cr
}
  \item{dirs}{
     Matrix with one row per MARS term, and with with ij-th element equal to\cr

\code{0}  if predictor j is not in term i\cr
\code{-1} if an expression of the form \code{pmax(const - xj)} is in term i\cr
\code{1}  if an expression of the form \code{pmax(xj - const)} is in term i\cr
\code{2}  if predictor j enters term i linearly.\cr

     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Example \code{dirs}:\preformatted{                       Girth Height
(Intercept)                0  0 #intercept
h(Girth-12.9)              1  0 #2nd term uses Girth
h(12.9-Girth)             -1  0 #3rd term uses Girth
h(Girth-12.9)*h(Height-76) 1  1 #4th term uses Girth and Height
...
}
% \cr
  }
  \item{cuts}{
     Matrix with ij-th element equal to the cut point
     for predictor j in term i.
     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that the terms may not be in pairs, because the forward pass
     deletes linearly dependent terms before handing control to the pruning pass.
     Note for programmers: the precedent is to use \code{dirs}
     for term names etc. and to only use \code{cuts} where cut information needed.
     Example \code{cuts}:\preformatted{                           Girth Height
(Intercept)                    0   0  #intercept, no cuts
h(Girth-12.9)               12.9   0  #2nd term has cut at 12.9
h(12.9-Girth)               12.9   0  #3rd term has cut at 12.9
h(Girth-12.9)*h(Height-76)  12.9  76  #4th term has two cuts
...}
% \cr
  }
  \item{selected.terms}{
     Vector of term numbers in the best model.
     Can be used as a row index vector into \code{cuts} and \code{dirs}.
     The first element \code{selected.terms[1]} is always 1, the intercept.
  }
  \item{prune.terms}{
     A matrix specifying which terms appear in which pruning pass subsets.
     The row index of \code{prune.terms} is the model size.
     (The model size is the number of terms in the model.
     The intercept is counted as a term.)
     Each row is a vector of term numbers for the best model of that size.
     An element is 0 if the term is not in the model, thus \code{prune.terms} is a
     lower triangular matrix, with dimensions \code{nprune x nprune}.
     The model selected by the pruning pass is at row number \code{length(selected.terms)}.
     Example \code{prune.terms}:\preformatted{[1,]    1  0  0  0  0  0  0 #intercept-only model
[2,]    1  2  0  0  0  0  0 #best 2 term model uses terms 1,2
[3,]    1  2  4  0  0  0  0 #best 3 term model uses terms 1,2,4
[4,]    1  2  6  9  0  0  0 #and so on
...}
% \cr
  }
  \item{rss.per.response}{
    A vector of the RSS for each response.
    Length is the number of responses, i.e., \code{ncol(y)} after factors in \code{y} have been expanded.
    The \code{rss} component above is  equal to \code{sum(rss.per.response)}.
  }
  \item{rsq.per.response}{
    A vector of the R-Squared for each response.
    Length is the number of responses.
  }
  \item{gcv.per.response}{
    A vector of the GCV for each response.
    Length is the number of responses.
    The \code{gcv} component above is equal to \code{sum(gcv.per.response)}.
  }
  \item{grsq.per.response}{
    A vector of the GRSq for each response.
    Length is the number of responses.
  }
  \item{rss.per.subset}{
     A vector of the RSS for each model subset generated by the pruning pass.
     Length is \code{nprune}.
     For multiple responses, the RSS is summed over all responses for each subset.
     The null RSS (i.e., the RSS of an intercept only-model) is \code{rss.per.subset[1]}.
     The \code{rss} above is\cr
     \code{rss.per.subset[length(selected.terms)]}.
    }
  \item{gcv.per.subset}{
     A vector of the GCV for each model in \code{prune.terms}.
     Length is \code{nprune}.
     For multiple responses, the GCV is summed over all responses for each subset.
     The null GCV (i.e., the GCV of an intercept-only model) is \code{gcv.per.subset[1]}.
     The \code{gcv} above is \code{gcv.per.subset[length(selected.terms)]}.
  }
  \item{fitted.values}{
    Fitted values.
    A matrix with dimensions \code{nrow(y)} x \code{ncol(y)}
    after factors in \code{y} have been expanded.
  }
  \item{residuals}{
    Residuals.
    A matrix with dimensions \code{nrow(y)} x \code{ncol(y)}
    after factors in \code{y} have been expanded.
  }
  \item{coefficients}{
    Regression coefficients.
    A matrix with dimensions \code{length(selected.terms)} x \code{ncol(y)}
    after factors in \code{y} have been expanded.
    Each column holds the least squares coefficients from regressing that
    column of \code{y} on \code{bx}.
    The first row holds the intercept coefficient(s).
  }
  \item{penalty}{
     The GCV penalty used during pruning.
     A copy of \code{earth}'s \code{penalty} argument.
  }
  \item{call}{
     The call used to invoke \code{earth}.
  }
  \item{terms}{
     Model frame terms.
     This component exists only if the model was built using \code{earth.formula}.
  }
  \item{namesx}{
     Column names of \code{x}, generated internally by earth when necessary
     so each column of \code{x} has a name.
     Used, for example, by \code{\link{predict.earth}} to name columns if necessary.
  }
  \item{namesx.org}{
     Original column names of \code{x}.
  }
  \item{levels}{
     Levels of \code{y} if \code{y} is a \code{\link{factor}}\cr
     \code{c(FALSE,TRUE)} if \code{y} is \code{\link{logical}}\cr
     Else NULL
  }
  \item{wp}{
     Copy of the \code{wp} argument to earth.\cr\cr
  \bold{The following fields appear only if earth's argument \code{keepxy} is TRUE.}
  }
  \item{x}{}
  \item{y}{}
  \item{data}{}
  \item{subset}{}
  \item{weights}{
     Copies of the corresponding arguments to earth.
     Only exist if \code{keepxy=TRUE}.\cr\cr
  \bold{The following fields appear only if earth's \code{glm} argument is used.}
  }
  \item{glm.list}{
     List of GLM models. Each element is the value returned by earth's
     internal call to \code{\link{glm}} for each response.\cr
     Thus if there is a single response (or a single binomial pair, see
     the "Binomial pairs" section below)
     this will be a one element list and you access the GLM model with
     \code{my.earth.model$glm.list[[1]]}.
  }
  \item{glm.coefficients}{
    GLM regression coefficients.
    Analogous to the \code{coefficients} field described above but for the GLM model(s).
    A matrix with dimensions \code{length(selected.terms)} x \code{ncol(y)}
    after factors in \code{y} have been expanded.
    Each column holds the coefficients from the GLM regression of that
    column of \code{y} on \code{bx}.
    This duplicates, for convenience, information buried in \code{glm.list}.
    }
   \item{glm.bpairs}{
       NULL unless there are paired binomial columns.
       A logical vector, derived internally by earth, or a copy
       the \code{bpairs} specified by the user in the \code{glm} list.
       See the "Binomial pairs" section below.\cr\cr
   \bold{The following fields appear only if the \code{nfold} argument is greater than 1.}
   }
   \item{cv.rsq.tab}{
      Matrix with \code{nfold+1} rows and \code{nresponse+1} columns,
      where \code{nresponse} is the number of responses,
      i.e., \code{ncol(y)} after factors in \code{y} have been expanded.
      The first \code{nresponse} elements of a row are the RSq's on
      the left-out data for each response of the model generated at that row's fold.
      The final column holds the row mean (a weighted mean if \code{wp} if specified).
      The final row of the table holds the column means.
      The values in this final row are the CV-RSqs printed by \link{summary.earth}.
\cr\cr
      Example for a single response model:
\preformatted{            y  mean
fold 1  0.909 0.909
fold 2  0.869 0.869
fold 3  0.952 0.952
fold 4  0.157 0.157
fold 5  0.961 0.961
mean    0.769 0.769
}
      Example for a multiple response model:
\preformatted{          y1   y2    y3  mean
fold 1 0.915 0.951 0.944 0.937
fold 2 0.962 0.970 0.970 0.968
fold 3 0.914 0.940 0.942 0.932
fold 4 0.907 0.929 0.925 0.920
fold 5 0.947 0.987 0.979 0.971
mean   0.929 0.955 0.952 0.946
}
   }
   \item{cv.maxerr.tab}{
      Like \code{cv.rsq.tab} but is the \code{MaxErr} at each fold.
      This is the signed max absolute value at each fold.
      Also, results are aggregrated for the final column and final row
      using the signed max absolute value instead of the mean.
      The \emph{signed max absolute value} is defined here as the maximum of the absolute difference
      between the predicted and observed response values, multiplied
      by \code{-1} if the sign of the difference is negative.
   }
   \item{cv.deviance.tab}{
      Like \code{cv.rsq.tab} but is the \code{MeanDev} at each fold.
      Binomial models only.
   }
   \item{cv.calib.int.tab}{
      Like \code{cv.rsq.tab} but is the \code{CalibInt} at each fold.
      Binomial models only.
   }
   \item{cv.calib.slope.tab}{
      Like \code{cv.rsq.tab} but is the \code{CalibSlope} at each fold.
      Binomial models only.
   }
   \item{cv.auc.tab}{
      Like \code{cv.rsq.tab} but is the \code{AUC} at each fold.
      Binomial models only.
   }
   \item{cv.cor.tab}{
      Like \code{cv.rsq.tab} but is the \code{cor} at each fold.
      Poisson models only.
   }
   \item{cv.nterms}{
      Vector of length \code{nfold}+1.
      Number of MARS terms in the model generated at each cross validation fold,
      with the final element being the mean of these.
   }
   \item{cv.nvars}{
      Vector of length \code{nfold}+1.
      Number of predictors in the model generated at each cross validation fold,
      with the final element being the mean of these.
   }
   \item{cv.groups}{
     Specifies which cases went into which folds.
     Vector of length equal to the number of cases, with elements taking values in \code{1:nfold}.
   }
   \item{cv.list}{
      List of earth models, one model for each fold.
      These fold models have extra fields
      \code{cv.rsq} and \code{cv.rsq.per.response}
      (and, if \code{keepxy} is set, also \code{cv.test.y} and \code{cv.test.fitted.values}).
      To save memory, lengthy fields
      in the fold models are removed, unless you use the \code{keepxy} argument.
      The "lengthy fields" are \code{$bx}, \code{$fitted.values}, and \code{$residuals}.
   }
}
\note{
Many users will find that it is unnecessary to read this entire section.
Just read the parts you need and skim the rest.
\cr\cr
\bold{Contents}\cr\cr
. Other implementations\cr
. Limitations\cr
. Multiple response models\cr
. Generalized linear models\cr
. Factors\cr
. Binomial pairs\cr
. The forward pass\cr
. The pruning pass\cr
. Execution time\cr
. Memory use\cr
. Cross validation\cr
. Cross validating binomial and poisson models\cr
. Using earth with fda and mda\cr
. Migrating from mda::mars\cr
. Standard model functions\cr
. Frequently asked questions

\bold{Other implementations}

The results are similar to but not identical to other
Multivariate Adaptive Regression Splines implementations.
The differences stem from the forward pass where very small
implementation differences (or perturbations of the input data) can cause
rather different selection of terms and knots (although similar GRSq's).
The backward passes give identical or near identical results,
given the same forward pass results.

The source code of \code{earth} is derived from \code{\link[mda]{mars}}
in the \code{mda} package written by
by Trevor Hastie and Robert Tibshirani.
See also \code{\link{mars.to.earth}}.

The term "MARS" is trademarked and licensed exclusively to
Salford Systems \url{http://www.salfordsystems.com}.
Their implementation uses an engine written by Friedman and
has some features not in earth.

StatSoft also have an implementation which they call MARSplines
\url{http://www.statsoft.com/textbook/stmars.html}.

\bold{Limitations}

The following aspects of MARS
are mentioned in Friedman's papers but not implemented in \code{earth}:

i) Piecewise cubic models\cr
ii) Model slicing (\code{\link{plotmo}} goes part way)\cr
iii) Handling missing values\cr
iv) Automatic grouping of categorical predictors into subsets\cr
v) Fast MARS h parameter

\bold{Multiple response models}

If \code{y} has \code{k} columns
then earth builds \code{k} simultaneous models.
(Note that \code{y} will have multiple columns
if a factor in \code{y} is expanded by earth;
see the "Factors" section below for details.)
Each model has the same set of basis functions
(the same \code{bx}, \code{selected.terms}, \code{dirs} and \code{cuts})
but different coefficients (the returned \code{coefficients} will have \code{k} columns).
The models are built and pruned as usual but with the GCVs
and RSSs averaged across all \code{k} responses.

Since earth attempts to optimize for all models simultaneously,
the results will not be as "good" as building the models independently,
i.e., the GCV of the combined model will usually not be as good as
the GCVs for independently built models.
However, the combined model may be a better model in other senses,
depending on what you are trying to achieve.
For example, it could be useful for earth to select
the set of MARS terms that is best across \emph{all} responses.
This would typically be the case in a multiple response logistic model
if some responses have a very small number of successes.

Note that automatic scaling of \code{y} (via the \code{scale.y} argument)
does not take place if \code{y} has multiple columns.
You may want to scale your \code{y} columns before calling earth
so each \code{y} column gets the appropriate weight during model building
(a \code{y} column with a big variance will influence the
model more than a column with a small variance).
You could do this by calling \code{\link{scale}} before invoking \code{earth},
or by setting the \code{scale.y} argument, or by using the
the \code{wp} argument.

Here are a couple of (artificial) examples to show some of the ways
multiple responses can be specified.
Note that \code{\link[=data.frame]{data.frames}} can't be used on the left side
of a formula, so \code{\link{cbind}} is used in the first example.
The examples use the standard technique of specifying
a tag \code{lvol=} to name a column.
\preformatted{
earth(cbind(Volume,lvol=log(Volume)) ~ ., data=trees)
attach(trees)
earth(data.frame(Girth,Height), data.frame(Volume,lvol=log(Volume)))
}
Don't use a plus sign on the left side of the tilde.
You might think that specifies a multiple response, but instead
it arithmetically adds the columns.

For more details on using residual errors averaged over multiple responses see section 4.1 of
Hastie, Tibshirani, and Buja \emph{Flexible Discriminant Analysis by Optimal Scoring},
JASA, December 1994 \url{http://www-stat.stanford.edu/~hastie/Papers/fda.pdf}.

\bold{Generalized linear models}

Earth builds a GLM model if the \code{glm} argument is specified.
Earth builds the model as usual and then invokes
\code{\link{glm}} on the MARS basis matrix \code{bx}.

In more detail, the model is built as follows.
Earth first builds a standard MARS model, including
the internal call to \code{\link{lm.fit}} on \code{bx} after the pruning pass.
(See "The forward pass" and "The pruning pass" sections below).
Thus knot positions and terms are determined as usual and
all the standard fields in earth's return value will be present.
Earth then invokes \code{glm} for the response on \code{bx}
with the parameters specified in the \code{glm} argument to earth.
For multiple response models
(when y has multiple columns), the call to \code{glm} is repeated independently for each response.
The results go into three extra fields in earth's return value:
\code{glm.list}, \code{glm.coefficients}, and \code{glm.bpairs}.

Earth's internal call to \code{glm} is made with
the \code{glm} arguments \code{x}, \code{y}, and \code{model} set TRUE
(see the documentation for \code{\link{glm}} for more information about those arguments).
% [TODO Would like to make this depend on \code{keepxy} but currently \emph{have} to
% keep \code{x} and \code{y} else later calls to \code{plot.glm} on the GLM models fail.]

Use \code{summary(my.model)} as usual to see the model.
Use \code{summary(my.model, details=TRUE)} to see more details, but
note that the printed P-values of the GLM coefficients are
meaningless.
This is because of the amount of preprocessing by earth ---
the mantra is "variable selection overstates significance of the selected variables".
Use \code{plot(my.model$glm.list[[1]])} to plot the (first) \code{glm} model.

The examples below show how to specify earth-glm models.
The examples are only to illustrate the syntax and not necessarily useful models.
In the examples \code{pmethod="none"}, otherwise with these artificial
models earth tends to prune away everything except the intercept term.
You wouldn't normally use \code{pmethod="none"}.
Also, \code{trace=1}, so if you run these examples you can
see how earth expands the input matrices
(as explained in the "Factors" and "Binomial pairs" sections below).

1. \emph{Two-level factor or logical response.}
The response is converted to a single column of 1s and 0s.
\preformatted{
a1 <- earth(survived ~ ., data=etitanic, # c.f. Harrell "Reg Mod Strat" ch. 12
             degree=2, trace=1,
             glm=list(family=binomial))

a1a <- earth(etitanic[,-2], etitanic[,2], # equivalent but using earth.default
             degree=2, trace=1,
             glm=list(family=binomial))
}
2. \emph{Factor response.}
This example is for a factor with more than two levels.
(For factors with just two levels, see the previous example.)
The factor \code{pclass} is expanded to three indicator columns
(whereas in a direct call to \code{\link{glm}}, \code{pclass} would be treated
as logical: the first level versus all other levels).
\preformatted{
a2 <- earth(pclass ~ ., data=etitanic, trace=1,
            glm=list(family=binomial))
}
3. \emph{Binomial model} specified with a column pair.
This is a single response model but specified with a pair of columns:
see the "Binomial pairs" section below.
For variety, this example uses a \code{probit} link and (unnecessarily) increases \code{maxit}.
\preformatted{
ldose <- rep(0:5, 2) - 2 # V&R 4th ed. p. 191
sex <- factor(rep(c("male", "female"), times=c(6,6)))
numdead <- c(1,4,9,13,18,20,0,2,6,10,12,16)
pair <- cbind(numdead, numalive=20 - numdead)
a3 <- earth(pair ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family=binomial(link=probit), maxit=100))
}
4. \emph{Double binomial response} (i.e., a multiple response model)
specified with two column pairs.
\preformatted{
numdead2 <- c(2,8,11,12,20,23,0,4,6,16,12,14) # bogus data
doublepair <- cbind(numdead, numalive=20-numdead,
                    numdead2=numdead2, numalive2=30-numdead2)
a4 <- earth(doublepair ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family="binomial"))
}
5. \emph{Poisson model.}
\preformatted{
counts <- c(18,17,15,20,10,20,25,13,12) # Dobson 1990 p. 93
outcome <- gl(3,1,9)
treatment <- gl(3,3)
a5 <- earth(counts ~ outcome + treatment, trace=1, pmethod="none",
            glm=list(family=poisson))
}
6. \emph{Standard earth model}, the long way.
\preformatted{
a6 <- earth(numdead ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family=gaussian(link=identity)))
print(a6$coefficients == a6$glm.coefficients)  # all TRUE
}

% [TODO investigate offsets in earth-glm models.]

\bold{Factors}

\bold{Factors in x}:
Earth treats factors in \code{x} in
the same way as standard R models such as \code{\link{lm}} (where
\code{x} is taken to mean the right hand side of the formula).
Thus factors are expanded using the current setting
of \code{\link{options}("\link{contrasts}")}.

\bold{Factors in y}:
Earth treats factors in the response in a non-standard way that makes use
of earth's ability to handle multiple responses.
A \emph{two level factor} (or logical) is converted to a single indicator column of 1s and 0s.
A \emph{factor with three or more levels}
is converted into \code{k} indicator columns of 1s and 0s, where \code{k} is the number of levels
(the \code{\link{contrasts}} matrix is diagonal, see \code{\link{contr.earth.response}}).
This happens regardless of the \code{options("contrasts")}
setting and regardless of whether the factors are ordered or unordered.
For example, if a column in \code{y} is a factor with levels
\code{A}, \code{B}, and \code{C},
the column will be expanded to three columns like this
(the actual data will vary but each row will have a single 1):
\preformatted{
A B C  # one column for each factor level
0 1 0  # each row has a single 1
1 0 0
0 0 1
0 0 1
0 0 1
...
}
In distinction, a standard \code{\link[=contrasts]{treatment contrast}}
on the rhs of a model with an intercept would have no first "A" column
(to prevent linear dependencies on the rhs of the model formula).

This expansion to multiple columns (which only happens for factors with more than two levels)
means that earth will build a multiple response model
as described in the "Multiple responses" section above.

% [TODO ordered factors should be treated differently?]\cr
% [TODO should be able to pass a y contrasts argument to earth --- overkill?]

Paired binomial response columns in \code{y} are treated specially
--- see the "Binomial pairs" section below.

Use \code{trace=1} or higher to see the column names of the \code{x} and \code{y}
matrices after factor processing.
Use \code{trace=4} to see the first few rows of \code{x} and \code{y}
after factor processing.

Here is an example which uses the \code{\link{etitanic}} data to
predict the passenger class (not necessarily a sensible thing to
do but provides a good example here):
\preformatted{
> data(etitanic)
> head(etitanic) # pclass and sex are unordered factors

  pclass survived    sex    age sibsp parch
1    1st        1 female 29.000     0     0
2    1st        1   male  0.917     1     2
3    1st        0 female  2.000     1     2

> earth(pclass ~ ., data=etitanic, trace=1) # note col names in x and y below

x is a 1046 by 5 matrix: 1=survived, 2=sexmale, 3=age, 4=sibsp, 5=parch
y is a 1046 by 3 matrix: 1=1st, 2=2nd, 3=3rd
rest not shown here...
}
% Note that earth does not automatically group categorical variables into subsets as
% described in Friedman 1991 \emph{Estimating Functions Of Mixed Ordinal and Categorical
% Variables Using Adaptive Splines}.

\bold{Binomial pairs}

This section is only relevant if you use earth's \code{glm} argument
with a binomial or quasibinomial \code{\link{family}}.

Users of the \code{\link{glm}} function will be familiar with
the technique of specifying a binomial response as a two-column matrix,
with a column for the number of successes and a column for the failures.
Earth automatically detects when such columns are present in \code{y}
(by looking for adjacent columns which both have entries greater than 1).
The first column only is used to build the standard earth model.
Both columns are then passed to earth's internal call to \code{glm}.
As always, use \code{trace=1} to see how the columns
of \code{x} and \code{y} are expanded.

You can override this automatic detection by including a \code{bpairs} parameter.
This is usually (always?) unnecessary. For example
\preformatted{
glm=list(family=binomial, bpairs=c(TRUE, FALSE))
}
specifies that there are two columns in the response with the second
paired with the first.
These examples
\preformatted{
glm=list(family=binomial, bpairs=c(TRUE, FALSE, TRUE, FALSE))
glm=list(family=binomial, bpairs=c(1,3)) # equivalent
}
specify that the 1st and 2nd columns are a binomial pair
and the 3rd and 4th columns another binomial pair.

\bold{The forward pass}

Understanding the details of the forward and pruning passes
will help you understand earth's return value and
the admittedly large number of arguments.
The result of the forward pass is the MARS basis matrix \code{bx} and
the set of terms defined by \code{dirs}
and \code{cuts} (these are all fields in earth's return value,
but the \code{bx} here includes all terms before trimming
back to \code{selected.terms}).

The forward pass adds terms in pairs until the first of the
following conditions is met:

i) reach maximum number of terms \code{(nterms >= nk)}\cr
ii) reach DeltaRSq threshold \code{(DeltaRSq < thresh)}, where
DeltaRSq is the difference in R-Squared caused by adding the current term pair,
and \code{thresh} is the argument to earth\cr
iii) reach max RSq \code{(RSq > 1-thresh)}\cr
iv) reach min GRSq \code{(GRSq < -10)} (-10 is a pathologically bad GRSq)\cr
v) no new term increases the RSq (reached numerical limits).

Set \code{trace>=1} to see the stopping condition and
\code{trace>=2} to trace the forward pass.

You can disable all termination conditions
except (i) and (v) by setting \code{thresh=0}.
See the FAQ below "Why do I get fewer terms than nk?".

Note that GCVs (via GRSq) are used during the forward pass only as one of the
(more unusual) stopping conditions and in \code{trace} prints.
Changing the \code{penalty} argument does not change the knot positions.

The various stopping conditions mean that the actual number of terms
created by the forward pass may be less than \code{nk}.
There are other
reasons why the actual number of terms may be less than \code{nk}:
(i) the forward pass discards one side of a term pair
if it adds nothing to the model ---
but the forward pass counts terms as if they were actually created in pairs,
and, (ii) as a final step, the forward pass deletes linearly dependent terms, if any,
so all terms in \code{dirs} and \code{cuts} are independent.
And remember that the pruning pass will further discard terms.

\bold{The pruning pass}

The pruning pass is handed the sets of terms created by the forward pass.
Its job is to find the subset of those terms that gives the lowest GCV.
The following description of the pruning pass
explains how various fields in earth's returned value are generated.

The pruning pass works like this:
it determines the subset of terms in \code{bx} (using \code{pmethod})
with the lowest RSS for each model size in \code{1:nprune}
(see the \code{Force.xtx.prune} argument above for some details).
It saves the RSS and term numbers for each such subset in \code{rss.per.subset}
and \code{prune.terms}.
It then applies the \code{Get.crit} function with \code{penalty}
to \code{rss.per.subset} to yield \code{gcv.per.subset}.
Finally it chooses the model with the lowest value in \code{gcv.per.subset},
puts its term numbers into \code{selected.terms},
and updates \code{bx} by keeping only the \code{selected.terms}.

After the pruning pass, earth runs \code{\link{lm.fit}} to determine the
\code{fitted.values}, \code{residuals}, and \code{coefficients}, by regressing
the response \code{y} on \code{bx}.
If \code{y} has multiple columns then \code{lm.fit} is called
for each column.

If a \code{glm} argument is passed to earth,
earth runs \code{\link{glm}} on (each column of) \code{y}
in addition to the above call to \code{lm.fit}.

Set \code{trace>=3} to trace the pruning pass.

By default \code{Get.crit} is \code{earth:::get.gcv}.
Alternative \code{Get.crit} functions can be defined.
See the source code of \code{get.gcv} for an example.

\bold{Execution time: "I wanna go fast"}

For a given set of input data,
the following can increase the speed of the forward pass:

i) decreasing \code{fast.k}\cr
ii) decreasing \code{nk} (because fewer forward pass terms)\cr
iii) decreasing \code{degree}\cr
iv) increasing \code{thresh} (because fewer forward pass terms)\cr
v) increasing \code{min.span}.

The backward pass is normally much faster than the forward pass,
unless \code{pmethod="exhaustive"}.
Reducing \code{nprune} reduces exhaustive search time.
One strategy is to first build a large model
and then adjust pruning parameters such as \code{nprune} using \code{\link{update.earth}}.

The following very rough rules of thumb apply for large models.
Using \code{minspan=1} instead of the default \code{0} will increase times by 20 to 50\%.
Using \code{fast.k=5} instead of the default \code{20} can give substantial speed gains
but will sometimes give a much smaller \code{GRSq}.
Using an \code{allowed} function slows down model building by about 10\%.

\bold{Memory use}

Earth does not impose specific limits on the model size.
Model size is limited only by the amount of memory on your system,
the maximum memory addressable by R, and your patience.
On a 32 bit machine with \code{x} and \code{y} of type
double (no factors),
the number of bytes of memory used by earth is about
\preformatted{
8 * (nk^2 * ncol(x) + (nrow(x) * (3 + 2*nk + ncol(x)/2))).
}
Earth prints the results of the above calculation if \code{trace>=4}.
Memory use peaks in the forward pass.
The bulk of the forward pass is implemented in C.
It allocates memory "outside of R" and so
\link{memory.size} will not report the memory it uses.

Before calling \code{earth}, R itself will of course allocate
memory over and above the amount calculated above.
To reduce total memory usage, it sometimes helps to \code{\link{remove}} variables
and call \link{gc} before invoking \code{earth}.

Earth uses more memory if any elements of the
\code{x} and \code{y} arguments are not \link{double},
because it must convert them to double internally.
The same applies if the \code{subset} argument is used.
Earth uses more memory if \code{trace>=2}
(because \code{DUP=TRUE} is required to pass predictor names
to earth's internal call to \code{\link{.C}}).
Increasing the \code{degree} does not change the memory requirement
but greatly increases the running time.

Here is an example of memory use:
the earth test suite builds a model using \code{earth.default}
with a 1e4 by 100 input matrix with \code{nk=21}.
The Windows XP task manager
reports that the peak memory use when building this model is 47 MBytes.
Using the formula interface to earth pushes memory to 62 MBytes.
Increasing the number of rows in the input matrix to 1e5 pushes memory to 240 MBytes.

\bold{Cross validation}

Use cross validation to get an estimate of RSq on independent data.
Example (note the \code{nfold} parameter):
\preformatted{
a <- earth(survived ~ ., data=etitanic, degree=2, nfold=10)
summary(a) # note the CV-RSq field
}
Cross validation is done if \code{nfold} is greater than 1 (typically 10).
Earth first builds a standard model with all the data as usual.
This means that all the standard fields in earth's return value
appear as usual.
Earth then builds \code{nfold} cross validated models.
It measures RSq on the test data (i.e., the left-out data) for each fold.
The final cross validation RSq is the mean of these RSq's.
Use \link{summary.earth} to see this final value and its standard deviation
across cross validation folds.

The cross validation results go into extra fields in earth's return value.
All of these have a \code{cv} prefix ---
see the "Value" section above for details.

For multiple response models, at each fold earth calculates the RSq for each
response independently, and combines these by taking their mean
(or weighted mean if the \code{wp} argument is used).

With \code{trace=.5} or higher,
earth prints out progress information as cross validation proceeds.
For example
\preformatted{
CV fold 3: CV-RSq 0.622  ntrain-nz 384 ntest-nz 43
}
shows that for cross validation fold number 3,
the RSq on the test set (i.e., the left-out data) is 0.622.
The printout also shows the number of non-zero values in the observed response
in the fold's training set and test set.
This is useful if you have a binary or factor response and want to check
that you have enough examples of each factor level in each fold.
With the \code{stratify} argument (which is enabled by default),
earth attempts to keep the numbers of each level
constant across folds.

For reproducibility, call \link{set.seed} before calling earth with \code{nfold}.

\bold{Cross validating binomial and poisson models}

If you cross validate a binomial or poisson model
(specified using earth's \code{nfold} and \code{glm} arguments),
earth returns the following additional statistics,
Each of these is measured on the test set for each fold,
and averaged across all folds
(except that the signed max absolute value instead of the average is used for \code{MaxErr}).
Use \link{summary.earth} to see these statistics
and their standard deviation across folds.

\bold{CV-RSq} cross validated R-Squared, identical to CV-RSq for non-glm models

\bold{MaxErr} signed max absolute difference between the predicted and observed response.
  This is the maximum of the absolute differences
  between the predicted and observed response values, multiplied
  by \code{-1} if the sign of the difference is negative.

\bold{MeanDev} deviation divided by the response length

\bold{CalibInt} \bold{CalibSlope}
calibration intercept and slope (from regressing the observed response
on the predicted response)

\bold{AUC} (binomial models only) area under the ROC curve

\bold{cor} (poisson models only) correlation between the predicted and observed response\cr

See the source code in \code{earth.cv.R} for details.

For multiple response models, at each fold earth calculates these statistics for each
response independently, and combines them by taking their mean,
or weighted mean if the \code{wp} argument is used
(but takes the signed max absolute value instead of the mean for \code{MaxErr})
[TODO should do the same for \code{CalibInt}, \code{CalibSlope}?].
Taking the mean is a rather dubious way of combining results from
what are essentially quite different models,
but can nevertheless be useful.

Explanations of the above statistics can be found in the following
(and many other) references:

T. Fawcett (2004)
\emph{ROC Graphs: Notes and Practical Considerations for Researchers}.
Revised version of Technical report HP Laboratories.
\url{http://home.comcast.net/~tom.fawcett/public_html/papers}

J. Pearce and S. Ferrier (2000)
\emph{Evaluating the predictive performance
of habitat models developed using logistic regression}

F. Harrell (2001)
\emph{Regression Modeling Strategies with Applications to Linear Models,
Logistic Regression, and Survival Analysis}
\url{http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/RmS}

% \bold{Balancing prediction error in multiple response models}
%
% If the number of occurrences of one class is much larger than the others,
% then earth, trying to minimize the overall error, will keep the
% error rate low on the large class while letting the error on
% the small classes increase.
% You can see when this happens by looking at the relative error rate for
% each class.
%
% To deal with this issue, use the \code{wp} argument to
% weight important but seldom occuring classes.
% TODO next sentence is misleading if y is a two level factor
% (Remember that after earth expands a factor with more than two levels,
% \code{y} will have have a column for for each class, and these
% are the columns weighted by \code{wp}.)
% Start off by making the weights inversely proportional to
% the class populations, build the model,
% look at the error rates, adjust the weights, and iterate.
%
% For more information (in a Random Forests setting) see
% Breiman and Cutler \emph{Random Forests}
% \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#balance}.

\bold{Using earth with fda and mda}

\code{Earth} can be used with \code{\link[mda]{fda}} and \code{\link[mda]{mda}}
in the \code{mda} package.  Earth will generate a multiple response model.
Use the \code{fda/mda} argument \code{keep.fitted=TRUE} if
you want to call \code{\link{plot.earth}} later
(actually only necessary for large datasets, see the description of
\code{keep.fitted} in \code{fda}).
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

\bold{Migrating from mda::mars}

Changing code from \code{mda::mars} to \code{earth}
is usually just a matter
of changing the call from "\code{mars}" to "\code{earth}".
But there are a few argument differences and
earth will issue a warning if you give it a \code{mars}-only argument.

The resulting model will be similar but not identical because
of small implementation differences which are magnified
by the inherent instability of the MARS forward pass.

If you are further processing the output of earth you will need to
consider differences in the returned value.  The header of the
source file \code{mars.to.earth.R} describes these.
Perhaps the most important is that \code{mars} returns the
MARS basis matrix in a field called "\code{x}"
whereas \code{earth} returns "\code{bx}".
Also, \code{earth} returns "\code{dirs}" rather than "\code{factors}",
and in \code{earth} this matrix can have entries of value 2 for linear predictors.

See also \link{mars.to.earth}.

\bold{Standard model functions}

Standard model functions such as \code{\link{case.names}}
are provided for \code{earth} objects and are not explicitly documented.
Many of these give warnings when the results are not what you may expect.
Pass \code{warn=FALSE} to these functions to turn of just these warnings.\cr\cr

\bold{FREQUENTLY ASKED QUESTIONS}\cr\cr
\bold{What are your plans for earth?}

We would like to add support of case weights (to allow boosting),
but that won't happen anytime soon.

\bold{How can I establish variable importance?}

Use the \code{\link{evimp}} function.
See its help page for more details.

The \code{\link{summary.earth}} function lists the predictors
in order of estimated importance
using the \code{nsubsets} criterion of \code{\link{evimp}}.

\bold{Which predictors were added to the model first?}

You can see the forward pass adding terms with \code{trace=2} or higher.
But remember, pruning will usually remove some of the terms.
You can also use
\preformatted{
summary(my.model, decomp="none")
}
which will list the basis functions remaining after pruning,
in the order they were added by the forward pass.

\bold{Which predictors are actually in the model?}

The following function will give a list of predictors in the model:
\preformatted{
get.used.pred.names <- function(obj) # obj is an earth object
{
  any1 <- function(x) any(x != 0)    # like any but no warning if x is double
  names(which(apply(obj$dirs[obj$selected.terms,,drop=FALSE],2,any1)))
}
}
\bold{How can I train on one set of data and test on another?}

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
times to see how the model changes as \code{\link{sample}}
splits the dataset differently on each run.
Also, remember that the test set should not be used for parameter tuning
because you will be optimizing for the test set ---
instead use GCVs, separate parameter selection sets, or techniques
such as cross-validation.

\bold{Why do I get fewer terms than \code{nk}, even with \code{prune="none"}?}

There are several conditions that can terminate the forward pass,
and reaching \code{nk} is just one of them.
See the "Forward pass" section above.

Setting earth's argument \code{thresh} to zero is treated as a special case:
\code{thresh=0} disables all termination conditions except \code{nk}
and conditions involving numerical limits.
With \code{thresh=0}, the measured GRSq
(and thus the efficacy of the pruning pass)
should be treated with skepticism,
especially if you get the warning\cr
\code{effective number of GCV parameters >= number of cases}.

\bold{Why do I get fewer terms than \code{nprune}?}

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

\bold{Is it best to hold down model size with \code{nk} or \code{nprune}?}

If you want the best possible small model, build a big model
(by specifying a big \code{nk})
and prune it back (by specifying a small \code{nprune}).
This is better than directly building a small model by
specifying a small \code{nk}, because the pruning pass can look at all the
terms whereas the forward pass can only see one term ahead.
However, it is much faster building a small model by specifying a small \code{nk}.

\bold{Can you give an example of the \code{linpreds} argument?}

With \code{linpreds} you can specify which predictors should enter linearly,
instead of in hinge functions.
The \code{linpreds} argument does not stipulate that a predictor \emph{must} enter the model,
only that if it enters it should enter linearly.
Starting with
\preformatted{
a1 <- earth(Volume ~ ., data = trees)
plotmo(a1)
}
we see in the \code{\link{plotmo}} graphs or by running \code{\link{evimp}}
that \code{Height} is not as important as \code{Girth}.
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
a3 <- earth(Volume ~ ., data = trees, linpreds = c(1,2))
}
or with (the single TRUE is recycled to the length of \code{linpreds})
\preformatted{
a4 <- earth(Volume ~ ., data = trees, linpreds = TRUE)
}
But specifying that all predictors should enter linearly is not really a useful thing to do.
In our simple example, the all-linear MARS model is the
same as a standard linear model
\preformatted{
a5 <- lm(Volume ~ ., data = trees)
}
(compare the \code{\link{summary}} for each) but in general that will not be true.
Earth will not include a linear predictor if that predictor does not improve the model.

\bold{Can you give an example of the \code{allowed} argument?}

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
a2 <- earth(Volume ~ . - Height, data = trees)
}
The example below is more useful.
It prevents the specified predictor
from being used in interaction terms.
(The example is artificial because it's unlikely you
would want to single out humidity from interactions in the ozone data.)

The \code{parents} argument is the candidate parent's row in the \code{dirs} matrix
(\code{dirs} is described in the "Value" section above).
Each entry of \code{parents} is 0, 1, -1, or 2, and you index
\code{parents} on the predictor index.
Thus \code{parents[pred]} is 0 if \code{pred} is not in the parent term.
\preformatted{
example2 <- function(degree, pred, parents)
{
    # disallow humidity in terms of degree > 1
    # 3 is the "humidity" column in the input matrix
    if (degree > 1 && (pred == 3 || parents[3]))
        return(FALSE)
    TRUE
}
a3 <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example2)
print(summary(a3))
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
a4 <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example3)
print(summary(a4))
}
\bold{Further notes.}
The basic MARS model building strategy is always
applied even when there is an \code{allowed} function.
For example, earth considers a term for addition only
if all factors of that term except the new one are already in a model term.
This means that an \code{allowed} function that inhibits, say, all degree 2
terms will also effectively inhibit higher degrees too, because
there will be no degree 2 terms for earth to extend to degree 3.

You can expect model building to be about 10\% slower with an \code{allowed} function
because of the time taken to invoke the \code{allowed} function.

\bold{Using predictor names instead of indices in the "allowed" function.}
You can use predictor names instead of indices using
the optional \code{namesx} argument.
If present, \code{namesx} is the column names of \code{x} after factors have been expanded.
The first example above (the one that disallows \code{Height}) can be rewritten as
\preformatted{
example1a  <- function(degree, pred, parents, namesx)
{
    namesx[pred] != "Height"
}
}
Comparing strings is inefficient and
the above example can be rewritten a little more efficiently
using the optional \code{first} argument.
If present, this is TRUE the first time your allowed function is called for
the current model and thereafter FALSE.
\preformatted{
iheight <- 0    # column index of "Height"

example1b  <- function(degree, pred, parents, namesx, first)
{
    if (first) {
        # first time this function is invoked, so
        # stash column index of "Height" in iheight
        iheight <<- which(namesx == "Height")   # note use of <<- not <-
        if (length(iheight) != 1) # sanity check
            stop("no Height in ", paste(namesx, collapse=" "))
    }
    pred != iheight
}
}
\bold{How does \code{summary.earth} order terms?}

With \code{decomp="none"},
the terms are ordered as created by the forward pass.

With the default \code{decomp="anova"},
the terms are ordered in increasing order of interaction.
In detail:\cr
(i) terms are sorted first on degree of interaction\cr
(ii) then terms with a \code{linpreds} linear factor before standard terms\cr
(iii) then on the predictors (in the order of the columns in the input matrix)\cr
(ii) and finally on increasing knot values.

It's actually \code{earth:::reorder.earth} that does the ordering.

\bold{\code{summary.earth} lists predictors with weird names that aren't in \code{x}.  What gives?}

You probably have factors in your \code{x} matrix,
and earth is applying \code{\link{contrasts}}.
See the "Factors" section above.

\bold{Why \code{pmax} and not \code{max} in the output from \code{summary.earth} (with \code{style="pmax"})?}

With \code{\link{pmax}} the earth equation is an R expression
that can handle multiple cases.
Thus the expression is consistent with the
way \code{\link{predict.earth}} works --- you can give \code{\link{predict}}
multiple cases (i.e., multiple rows in the input matrix)
and it will return a vector of predicted values.

\bold{What about boosting MARS?}

If you want to boost, use boosted trees rather than boosted MARS --- you
will get better results.

More precisely, although gradient boosted MARS gives
better results than plain MARS,
if you would like to improve prediction performance (at the cost
of a more complicated and less interpretable model)
you will usually get better results with
boosted trees (via, say, the \code{\link[gbm]{gbm}} package) than with boosted MARS.
See Gillian Ward (2007) \emph{Statistics in Ecological Modeling:
Presence-Only Data and Boosted Mars (Doctoral Thesis)}
\url{http://www-stat.stanford.edu/~hastie/THESES/Gill\_Ward.pdf}.

This could change as the state of the art advances.

\bold{What about bagging MARS?}

The \code{caret} package provides functions for bagging earth
(and for parameter selection).
Our personal experience has been that bagging earth does not give models
with better predictive ability.
Your mileage may vary.

\bold{What is a GCV, in simple terms?}

GCVs are important for MARS because the pruning pass
uses GCVs to evaluate model subsets.

In general terms,
when testing a model
(not necessarily a MARS model)
we want to test \emph{generalization} performance and so
want to measure error on independent data, i.e., not on the training data.
Often a decent set of independent data is unavailable and so we
resort to cross validation or leave-one-out methods.
But that can be painfully slow.
As an alternative, for certain forms of model we can use a formula to
approximate the error that would be determined by leave-one-out
validation --- that approximation is the GCV.
The formula adjusts (i.e., increases) the training RSS to take into
account the flexibility of the model.
Summarizing, the GCV approximates the RSS (divided by the number of cases)
that would be measured on independent data.
Even when the approximation is not that good,
it is usually good enough for comparing models during pruning.

GCVs were introduced by Craven and Wahba, and extended by Friedman for MARS.
See Hastie et al. p216 and the Friedman MARS paper.
GCV stands for "Generalized Cross Validation", a perhaps misleading term.
% For example, the terms selected by the pruning pass for each subset size
% are the same whether the GCV or RSS is used to select terms.
% Actual cross validation during pruning would choose terms for each subset that
% are different in general from those selected by the RSS on the full training set.

The \code{GRSq} measure used in the earth package standardizes the raw GCV,
in the same way that R-Squared standardizes the RSS.

\bold{If GCVs are so important, why don't linear models use them?}

First a few words about overfitting.
An overfit model fits the training data well but will not give good predictions on new data.
The idea is that the training data captures the underlying structure
in the system being modeled, plus noise.
We want to model the underlying structure and ignore the noise.
An overfit model models the specific realization of noise in the training data
and thus is too specific to the training data.

The more flexible a model, the more its propensity to overfit the training data.
Linear models are constrained, with usually only a few parameters,
and don't have the tendency to overfit like more flexible models such as MARS.
This means that for linear models, the RSS on the data used to build the model is usually an
adequate measure of generalization ability.

This is no longer true if you do automatic variable selection on linear models,
because the process of selecting variables increases the flexibility
of the model. Hence the AIC --- as used in, say, \code{\link{drop1}}.
The GCV, AIC, and friends are means to the same end.
Depending on what information is available during model building.
we use one of these statistics to estimate model generalization performance
for the purpose of selecting a model.

\bold{What happened to \code{get.nterms.per.degree},
\code{get.nused.preds.per.subset}, and \code{reorder.earth}?}

From release 1.3.0, some earth functions are no longer public,
to help simplify the user interface.
The functions are still available (and stable) if you need them ---
use for example \code{earth:::reorder.earth()}.

\bold{What happened to the \code{ppenalty} argument?}

This was removed (release 1.3.1) because it is no longer needed.
The \code{ponly} argument of \code{\link{update.earth}}
is a more flexible way of achieving the same end.
}
\author{
  Stephen Milborrow, derived from \code{mda::\link[mda]{mars}}
  by Trevor Hastie and Robert Tibshirani.

  The approach used for GLMs was motivated by work done by
  Jane Elith and John Leathwick
  (a representative paper is listed in the references below).

  The \code{\link{evimp}} function uses ideas from Max Kuhn's \code{caret} package
  \url{http://cran.r-project.org/web/packages/caret/index.html}.

  Users are encouraged to send feedback --- use
  milbo AT sonic PERIOD net \url{http://www.milbo.users.sonic.net}.
}
\references{
  The primary references are the Friedman papers.
  Readers may find the MARS section in Hastie, Tibshirani,
  and Friedman a more accessible introduction.
  The Wikipedia article is recommended for an elementary introduction.
  Faraway takes a hands-on approach,
  using the \code{\link[=ozone1]{ozone}} data to compare \code{mda::mars} with other techniques.
  (If you use Faraway's examples with \code{earth} instead of \code{mars}, use \code{$bx}
  instead of \code{$x}.)
  Friedman and Silverman is recommended background reading for the MARS paper.
  Earth's pruning pass uses the \code{\link[leaps]{leaps}} package which is based on
  techniques in Miller.

  Faraway (2005) \emph{Extending the Linear Model with R}
  \url{http://www.maths.bath.ac.uk/~jjf23}

  Friedman (1991) \emph{Multivariate Adaptive Regression Splines (with discussion)}
  Annals of Statistics 19/1, 1--141
  \url{http://www.salfordsystems.com/doc/MARS.pdf}

  Friedman (1993) \emph{Fast MARS}
  Stanford University Department of Statistics, Technical Report 110
  \url{http://www.milbo.users.sonic.net/earth/Friedman-FastMars.pdf},
  \url{http://www-stat.stanford.edu/research/index.html}

  Friedman and Silverman (1989)
  \emph{Flexible Parsimonious Smoothing and Additive Modeling}
  Technometrics, Vol. 31, No. 1.
  \url{http://links.jstor.org/sici?sici=0040-1706\%28198902\%2931\%3A1\%3C3\%3AFPSAAM\%3E2.0.CO\%3B2-Z}

  Hastie, Tibshirani, and Friedman (2001) \emph{The Elements of Statistical Learning}
  \url{http://www-stat.stanford.edu/~hastie/pub.htm}

  Leathwick, J.R., Rowe, D., Richardson, J., Elith, J., & Hastie, T. (2005)
  \emph{Using multivariate adaptive regression splines to predict the distributions
  of New Zealand's freshwater diadromous fish} Freshwater Biology, 50, 2034-2052
  \url{http://www-stat.stanford.edu/~hastie/pub.htm},
  \url{http://www.botany.unimelb.edu.au/envisci/about/staff/elith.html}

  Miller, Alan (1990, 2nd ed. 2002) \emph{Subset Selection in Regression}
  \url{http://users.bigpond.net.au/amiller}

  Wikipedia article on MARS
  \url{http://en.wikipedia.org/wiki/Multivariate_adaptive_regression_splines}
}
\seealso{
  Start with \code{\link{summary.earth}}, \code{\link{plot.earth}},
  \code{\link{plotmo}}, and \code{\link{evimp}}.

  \code{\link{etitanic}}\cr
  \code{\link{evimp}}\cr
  \code{\link{format.earth}}\cr
  \code{\link{mars.to.earth}}\cr
  \code{\link{model.matrix.earth}}\cr
  \code{\link{ozone1}}\cr
  \code{\link{plot.earth.models}}\cr
  \code{\link{plot.earth}}\cr
  \code{\link{plotd}}\cr
  \code{\link{plotmo}}\cr
  \code{\link{predict.earth}}\cr
  \code{\link{residuals.earth}}\cr
  \code{\link{summary.earth}}\cr
  \code{\link{update.earth}}
}
\examples{
a <- earth(Volume ~ ., data = trees)
summary(a, digits = 2, style = "pmax")

# yields:
#    Call: earth(formula=Volume~., data=trees)
#
#    Volume =
#      23
#      +  5.7 * pmax(0,  Girth -     13)
#      -  2.9 * pmax(0,     13 -  Girth)
#      + 0.72 * pmax(0, Height -     76)
#
#    Selected 4 of 5 terms, and 2 of 2 predictors
#    Estimated importance: Girth Height
#    Number of terms at each degree of interaction: 1 3 (additive model)
#    GCV 11    RSS 213    GRSq 0.96    RSq 0.97
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

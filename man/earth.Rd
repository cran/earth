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

See the package vignette "Notes on the \code{earth} package"
\url{../doc/earth-notes.pdf}.
}
\usage{
\method{earth}{formula}(formula = stop("no 'formula' arg"),
   data, weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
   na.action = na.fail, glm = NULL, trace = 0,
   keepxy = FALSE, ncross=1, nfold=0, stratify=TRUE, \dots)

\method{earth}{default}(x = stop("no 'x' arg"), y = stop("no 'y' arg"),
    weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
    na.action = na.fail, glm = NULL, trace = 0,
    keepxy = FALSE, ncross=1, nfold=0, stratify=TRUE, \dots)

\method{earth}{fit}(x = stop("no 'x' arg"), y = stop("no 'y' arg"),
    weights = NULL, wp = NULL, scale.y = (NCOL(y)==1), subset = NULL,
    na.action = na.fail, glm = NULL, trace = 0,
    nk = max(21, 2 * ncol(x) + 1), degree = 1,
    penalty = if(degree > 1) 3 else 2, thresh = 0.001,
    minspan = 0, newvar.penalty = 0, fast.k = 20, fast.beta = 1,
    linpreds = FALSE, allowed = NULL,
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep"),
    nprune = NULL,
    Object = NULL, Get.crit = get.gcv, Eval.model.subsets = eval.model.subsets,
    Force.xtx.prune = FALSE, Use.beta.cache = TRUE, Exhaustive.tol = 1e-10, \dots)
}
\arguments{
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
% %      You do not need to normalize \code{weights} --- \code{earth} does that.
% %      \code{wp} is normalized so that its elements sum to
% %      \code{nrow(x)} and all residuals calculated in the forward and pruning
% %      pass are weighted by the corresponding element of \code{weights}.
% %      Zero values are converted to a small non zero value.
  }
  \item{wp}{
     Vector of response weights.
     Default is NULL, meaning no response weights.
     If specified, \code{wp} must have an element for each column of
     \code{y} (after \code{\link[=factor]{factors}}, if any, have been expanded).\cr
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
    \code{\link[=scale]{Scale}} \code{y} in the forward pass for better numeric stability.
    Scaling here means subtract the mean and divide by the standard deviation.
    Default is \code{NCOL(y)==1},
    i.e., scale \code{y} unless \code{y} has multiple columns.
  }
  \item{na.action}{
     NA action. Default is \code{na.fail}, and only \code{na.fail} is supported.
  }
  \item{glm}{
     NULL (default) or a list of arguments to pass on to \code{\link{glm}}.
     See the documentation of \code{\link{glm}} for a description of these arguments
     (but not all \code{glm}'s arguments are supported, \code{earth} will give an error message
     if you use a glm argument it does not support).\cr
     See also ``Generalized linear models'' in the package vignette.
     Example:\cr
     \code{earth(y~x, glm=list(family=binomial))}
  }
  \item{trace}{
     Trace \code{earth}'s execution.  Default is 0.  Values:\cr
        0 no tracing\cr
        0.5 cross validation\cr
        1 overview\cr
        2 forward pass\cr
        3 pruning\cr
        4 model mats, pruning details\cr
        5 internal details of operation\cr

  }
  \item{keepxy}{
     Set to TRUE to retain the following in the returned value: \code{x} and \code{y} (or \code{data}),
     \code{subset}, and \code{weights}.
     Default is FALSE.
     The function \code{\link{update.earth}} and friends will use these
     if present instead of searching for them
     in the environment at the time \code{\link{update.earth}} is invoked.
     This argument also affects the amount of data kept when
     the \code{nfold} argument is used
     (see \code{cv.list} in the ``Value'' section below.\cr

     \bold{The following arguments are for the forward pass.}
  }
  \item{nk}{
     Maximum number of model terms before pruning, i.e., the
     maximum number of terms created by the forward pass.
     Includes the intercept.
     Default is \code{max(21,2*NCOL(x)+1)}.\cr
     Use \code{trace=1} or more to see why the forward pass terminated --- if it
     is "Reached max number of terms" then you could consider increasing \code{nk}.
     The number of terms created by the forward pass will be
     less than \code{nk} if a forward stopping condition is reached before \code{nk} terms,
     or if the forward pass drops one side of a hinge pair to prevent linear dependencies.
     See also ``The forward pass'' in the package vignette.
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
    Simulation studies have suggested values in the range of about 2 to 4.
    In practice, for big datasets larger values can be useful to force a smaller model.
    The FAQ section in the package vignette has some information on GCVs.
  }
  \item{thresh}{
     Forward stepping threshold.
     Default is 0.001.
     This is one of the arguments used to decide when forward stepping
     should terminate.
     See also ``The forward pass'' in the package vignette.
  }
  \item{minspan}{
    Minimum distance between knots.
    The \code{minspan} argument is intended to increase resistance to runs of noise in the input data.
    (Note: predictor value extremes are ineligible for knots
    regardless of the \code{minspan} setting, as per the MARS paper equation 45.)\cr
    Use \code{minspan=1} to consider all x values (which is good if the
    data are not noisy).\cr
    The default is \code{minspan=0}.
    This value 0 is treated specially and
    means calculate the \code{minspan} internally as per
    Friedman's MARS paper section 3.8 with \eqn{alpha} = 0.05.
    Set \code{trace>=2} to see the calculated value.\cr
    \code{minspan=-1} is also treated specially. It is for back compatibility, and means
    calculate minspan using the (incorrect) method used in
    versions of \code{earth} prior to 2.4-0.
    Using \code{minspan=-1} instead of the default \code{minspan=0}
    will usually build a model with a very similar GCV,
    although with slightly different knots and terms.\cr
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
     would do better to delete one of them before calling \code{earth}.
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
     However it is not unusual to get a slightly better model with a lower \code{fast.k}.
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
     Examples are given in the package vignette.\cr
  }
  \item{allowed}{
     Function specifying which predictors can interact and how.
     Default is NULL, meaning all standard MARS terms are allowed.\cr
     During the forward pass, \code{earth} calls the \code{allowed} function
     before considering a term for inclusion; the term can go into the
     model only if the \code{allowed} function returns TRUE.
     See the package vignette for details.
\cr\cr
     \bold{The following arguments are for cross validation.}
  }
  \item{ncross}{
    Only applies if \code{nfold>1}.
    Number of cross-validations.  Each cross-validation has \code{nfold} folds.
    Default 1.
  }
  \item{nfold}{
     Number of cross-validation folds.
     Default is 0, no cross validation.
     If greater than 1, \code{earth} first builds a standard model as usual with all the data.
     It then builds \code{nfold} cross-validated models,
     measuring R-Squared on the out-of-fold (left out) data each time.
     The final cross validation R-Squared (\code{cv.rsq}) is the mean of these fold R-Squareds.\cr
     The above process of building \code{nfold} models is repeated
     \code{ncross} times (by default, once).\cr
     Further statistics are calculated if \code{keepxy=TRUE} or
     if a binomial or poisson model (the \code{glm} argument).
     For details, see ``Cross validation'' in the package vignette.
  }
  \item{stratify}{
     Only applies if \code{nfold>1}.
     Default is TRUE.
     Stratify the cross-validation samples so that
     an approximately equal number of cases with a non-zero response
     occur in each cross validation subset.
     So if the response \code{y} is logical, the \code{TRUE}s will be spread
     evenly across folds.
     If the response is a factor, there will be an
     approximately equal number of each factor level in each fold
     (because a multilevel factor response gets expanded to columns of zeros and ones,
     see ``Factors'' in  the package vignette).
     We say ``approximately equal'' because the number of occurrences of a factor
     level may not be exactly divisible by the number of folds.
\cr\cr
     \bold{The following arguments are for the pruning pass.}
  }
  \item{pmethod}{
     Pruning method.
     One of: \code{backward none exhaustive forward seqrep}.\cr
     Default is \code{"backward"}.\cr
     Use \code{none} to retain all the terms created by the forward pass.\cr
     If \code{y} has multiple columns, then only \code{backward} or \code{none}
     is allowed.\cr
     Pruning can take a while if \code{exhaustive} is chosen and
     the model is big (more than about 30 terms).
     The current version of the \code{\link[leaps]{leaps}} package
     used during pruning does not allow user interrupts
     (i.e., you have to kill your R session to interrupt;
     in Windows use the Task Manager or from the command line use \code{taskkill}).
  }
  \item{nprune}{
     Maximum number of terms (including intercept) in the pruned model.
     Default is NULL, meaning all terms created by the forward pass
     (but typically not all terms will remain after pruning).
     Use this to enforce an upper bound on the model size (that is less than \code{nk}),
     or to reduce exhaustive search time with \code{pmethod="exhaustive"}.
\cr\cr
    \bold{The following arguments are for internal or advanced use.}
  }
  \item{Object}{
     Earth object to be updated, for use by \code{\link{update.earth}}.
  }
  \item{Get.crit}{
     Criterion function for model selection during pruning.
     By default a function that returns the GCV.
     See ``The pruning pass'' in the package vignette.
  }
  \item{Eval.model.subsets}{
     Function to evaluate model subsets
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
    Using the ``beta cache'' takes more memory but is faster
    (by 20\% and often much more for large models).
    The beta cache uses \code{nk * nk * ncol(x) * sizeof(double)} bytes.
    Set \code{Use.beta.cache=FALSE} to save a little memory.
    (The beta cache is an innovation in this implementation of MARS
    and does not appear in Friedman's papers. It is not related to
    the \code{fast.beta} argument.  Certain regression coefficients
    in the forward pass can be re-used, and are saved, thus
    saving recalculation time.)
  }
  \item{Exhaustive.tol}{
     Default \code{1e-10}.
     Applies only when \code{pmethod="exhaustive"}.
     If the reciprocal of the condition number of \code{bx}
     is less than \code{Exhaustive.tol}, \code{earth} forces \code{pmethod="backward"}.
     See ``XHAUST returned error code -999'' in the package vignette.
  }
  \item{\dots}{
     Dots are passed on to \code{earth.fit}.
 }
}
\value{
  An object of class \code{"earth"} which is a list with the components
  listed below.
  \emph{Term} refers to a term created during the forward pass
  (each line of the output from \code{\link{format.earth}} is a term).
  Term number 1 is always the intercept.
  \item{rss}{
     Residual sum-of-squares (RSS) of the model (summed over all responses
     if \code{y} has multiple columns).
  }
  \item{rsq}{
    \code{1-rss/tss}.
     R-Squared of the model (calculated over all responses).
     A measure of how well the model fits the training data.
     Note that \code{tss} is the total sum-of-squares, \code{sum((y - mean(y))^2)}.
  }
  \item{gcv}{
     Generalized Cross Validation (GCV) of the model (summed over all responses)
     The GCV is calculated using the \code{penalty} argument.
     For details of the GCV calculation, see
     equation 30 in Friedman's MARS paper and \code{earth:::get.gcv}.
  }
  \item{grsq}{
     \code{1-gcv/gcv.null}.
     An estimate of the predictive power of the model (calculated over all responses,
     \code{gcv.null} is the GCV of an intercept-only model).
     See ``Can \code{GRSq} be negative?'' in the
     package vignette for some discussion.
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
\code{-1} if an expression of the form \code{h(const - xj)} is in term i\cr
\code{1}  if an expression of the form \code{h(xj - const)} is in term i\cr
\code{2}  if predictor j enters term i linearly (the \code{linpreds} argument).\cr

     This matrix includes all terms generated by the forward.pass,
     including those not in \code{selected.terms}.
     Note that here the terms may not all be in pairs, because
     although the forward pass add terms as hinged pairs (so both sides of
     the hinge are available as building blocks for further terms), it also
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
     Vector of term numbers in the selected model.
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
     A vector of the RSS (measured on the training data)
     for each model subset generated by the pruning pass.
     Length is \code{nprune}.
     For multiple responses, the RSS is summed over all responses for each subset.
     The \code{rss} above is\cr
     \code{rss.per.subset[length(selected.terms)]}.
     The RSS of an intercept only-model is \code{rss.per.subset[1]}.
    }
  \item{gcv.per.subset}{
     A vector of the GCV for each model in \code{prune.terms}.
     Length is \code{nprune}.
     For multiple responses, the GCV is summed over all responses for each subset.
     The \code{gcv} above is \code{gcv.per.subset[length(selected.terms)]}.
     The GCV of an intercept-only model is \code{gcv.per.subset[1]}.
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
     Column names of \code{x}, generated internally by \code{earth} when necessary
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
     Copy of the \code{wp} argument to \code{earth}.\cr\cr
  \bold{The following fields appear only if \code{earth}'s argument \code{keepxy} is TRUE.}
  }
  \item{x}{}
  \item{y}{}
  \item{data}{}
  \item{subset}{}
  \item{weights}{
     Copies of the corresponding arguments to \code{earth}.
     Only exist if \code{keepxy=TRUE}.
  \cr\cr
  \bold{The following fields appear only if \code{earth}'s \code{glm} argument is used.}
  }
  \item{glm.list}{
     List of GLM models. Each element is the value returned by \code{earth}'s
     internal call to \code{\link{glm}} for each response.\cr
     Thus if there is a single response (or a single binomial pair, see
     ``Binomial pairs'' in the package vignette)
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
       A logical vector, derived internally by \code{earth}, or a copy
       the \code{bpairs} specified by the user in the \code{glm} list.
       See ``Binomial pairs'' in the package vignette.
   \cr\cr
   \bold{The following fields appear only if the \code{nfold} argument is greater than 1.}
   }
   \item{cv.list}{
      List of \code{earth} models, one model for each fold (\code{ncross * nfold} models).
      To save memory, lengthy fields
      in the fold models are removed unless you use \code{keepxy=TRUE}.
      The ``lengthy fields'' are \code{$bx}, \code{$fitted.values}, and \code{$residuals}.
      The fold models have two extra fields,
      \code{icross} (the cross-validation index, \code{1:ncross})
      and \code{ifold} (the fold index, \code{1:nfold}).
   }
   \item{cv.nterms}{
      Vector of length \code{ncross * nfold + 1}.
      Number of MARS terms in the model generated at each cross-validation fold,
      with the final element being the mean of these.
   }
   \item{cv.nvars}{
      Vector of length \code{ncross * nfold + 1}.
      Number of predictors in the model generated at each cross-validation fold,
      with the final element being the mean of these.
   }
   \item{cv.groups}{
     Specifies which cases went into which folds.
     Matrix with two columns and number of rows equal to the the number of cases \code{nrow(x)}
     Elements of the first column specify the cross-validation number, \code{1:ncross}.
     Elements of the second column specify the fold number, \code{1:nfold}.
   }
   \item{cv.rsq.tab}{
      Matrix with \code{ncross * nfold + 1} rows and \code{nresponse+1} columns,
      where \code{nresponse} is the number of responses,
      i.e., \code{ncol(y)} after factors in \code{y} have been expanded.
      The first \code{nresponse} elements of a row are the \code{cv.rsq}'s on
      the out-of-fold data for each response of the model generated at that row's fold.
      (A \code{cv.rsq} is calculated from predictions on the out-of-fold data
      using the best model built from the in-fold data;
      where ``best'' means the model was selected using the in-fold GCV.)
      The final column holds the row mean (a weighted mean if \code{wp} if specified).
      The final row holds the column means.
      The values in this final row is the mean \code{cv.rsq}
      printed by \code{\link{summary.earth}}.
\cr\cr
      Example for a single response model (where the \code{mean} column
is redundant but included for uniformity with multiple response models):
\preformatted{            y  mean
fold1  0.909 0.909
fold2  0.869 0.869
fold3  0.952 0.952
fold4  0.157 0.157
fold5  0.961 0.961
mean    0.769 0.769
}
      Example for a multiple response model:
\preformatted{          y1   y2    y3  mean
fold1 0.915 0.951 0.944 0.937
fold2 0.962 0.970 0.970 0.968
fold3 0.914 0.940 0.942 0.932
fold4 0.907 0.929 0.925 0.920
fold5 0.947 0.987 0.979 0.971
mean   0.929 0.955 0.952 0.946

}
   }
   \item{cv.oof.rsq.tab}{
     Matrix with \code{ncross * nfold + 1} rows and \code{max.nterms} columns,
     Only calculated and kept if \code{keepxy=TRUE}.
     Each element holds an out-of-fold RSq (\code{oof.rsq}),
     calculated from predictions from the out-of-fold observations using
     the model built with the in-fold data.  The final row is the mean over
     all folds.
   }
   \item{cv.infold.rsq.tab}{
     Like \code{cv.oof.rsq.tab} but from predictions made on the in-fold observations.
   }
   \item{cv.class.rate.tab}{
      Like \code{cv.rsq.tab} but is the classification rate at each fold
      i.e. the fraction of classes correctly predicted.
      Models with discrete response only.
      Calculated with \code{thresh=.5} for binary responses.
      For responses with more than two
      levels, the final row is the overall classification rate.  The other
      rows are the classification rates for each level (the level versus
      versus not-the-level), which are usually higher than the overall
      classification rate (predicting the level versus not-the-level is
      easier than correctly predicting one of many levels).
   }
   \item{cv.maxerr.tab}{
      Like \code{cv.rsq.tab} but is the \code{MaxErr} at each fold.
      This is the signed max absolute value at each fold.
      Results are aggregated for the final column and final row
      using the signed max absolute value.
      The \emph{signed max absolute value} is defined
      as the maximum of the absolute difference
      between the predicted and observed response values, multiplied
      by \code{-1} if the sign of that difference is negative.
   }
   \item{cv.auc.tab}{
      Like \code{cv.rsq.tab} but is the \code{AUC} at each fold.
      Binomial models only.
   }
   \item{cv.cor.tab}{
      Like \code{cv.rsq.tab} but is the \code{cor} at each fold.
      Poisson models only.
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
}
\note{
Please see the package vignette "Notes on the \code{earth} package"
\url{../doc/earth-notes.pdf}.\cr
The vignette can also be downloaded from
\url{www.milbo.org/doc/earth-notes.pdf}.
}
\author{
  Stephen Milborrow, derived from \code{mda::\link[mda]{mars}}
  by Trevor Hastie and Robert Tibshirani.

  The approach used for GLMs was motivated by work done by
  Jane Elith and John Leathwick
  (a representative paper is given below).

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
  instead of \code{$x}, and check out the book's errata.)
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

  Hastie, Tibshirani, and Friedman (2009) \emph{The Elements of Statistical Learning (2nd ed.)}
  \url{http://www-stat.stanford.edu/~hastie/pub.htm}

  Leathwick, J.R., Rowe, D., Richardson, J., Elith, J., & Hastie, T. (2005)
  \emph{Using multivariate adaptive regression splines to predict the distributions
  of New Zealand's freshwater diadromous fish} Freshwater Biology, 50, 2034-2052
  \url{http://www-stat.stanford.edu/~hastie/pub.htm},
  \url{http://www.botany.unimelb.edu.au/envisci/about/staff/elith.html}

  Miller, Alan (1990, 2nd ed. 2002) \emph{Subset Selection in Regression}
  \url{http://www.cmis.csiro.au/Alan_Miller/index.html}

  Wikipedia article on MARS
  \url{http://en.wikipedia.org/wiki/Multivariate_adaptive_regression_splines}
}
\seealso{
  Start with \code{\link{summary.earth}}, \code{\link{plot.earth}},
  \code{\link{evimp}}, and \code{\link[plotmo]{plotmo}}.

  See the package vignette \url{../doc/earth-notes.pdf}.
  The vignette can also be downloaded from
  \url{www.milbo.org/doc/earth-notes.pdf}.
}
\examples{
a <- earth(Volume ~ ., data = trees)
plotmo(a)
summary(a, digits = 2, style = "pmax")
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

\name{expand.bpairs}
\alias{expand.bpairs}
\alias{expand.bpairs.default}
\alias{expand.bpairs.formula}
\title{Expand binomial-pair data from short to long form}
\description{
Expand binomial-pair data from ``short'' to ``long'' form.

The short form specifies the response with two columns
giving the numbers of successes and failures.
Example short form:\preformatted{
    survived died dose    sex
           3    0   10   male
           2    1   10 female
           1    2   20   male
           1    2   20 female
}
The long form specifies the response as single column of
\code{TRUE}s and \code{FALSE}s.
Long form of the above data (spaces and comments added):\preformatted{
    survived dose    sex
        TRUE   10   male     # row 1 of short data: 0 died, 3 survived
        TRUE   10   male
        TRUE   10   male

       FALSE   10 female     # row 2 of short data: 1 died, 2 survived
        TRUE   10 female
        TRUE   10 female

       FALSE   20   male     # row 3 of short data: 2 died, 1 survived
       FALSE   20   male
        TRUE   20   male

       FALSE   20 female     # row 4 of short data: 2 died, 1 survived
       FALSE   20 female
        TRUE   20 female
}
In this example the total number of survived and died for each row in
the short data is the same, but in general that need not be true.
}
\usage{
\method{expand.bpairs}{formula}(formula = stop("no 'formula' argument"), data = NULL, sort = FALSE, \dots)

\method{expand.bpairs}{default}(data = stop("no 'data' argument"), y = NULL, sort = FALSE, \dots)
}
\arguments{
\item{formula}{
   Model formula such as \code{survived + died ~ dose + temp}.
}
\item{data}{
     Matrix or dataframe containing the data.
}
\item{y}{
     Model response.  One of:
\cr
     o Two column matrix or dataframe of binomial pairs e.g. \code{cbind(survived, died=20-survived)}
\cr
     o Two-element numeric vector specifying the response columns in \code{data} e.g. \code{c(1,2)}
\cr
     o Two-element character vector specifying the response column names in \code{data}
       e.g. \code{c("survived", "died")}.
       The full names must be used (partial matching isn't supported).
}
\item{sort}{
     Default \code{FALSE}.
     Use \code{TRUE} to sort the rows of the long data so it is returned
     in canonical form, independent of the row order of the short data.
     The long data is sorted on predictor values;
     predictors on the left take precedence in the sort order.
}
\item{\dots}{
     Unused, but provided for generic/method consistency.
}
}
\value{
A dataframe of the data in the long form, with expanded binomial pairs.
The first column of the data will be the response column
(a column of \code{TRUE}s and \code{FALSE}s).

Additionally, the returned value has two attached attributes:

\code{bpairs.index} A vector of row indices into the returned data.
Can be used to reconstruct the short data from the long data
(although this package does not yet provide a function to do so).

\code{ynames} Column names of the original response (a two-element character vector).
}
\examples{
survived <- c(3,2,1,1) # short data for demo (too short to build a real model)
died     <- c(0,1,2,2)
dose <- c(10,10,20,20)
sex  <- factor(c("male", "female", "male", "female"))

short.data <- data.frame(survived, died, dose, sex)

expand.bpairs(survived + died ~ ., short.data) # returns long form of the data

# expand.bpairs(data=short.data, y=cbind(survived, died)) # equivalent
# expand.bpairs(short.data, c(1,2))                       # equivalent
# expand.bpairs(short.data, c("survived", "died"))        # equivalent

# For example models, see the earth vignette
# section "Short versus long binomial data".
}

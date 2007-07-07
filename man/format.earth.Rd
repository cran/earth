\name{format.earth}
\alias{format.earth}
\title{Format 'earth' objects}
\description{
Return a string representing the \sQuote{earth} expression.
}
\usage{
\method{format}{earth}(x = stop("no 'x' arg"),
       digits = getOption("digits"), use.names = TRUE,
       add.labels = FALSE,  decomp = "anova",  \dots)
}
\arguments{
  \item{x}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{digits}{
    Number of significant digits.
    The default is \code{getOption(digits)}.
  }
  \item{use.names}{
    If \code{TRUE} (default), use variable names.  Else use names of the form \code{x[,1]}.
  }
  \item{add.labels}{
    If \code{TRUE} add comments numbering each term. The default is \code{FALSE}.
  }
  \item{decomp}{One of\cr
\code{"anova"}  (the default) order the terms using the \sQuote{anova decomposition}
          i.e. in increasing order of interaction\cr
\code{"none"}   order the terms as created during the earth forward pass.\cr
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
A character representation of the earth object.
}
\seealso{
  \code{\link{earth}},
  \code{\link{pmax}}
}
\note{
Using \code{format.earth}, perhaps after hand editing the returned string,
you can create an alternative to \code{predict.earth}.
For example:
\preformatted{
as.func <- function(
               object, digits = 8, use.names = TRUE, ...)
  eval(parse(text=paste(
    "function(x)\n",
    "{\n",
    "if(is.vector(x))\n",
    "  x <- matrix(x, nrow = 1, ncol = length(x))\n",
    "with(as.data.frame(x),\n",
    format(object, digits = digits, use.names = use.names, ...),
    ")\n",
    "}\n", sep = "")))
a <- earth(Volume ~ ., data = trees)
my.func <- as.func(a, use.names = FALSE)
my.func(c(10,80))     # yields 17.82973
predict(a, c(10,80))  # yields 17.82973, but is slower
}
}
\examples{
a <- earth(Volume ~ ., data = trees)
cat(format(a))

# yields:
#    26.6 
#    +  6.05 * pmax(0,  Girth -   13.7) 
#    -   3.2 * pmax(0,   13.7 -  Girth) 
#    + 0.608 * pmax(0, Height -     75) 

cat(format(a, use.names = FALSE, add.labels = TRUE))

# yields:
#    26.6  # 1
#    +  6.05 * pmax(0, x[,1] -  13.7)  # 2
#    -   3.2 * pmax(0,  13.7 - x[,1])  # 3
#    + 0.608 * pmax(0, x[,2] -    75)  # 4
}
\keyword{models}

\name{format.earth}
\alias{format.earth}
\title{Format "earth" objects}
\description{
Return a string representing an \code{\link{earth}} expression.
}
\usage{
\method{format}{earth}(x = stop("no 'x' arg"),
       digits = getOption("digits"), use.names = TRUE,
       decomp = "anova",  style = "h", colon.char = ":", \dots)
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
  \item{decomp}{One of\cr
\code{"anova"}  (default) order the terms using the "anova decomposition"
          i.e. in increasing order of interaction\cr
\code{"none"}   order the terms as created during the earth forward pass.\cr
  }
  \item{style}{
    Formatting style. One of\cr
\code{"h"} (default) more compact\cr
\code{"pmax"} for those who prefer it and for compatibility with old versions of earth\cr
\code{"bf"} basis function format.
  }
  \item{colon.char}{
    Change colons in the returned string to colon.char.
    Default is ":" i.e. no change.
    Specifying \code{colon.char="*"} can be useful in some contexts to change
    names of the form \code{x1:x2} to \code{x1*x2}.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
A character representation of the earth object.

If there are multiple responses, \code{format.earth} will return multiple strings.

If there are embedded GLM model(s), the string(s) for the GLM model(s)
come after the strings for the standard earth model(s).
}
\seealso{
  \code{\link{earth}},
  \code{\link{pmax}},
}
\note{
The FAQ section for \code{\link{earth}} has some comments on the \code{"anova"} ordering.

Using \code{format.earth}, perhaps after hand editing the returned string,
you can create an alternative to \code{predict.earth}.
For example:
\preformatted{
as.func <- function(
               object, digits = 8, use.names = FALSE, ...)
  eval(parse(text=paste(
    "function(x)\n",
    "{\n",
    "if(is.vector(x))\n",
    "  x <- matrix(x, nrow = 1, ncol = length(x))\n",
    "with(as.data.frame(x),\n",
    format(object, digits = digits, use.names = use.names, style = "pmax", ...),
    ")\n",
    "}\n", sep = "")))
a <- earth(Volume ~ ., data = trees)
my.func <- as.func(a, use.names = FALSE)
my.func(c(10,80))     # yields 17.82973
predict(a, c(10,80))  # yields 17.82973, but is slower
}
The  earth package also provides a function \code{format.lm}.
It has arguments as follows\cr
\code{format.lm(x, digits=getOption("digits"), use.names=TRUE, colon.char=":")}\cr
(Strictly speaking, \code{format.lm} doesn't belong in the
earth package.) Example:
\preformatted{
a <- lm(Volume ~ Height*Girth, data = trees)
cat(format(a, colon.char="*"))

# yields:
#    69.4
#    -  1.30 * Height
#    -  5.86 * Girth
#    + 0.135 * Height*Girth
}
}
\examples{
a <- earth(Volume ~ ., data = trees)
cat(format(a))

# yields:
#    23.2
#    +  5.75 * h(Girth-12.9)
#    -  2.87 * h(12.9-Girth)
#    + 0.718 * h(Height-76)

cat(format(a, style="pmax")) # default formatting style prior to earth version 1.4

# yields:
#    23.2
#    +  5.75 * pmax(0,  Girth -   12.9)
#    -  2.87 * pmax(0,   12.9 -  Girth)
#    + 0.718 * pmax(0, Height -     76)

cat(format(a, style="bf"))

# yields:
#    23.2
#    +  5.75 * bf1
#    -  2.87 * bf2
#    + 0.718 * bf3
#  
#     bf1  h(Girth-12.9)
#     bf2  h(12.9-Girth)
#     bf3  h(Height-76)
}
\keyword{models}

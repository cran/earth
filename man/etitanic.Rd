\name{etitanic}
\alias{etitanic}
\title{Titanic data with incomplete cases removed}
\description{
Titanic data with incomplete cases,
passenger names, and other details removed.
}
\format{
    \item{pclass}{passenger class, unordered factor: 1st 2nd 3rd
    }
    \item{survived}{integer: 0 or 1
    }
    \item{sex}{unordered factor: male female
    }
    \item{age}{age in years, double: min 0.167 max 80.0
    }
    \item{sibsp}{number of siblings or spouses aboard, integer: 0...8
    }
    \item{parch}{number of parents or children aboard, integer: 0...6
    }
}
\source{
This dataset is included in the earth package because it is a convenient
vehicle for illustrating earth's GLM and factor handling.

The dataset was compiled by Frank Harrell and Robert Dawson:
\url{http://biostat.mc.vanderbilt.edu/twiki/pub/Main/DataSets/titanic.html}\cr
See also:\cr
\url{http://biostat.mc.vanderbilt.edu/twiki/pub/Main/DataSets/titanic3info.txt}.

For this version of the Titanic data,
passenger details and incomplete cases were deleted
and the name changed to \code{etitanic} to minimise confusion with other
versions ("e" because it is part of the earth package).

Contents of etitanic:
    \preformatted{
         pclass survived    sex    age sibsp parch
    1       1st        1 female 29.000     0     0
    2       1st        1   male  0.917     1     2
    3       1st        0 female  2.000     1     2
    4       1st        0   male 30.000     1     2
    5       1st        0 female 25.000     1     2
    ...
    1309    3rd        0   male 29.000     0     0
    }
How \code{etitanic} was built:
    \preformatted{
    load("titanic3") # from Harrell's web site
    # discard name, ticket, fare, cabin, embarked, body, home.dest
    etitanic = titanic3[,c(1,2,4,5,6,7)]
    etitanic = etitanic[!is.na(etitanic$age),]
    save(etitanic, file="etitanic.rda")
    }
}
% \references{
% Harrell (2001) \emph{
% Regression Modeling Strategies with Applications to Linear Models, Logistic Regression, and Survival Analysis }
% \url{http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/RmS}
% }
\seealso{
  \code{\link{earth}}
}
\keyword{datasets}
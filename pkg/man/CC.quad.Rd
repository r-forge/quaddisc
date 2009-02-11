\name{CC.quad}
\alias{CC.quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CC.quad(iN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iN}{ ~~Describe \code{iN} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (iN) 
{
    vI <- matrix(0:(iN - 1), iN, 1)
    vZ <- as.vector(0.5 - 0.5 * cos(((2 * vI + 1) * pi)/(2 * 
        iN)))
    vD <- c(1, rep(2, iN - 1))/iN
    vK <- matrix(2 * (1:iN) - 1, 1, iN)
    mCt <- cos(pi * (vI \%*\% vK)/(2 * iN))
    vS <- as.matrix(1:iN)
    vS <- -2/(vS * (vS - 2))
    vS[seq(2, length(vS), by = 2)] <- 0
    vW <- as.vector(t(0.5 * crossprod(vS * vD, mCt)))
    return(list(nodes = vZ, weights = vW))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

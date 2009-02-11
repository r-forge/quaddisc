\name{hDiaphony}
\alias{hDiaphony}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hDiaphony(vX, vY, lExtraArgsKernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vX}{ ~~Describe \code{vX} here~~ }
  \item{vY}{ ~~Describe \code{vY} here~~ }
  \item{lExtraArgsKernel}{ ~~Describe \code{lExtraArgsKernel} here~~ }
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
function (vX, vY, lExtraArgsKernel = NULL) 
{
    f <- function(vX) -1 + prod(1 + 2 * pi^2 * (vX^2 - vX + 1/6))
    iD <- length(vX)
    vZ <- vX - vY
    for (i in 1:iD) vZ[i] <- mod(vZ[i], 1)
    dResult <- (1/((1 + pi^2/3)^iD - 1)) * f(vZ)
    return(dResult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

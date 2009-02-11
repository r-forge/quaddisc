\name{qgausschaos}
\alias{qgausschaos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qgausschaos(dP, vLambda, lArgsUniroot = list(interval = c(0, 100)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dP}{ ~~Describe \code{dP} here~~ }
  \item{vLambda}{ ~~Describe \code{vLambda} here~~ }
  \item{lArgsUniroot}{ ~~Describe \code{lArgsUniroot} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function (dP, vLambda, lArgsUniroot = list(interval = c(0, 100)), 
    ...) 
{
    if ((dP <= 0) | (dP >= 1)) 
        return("wrong argument (should be in (0, 1))")
    fnObj <- function(dC) {
        lP <- pgausschaos(dC, vLambda)
        if (lP$iFault != 0) 
            print("warning: error in computation of pgausschaos")
        return(lP$dResult - dP)
    }
    lArgsUniroot <- c(f = fnObj, lArgsUniroot)
    return(do.call("uniroot", lArgsUniroot))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{intrep}
\alias{intrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intrep(iK, iM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iK}{ ~~Describe \code{iK} here~~ }
  \item{iM}{ ~~Describe \code{iM} here~~ }
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
function (iK, iM) 
{
    vB <- 0
    i <- 1
    if (iK < iM) 
        vB[i] <- iK
    else {
        while (iK >= iM) {
            vB[i] <- iK\%\%iM
            iK <- iK\%/\%iM
            vB[i + 1] <- 1
            i <- i + 1
        }
    }
    return(vB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

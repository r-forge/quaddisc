\name{gGD}
\alias{gGD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gGD(vXk, sType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vXk}{ ~~Describe \code{vXk} here~~ }
  \item{sType}{ ~~Describe \code{sType} here~~ }
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
function (vXk, sType) 
{
    dResult <- 1
    if (sType == "star") 
        dResult <- gGDStar(vXk)
    else if (sType == "symm") 
        dResult <- gGDSymm(vXk)
    else if (sType == "cent") 
        dResult <- gGDCent(vXk)
    else return("Error: wrong type")
    return(dResult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

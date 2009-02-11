\name{haltset}
\alias{haltset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
haltset(iNmin, iNmax, vP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iNmin}{ ~~Describe \code{iNmin} here~~ }
  \item{iNmax}{ ~~Describe \code{iNmax} here~~ }
  \item{vP}{ ~~Describe \code{vP} here~~ }
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
function (iNmin, iNmax, vP) 
{
    iD <- length(vP)
    iN <- iNmax - iNmin + 1
    mX <- matrix(nrow = iN, ncol = iD)
    for (k in 1:iN) {
        for (s in 1:iD) {
            mX[k, s] <- radinv(k + iNmin, vP[s])
        }
    }
    return(mX)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

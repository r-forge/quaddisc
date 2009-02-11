\name{hammset}
\alias{hammset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hammset(iN, vP = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iN}{ ~~Describe \code{iN} here~~ }
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
function (iN, vP = 0) 
{
    if (vP == 0) {
        mX <- matrix(nrow = iN, ncol = 1)
        for (k in 1:iN) mX[k, 1] <- (2 * k - 1)/(2 * iN)
        return(mX)
    }
    else {
        iD <- length(vP) + 1
        mX <- matrix(nrow = iN, ncol = iD)
        for (k in 1:iN) {
            mX[k, 1] <- (2 * k - 1)/(2 * iN)
            for (s in 2:iD) {
                dY <- radinv(k, vP[s - 1])
                mX[k, s] <- dY
            }
        }
        return(mX)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

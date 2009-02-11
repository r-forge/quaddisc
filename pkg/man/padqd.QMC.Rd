\name{padqd.QMC}
\alias{padqd.QMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
padqd.QMC(dC, fnH, mZ, sForm = "v", lExtraArgsKernel, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dC}{ ~~Describe \code{dC} here~~ }
  \item{fnH}{ ~~Describe \code{fnH} here~~ }
  \item{mZ}{ ~~Describe \code{mZ} here~~ }
  \item{sForm}{ ~~Describe \code{sForm} here~~ }
  \item{lExtraArgsKernel}{ ~~Describe \code{lExtraArgsKernel} here~~ }
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
function (dC, fnH, mZ, sForm = "v", lExtraArgsKernel, ...) 
{
    if (!is.null(lExtraArgsKernel)) 
        lArgs <- c(list(mZ = mZ, fnH = fnH), lExtraArgsKernel)
    else lArgs <- list(mZ = mZ, fnH = fnH)
    if (sForm == "u") 
        mH <- do.call("matHu.QMC", lArgs)
    if (sForm == "v") 
        mH <- do.call("matH.QMC", lArgs)
    vLambda <- eigen(mH, symmetric = TRUE)$values
    return(pgausschaos(dC = dC, vLambda = vLambda, ...))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

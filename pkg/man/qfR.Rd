\name{qfR}
\alias{qfR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qfR(vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt = 1e+05, dAcc = 1e-05, vTrace = rep(0, 7), iFault = 0, dValue = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vLambda}{ ~~Describe \code{vLambda} here~~ }
  \item{vNoncent}{ ~~Describe \code{vNoncent} here~~ }
  \item{vDf}{ ~~Describe \code{vDf} here~~ }
  \item{iTermsSum}{ ~~Describe \code{iTermsSum} here~~ }
  \item{dSigma}{ ~~Describe \code{dSigma} here~~ }
  \item{dC}{ ~~Describe \code{dC} here~~ }
  \item{iMaxTermsInt}{ ~~Describe \code{iMaxTermsInt} here~~ }
  \item{dAcc}{ ~~Describe \code{dAcc} here~~ }
  \item{vTrace}{ ~~Describe \code{vTrace} here~~ }
  \item{iFault}{ ~~Describe \code{iFault} here~~ }
  \item{dValue}{ ~~Describe \code{dValue} here~~ }
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
function (vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt = 1e+05, 
    dAcc = 1e-05, vTrace = rep(0, 7), iFault = 0, dValue = 0) 
{
    .C("qfR", vLambda = as.double(vLambda), vNoncent = as.double(vNoncent), 
        vDf = as.integer(vDf), iTermsSum = as.integer(iTermsSum), 
        dSigma = as.double(dSigma), dC = as.double(dC), iMaxTermsInt = as.integer(iMaxTermsInt), 
        dAcc = as.double(dAcc), vTrace = as.double(vTrace), iFault = as.integer(iFault), 
        dResult = as.double(iFault))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

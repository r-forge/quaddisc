\name{eigfun}
\alias{eigfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eigfun(vX, fnH, iN, sMethod, bPoints = FALSE, vZ = NULL, vW = NULL, lExtraArgsKernel = NULL, vJ = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vX}{ ~~Describe \code{vX} here~~ }
  \item{fnH}{ ~~Describe \code{fnH} here~~ }
  \item{iN}{ ~~Describe \code{iN} here~~ }
  \item{sMethod}{ ~~Describe \code{sMethod} here~~ }
  \item{bPoints}{ ~~Describe \code{bPoints} here~~ }
  \item{vZ}{ ~~Describe \code{vZ} here~~ }
  \item{vW}{ ~~Describe \code{vW} here~~ }
  \item{lExtraArgsKernel}{ ~~Describe \code{lExtraArgsKernel} here~~ }
  \item{vJ}{ ~~Describe \code{vJ} here~~ }
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
function (vX, fnH, iN, sMethod, bPoints = FALSE, vZ = NULL, vW = NULL, 
    lExtraArgsKernel = NULL, vJ = 1) 
{
    lH <- matH(fnH, iN, sMethod, bPoints = FALSE, vZ = NULL, 
        vW = NULL)
    vZ <- as.matrix(lH$nodes)
    mResult <- matrix(nrow = length(vX), ncol = length(vJ))
    vW <- as.matrix(lH$weights)
    mEigVec <- as.matrix(lH$eigvecnorm[, vJ])
    if (sMethod \%in\% c("MC", "QMC", "TR", "GL")) {
        mH <- matrix(nrow = length(vX), ncol = length(vZ))
        for (j in 1:length(vZ)) {
            for (i in 1:length(vX)) {
                mH[i, j] <- fnH(vX[i], vZ[j], lExtraArgsKernel)
            }
        }
        vEigVal <- as.matrix(lH$eigval[vJ])
        mResult <- mH \%*\% ((vW \%*\% t(1/vEigVal)) * mEigVec)
    }
    else if (sMethod == "CC") {
        vI1 <- matrix(1:iN, iN, 1)
        vI0 <- vI1 - 1
        vD <- c(1, rep(2, iN - 1))/iN
        vK <- matrix(2 * vI0 + 1, 1, iN)
        mCt <- cos(pi * (vI0 \%*\% vK)/(2 * iN))
        mTt <- matrix(nrow = length(vX), ncol = iN)
        mTt <- cos(acos(2 * vX - 1) \%*\% matrix(vI0, 1, iN))
        mResult <- mTt \%*\% diag(as.vector(vD)) \%*\% mCt \%*\% mEigVec
    }
    return(mResult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

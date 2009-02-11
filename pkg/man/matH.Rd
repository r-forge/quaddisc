\name{matH}
\alias{matH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matH(fnH, iN, sMethod, bPoints = FALSE, vZ = NULL, vW = NULL, lExtraArgsKernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnH}{ ~~Describe \code{fnH} here~~ }
  \item{iN}{ ~~Describe \code{iN} here~~ }
  \item{sMethod}{ ~~Describe \code{sMethod} here~~ }
  \item{bPoints}{ ~~Describe \code{bPoints} here~~ }
  \item{vZ}{ ~~Describe \code{vZ} here~~ }
  \item{vW}{ ~~Describe \code{vW} here~~ }
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
function (fnH, iN, sMethod, bPoints = FALSE, vZ = NULL, vW = NULL, 
    lExtraArgsKernel = NULL) 
{
    if (!bPoints) {
        if (sMethod == "MC") {
            vZ <- runif(iN)
            vW <- rep(1/iN, iN)
        }
        else if (sMethod == "QMC") {
            vZ <- hammset(iN)
            vW <- rep(1/iN, iN)
        }
        else if (sMethod == "TR") {
            vZ <- as.matrix(0:(iN - 1))/(iN - 1)
            vW <- c(0.5, matrix(1, iN - 2, 1), 0.5)/(iN - 1)
        }
        else if (sMethod == "GL") {
            lGL <- gauss.quad(iN)
            vZ <- (lGL$nodes + 1)/2
            vW <- (lGL$weights)/2
        }
        else if (sMethod == "CC") {
            lCC <- CC.quad(iN)
            vZ <- lCC$nodes
            vW <- lCC$weights
        }
        mW <- diag(as.vector(vW))
        mHtemp <- matrix(0, iN, iN)
        for (k in 1:iN) {
            for (l in 1:iN) {
                mHtemp[k, l] <- fnH(vZ[k], vZ[l], lExtraArgsKernel)
            }
        }
        if (sMethod \%in\% c("MC", "QMC")) {
            mH <- mHtemp/iN
            mHsymm <- mH
        }
        else if (sMethod \%in\% c("TR", "GL", "CC")) {
            mH <- mHtemp \%*\% mW
            mHsymm <- sqrt(mW) \%*\% mHtemp \%*\% sqrt(mW)
        }
    }
    else {
        iN <- length(vZ)
        mHtemp <- matrix(0, iN, iN)
        for (k in 1:iN) {
            for (l in 1:iN) {
                mHtemp[k, l] <- fnH(vZ[k], vZ[l], lExtraArgsKernel)
            }
        }
        mH <- mHtemp \%*\% mW
        mHsymm <- sqrt(mW) \%*\% mHtemp \%*\% sqrt(mW)
    }
    lEigen <- eigen(mHsymm, symmetric = TRUE)
    vEigVal <- lEigen$values
    lEigen <- eigen(mH, symmetric = FALSE)
    mEigVec <- lEigen$vectors
    if (sMethod \%in\% c("MC", "QMC")) {
        mEigVecNorm <- mEigVec * sqrt(iN)
    }
    else if (sMethod \%in\% c("TR", "GL", "CC")) {
        mEigVecNorm <- mEigVec
        for (k in 1:length(vZ)) {
            dTemp <- sqrt(sum(t(vW) * mEigVec[, k] * mEigVec[, 
                k]))
            mEigVecNorm[, k] <- mEigVec[, k]/as.numeric(dTemp)
        }
    }
    return(list(nodes = vZ, weights = vW, eigval = vEigVal, eigvec = mEigVec, 
        eigvecnorm = mEigVecNorm))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

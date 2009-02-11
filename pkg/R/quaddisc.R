###-------------------------------------------------------------------
### quaddisc.R - Statistical analysis of discrepancies
###
###        code (C)  Christine Choirat and Raffaello Seri
###
###  This file:
###    created 2004-09-08
###    last modified 2009-02-11
###-------------------------------------------------------------------


##----------------------------Quasi-Monte Carlo-----------------------
## This part of the code implements the methods described in Chapter 1
## of:
##
##     K.-T. Fang and Y. Wang (1994), "Number-theoretic Methods in
##     Statistics", Monographs on Statistics and Probability 51,
##     Chapman & Hall, London.
##--------------------------------------------------------------------

intrep <- function(iK, iM) {
  ## SYNTAX: intrep(iK, iM)
  ## DESCRIPTION: Computes the integer representation of an integer iK
  ##              in a base iM
  ## INPUT: iK integer to decompose.
  ##        iM integer, number of digit of representation
  ## OUTPUT: Returns a vector of size iR containing the coefficient b_i of
  ##         the representation:
  ##                iK = b_0 + b_1*iM + b_2*iM^2 + ...+ b_iR*iM^iR
  ##         where iR is an integer computed within the function
  ##               0 <= b_i <iM.
  vB <- 0
  i <- 1
  if (iK < iM)
    vB[i] <- iK
  else {
    while (iK >= iM) {
      vB[i] <- iK %% iM
      iK <- iK %/% iM
      vB[i+1] <- 1
      i <- i+1
    }
  }
  return(vB)
}

radinv <- function(iK, iM) {
  ## SYNTAX: radinv(iK, iM)
  ## DESCRIPTION: Computes the radical inverse of an integer iK an a base iM.
  ## INPUT: iK integer to take its radical inverse
  ##        iM integer, base of radical inverse
  ## OUTPUT: a double 'dY(iM)' :
  ##              dY(iM) = b_0*iM^(-1) + b_1*iM^(-2)+ ... + b_iR*iM^(-iR-1)
  ##         where the b-i's and iR are obtained through the function
  ##         intrep(iK, iM).
  dY <- 0
  vB <- intrep(iK, iM)
  iR <- length(vB)
  for (i in 1:iR)
    dY <- dY + vB[i] * iM^(-i)
  return(dY)
}

hammset <- function(iN, vP=0) {
  ## SYNTAX: This function has two different syntaxes, hammset(iN) and hammset(iN, vP)
  ## DESCRIPTION: Computes a matrix  of Hammersley points.
  ##   1) hammset(iN)
  ## INPUT: iN number of points to be comptuted
  ## OUTPUT: a (iN)-vector
  ##
  ##   2) hammset(iN, vP)
  ## INPUT: vP vector of distinct prime numbers, e.g. <2; 3; 5; 7>
  ##        iN numbers of vectors with length(vP) components
  ## OUTPUT: a (iN) x (lenght(vP)+1) matrix: the Hammersley set.
  if (vP == 0) {
    mX <- matrix(nrow=iN, ncol=1)
    for (k in 1:iN)
      mX[k, 1] <-  (2*k - 1) / (2*iN)
    return(mX)
  }
  else {
    iD <- length(vP) + 1
    mX <- matrix(nrow=iN, ncol=iD)
    for (k in 1:iN) {
      mX[k, 1] <- (2*k - 1) / (2*iN)
      for (s in 2:iD) {
        dY <- radinv(k, vP[s-1])
        mX[k, s] <- dY
      }
    }
    return(mX)
  }
}

haltset <- function(iNmin, iNmax, vP) {
  ## SYNTAX: haltset(iNmin, iNmax, vP)
  ## DESCRIPTION: Computes a matrix  of Halton points.
  ## INPUT: vP vector of distinct prime numbers, e.g. c(2, 3, 5, 7)
  ##        iNmin integer, starting line satisfying 0 < iNmin
  ##        iNmax integer, ending line satisfying iNmin <= iNmax
  ## OUPUT: a (iNmax - iNmin + 1) x (length(vP)) matrix: the H-set.
  ## REMARK: the iN initial Halton vectors with sizec(vP) components
  ##         are computed with: haltset(vP, 1, iN).
  ##         The subesequent iN vectors are computed with:
  ##                   haltset(vP, iN+1, 2*iN).
  iD <- length(vP)
  iN <- iNmax - iNmin + 1
  mX <- matrix(nrow=iN, ncol=iD)
  for (k in 1:iN) {
    for (s in 1:iD) {
      mX[k, s] <- radinv(k+iNmin, vP[s])
    }
  }
  return(mX)
}

##-----------------------------END Quasi-Monte Carlo------------------

##------------------------------SAMPLE CODE --------------------------
##  iN <- 10; iNmin <- 10; iNmax <- 20; vP <- c(2, 3, 7)
##  mXham <- hammset(iN, vP)
##  mXhal <- haltset(iNmin, iNmax, vP)
##
##  print("Hammersley points "); print(mXham)
##  print("Halton points: "); print(mXhal)
##-----------------------------END SAMPLE CODE------------------------

##-----------------------Generalized Discrepancies--------------------
## This part of the code implements the methods described in the
##  article:
##
##    "Statistical Properties of Generalized Discrepancies",
##            by Christine Choirat and Raffaello Seri
##--------------------------------------------------------------------

##--------------------------------------------------------------------
## fGDStar, fGDSymm, fGDCent calculate the function f
## INPUT: iD integer
##        vXk, vXl vectors
## OUTPUT: a scalar

fGDStar <- function(vXk, vXl)
  .C("fGDStar",
     as.integer(length(vXk)),
     as.double(vXk),
     as.double(vXl),
     dResult=as.double(0),
     PACKAGE="quaddisc"
    )$dResult

fGDSymm <- function(vXk, vXl)
  .C("fGDSymm",
     as.integer(length(vXk)),
     as.double(vXk),
     as.double(vXl),
     dResult=as.double(0),
     PACKAGE="quaddisc"
   )$dResult

fGDCent <- function(vXk, vXl)
  .C("fGDCent",
     as.integer(length(vXk)),
     as.double(vXk),
     as.double(vXl),
     dResult=as.double(0),
     PACKAGE="quaddisc"
   )$dResult

fGD <- function(vXk, vXl, sType) {
  dResult <- 1
  if (sType == "star")
    dResult <- fGDStar(vXk, vXl)
  else if (sType == "symm")
    dResult <- fGDSymm(vXk, vXl)
  else if (sType == "cent")
    dResult <- fGDCent(vXk, vXl)
  else
    return("Error: wrong type")
  return(dResult)
}

##--------------------------------------------------------------------
## gGDStar, gGDSymm, gGDCent calculate the function g_1
## INPUT: iD integer
##        vXk vector
## OUTPUT: a scalar

gGDStar <- function(vXk)
  .C("gGDStar",
     as.integer(length(vXk)),
     as.double(vXk),
     dResult=as.double(0),
     PACKAGE="quaddisc"
   )$dResult

gGDSymm <- function(vXk)
  .C("gGDSymm",
     as.integer(length(vXk)),
     as.double(vXk),
     dResult=as.double(0),
     PACKAGE="quaddisc"
   )$dResult

gGDCent <- function(vXk)
  .C("gGDCent",
     as.integer(length(vXk)),
     as.double(vXk),
     dResult=as.double(0),
     PACKAGE="quaddisc"
   )$dResult

gGD <- function(vXk, sType) {
  dResult <- 1
  if (sType == "star")
    dResult <- gGDStar(vXk)
  else if (sType == "symm")
    dResult <- gGDSymm(vXk)
  else if (sType == "cent")
    dResult <- gGDCent(vXk)
  else
    return("Error: wrong type")
  return(dResult)
}

##--------------------------------------------------------------------
## Definition of hGD: kernel
hGD <- function(vX, vY, sType) {
  fnK <- function(sType) {
    if (sType == "cent")
      return(13/12)
    else
      return(4/3)
  }
  return(fGD(vX, vY, sType) - gGD(vX, sType) - gGD(vY, sType) + fnK(sType)^length(vX))
}

##--------------------------------------------------------------------
## Definition of hAD, Anderson and Darling: kernel
hCvM <- function(dX, dY, lExtraArgsKernel=NULL) {
  return( 0.5*(dX^2+dY^2)-max(dX,dY)+1/3 )
  ##  return((max(dX,dY)-dX*dY))##/sqrt(dX*dY*(1-dX)*(1-dY)))
}

##--------------------------------------------------------------------
## Definition of hAD, Anderson and Darling: kernel
hAD <- function(dX, dY, lExtraArgsKernel=NULL) {
  return(-log(max(dX,dY)-dX*dY)-1)
##  return((max(dX,dY)-dX*dY))##/sqrt(dX*dY*(1-dX)*(1-dY)))
}

##--------------------------------------------------------------------
## Definition of hW, Watson: kernel
hW <- function(dX, dY, lExtraArgsKernel=NULL) {
  return(.5*(dX - dY)^2 - .5*(dX + dY) + min(dX, dY) + 1/12)
}

##--------------------------------------------------------------------
## Definition of hChi, ChiSquared: kernel
hChi <- function(dX, dY, lExtraArgsKernel=NULL) {
  iJ <- 5
  if ((dX %/% (1/iJ) - dY %/% (1/iJ)) == 0)
   return(1)
  else return(0)
}

##--------------------------------------------------------------------
## Definition of hDiaphony, Diaphony: kernel
mod <- function(dX, dY) {
  if (dY == 0)
    return("Error: second argument must be != 0 in function mod")
  else {
    iResult <- dX %% dY
    if (iResult < 0)
      iResult <- iResult + abs(dY)
    return(iResult)
  }
}
hDiaphony <- function(vX, vY, lExtraArgsKernel=NULL) {
  f <- function(vX)
    -1 + prod(1 + 2 * pi^2 * (vX^2 - vX + 1/6))
  iD <- length(vX)
  vZ <- vX - vY
  for (i in 1:iD)
    vZ[i] <- mod(vZ[i], 1)
  dResult <- (1 / ((1 + pi^2 / 3)^iD - 1)) * f(vZ)
  return(dResult)
}

##------------------------------SAMPLE CODE --------------------------
## mZ <- haltset(1, 100, 2)
## mA <- matH.QMC(mZ, hGD, "star")
## lEig <- eigen(mA, symmetric=TRUE)
## vX <- sapply(1:100, function(i) 1/(pi^2 * i^2))
## plot(vX, lEig$values, type="b", col="red")
##
## plot(1:100, lEig$values, type="b", col="red")
## lines(1:100, vX, col="blue")
##-----------------------------END SAMPLE CODE------------------------

##-----------------------------Gaussian chaos-------------------------
## This part of the code implements a gaussian chaos. It includes, as a
## dynamic library, an unaltered version of 'qf.C' by R. Davies.
##--------------------------------------------------------------------

qfR <- function(vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt=100000, dAcc=1e-5, vTrace=rep(0,7), iFault=0, dValue=0) {
  ## SYNTAX: qfR(vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt, dAcc, vTrace, iFault, dValue)
  ## DESCRIPTION: Computes the cdf of quadratric forms in normal random variables.
  ## INPUT: vLambda vector of weights of the chi-square terms
  ##        vNoncent vector of noncentrality parameters of the chi-square terms
  ##        vDf vector of degrees of freedom of the chi-square terms
  ##        iTermsSum common length of the previous vectors
  ##        dSigma weight of the normal term
  ##        dC number in which the cdf is calculated
  ##        iMaxTermsInt maximal number of terms in the C integration routine
  ##        dAcc accuracy of the computation
  ##        vTrace vector of intermediate computations
  ##        iFault error code
  ##               1 "no error"
  ##               2 "requested accuracy could not be obtained"
  ##               3 "round-off error possibly significant"
  ##               4 "invalid parameters"
  ##               5 "unable to locate integration parameters"
  ##        dValues returned value of the cdf
  ## OUPUT: none (procedure and not function). Side-effects: vTrace, iFault, dValue.
  ## REMARK: R interface to R.B. Davies' code.
  ## REFERENCES: ??Davies 1973??
  .C("qfR",
     vLambda=as.double(vLambda),
     vNoncent=as.double(vNoncent),
     vDf=as.integer(vDf),
     iTermsSum=as.integer(iTermsSum),
     dSigma=as.double(dSigma),
     dC=as.double(dC),
     iMaxTermsInt=as.integer(iMaxTermsInt),
     dAcc=as.double(dAcc),
     vTrace=as.double(vTrace),
     iFault=as.integer(iFault),
     dResult=as.double(iFault),
     PACKAGE="quaddisc"
    )
}

pquadform <- function(dC, vLambda, vNoncent, vDf, dSigma, ...) {
  ## SYNTAX: pquadform(dC, vLambda, vNoncent, vDf, dSigma, ...)
  ##         where '...' are additional parameters passed to 'qfR'
  ## DESCRIPTION: Computes the cdf of quadratric forms in normal random variables.
  ## INPUT: dC number in which the cdf is calculated
  ##        vLambda vector of weights of the chi-square terms
  ##        vNoncent vector of noncentrality parameters of the chi-square terms
  ##        vDf vector of degrees of freedom of the chi-square terms
  ##        dSigma weight of the normal term
  ##        vTrace vector of intermediate computations
  ## OUPUT: a list of class 'quadform' with elements:
  ##        vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt0, dAcc, vTrace, iFault, dValue
  ##        as described in 'qfR'
  ## REMARK: User-friendly version of qfR. Computes the cdf of quadratric forms in normal random variables.
  iTermsSum <- length(vLambda)
  if (length(vNoncent) != iTermsSum)
    return("Error: wrong number of noncentrality coefficients vNoncent")
  if (length(vDf) != iTermsSum)
    return("Error: wrong number of degrees-of-freedom coefficients vDf")
  lResult <- qfR(vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, ...)
  class(lResult) <- "quadform"
  return(lResult)
}

print.quadform <- function(lQuadform) {
  ## SYNTAX: print.quadform(lQuadform)
  ##         print(lQuadform)
  ## DESCRIPTION: 'print' method for objects of class 'quadform': pretty-prints the output of 'pquadform'
  ## INPUT: lQuadform a list of class 'quadform' as obtained by the 'pquadform' function.
  ## OUPUT: Displays the value of the cdf and the associated error code.
  ## REMARK: ??.
  vError <- character(5)
  vError[1] <- "no error"
  vError[2] <- "requested accuracy could not be obtained"
  vError[3] <- "round-off error possibly significant"
  vError[4] <- "invalid parameters"
  vError[5] <- "unable to locate integration parameters"
  cat(sep="", "value:\t", lQuadform$dResult, "\n")
  cat(sep="", "fault code:\t", lQuadform$iFault, " (", vError[lQuadform$iFault + 1], ")\n")
}

pgausschaos <- function(dC, vLambda, ...) {
  ## SYNTAX: pgausschaos(dC, vLambda, ...)
  ##         where '...' are additional parameters passed to 'qfR'
  ## DESCRIPTION: Computes the cdf of a second order Gaussian chaos in normal random variables.
  ## INPUT: dC number in which the cdf is calculated
  ##        vLambda vector of weights of the chi-square terms
  ## OUPUT: a list of class 'quadform' with elements:
  ##        vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt0, dAcc, vTrace=, iFault, dValue
  ##        as described in 'qfR'
  ## REMARK: Special case of 'pquadform'.
  ## REFERENCES: See ?? for definitions and main properties.
  iLen <- length(vLambda)
  vNoncent <- rep(0, iLen)
  vDf <- rep(1, iLen)
  dSigma <- 0
  return(pquadform(dC, vLambda, vNoncent, vDf, dSigma,  ...))
}

qgausschaos <- function(dP, vLambda, lArgsUniroot=list(interval=c(0, 100)), ...) {
  ## SYNTAX: qgausschaos(dP, vLambda, ...)
  ##         where '...' are additional parameters passed to 'qfR'
  ## DESCRIPTION: Computes the quantiles of a second order Gaussian chaos in normal random variables.
  ## INPUT: dP number between 0 and 1 (excluded) to get the quantile from
  ##        vLambda vector of weights of the chi-square terms
  ## OUPUT: a list, output of R build-in function 'uniroot'
  if ( (dP <= 0) | (dP >= 1) )
    return("wrong argument (should be in (0, 1))")
  fnObj <- function(dC) {
    lP <- pgausschaos(dC, vLambda)
    if (lP$iFault !=0)
      print("warning: error in computation of pgausschaos")
    return(lP$dResult - dP)
  }
  lArgsUniroot <- c(f=fnObj, lArgsUniroot)
  return(do.call("uniroot", lArgsUniroot))
}

padqd.QMC <- function(dC, fnH, mZ, sForm="v", lExtraArgsKernel,  ...) {
  ## SYNTAX: padqd(dC, fnH, mZ, sForm, lExtraArgsKernel)
  ##         where '...' are additional parameters passed to 'qfR'
  ## DESCRIPTION: Computes the cdf at dC of the asymptotic distribution of quadratic discrepancies
  ##              with kernel fnH based on points mZ.
  ## INPUT: dC number in which the cdf is calculated
  ##        fnH kernel function of the degenerate V-statistic
  ##        mZ integration nodes in which the kernel is evaluated
  ##        sForm string stating whether the approximate integral operator should include the constant
  ##              ("v") or not ("u").
  ##        lExtraArgsKernel list of extra arguments passed to fnH
  ## OUPUT: a list of class 'quadform' with elements:
  ##        vLambda, vNoncent, vDf, iTermsSum, dSigma, dC, iMaxTermsInt0, dAcc, vTrace, iFault, dValue
  ##        as described in 'qfR'
  ## REMARK: Special case of 'pgausschaos'.
  if (!is.null(lExtraArgsKernel))
    lArgs <- c(list(mZ=mZ, fnH=fnH), lExtraArgsKernel)
  else
    lArgs <- list(mZ=mZ, fnH=fnH)
  if (sForm == "u")
    mH <- do.call("matHu.QMC", lArgs)
  if (sForm == "v")
    mH <- do.call("matH.QMC", lArgs)
  vLambda <- eigen(mH, symmetric=TRUE)$values
  return(pgausschaos(dC=dC, vLambda=vLambda, ...))
}

qadqd.QMC <- function(dP, fnH, mZ, sForm="v", lExtraArgsKernel, lArgsUniroot=list(interval=c(0, 100)),  ...) {
  ## SYNTAX: qadqd(dP, fnH, mZ, sForm="v", lExtraArgsKernel,  ...)
  ##         where '...' are additional parameters passed to 'qfR'
  ## DESCRIPTION: Computes the cdf at dC of the asymptotic distribution of quadratic discrepancies
  ##              with kernel fnH based on points mZ.
  ## INPUT: dP number between 0 and 1 (excluded) to get the quantile from
  ##        fnH kernel function of the degenerate V-statistic
  ##        mZ integration nodes in which the kernel is evaluated
  ##        dC number in which the cdf is calculated
  ##        sForm string stating whether the approximate integral operator should include the constant
  ##              ("v") or not ("u").
  ##        lExtraArgsKernel list of extra arguments passed to fnH
  ## OUPUT: a list, output of R build-in function 'uniroot'
  if ( (dP <= 0) | (dP >= 1) )
    return("wrong argument (should be in (0, 1))")
  if (!is.null(lExtraArgsKernel))
    lArgs <- c(list(mZ=mZ, fnH=fnH), lExtraArgsKernel)
  else
    lArgs <- list(mZ=mZ, fnH=fnH)
  if (sForm == "u")
    mH <- do.call("matHu.QMC", lArgs)
  if (sForm == "v")
    mH <- do.call("matH.QMC", lArgs)
  vLambda <- eigen(mH, symmetric=TRUE)$values
  fnObj <- function(dC) {
    lP <- pgausschaos(dC, vLambda, ...)
    if (lP$iFault !=0)
      print("warning: error in computation of pgausschaos")
    return(lP$dResult - dP)
  }
  lArgsUniroot <- c(f=fnObj, lArgsUniroot)
  return(do.call("uniroot", lArgsUniroot))
}


##------------------------------SAMPLE CODE --------------------------
## pquadform(1.0, c(6.0, 3.0, 1.0), c(0.0, 0.0, 0.0), c(1, 1 ,1 ), 0.0)
## mZ <- haltset(1, 100, 2)
## mA <- matH.QMC(mZ, hGD, "star")
## vLambda <- eigen(mA, symmetric=TRUE)$values
## vX <- sapply(1:100, function(i) 1/(pi^2 * i^2))
## plot(vX, lEig$values, type="b", col="red")
##
## plot(1:100, lEig$values, type="b", col="red")
## lines(1:100, vX, col="blue")
##-----------------------------END SAMPLE CODE------------------------


##----------------------Clenshaw-Curtis Integration-------------------
## This part of the code computes the Clenshaw-Curtis integration
##--------------------------------------------------------------------
## Definition of CC.quad: matrix approximation of the integral operator
CC.quad <- function(iN) {
  vI <- matrix(0:(iN - 1), iN, 1)
  vZ <- as.vector( 0.5 - 0.5 * cos(((2 * vI + 1) * pi) / (2 * iN)) )
  vD <- c(1, rep(2, iN-1)) / iN
  vK <- matrix(2*(1:iN)-1, 1, iN)
  ## mCt is the transpose of the matrix C
  mCt <- cos( pi * (vI %*% vK) / (2*iN) )
  vS <- as.matrix(1:iN)
  vS <- -2/(vS * (vS - 2))
  vS[seq(2, length(vS), by=2)] <- 0
  vW <- as.vector(t(0.5 * crossprod(vS*vD, mCt)))
  return(list(nodes=vZ, weights=vW))
}


##------------Matrix Approximation to the Integral Operator-----------
## This part of the code computes the matrix approximation of an
## integral operator.
##--------------------------------------------------------------------
## Definition of matH: matrix approximation of the integral operator
matH <- function(fnH, iN, sMethod, bPoints=FALSE, vZ=NULL, vW=NULL, lExtraArgsKernel=NULL) {
  if (!bPoints) {
    if (sMethod=="MC") {
      vZ <- runif(iN)
      vW <- rep(1 / iN, iN)
    } else if (sMethod=="QMC") {
      vZ <- hammset(iN)
      vW <- rep(1 / iN, iN)
    } else if (sMethod=="TR") {
      vZ <- as.matrix(0:(iN-1))/(iN-1)
      vW <- c(0.5, matrix(1, iN - 2, 1), 0.5) / (iN-1)
    } else if (sMethod=="GL") {
      lGL <- gauss.quad(iN)
      vZ <- (lGL$nodes + 1) / 2
      vW <- (lGL$weights) / 2
    } else if (sMethod=="CC") {
      lCC <- CC.quad(iN)
      vZ <- lCC$nodes
      vW <- lCC$weights
    }
    mW <- diag(as.vector(vW))
    mHtemp <- matrix(0, iN, iN)
    for (k in 1:iN) {
      for (l in 1:iN) {
        mHtemp[k,l] <- fnH(vZ[k], vZ[l], lExtraArgsKernel)
      }
    }
    if (sMethod %in% c("MC","QMC")) {
      mH <- mHtemp / iN
      mHsymm <- mH
    } else if (sMethod %in% c("TR","GL","CC")) {
      mH <- mHtemp %*% mW
      mHsymm <- sqrt(mW) %*% mHtemp %*% sqrt(mW)
    }
  }
  else {
    iN <- length(vZ)
    ## Check that length(vZ) == length(vW) ??
    mHtemp <- matrix(0, iN, iN)
    for (k in 1:iN) {
      for (l in 1:iN) {
        mHtemp[k, l] <- fnH(vZ[k], vZ[l], lExtraArgsKernel)
      }
    }
    mH <- mHtemp %*% mW
    mHsymm <- sqrt(mW) %*% mHtemp %*% sqrt(mW)
  }
  ## We use symmetric=TRUE since the eigenvalues are computed with greater
  ## precision using a symmetric version of the matrix
  lEigen <- eigen(mHsymm, symmetric=TRUE)
  vEigVal <- lEigen$values
  ## We use symmetric=FALSE since the eigenfunctions have to be computed
  ## using an unsymmetric version of the matrix
  lEigen <- eigen(mH, symmetric=FALSE)
  mEigVec <- lEigen$vectors
  ## We normalize the eigenvectors as in the paper
  if (sMethod %in% c("MC","QMC")) {
    mEigVecNorm <- mEigVec * sqrt(iN)
  } else if (sMethod %in% c("TR","GL","CC")) {
  ##    vNorm <- sqrt(crossprod(vW, mEigVec^2))
  ##    mEigVecNorm <- mEigVec / (matrix(1, iN, 1) %*% vNorm)
  mEigVecNorm <- mEigVec
  for (k in 1:length(vZ)) {
    dTemp <- sqrt( sum( t(vW) * mEigVec[, k] * mEigVec[, k] ) )
    mEigVecNorm[, k] <- mEigVec[, k] / as.numeric(dTemp)
  }
  }
  ##  return(list(symm=mHsymm, unsymm=mH, eigval=vEigVal, eigvec=mEigVec, eigvecnorm=mEigVecNorm))
  return(list(nodes=vZ, weights=vW, eigval=vEigVal, eigvec=mEigVec, eigvecnorm=mEigVecNorm))
}


eigfun <- function(vX, fnH, iN, sMethod, bPoints=FALSE, vZ=NULL, vW=NULL, lExtraArgsKernel=NULL, vJ=1) {
  lH <- matH(fnH, iN, sMethod, bPoints=FALSE, vZ=NULL, vW=NULL)
  vZ <- as.matrix(lH$nodes)
  mResult <- matrix(nrow=length(vX), ncol=length(vJ))
  vW <- as.matrix(lH$weights)
  mEigVec <- as.matrix(lH$eigvecnorm[,vJ])

  if (sMethod %in% c("MC","QMC","TR","GL")) {
    mH <- matrix(nrow=length(vX), ncol=length(vZ))
    for (j in 1:length(vZ)) {
      for (i in 1:length(vX)) {
        mH[i,j] <- fnH(vX[i], vZ[j], lExtraArgsKernel)
      }
    }
    vEigVal <- as.matrix(lH$eigval[vJ])
    mResult <- mH %*% ((vW %*% t(1 / vEigVal)) * mEigVec)
  }

  else if (sMethod == "CC") {
    vI1 <- matrix(1:iN, iN, 1)
    vI0 <- vI1 - 1
    vD <- c(1, rep(2, iN - 1)) / iN
    vK <- matrix(2 * vI0 + 1, 1, iN)
    ## mCt is the transpose of the matrix C
    mCt <- cos( pi * (vI0 %*% vK) / (2*iN) )
    ## mTt is the transpose of the matrix T
    mTt <- matrix(nrow=length(vX), ncol=iN)
    mTt <- cos( acos(2 * vX - 1) %*% matrix(vI0, 1, iN))
    mResult <- mTt %*% diag(as.vector(vD)) %*% mCt %*% mEigVec
  }
  return(mResult)
}

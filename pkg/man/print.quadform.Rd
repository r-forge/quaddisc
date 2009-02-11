\name{print.quadform}
\alias{print.quadform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.quadform(lQuadform)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lQuadform}{ ~~Describe \code{lQuadform} here~~ }
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
function (lQuadform) 
{
    vError <- character(5)
    vError[1] <- "no error"
    vError[2] <- "requested accuracy could not be obtained"
    vError[3] <- "round-off error possibly significant"
    vError[4] <- "invalid parameters"
    vError[5] <- "unable to locate integration parameters"
    cat(sep = "", "value:\t", lQuadform$dResult, "\n")
    cat(sep = "", "fault code:\t", lQuadform$iFault, " (", vError[lQuadform$iFault + 
        1], ")\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

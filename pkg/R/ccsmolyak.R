## The R interface to Knut Petras' SMOLPACK is mostly taken from the
## R package 'gss' by Chong Gu (file "smolyak.R").

ccsmolyak.quad <- function(d, k) {
  size <- .C("size_smolyak",
             as.integer(d),
             as.integer(d + k),
             size=integer(1),
	     PACKAGE="quaddisc"
             )$size
  z <- .C("quad_smolyak",
          as.integer(d),
          as.integer(d + k),
          pt=double(d * size),
          wt=as.double(1:size),
          PACKAGE="quaddisc"
          )
  list(pt=t(matrix(z$pt, d, size)), wt=z$wt)
}

ccsmolyak.size <- function(d, k) {
  .C("size_smolyak",
     as.integer(d),
     as.integer(d + k),
     size=integer(1),
     PACKAGE="quaddisc"
     )$size
}

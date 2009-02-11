## .First.lib <- function(lib, pkg) {
##     require(statmod, quietly = TRUE)
##     library.dynam("quaddisc", pkg, lib)
## }

.First.lib <-
  function (libname, pkgname) {
    ## figure out this year automatically
    this.year <- substr(as.character(Sys.Date( )), 1, 4)
    ## echo output to screen
    cat("##\n## quaddisc - Quadratic discrepancies and related quantities\n")
    cat("## Copyright (C) 2009-", this.year,
        " C. Choirat and R. Seri\n", sep="")
    cat("##\n")
    ## assuming you would need the package MASS for your package
    require(statmod, quietly=TRUE)
    library.dynam(pkgname, pkgname, lib.loc=libname)
  }

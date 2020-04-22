.packageName <- "rcellminer"

#' Make sure that rcellminerData is loaded
#' @param libname a character string giving the library directory where the package defining the namespace was found.
#' @param pkgname a character string giving the name of the package. 
#' 
#' @import rcellminerData
.onLoad <- function(libname, pkgname) {
}

#' Display citation message
#' @param libname a character string giving the library directory where the package defining the namespace was found.
#' @param pkgname a character string giving the name of the package. 
.onAttach <- function(libname, pkgname) {
    packageStartupMessage('Consider citing this package: Luna A, et al. rcellminer: exploring molecular profiles and drug response of the NCI-60 cell lines in R. PMID: 26635141; citation("rcellminer")')
}

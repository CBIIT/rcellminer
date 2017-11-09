.packageName <- "rcellminer"

#' Make sure that rcellminerData is loaded
#' 
#' @import rcellminerData
.onLoad <- function(lib, pkg){
}

.onAttach <- function(libname, pkgname){
    packageStartupMessage('Consider citing this package: Luna A, et al. rcellminer: exploring molecular profiles and drug response of the NCI-60 cell lines in R. PMID: 26635141; citation("rcellminer")')
}

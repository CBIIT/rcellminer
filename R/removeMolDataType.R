#' Remove molecular data type prefixes from features. 
#' 
#' @param features A vector of features. 
#' @param prefixLen The length of the molecular data type prefix.
#' 
#' @details This function is primarily used to remove prefixes from elastic 
#'   net features.
#'   
#' @return A named vector of features without molecular data type prefixes.
#' 
#' @examples 
#' removeMolDataType(c("expTP53", "copMDM2", "mutCHEK2", "mutBRAF"))
#' 
#' @concept rcellminer
#' @export
removeMolDataType <- function(features, prefixLen=3){
    trimmedFeatures <- vapply(features, function(s) str_sub(s, start = (prefixLen+1), end = nchar(s)),
                             character(1))
    return(trimmedFeatures)
}

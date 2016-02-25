#' Check if NSC has Mechanism of Action (MOA) Annotation
#' 
#' @param nsc a string, an NSC identifier
#' 
#' @return a boolean whether the NSC has an MOA
#' 
#' @examples
#' hasMoa("754365")
#' 
#' @concept rcellminer
#' @export
hasMoa <- function(nsc) {
  if(!is.null(getDrugMoaList(nsc))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

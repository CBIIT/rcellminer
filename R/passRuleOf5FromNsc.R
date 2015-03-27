#' Checks if NSC passes Lipinski's Rule of 5 
#' 
#' @param nsc a string, the NSC identifier to be checked 
#' @param acceptableViolations, a number, the number of acceptable rule violations (default: 0)
#' @param verbose a boolean, whether to write out the failing criteria (default: FALSE)
#' @return a boolean or NA, whether the NSC passes the criteria or NA if the NSC had no structure
#' 
#' @details Uses RCDK: org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor
#' 
#' @examples
#' passRuleOf5FromNsc("94600", verbose=TRUE)
#' 
#' @concept rcellminer
#' @export
passRuleOf5FromNsc <- function(nsc, acceptableViolations=0, verbose=FALSE) {
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
  smiles <- drugAnnot[nsc, "SMILES"]

  if(!is.na(smiles)) {
    passRuleOf5(smiles, acceptableViolations, verbose)
  } else {
    return(NA)
  }
}



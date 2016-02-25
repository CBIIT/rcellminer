#' Checks if SMILES passes Lipinski's Rule of 5 
#' 
#' @param smiles a string, the SMILES structure to be checked 
#' @param acceptableViolations, a number, the number of acceptable rule violations (default: 0)
#' @param verbose a boolean, whether to write out the failing criteria (default: FALSE)
#' @return a boolean, whether the SMILES passes the criteria 
#' 
#' @details Uses RCDK: org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor
#' 
#' @examples
#' # Docetaxel
#' passRuleOf5("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC
#' (=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O", verbose=TRUE)
#' 
#' # Gemcitabine 
#' passRuleOf5("C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)(F)F", verbose=TRUE)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom rcdk eval.desc parse.smiles
passRuleOf5 <- function(smiles, acceptableViolations=0, verbose=FALSE) {
	tryCatch({
		tmp <- parse.smiles(smiles)[[1]]
		failures <- eval.desc(tmp,"org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor")
		
		if(verbose) {
			cat("Rule of 5 Violations:", as.numeric(failures), "\n")
		}
		
		if(failures-acceptableViolations <= 0) {
			return(TRUE) 
		} else {
			return(FALSE)
		}		
	}, error = function(e) {
    return(NA)
  })
}



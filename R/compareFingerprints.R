#' Compare Structure Fingerprints to NCI DTP Compounds
#' 
#' This function compares the first SMILES structure to all other structures given.
#' 
#' @param ids a vector of IDs corresponding to structures
#' @param smiles a vector of strings SMILES structures 
#' @param fpType the type of fingerprint to be used; uses the RCDK 
#'   get.fingerprint() (default: standard)
#' @param verbose a boolean whether to display debugging information
#' @param fingerprint.list a list of fingerprints generated with \code{getFingerprintList}
#' 
#' @return a sorted named vector of Tanimoto distances  
#' 
#' @seealso rcdk::get.fingerprint
#' 
#' @examples
#' drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
#' ids <- head(drugAnnot$NSC)
#' smiles <- head(drugAnnot$SMILES)
#' fingerprintList <- getFingerprintList(ids, smiles)
#' compareFingerprints(fingerprint.list=fingerprintList)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom fingerprint distance
compareFingerprints <- function(ids=NULL, smiles=NULL, fpType="standard", 
																verbose=TRUE, fingerprint.list=NULL) {
	found <- NULL 
	
	if(is.null(fingerprint.list)) {
		fingerprint.list <- getFingerprintList(ids, smiles, fpType=fpType, verbose)
	}
	
	tmp <- sapply(names(fingerprint.list), function(x) {!is.null(fingerprint.list[[x]])})
	found <- which(tmp)
	
	#fp.sim <- fp.sim.matrix(fingerprint.list, method = "tanimoto")
	#fp.dist <- 1 - fp.sim
	
	results <- NULL 
	
	if(verbose) {
  	cat("\nGet Distances\n")
	  pb <- txtProgressBar(min=1, max=length(fingerprint.list), style=3)
	}
	
	for(i in seq_along(fingerprint.list)) {
		tmp <- distance(fingerprint.list[[i]], fingerprint.list[[1]])
		results <- c(results, tmp)
		
		if(verbose) {
  		setTxtProgressBar(pb, i)
		}
	}
	
	names(results) <- names(found)
	results <- sort(results, decreasing=TRUE)
	
	return(results)
}


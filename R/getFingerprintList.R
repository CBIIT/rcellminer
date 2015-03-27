#' Get a list of fingerprints for a set of compounds
#'
#' @param ids a vector of IDs corresponding to structures
#' @param smiles a vector of strings SMILES structures
#' @param fpType the type of fingerprint to be used; uses the RCDK
#'   get.fingerprint() (default: standard)
#' @param verbose show debugging output
#' @return a list of fingerprints returned from RCDK where the list names are
#'   the IDs
#'
#' @seealso rcdk::get.fingerprint
#'
#' @examples
#' drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
#' ids <- head(drugAnnot$NSC)
#' smiles <- head(drugAnnot$SMILES)
#' fingerprintList <- getFingerprintList(ids, smiles)
#'
#' \dontrun{
#' # All fingerprints
#' ids <- drugAnnot$NSC
#' smiles <- drugAnnot$SMILES
#' fingerprintList <- getFingerprintList(ids, smiles)
#' save(fingerprintList, file="fingerprintList.RData")
#' }
#'
#' @concept rcellminer
#' @export
#'
#' @importFrom rcdk get.fingerprint
getFingerprintList <- function(ids, smiles, fpType="standard", verbose=TRUE) {
	fingerprint.list <- list()

	if(verbose && length(ids) > 1) {
		cat("Get Fingerprints\n")
		pb <- txtProgressBar(min=1, max=length(ids), style=3)
	}

	for(i in seq_along(ids)) {
		id <- ids[i]
		smile <- smiles[i]

		fp <- tryCatch({
			mol <- parse.smiles(smile)[[1]]
			get.fingerprint(mol, type=fpType)
		}, error=function(e) { NULL })

		if(!is.null(fp)) {
			fingerprint.list[[as.character(id)]] <- fp
		}

		if(verbose && length(ids) > 1) {
			setTxtProgressBar(pb, i)
		}
	}

	return(fingerprint.list)
}

#' A list of pre-computed fingerprints using getFingeprintList() 
#'
#' @name fingerprintList
#' @docType data
#' 
#' @concept rcellminer
#' @keywords data
NULL

#' Plot the structures for NSCs
#' 
#' @param nscs a vector of strings of NSCs used as structure titles
#' 
#' @details This is a wrapper for the plotStructures() function that takes only NSCs
#' 
#' @return the function does not return anything 
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @examples
#' plotStructuresFromNscs("94600")
#' plotStructuresFromNscs(c("609699", "94600"))
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom rcdk parse.smiles
plotStructuresFromNscs <- function(nscs) {
	drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
	plotStructures(nscs, drugAnnot[nscs, "SMILES"], mainLabel=nscs)
}

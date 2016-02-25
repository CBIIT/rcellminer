#' Plot Structures
#' 
#' @param ids a vector of strings of IDs used as structure titles
#' @param smiles a vector of strings where the strings are SMILES structures
#' @param titleCex a number, the scaling factor for the title (default: 1)
#' @param structSize a number, the size of the structure image (default: 200)
#' @param structAnnotPos a number, how far above the structure to display the title (default: 50)
#' @param mainLabel a string, the main plot label
#' @param rows number of rows in figure (default: 1)
#' @param cols number of columns in figure (default: input structures number)
#' 
#' @details The parameter ids is used as a title, this function does not search for IDs, but works
#'   based off the smiles given. This is a wrapper around rcdkplot for plotting multiple structures.
#' 
#' @return the function does not return anything 
#' 
#' @author Augustin Luna <augustin AT mail.nih.gov>
#' 
#' @examples
#' drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
#' plotStructures("94600", drugAnnot["94600","SMILES"])
#' plotStructures(c("609699", "94600"), drugAnnot[c("609699", "94600"),"SMILES"], 
#'   mainLabel=c("609699", "94600"))
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom rcdk parse.smiles
plotStructures <- function(ids, smiles, titleCex=1, 
													 structSize=300, structAnnotPos=50, 
													 mainLabel="", rows=1, cols=length(ids)) {  
	#   # Pre-process SMILES 
	#   # Remove Cl groups
	#   #smiles <- sub("^\\[Cl-\\]\\.", "", smiles)
	#   
	# Keep only the largest SMILES component
	tmp <- strsplit(smiles, "\\.")
	smilesTmp <- NULL
	
	for(i in seq_along(tmp)) {
		if(length(tmp[[i]]) == 1) {
			smilesTmp <- c(smilesTmp, smiles[i])	
		} else {
			tmp2 <- tmp[[i]][nchar(tmp[[i]]) == max(nchar(tmp[[i]]))]
			smilesTmp <- c(smilesTmp, tmp2)
		}
	}	
	
	smiles <- smilesTmp 
	
	tmp <- parse.smiles(smiles)
	
	par(mfrow=c(rows,cols))
	
	for(i in seq_along(tmp)) {
		rcdkplot(tmp[[i]], width=structSize, height=structSize, main=ids[i])			
	}
}

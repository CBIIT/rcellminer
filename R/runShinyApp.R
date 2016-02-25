#' Run Shiny App
#' 
#' @param app string Shiny app to run (See Details, OPTIONS: shinyComparePlots, shinyCompareStructures, 
#'   shinyReprodPlots, DEFAULT: shinyComparePlots)
#' @param launch.browser launch in browser? (default: TRUE)
#' @param port port number to use
#' @return None
#' 
#' @details 
#'   \itemize{
#'     \item shinyComparePlots: The "Comparison" application allows users to plot 
#'       any two variables from the CellMiner data against each other. 
#'       It additionally allows users to search for compound NSC IDs using 
#'       names and mechanisms of action.
#'     \item shinyCompareStructures: The "Compound Browser" application allows 
#'       users to see information about each compound, including structures 
#'       and any repeat assay information.
#'     \item shinyReprodPlots: The "Structure Comparison" application allows 
#'       users to identify similar compounds within the dataset either by NSC ID or SMILES string.
#'   }
#'
#' @examples
#' # Uncomment first
#' #runShinyApp()
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom shiny runApp
runShinyApp <- function(app="shinyComparePlots", launch.browser=TRUE, port=3838) {
	apps <- c("shinyComparePlots", "shinyCompareStructures", "shinyReprodPlots")
	if(app %in% apps)	{
		runApp(system.file(app, package='rcellminer'), launch.browser=launch.browser, port=port)	
	} else {
		stop(paste('"app" must be one of the following:', paste(apps, collapse=" ")))
	}
}

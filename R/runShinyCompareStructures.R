#' Run the Compare Structures Shiny App
#' 
#' @param launch.browser launch in browser? (default: TRUE)
#' @param port port number to use
#' @return None
#'
#' @examples
#' port <- 3838
#' # Uncomment first
#' #runShinyCompareStructures(port=port)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom shiny runApp
runShinyCompareStructures <- function(launch.browser=TRUE, port=3838) {
	runApp(system.file('shinyCompareStructures', package='rcellminer'), 
								launch.browser=launch.browser, port=port)	
}

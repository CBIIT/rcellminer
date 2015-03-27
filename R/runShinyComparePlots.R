#' Run the Compare Plots Shiny App
#' 
#' @param launch.browser launch in browser? (default: TRUE)
#' @param port port number to use
#' @return None
#'
#' @examples
#' port <-3838
#' # Uncomment first
#' #runShinyComparePlots(port=port)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom shiny runApp
runShinyComparePlots <- function(launch.browser=TRUE, port=3838) {
	runApp(system.file('shinyComparePlots', package='rcellminer'), 
								launch.browser=launch.browser, port=port)	
}

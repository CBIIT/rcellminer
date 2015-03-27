#' Run the Compound Browser
#' 
#' @param launch.browser launch in browser? (default: TRUE)
#' @param port port number to use
#' @return None
#' 
#' @examples
#' port <- 3838
#' # Uncomment first
#' #runShinyCompoundBrowser(port=port)
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom shiny runApp
runShinyCompoundBrowser <- function(launch.browser=TRUE, port=3838) {
	runApp(system.file('shinyReprodPlots', package='rcellminer'), 
								launch.browser=launch.browser, port=port)	
}

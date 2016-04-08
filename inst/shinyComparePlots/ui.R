library(shiny)
library(rcellminer)

#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON(system.file("shinyComparePlots", "config.json", package="rcellminer"))
appConfig <- jsonlite::fromJSON(system.file("shinyComparePlots", "appConfig.json", package="rcellminer"))
source(system.file("shinyComparePlots", "modal.R", package="rcellminer"))

if (!is.null(appConfig$appName)){
	appTitle <- appConfig$appName
} else{
	appTitle <- "CellMiner"
}

# Construct named character vector mapping displayed data source names to
# internally used source identifiers.
dataSourceChoices <- setNames(names(config),
														  vapply(config, function(x) { x[["displayName"]] }, 
														  			 character(1)))

for (configSrcId in names(config)){
	srcName <- config[[configSrcId]][["displayName"]]
	srcPackages <- names(config[[configSrcId]][["packages"]])
	
	for (pkgName in srcPackages){
		if (!require(pkgName, character.only = TRUE)){
			dataSourceChoices[srcName] <- NA
			break
		}
	}
}

if (any(is.na(dataSourceChoices))){
	stop("Check configuration file: one or more required data source packages must be installed.")
} 
#--------------------------------------------------------------------------------------------------

if("rCharts" %in% installed.packages()) {
	options(RCHART_LIB='highcharts')	
	library(rCharts)
	hasRCharts <- TRUE
} else {
	hasRCharts <- FALSE
}

colorSet <- loadNciColorSet(returnDf=TRUE)

# NOTE: Size is not automatically set for rChartsAlternative output
plotHeight <- 800
plotWidth <- 800

shinyUI(
	navbarPage(appTitle,
		inverse=FALSE,
		tabPanel("Analysis",
			fluidPage(
    		tags$head(includeScript("www/js/google-analytics.js")),
    		loadingModal(),
    		# Add/run startup Javascript
    		tags$head(
    			tags$script(onloadJs)
    		),
    		# Use JQuery (built into Shiny) calls to show/hide modal based on message
    		tags$head(
					tags$script('
		        Shiny.addCustomMessageHandler("showLoading",
		          function(message) {
		            //console.log(message);
		
		            if(message.show) {
		              $("#loadingModal").modal("show");
		            } else {
		              $("#loadingModal").modal("hide");
		            }
		          }
		        );
    			')
				),
				tags$head(
					tags$style(HTML(paste0("
						.rChart {
			  		display: block;
			  		margin-left: auto; 
			  		margin-right: auto;
			  		width: ", plotWidth, "px;
			  		height: ", plotHeight, "px;
					}")))
				),
    
	    	sidebarLayout(
	        sidebarPanel(
	        	width=3, 
	        	tags$div(
	          	id="input_container", 
	            selectInput("xDataset", "x-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
	            uiOutput("xPrefixUi"),
	            textInput("xId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "SLFN11"), 
	          	#uiOutput("xIdUi"),
	          	
	            selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
	            uiOutput("yPrefixUi"),
	          	textInput("yId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "94600"), 
	          	#uiOutput("yIdUi"),
	                
	            checkboxInput("showColor", "Show Color?", value=TRUE),
	        	  checkboxInput("selectedTissuesOnly", "Selected Tissues Only?", value=FALSE),
	                
	            uiOutput("showColorTissuesUi"),
	            
	            # Generate a hidden input with TRUE or FALSE if rCharts is installed
	        		tags$input(id="hasRCharts", type="text", value=hasRCharts, style="display:none"),
	        					 
	            # load Javascript snippet to parse the query string.
	            singleton(tags$script(type="text/javascript", src="js/parse_input.js"))
	        	)
	        ),
	        
	        mainPanel(
	        	uiOutput('tabsetPanel') 
	        )
	    	)
			)
		),
    tabPanel("About",
    	includeMarkdown("www/files/about.md") 
    )
	)
)

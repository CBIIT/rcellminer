library(shiny)
library(rcellminer)

#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON("config.json")
appConfig <- jsonlite::fromJSON("appConfig.json")
metaConfig <- jsonlite::fromJSON("configMeta.json")
source("modal.R")
source("appUtils.R")

if (!is.null(appConfig$appName)){
	appTitle <- appConfig$appName
} else{
	appTitle <- "CellMiner"
}

dataSourceChoices <- setNames(names(config),
															vapply(config, function(x) { x[["displayName"]] }, 
																		 character(1)))

metaChoices <- setNames(names(metaConfig),
												vapply(metaConfig, function(x) { x[["displayName"]] }, 
															 character(1)))

if("rCharts" %in% installed.packages()) {
	options(RCHART_LIB='highcharts')	
	library(rCharts)
	hasRCharts <- TRUE
} else {
	hasRCharts <- FALSE
}
## ---

shinyUI(
	navbarPage(appTitle, 
						 inverse=FALSE,
						 header = list(tags$head(includeCSS("www/css/hacks.css")),
						 							 #tags$head(includeCSS("www/css/tooltip.css")),
						 							 # Add/run startup Javascript
						 							 tags$head(tags$script(onloadJs)),
						 							 # Use JQuery (built into Shiny) calls to show/hide modal based on message
						 							 tags$head(includeScript("www/js/showLoading.js")),
						 							 # load Javascript snippet to parse the query string.
						 							 #tags$script(includeScript("www/js/parse_input.js")),
						 							 tags$head(includeScript("www/js/google-analytics.js")),
						 							 tags$head(
						 								 tags$style(HTML(paste0("
															 .rChart {
												  		 display: block;
												  		 margin-left: auto; 
												  		 margin-right: auto;
												  		 width: ", plotWidth, "px;
												  		 height: ", plotHeight, "px;
														 }")))
													 )),
		#------[NavBar Tab: Univariate Analyses]---------------------------------------------------------
		tabPanel("Univariate Analyses",
			fluidPage(
    		loadingModal(),
	    	sidebarLayout(
	        sidebarPanel(
	        	width=3, 
	        	tags$div(
	          	id="input_container", 
	            selectInput("xDataset", "x-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
	            uiOutput("xPrefixUi"),
	            textInput("xId", "ID: (e.g. topotecan or SLFN11)", "SLFN11"),
	          	uiOutput("xAxisRangeUi"),
	          	
	            selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
	            uiOutput("yPrefixUi"),
	          	textInput("yId", "ID: (e.g. topotecan or SLFN11)", "topotecan"),
	          	uiOutput("yAxisRangeUi"),
	          	
	            checkboxInput("showColor", "Show Color?", value=TRUE),

	          	radioButtons("tissueSelectionMode", "Select Tissues", c("Include", "Exclude")),
	          	uiOutput("selectTissuesUi"),
	            uiOutput("showColorTissuesUi"),
	            
	            # Generate a hidden input with TRUE or FALSE if rCharts is installed
	        		tags$input(id="hasRCharts", type="text", value=hasRCharts, style="display:none")
	        	)
	        ),
        mainPanel(
        	uiOutput('tabsetPanel')
        )
    	 )
			)
		),
		#-----[NavBar Tab: Regression Models]------------------------------------------------------------
		regressionModelsInput("rm", dataSourceChoices),
		#-----[NavBar Tab: Metadata]---------------------------------------------------------------------
		tabPanel("Metadata", 
						 fluidPage(	
						 	sidebarLayout(
						 		sidebarPanel(
						 			width=3, 
						 			tags$div(
						 				id="input_container", 
						 				selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
						 				#uiOutput(""),
						 			)
						 		), #end sidebarPanel
						 		mainPanel(
						 			uiOutput('metadataPanel'),
						 			h4(htmlOutput('sourceLink'))
						 		)
						 	)
						 ) #end fluidPage
		), #end tabPane 
		#-----[NavBar Tab: About]------------------------------------------------------------------------
    tabPanel("About",
    	includeMarkdown("www/files/about.md")
    	#h1("For testing"),
    	#textOutput("ipAddress")
    )
	)
)

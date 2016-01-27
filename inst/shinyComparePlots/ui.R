library(shiny)
library(rcellminer)

#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
source(system.file("shinyComparePlots", "functions.R", package="rcellminer"))
config <- jsonlite::fromJSON(system.file("shinyComparePlots", "config.json", package="rcellminer"))

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

shinyUI(fluidPage(
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
	),
	
    # Application title
    titlePanel(textOutput("public"), windowTitle="CellMiner"),
    
    sidebarLayout(
        sidebarPanel(
        	width=3, 
        	tags$div(id="input_container", 
                selectInput("xDataset", "x-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
                #selectInput("xPrefix", "x-Axis Type", choices=prefixChoices, selected="exp"), 
                uiOutput("xPrefixUi"),
                textInput("xId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "SLFN11"), 
    
                selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
                #selectInput("yPrefix", "y-Axis Type", choices=prefixChoices, selected="act"), 
                uiOutput("yPrefixUi"),
                textInput("yId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "94600"), 
                
                checkboxInput("showColor", "Show Color?", value=TRUE),
        				checkboxInput("selectedTissuesOnly", "Selected Tissues Only?", value=FALSE),
                
                #selectInput("showColorTissues", "Color Specific Tissues?", 
                #						choices=c("all", unique(colorSet$tissues)), multiple=TRUE, selected="all"),
                uiOutput("showColorTissuesUi"),
                
                helpText(a("CellMiner Website", href="http://discover.nci.nih.gov/", target="_blank")),
        				
            	# Generate a hidden input with TRUE or FALSE if rCharts is installed
        		tags$input(id="hasRCharts", type="text", value=hasRCharts, style="display:none"),
        					 
                #submitButton("Update"),
                
                # load Javascript snippet to parse the query string.
                singleton(tags$script(type="text/javascript", src="js/parse_input.js"))
        	)
        ),
        
        mainPanel(
        	uiOutput('tabsetPanel') 
        )
    )
))

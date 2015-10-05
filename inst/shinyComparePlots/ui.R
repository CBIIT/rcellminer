library(shiny)
library(rcellminer)

# Emulate this to check for SWATH package
if("rCharts" %in% installed.packages()) {
	options(RCHART_LIB='highcharts')	
	library(rCharts)
	hasRCharts <- TRUE
} else {
	hasRCharts <- FALSE
}

colorSet <- loadNciColorSet(returnDf=TRUE)

dataSourceChoices <- c("NCI-60"="nci60")
if (require(ccleData)){
  dataSourceChoices <- c(dataSourceChoices, "CCLE" = "ccle")
}
if (require(cgpData)){
  dataSourceChoices <- c(dataSourceChoices, "CGP" = "cgp")
}
if (require(gdscData)){
	dataSourceChoices <- c(dataSourceChoices, "GDSC" = "gdsc")
}

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

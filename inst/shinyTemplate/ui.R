library(shiny)
library(rcellminer)

colorSet <- loadNciColorSet(returnDf=TRUE)

dataSourceChoices <- setNames(names(config), sapply(names(config), function(x) { config[[x]]$name }))

# NOTE: Size is not automatically set for rChartsAlternative output
plotHeight <- 600
plotWidth <- 960

shinyUI(fluidPage(
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
                uiOutput("xPrefixUi"),
                textInput("xId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "SLFN11"), 
    
                selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
                uiOutput("yPrefixUi"),
                textInput("yId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "94600"),
                
                checkboxInput("showColor", "Show Color?", value=TRUE),
                uiOutput("showColorTissuesUi"),
                
                helpText(a("CellMiner Website", href="http://discover.nci.nih.gov/", target="_blank"))
        	)
        ),
        
        mainPanel(
        	uiOutput('tabsetPanel') 
        )
    )
))

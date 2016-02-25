library(shiny)
library(rcellminer)

colorSet <- loadNciColorSet(returnDf=TRUE)

config <- jsonlite::fromJSON(system.file("shinyTemplate", "config.json", package="rcellminer"))
dataSourceChoices <- setNames(names(config), sapply(names(config), function(x) { config[[x]]$name }))

shinyUI(
	navbarPage("CellMiner",
	  inverse=TRUE,
		tabPanel("Analysis",
			sidebarLayout(
				sidebarPanel(
					width=3, 
					tags$div(id="input_container", 
						selectInput("xDataset", "x-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
						uiOutput("xPrefixUi"),
						#textInput("xId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "SLFN11"), 
						selectizeInput('xId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5)),
						
						selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
						uiOutput("yPrefixUi"),
						#textInput("yId", "ID: (e.g. 94600 or SLFN11); Case-Sensitive", "94600"),
						selectizeInput('yId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5)),
							 						 
						checkboxInput("showColor", "Show Color?", value=TRUE),
						uiOutput("showColorTissuesUi"),
							 						 
						helpText(a("CellMiner Website", href="http://discover.nci.nih.gov/", target="_blank"))
					)
				),
				mainPanel(
					uiOutput('tabsetPanel') 
				)
			)
		),
		tabPanel("About",
			includeMarkdown("www/about.md") 
		)
	)
)

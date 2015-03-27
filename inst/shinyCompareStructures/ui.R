library(shiny)

shinyUI(fluidPage(
	  # Application title
    titlePanel(textOutput("public"), windowTitle="CellMiner"),
    
    sidebarLayout(
        sidebarPanel(width=3, 
        		selectInput("structure", "Structure Type", choices=c("NSC", "SMILES")),            
        		uiOutput("ui")
        ),
        
        mainPanel(
        	plotOutput("structurePlot", height=250, width=1250),
        	plotOutput("zScorePlot", height=800, width=1250)
        )
    	)
		)
)

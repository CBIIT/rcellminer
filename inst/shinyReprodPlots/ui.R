library(shiny)

shinyUI(fluidPage(
    # Application title
    titlePanel(textOutput("public"), windowTitle="CellMiner"),
    
    sidebarLayout(
        sidebarPanel(width=3, 
            textInput("nsc", "NSC: (e.g. 94600)", "94600"), 
            checkboxInput("zscoreRepeats", "Repeats as Z-scores?", value=FALSE)
        ),
        
        mainPanel(
            # Aligns tab content to center
            #div(align="center",
                tabsetPanel(type="tabs", 
                            tabPanel("Summary", dataTableOutput("summary")), 
                            tabPanel("Activity", plotOutput("zScorePlot", height="100%", width="100%")), 
                            tabPanel("Structure", plotOutput("structurePlot", height="100%", width="100%")), 
                            tabPanel("Repeats", htmlOutput("repeatCorResults"), plotOutput("repeatResults", height="100%", width="100%"))
                )
            #)
        )
    )
))

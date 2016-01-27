library(shiny)
library(rcellminer)

shinyUI(fluidPage(
    titlePanel("CellMiner Formulas"),
    sidebarLayout(
        sidebarPanel(
        	textInput("formula", "Formula:", "act94600~expSLFN11"),
        	includeMarkdown("www/help.md")
        ),
        mainPanel(
        	plotOutput("formulaCorPlot"),
        	verbatimTextOutput("formulaCorSummary")
        )
    )
))
library(rcellminer)
library(shiny)

shinyServer(function(input, output, session) {
    drugAct <- exprs(getAct(rcellminerData::drugData))
    drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
    molData <- getMolDataMatrices()
    
    isPublicData <- all(isPublic(rownames(drugAct)))

    if(isPublicData) {
        output$public <- renderText("CellMiner: Public")        
    } else {
        output$public <- renderText("CellMiner: Private")      
    }

    output$zScorePlot <- renderPlot({  
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }

        validate(
            need(nsc %in% rownames(drugAct), "ERROR: Invalid NSC")
        )

        plotCellMiner(drugAct, molData, "drug", nsc, NULL)
    }, height=1024, width=600, res=150)
    
    output$structurePlot <- renderPlot({        
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }
        
        validate(
            need(nsc %in% rownames(drugAct), "ERROR: Invalid NSC")
        )
        
        publicIdx <- which(drugAnnot$NSC %in% nsc)
        
        nscs <- rownames(drugAnnot[publicIdx, ])
        smiles <- drugAnnot[publicIdx, "SMILES"]
        
        plotStructuresFromNscs(nscs)
    }, height=600, width=600, res=150)
    
    output$repeatResults <- renderPlot({
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }
        
        validate(
            need(nsc %in% rownames(drugAct), "ERROR: Invalid NSC")
        )

        plotDrugActivityRepeats(nsc, useZScore=input$zscoreRepeats)
    }, height=1024, width=1024, res=150)
    
    output$repeatCorResults <- renderText({
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }
        
        validate(
            need(nsc %in% rownames(drugAct), "ERROR: Invalid NSC")
        )
        
        results <- getMinDrugActivityRepeatCor(nsc)
        
        resultsStr <- paste("Least Signicant Pairwise Correlation Between Repeats:<br/>NSC:", 
                        results$NSC, "<br/>Correlation:", results$cor, "<br/>P-value:", results$pval, collapse=" ")
    })

    output$enResults <- renderText({
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }
        
        fileName <- file.path(.db, "CGP_DRUG_SELECTION_2014", "pub_en_results", paste0("nsc_", nsc, ".html"))
        
        validate(
            need(file.exists(fileName), "ERROR: EN results do not exist for this NSC")
        )

        readChar(fileName, file.info(fileName)$size)
    })
    
    # Generate an HTML table view of the data
    output$summary <- renderDataTable({
        query <- parseQueryString(session$clientData$url_search)
        
        if("nsc" %in% names(query)) {
            nsc <- query[[1]]
        } else {
            nsc <- input$nsc
        }
        
        validate(
            need(nsc %in% rownames(drugAct), "ERROR: Invalid NSC")
        )
        
        columns <- c("NSC", "NAME", "FDA_STATUS","MOA", "SMILES")
        tmp <- drugAnnot[which(drugAnnot$NSC == nsc), columns]
        
        # Use updated MOAs
        tmp[, "mechanism"] <- getMoaStr(nsc)

        tmp 
    }, options = list(bPaginate=FALSE, bFilter=FALSE))
})

library(rcellminer)
library(sqldf)
library(shiny)
library(fingerprint)

shinyServer(function(input, output, session) {
		drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")
		
		getSimilarCompounds <- reactive({       
			if(input$structure == "NSC") {
				validate(need(input$nsc %in% names(fingerprintList), "ERROR: Invalid NSC"))
				
				drugOfInterest <- paste0("QueryNSC", input$nsc)
				smilesOfInterest <- drugAnnot[input$nsc, "SMILES"] 	
				
				tmpFingerprintList <- c(list(drugOfInterest=fingerprintList[[input$nsc]]), fingerprintList) 
			} else {
				validate(need(input$smiles != "", "ERROR: Invalid SMILES"))
				
				drugOfInterest <- "QuerySMILES"
				smilesOfInterest <- input$smiles 	      
				
				tmp <- getFingerprintList(drugOfInterest, smilesOfInterest, fpType="standard", verbose=TRUE)
				tmpFingerprintList <- c(tmp, fingerprintList) 
			}
			
			# Run fingerprint comparison 
			#results <- compareFingerprints(fingerprint.list=tmpFingerprintList)

			#DEBUG
			#tmpFingerprintList <- tmpFingerprintList[1:100]
			
			tmp <- sapply(names(tmpFingerprintList), function(x) {!is.null(tmpFingerprintList[[x]])})
			found <- which(tmp)
			
			results <- NULL 
			
			withProgress(message="Calculating Similarities...", detail="0%", 
									 min=1, max=length(tmpFingerprintList), {
			 	for(i in seq_along(tmpFingerprintList)) {
			 		setProgress(value=i, detail=paste0(floor(i*100/length(tmpFingerprintList)), "%"))
			 		tmp <- fingerprint::distance(tmpFingerprintList[[i]], tmpFingerprintList[[1]])
			 		results <- c(results, tmp)
			 	}
			 	
			 	names(results) <- names(found)
			 	results <- sort(results, decreasing=TRUE)
			})
			
			# Set up necessary data
			## Compound annotations
			df <- drugAnnot
			## Drug activities 
			drugAct <- exprs(getAct(rcellminerData::drugData))
			## Molecular profiling data
			molData <- getMolDataMatrices() 
			
			# Example filter on particular properties of the compounds
			tmpDf <- sqldf("SELECT nsc, smiles 
										 FROM df 
										 WHERE smiles != ''")
	    		
  		# Compare against the 100 NSCs for demonstration
  		#ids <- head(tmpDf$nsc, 100)
  		#smiles <- head(tmpDf$smiles, 100)
  		
  		# All public
  		ids <- tmpDf$nsc
  		smiles <- tmpDf$smiles
  		
  		# Make a vector of all the compounds to be pairwise compared 
  		ids <- c(drugOfInterest, ids)
  		smiles <- c(smilesOfInterest, smiles)
  	
      # Plot top 5 results 
      resultsIdx <- sapply(names(results)[2:6], function(x) { which(tmpDf$nsc == x) })
      resultsIds <- paste0(names(results)[2:6], "\nSim: ", round(results[2:6], 2))
      resultsSmiles <- tmpDf$smiles[resultsIdx]
      
      resultsIds <- c(drugOfInterest, names(results)[2:6])
      resultsSmiles <- c(smilesOfInterest, resultsSmiles)
			resultsCoef <- round(results[1:6], 2)
			
  		similarCompounds <- data.frame(ids=resultsIds, smiles=resultsSmiles, coef=resultsCoef, 
  																	 stringsAsFactors=FALSE)

			#DEBUG
			#cat(paste(similarCompounds$ids, collapse=", "), "\n")
			#cat(paste(similarCompounds$smiles, collapse=", "), "\n")
			#cat(paste(similarCompounds$coef, collapse=", "), "\n")
			
			return(similarCompounds)
		}) 
	
		drugAct <- exprs(getAct(rcellminerData::drugData))
		molData <- getMolDataMatrices()
		
		similarNscs <- NULL
	
    isPublicData <- all(isPublic(rownames(drugAct)))

    if(isPublicData) {
        output$public <- renderText("CellMiner: Public")        
    } else {
        output$public <- renderText("CellMiner: Private")      
    }

    output$zScorePlot <- renderPlot({
    	tmp <- getSimilarCompounds()
    	
    	# Remove query
    	tmp <- tmp[2:nrow(tmp), ]
      plotCellMiner(drugAct, molData, rep("drug", length(tmp$ids)), tmp$ids, NULL)
    }, res=150)
    
    output$structurePlot <- renderPlot({         
    	tmp <- getSimilarCompounds()
    	labels <- paste0(tmp$ids, "\n", tmp$coef)
      plotStructures(labels, tmp$smiles, titleCex=0.5, structSize=400)
    }, res=150)

    output$ui <- renderUI({
    	switch(input$structure,
    				 "NSC" = textInput("nsc", "NSC (e.g. 94600)"),
    				 "SMILES" = textInput("smiles", "SMILES (e.g. c1cncnc1)")
    	)
    })
})

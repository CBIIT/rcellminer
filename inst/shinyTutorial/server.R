library(shiny)
library(rcellminer)
library(rcellminerData)

# GET DATA ----
nci60 <- list()

# Get molecular profiles 
nci60MolDataTypes <- c("exp", "xai", "mut", "cop", "mir", "mda", "pro")
nci60$molPharmData <- getMolDataMatrices(getAllFeatureData(molData))[nci60MolDataTypes]

# Get drug activity 
nci60$molPharmData[["act"]] <- exprs(getAct(drugData))
rownames(nci60$molPharmData$act) <- paste0("act", rownames(nci60$molPharmData$act))

# SERVER ----
shinyServer(function(input, output) {
	
	# LINEAR REGRESSION ----
	getTerms <- reactive({
		# Split formula string and convert to vector
		terms <- unlist(strsplit(input$formula, '[~+]'))
		
		return(terms)
	})
	
	# Use terms to generate a dataset for lm() and run lm()
 	getFit <- reactive({
 		# Get model terms
 		terms <- getTerms() 
 		
		tmp <- NULL 
		
		for(term in terms) {
			molDataType <- getMolDataType(term)
			tmp <- rbind(tmp, nci60$molPharmData[[molDataType]][term,])
		}
		
		# NOTE: Use the terms from the input or they will not be found by lm()
		rownames(tmp) <- terms
		df <- as.data.frame(t(tmp))
		
		# Do regression 
		fit <- lm(as.formula(input$formula), data=df)
 	})

	# PLOT ----
  output$formulaCorPlot <- renderPlot({
  	# Get model terms
  	terms <- getTerms()
  	
    # Get prediction
  	fit <- getFit()
    pred <- predict(fit)
    
    # Get response variable 
    molDataType <- getMolDataType(terms[1])
    y <- nci60$molPharmData[[molDataType]][terms[1],]
    
    # Plot only matching entries (avoid samples that had missing data)
    samples <- intersect(names(pred), names(y))
    y <- y[samples]
    
    # Add regression between predicted and observed 
    fit2 <- lm(pred ~ y)
    
    # Plot predicted versus observed
    title <- paste0("Formula: ", as.character(input$formula), "; R-squared: ", round(summary(fit2)$r.squared, 3))
    plot(y, pred, xlim=range(c(y, pred)), ylim=range(c(y,pred)), xlab="observed", ylab="predicted",  main=title)
    
    # Add regression line
    abline(fit2, lwd=2)
  })
  
  # SUMMARY ----
	output$formulaCorSummary <- renderPrint({
		fit <- getFit()
    summary(fit)
  })
})

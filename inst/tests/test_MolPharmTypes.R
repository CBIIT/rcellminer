test_that("MolData and DrugData class behavior is correct.", {
  skipTest <- FALSE

  if(skipTest) {
    skip("Assuming user does not have rcellminer installed to run tests.")
  }
  
  molDataMats <- getMolDataMatrices()
  drugAct <- exprs(getAct(rcellminerData::drugData))
  drugRepeatAct <- exprs(getRepeatAct(rcellminerData::drugData))
  
	#----[Check validation code]------------------------------------------------------------
	tmp <- loadNciColorSet(returnDf = TRUE)

	# BAD: must have list of eSet objects.
	nci60ESetList <- list(cop = as.matrix(molDataMats$cop))

	# BAD: first entry in samples must match sample names associated with eSet list members.
	nci60Miame <- new("MIAME", name="CellMiner", lab="NCILMP",
										samples=list(Name=as.character(1:60), TissueType=tmp$tissues))

	expect_error(new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame),
							 "All added elements must be subclasses of eSet.")
  
  expect_error(MolData(eSetList = nci60ESetList, sampleData = nci60Miame),
  						 "All added elements must be subclasses of eSet.")

	nci60ESetList <- list(cop = ExpressionSet(as.matrix(molDataMats$cop)),
												exp = ExpressionSet(as.matrix(molDataMats$exp)))

	expect_error(new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame),
							 "All eSet samples must match samples in sampleData.")
  
  expect_error(MolData(eSetList = nci60ESetList, sampleData = nci60Miame),
  						 "All eSet samples must match samples in sampleData.")

	nci60Miame <- new("MIAME", name="CellMiner", lab="NCILMP",
										samples=list(Name=colnames(drugAct), TissueType=tmp$tissues))

	# TO DO: Address how to make construction via MolData() available upon package load.
	nci60MolData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)
  
	expect_error({ nci60MolData[["pro"]] <- as.matrix(molDataMats$pro) },
							 "All added elements must be subclasses of eSet.")
  
  nci60MolData <- MolData(eSetList = nci60ESetList, sampleData = nci60Miame)
  
  expect_error({ nci60MolData[["pro"]] <- as.matrix(molDataMats$pro) },
  						 "All added elements must be subclasses of eSet.")
	#---------------------------------------------------------------------------------------
  nci60MolData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)
	
	nci60MolData[["pro"]] <- ExpressionSet(as.matrix(molDataMats$pro))

	expect_identical(names(getAllFeatureData(nci60MolData)), c("cop", "exp", "pro"))

	expect_identical(getAllFeatureData(nci60MolData)[["cop"]], molDataMats$cop)
	expect_identical(getAllFeatureData(nci60MolData)[["exp"]], molDataMats$exp)
	expect_identical(getAllFeatureData(nci60MolData)[["pro"]], molDataMats$pro)

	expect_identical(exprs(nci60MolData[["exp"]]), molDataMats$exp)
	expect_identical(exprs(nci60MolData[["cop"]]), molDataMats$cop)
	expect_identical(exprs(nci60MolData[["pro"]]), molDataMats$pro)

	expect_identical(tmp$abbrCellLines, getSampleData(nci60MolData)$Name)
	expect_identical(tmp$tissues, getSampleData(nci60MolData)$TissueType)
  #---------------------------------------------------------------------------------------
  nci60MolData <- MolData(eSetList = nci60ESetList, sampleData = nci60Miame)
  
  nci60MolData[["pro"]] <- ExpressionSet(as.matrix(molDataMats$pro))
  
  expect_identical(names(getAllFeatureData(nci60MolData)), c("cop", "exp", "pro"))
  
  expect_identical(getAllFeatureData(nci60MolData)[["cop"]], molDataMats$cop)
  expect_identical(getAllFeatureData(nci60MolData)[["exp"]], molDataMats$exp)
  expect_identical(getAllFeatureData(nci60MolData)[["pro"]], molDataMats$pro)
  
  expect_identical(exprs(nci60MolData[["exp"]]), molDataMats$exp)
  expect_identical(exprs(nci60MolData[["cop"]]), molDataMats$cop)
  expect_identical(exprs(nci60MolData[["pro"]]), molDataMats$pro)
  
  expect_identical(tmp$abbrCellLines, getSampleData(nci60MolData)$Name)
  expect_identical(tmp$tissues, getSampleData(nci60MolData)$TissueType)
	#---------------------------------------------------------------------------------------

  unmatchedSamplesAct <- as.matrix(drugAct)
  colnames(unmatchedSamplesAct) <- 1:60

	unmatchedSamplesRepeatAct <- as.matrix(drugRepeatAct)
	colnames(unmatchedSamplesRepeatAct) <- 1:60

	expect_error(new("DrugData", act = ExpressionSet(unmatchedSamplesAct),
	                   repeatAct = ExpressionSet(unmatchedSamplesRepeatAct),
	                   sampleData = nci60Miame),
	             "eSet act samples must match samples in sampleData.")
  
  expect_error(DrugData(act = ExpressionSet(unmatchedSamplesAct),
  								 repeatAct = ExpressionSet(unmatchedSamplesRepeatAct),
  								 sampleData = nci60Miame),
  						 "eSet act samples must match samples in sampleData.")

	expect_error(new("DrugData", act = ExpressionSet(as.matrix(drugAct)),
	                 repeatAct = ExpressionSet(unmatchedSamplesRepeatAct),
	                 sampleData = nci60Miame),
	             "eSet repeatAct samples must match samples in sampleData.")
  
  expect_error(DrugData(act = ExpressionSet(as.matrix(drugAct)),
  								 repeatAct = ExpressionSet(unmatchedSamplesRepeatAct),
  								 sampleData = nci60Miame),
  						 "eSet repeatAct samples must match samples in sampleData.")

  nci60RepeatAct <- drugRepeatAct
  rownames(nci60RepeatAct) <- as.character(1:nrow(nci60RepeatAct))

	nci60DrugData <- new("DrugData", act = ExpressionSet(as.matrix(drugAct)),
	                     repeatAct = ExpressionSet(nci60RepeatAct),
	                     sampleData = nci60Miame)

  expect_identical(as.matrix(drugAct), exprs(getAct(nci60DrugData)))

	expect_identical(nci60RepeatAct, exprs(getRepeatAct(nci60DrugData)))
  
  nci60DrugData <- DrugData(act = ExpressionSet(as.matrix(drugAct)),
  										 repeatAct = ExpressionSet(nci60RepeatAct),
  										 sampleData = nci60Miame)
  
  expect_identical(as.matrix(drugAct), exprs(getAct(nci60DrugData)))
  
  expect_identical(nci60RepeatAct, exprs(getRepeatAct(nci60DrugData)))

	#---------------------------------------------------------------------------------------

})

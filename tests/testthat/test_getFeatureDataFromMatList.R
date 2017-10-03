test_that("getFeatureDataFromMatList() returns correct information", {
	featureSet <- c("expSLFN11", "mutSLX4", "copTOP1")
	molDataMats <- getMolDataMatrices()
	
	featureMat <- getFeatureDataFromMatList(featureSet, molDataMats)
	expect_identical(featureSet, rownames(featureMat))
	
	for (featureName in featureSet){
		expect_identical(molDataMats[[getMolDataType(featureName)]][featureName, , drop=FALSE],
										 featureMat[featureName, , drop=FALSE])
	}
	
	singleFeature <- featureSet[1]
	featureMat <- getFeatureDataFromMatList(singleFeature, molDataMats)
	expect_identical(molDataMats[[getMolDataType(singleFeature)]][singleFeature, , drop=FALSE],
									 featureMat)
})
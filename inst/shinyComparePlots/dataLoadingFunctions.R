getTissue2SamplesMap <- function(sampleData, typeLevelSeparator = ":"){
	stopifnot(all(!duplicated(sampleData$Name)))
	rownames(sampleData) <- sampleData$Name
	tissueToSamples <- list()
	ocLevels <- paste0("OncoTree", 1:4)
	
	for (sample in rownames(sampleData)){
		sampleOcTypes <- as.character(sampleData[sample, ocLevels])
		typeName <- sampleOcTypes[1]
		if (is.na(typeName)){
			next
		}
		
		tissueToSamples[[typeName]] <- c(tissueToSamples[[typeName]], sample)
		for (i in (2:4)){
			if (is.na(sampleOcTypes[i])){
				break
			}
			typeName <- paste0(typeName, typeLevelSeparator, sampleOcTypes[i])
			tissueToSamples[[typeName]] <- c(tissueToSamples[[typeName]], sample)
		}
	}
	
	# Additional phenotype-based sample sets -------------------------------------
	if ("EMT" %in% colnames(sampleData)) {
		emtSampleData <- sampleData[!is.na(sampleData$EMT), , drop = FALSE]
		tmp <- split(emtSampleData$Name, emtSampleData$EMT)
		stopifnot(length(intersect(names(tmp), names(tissueToSamples))) == 0)
		tissueToSamples <- c(tissueToSamples, tmp)
	}
	# ----------------------------------------------------------------------------
	
	return(tissueToSamples)
}



loadSourceContent <- function(srcConfig){
	# Load necessary packages; return NULL otherwise. ----------------------------
	packageNames <- names(srcConfig$packages)
	pkgsInstalled <- setNames(logical(length(packageNames)), packageNames)
	for (pkgName in packageNames){
		if (!require(pkgName, character.only = TRUE)){
			warning(paste0("Missing package required for ", srcConfig$displayName, 
										 ": ", pkgName, "."))
		} else{
			pkgsInstalled[pkgName] <- TRUE
		}
	}
	if (!all(pkgsInstalled)){
		return(NULL)
	}
	#-----------------------------------------------------------------------------
	
	# Load molecular profiling data and drug activity data -----------------------
	# The data associated with a source may be distributed over multiple packages,
	# each containing a molData object and/or a drugData object.
	molDataMats <- list()
	drugFeaturePrefix <- NULL
	drugAct <- NULL
	drugAnnot <- NULL
	sampleData <- NULL
	
	# This will be a two column table [featurePrefix, displayName], aggregating
	# all feature information for the source.
	featureInfoTab <- NULL
	
	for (pkgName in packageNames){
		pkgEnv <- new.env()
		
		if ("MolData" %in% names(srcConfig$packages[[pkgName]])){
			data("molData", package=pkgName, envir=pkgEnv)
			pkgMolDataInfo <- srcConfig$packages[[pkgName]][["MolData"]]
			rownames(pkgMolDataInfo) <- pkgMolDataInfo$eSetListName
			
			if (!("dataMin" %in% colnames(pkgMolDataInfo))){
				pkgMolDataInfo$dataMin <- NA
			}
			if (!("dataMax" %in% colnames(pkgMolDataInfo))){
				pkgMolDataInfo$dataMax <- NA
			}
			
			featureInfoTab <- rbind(featureInfoTab, 
				pkgMolDataInfo[, c("featurePrefix", "displayName", "dataMin", "dataMax")])
			
			pkgMolDataMats <- getAllFeatureData(pkgEnv$molData)
			if (is.null(sampleData)){ 
				sampleData <- getSampleData(pkgEnv$molData) 
			}
			
			for (dataType in rownames(pkgMolDataInfo)){
				# validity checks ---------------------------------------------------------
				if (!(dataType %in% names(pkgMolDataMats))){
					stop("Check config file: indicated eSetListName '", dataType, 
							 "' is not found in ", pkgName, " package molData object.")
				}
				
				dat <- pkgMolDataMats[[dataType]]
				# restrict features ---------------------------------------------------------------------------
				if (!is.null(pkgMolDataInfo[dataType, "includedFeaturesFile"])){
					if (!is.na(pkgMolDataInfo[dataType, "includedFeaturesFile"])){
						inputFile <- pkgMolDataInfo[dataType, "includedFeaturesFile"]
						includedFeatures <- read.table(file=inputFile, header=FALSE, sep = "\t",
																					 stringsAsFactors=FALSE)[,1]
						if (!all(includedFeatures %in% rownames(dat))){
							stop("Check ", pkgName, " MolData included features file for data type ",
									 dataType, " because some specified features are not in data matrix.")
						}
						dat <- dat[includedFeatures, ]
					}
				}
				# ---------------------------------------------------------------------------------------------
				molFeaturePrefix <- pkgMolDataInfo[dataType, "featurePrefix"]
				rownames(dat) <- paste0(molFeaturePrefix, rownames(dat))
				
				if (!identical(colnames(dat), sampleData$Name)){
					stop("Sample names for ", pkgName, " molData data type '",
							 dataType, "' are inconsistent with loaded data source sample names.")
				}
				
				if (!is.null(molDataMats[[molFeaturePrefix]])){
					if (any(rownames(dat) %in% rownames(molDataMats[[molFeaturePrefix]]))){
						stop("Molecular features for ", pkgName, " molData type '", dataType,
								 "' duplicate loaded features of this type.")
					}
				}
				# -------------------------------------------------------------------------
				
				molDataMats[[molFeaturePrefix]] <- rbind(molDataMats[[molFeaturePrefix]], dat)
			}
		}
		
		if ("DrugData" %in% names(srcConfig$packages[[pkgName]])){
			data("drugData", package=pkgName, envir=pkgEnv)
			pkgDrugDataInfo <- srcConfig$packages[[pkgName]][["DrugData"]]
			
			if (!("dataMin" %in% colnames(pkgDrugDataInfo))){
				pkgDrugDataInfo$dataMin <- NA
			}
			if (!("dataMax" %in% colnames(pkgDrugDataInfo))){
				pkgDrugDataInfo$dataMax <- NA
			}
			
			featureInfoTab <- rbind(featureInfoTab, 
				pkgDrugDataInfo[, c("featurePrefix", "displayName", "dataMin", "dataMax")])
			
			dat <- exprs(getAct(pkgEnv$drugData))
			if (is.null(sampleData)){ 
				sampleData <- getSampleData(pkgEnv$drugData) 
			}
			if (is.null(drugFeaturePrefix)) { 
				drugFeaturePrefix <- pkgDrugDataInfo$featurePrefix 
			}
			
			annot <- getFeatureAnnot(pkgEnv$drugData)[["drug"]]
			
			# Packages should provide at least minimal drug annotation data,
			# but if this is not available, the code below constructs it 
			# from the activity data matrix feature (row) names.
			if (ncol(annot) == 0){
				if (is.na(pkgDrugDataInfo[, "drugAnnotIdCol"])){
					pkgDrugDataInfo[, "drugAnnotIdCol"] <- "ID"
					annot[["ID"]] <- rownames(dat)
				}
				if (is.na(pkgDrugDataInfo[, "drugAnnotNameCol"])){
					pkgDrugDataInfo[, "drugAnnotNameCol"] <- "NAME"
					annot[["NAME"]] <- rownames(dat)
				}
				if (is.na(pkgDrugDataInfo[, "drugAnnotMoaCol"])){
					pkgDrugDataInfo[, "drugAnnotMoaCol"] <- "MOA"
					annot[["MOA"]] <- ""
				}
			} 

			# validity checks ---------------------------------------------------------
			expectedAnnotCols <- as.character(pkgDrugDataInfo[, c("drugAnnotIdCol", 
																														"drugAnnotNameCol", 
																														"drugAnnotMoaCol")])
			
			if (any(is.na(expectedAnnotCols))){
				stop("Check config file: drug annotation table columns ",  
						 "for ", pkgName, " are innappropriately set to NA values.")
			}
			
			if (!all(expectedAnnotCols %in% colnames(annot))){
				stop("Check config file: indicated drug annotation table columns ",  
						 "are not found in ", pkgName, " package drugData object.")
			} else{
				annot <- annot[expectedAnnotCols]
				colnames(annot) <- c("ID", "NAME", "MOA")
				# Handle uncommon characters
				annot$ID <- iconv(enc2utf8(as.character(annot$ID)), sub="byte")
				annot$NAME <- iconv(enc2utf8(annot$NAME), sub="byte")
				annot$MOA <- iconv(enc2utf8(annot$MOA), sub="byte")
				rownames(annot) <- annot$ID
			}
			
			if (!identical(rownames(annot), rownames(dat))){
				stop("Check drug data for ", pkgName, 
						 ": annotation data table and activity data matrix row orders do not match.")
			} else{
				if (!identical(drugFeaturePrefix, pkgDrugDataInfo$featurePrefix)){
					stop("Check config file: use same drug activity data featurePrefix with all ",
							 "packages providing drug data for ", srcConfig$displayName, ".")
				}
				# restrict features ---------------------------------------------------------------------------
				if (!is.null(pkgDrugDataInfo[1, "includedFeaturesFile"])){
					inputFile <- pkgDrugDataInfo[1, "includedFeaturesFile"]
					includedFeatures <- read.table(file=inputFile, header=FALSE, sep = "\t",
																				 stringsAsFactors=FALSE)[,1]
					includedFeatures <- as.character(includedFeatures)
					if (!all(includedFeatures %in% rownames(dat))){
						stop("Check ", pkgName, " DrugData included features file for data type ",
								 drugFeaturePrefix, " because some specified features are not in data matrix.")
					}
					dat <- dat[includedFeatures, ]
					annot <- annot[rownames(dat), ]
				}
				# ---------------------------------------------------------------------------------------------
				rownames(dat) <- paste0(drugFeaturePrefix, rownames(dat))
				rownames(annot) <- rownames(dat)
			}
			
			if (!identical(colnames(dat), sampleData$Name)){
				stop("Sample names for ", pkgName, 
						 " drug activity data are inconsistent with loaded data source sample names.")
			}
			
			if (!is.null(drugAct)){
				if (any(rownames(dat) %in% rownames(drugAct))){
					stop("Drug names for ", pkgName, 
							 " drugData clash with already loaded drug names." )
				}
			}
			# -------------------------------------------------------------------------
			
			drugAct <- rbind(drugAct, dat)
			drugAnnot <- rbind(drugAnnot, annot)
		}
	}
	#-----------------------------------------------------------------------------
	
	# Assemble app data object to be returned. ---------------------------------
	src <- list()
	src$molPharmData <- molDataMats
	if (drugFeaturePrefix %in% names(src$molPharmData)){
		stop("Check config file: ", srcConfig$displayName, 
				 " drug featurePrefix must be different from all molecular data feature prefixes.")
	}
	src$molPharmData[[drugFeaturePrefix]] <- drugAct
	src$drugInfo <- drugAnnot
	
	src$sampleData <- sampleData
	rownames(src$sampleData) <- src$sampleData$Name
	
	# TO DO: Check whether spaces in tissue sample names creates any problems.
	src$tissueToSamplesMap <- getTissue2SamplesMap(src$sampleData)
	
	# TO DO: READ COLOR MAP INFORMATION FROM CONFIG FILE?
	src$tissueColorMap <- rep("rgba(0,0,255,0.5)", length(src$tissueToSamplesMap))
	names(src$tissueColorMap) <- names(src$tissueToSamplesMap)
	# ---------------------------------------------------------------------------
	
	if (any(duplicated(featureInfoTab$featurePrefix))){
		stop("Check configuration file: unexpected duplicate feature prefix entries.")
	} else{
		rownames(featureInfoTab) <- featureInfoTab$featurePrefix
		
		if (!is.null(drugFeaturePrefix)){
			rownamesDrugFirst <- c(drugFeaturePrefix, 
														 setdiff(rownames(featureInfoTab), drugFeaturePrefix))
			featureInfoTab <- featureInfoTab[rownamesDrugFirst, ]
		}
	}
	
	src$drugFeaturePrefix <- drugFeaturePrefix
	src$featurePrefixes <- setNames(featureInfoTab$featurePrefix, featureInfoTab$displayName)
	
	if (is.null(drugFeaturePrefix)){
		src$defaultFeatureX <- src$featurePrefixes[1]
		src$defaultFeatureY <- src$featurePrefixes[1]
	} else{
		src$defaultFeatureX <- src$featurePrefixes[2]
		src$defaultFeatureY <- src$featurePrefixes[1]
	}
	
	# Check/set feature range information ----------------------------------------
	stopifnot(all(rownames(featureInfoTab) %in% names(src$molPharmData)))
	stopifnot(c("dataMin", "dataMax") %in% colnames(featureInfoTab))
	
	src$featureValRanges <- list()
	for (dataType in rownames(featureInfoTab)){
		datRange <- range(as.numeric(src$molPharmData[[dataType]]), na.rm = TRUE)
		datRange[1] <- floor(datRange[1])
		datRange[2] <- ceiling(datRange[2])
		if (!is.na(featureInfoTab[dataType, "dataMin"])){
			datRange[1] <- min(datRange[1], featureInfoTab[dataType, "dataMin"])
		}
		if (!is.na(featureInfoTab[dataType, "dataMax"])){
			datRange[2] <- max(datRange[2], featureInfoTab[dataType, "dataMax"])
		}
		src$featureValRanges[[dataType]] <- datRange
	}
	#-----------------------------------------------------------------------------
	
	return(src)
}
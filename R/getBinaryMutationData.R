#' Compute a binary gene mutation data matrix from SNP and other mutation
#' event-level data.
#' 
#' @param mutInfo A data frame with the following named columns:
#'  Gene, the name of the gene associated with the mutation event;
#'  probe.ids, a unique identifier specifying the mutation event;
#'  SNP_1000_genome, the frequency of the mutation event in SNP 1000;
#'  ESP5400, the frequency of the mutation event in ESP5400;
#'  SNP_type, the type of mutation event, chosen from "MISSENSE", "FRAMESHIFT", 
#'  "NONFRAMESHIFT", "NONSENSE", "SPLICING";
#'  SIFT_score, the SIFT score;
#'  Polyphen_score, the POLYPHEN score.
#'  Rownames of mutInfo should be set to probe.ids, i.e., the unique mutation
#'  event specifier.
#' @param mutData A matrix with event level mutation information, with SNPs, etc.
#' along rows and samples along columns.  Rownames of mutData should exactly
#' match those of mutInfo.  The i-th row of mutInfo should thus give detailed
#' information for the mutation event with data specified in the i-th row of
#' mutData.
#' @param maxVariantFreq The maximum proportion of mutant samples (used to
#' exclude frequently occuring events); default value = 0.2.
#' @param maxNormalPopulationFreq  The maximum freqency of a mutation in the normal
#' population (used to exclude likely germline variants); default value = 0.005.
#' @param maxSiftScore The maximum accepted SIFT score (used to exclude 
#' presumed non-deleterious mutations); default value = 0.05.
#' @param minPolyPhenScore  The minimum accepted POLYPHEN score (used to
#' exclude presumed non-deleterious mutations); default value = 0.85.
#' @return A binary gene mutation matrix, with genes along rows, samples
#' along columns, and 1s indicating deleterious mutations.
#' 
#' @concept rcellminer
#' @export
#' 
#' @importFrom stringr str_match str_split str_trim
getBinaryMutationData <- function(mutInfo, mutData, maxVariantFreq = 0.2,
                                  maxNormalPopulationFreq = 0.005, maxSiftScore = 0.05,
                                  minPolyPhenScore = 0.85){
  if (!identical(rownames(mutInfo), rownames(mutData))){
    stop("rownames(mutInfo) must be identical to rownames(mutData).")
  }
  
  #----[filter: remove SNPs between two genes]-------------------------------------------
  # SNPs with gene names containing parentheses are between genes, with detailed
  # information between parentheses.
  hasParen <- !is.na(str_match(mutInfo$Gene, pattern = "[(]")[, 1])
  mutInfo <- mutInfo[!hasParen, ]
  mutData <- mutData[!hasParen, ]
  #--------------------------------------------------------------------------------------
  
  #----[filter: remove rows for SNPs mapped to >1 gene; add replicate rows by gene]------
  isMultiGene <- !is.na(str_match(mutInfo$Gene, pattern = "[,;]")[, 1])
  mutInfoMultiGene <- mutInfo[isMultiGene, ]
  mutDataMultiGene <- mutData[isMultiGene, ]
  
  mutInfo <- mutInfo[!isMultiGene, ]
  mutData <- mutData[!isMultiGene, ]
  
  # Initialize added mutation info matrix.
  numAddedRows <- length(c(str_split(mutInfoMultiGene$Gene, pattern = "[,;]"), recursive = TRUE))
  addedRowNames <- paste0("added_", 1:numAddedRows)
  mutInfoAdded <- data.frame(row.names = addedRowNames, stringsAsFactors = FALSE)
  for (colName in colnames(mutInfo)){
    mutInfoAdded[, colName] <- vector(class(mutInfo[, colName]), numAddedRows)
  }
  
  # Initialize added mutation data matrix.
  mutDataAdded <- matrix(NA_real_, nrow=numAddedRows, ncol=ncol(mutData))
  rownames(mutDataAdded) <- addedRowNames
  colnames(mutDataAdded) <- colnames(mutData)
  
  i <- 0
  geneColIndex <- which(colnames(mutInfoMultiGene) == "Gene")
  for (probeName in rownames(mutInfoMultiGene)){
    geneSet <- str_trim(str_split(mutInfoMultiGene[probeName, "Gene"], pattern = "[,;]")[[1]])
    
    for (gene in geneSet){
      i <- i + 1
      mutInfoAdded[i,  geneColIndex]  <- gene
      mutInfoAdded[i, -geneColIndex]  <- mutInfoMultiGene[probeName, -geneColIndex]
      mutDataAdded[i, ] <- mutDataMultiGene[probeName, ]
    }
  }
  if (!(i == numAddedRows)){
    stop("Unexpected state (with duplication of data rows for probes mapped to multiple genes).")
  }
  
  mutInfo <- rbind(mutInfo, mutInfoAdded)
  mutData <- rbind(mutData, mutDataAdded)
  mutData[which(is.na(mutData))] <- 0
  #--------------------------------------------------------------------------------------
  
  #----[filter: remove frequent variants]------------------------------------------------
  numMutatedLines <- rowSums(mutData > 0)
  isInfreqMutation <- numMutatedLines < (ncol(mutData) * maxVariantFreq)
  
  mutInfo <- mutInfo[isInfreqMutation, ]
  mutData <- mutData[isInfreqMutation, ]
  #--------------------------------------------------------------------------------------
  
  #----[filter: remove non-deleterious variants]-----------------------------------------
  isDeleteriousMisense <- (mutInfo$SNP_type == "MISSENSE") & 
                          ((mutInfo$SIFT_score < maxSiftScore) | (mutInfo$Polyphen_score >= minPolyPhenScore))
  isDeleterious <- (mutInfo$SNP_type %in% c("FRAMESHIFT", "NONFRAMESHIFT", "NONSENSE", "SPLICING")) | 
                    isDeleteriousMisense
  
  mutInfo <- mutInfo[isDeleterious, ]
  mutData <- mutData[isDeleterious, ]
  #--------------------------------------------------------------------------------------
  
  #----[filter: remove (putative germline) variants common in the normal population]-----
  isSomatic <- (mutInfo$SNP_1000_genome <= maxNormalPopulationFreq) |
    (mutInfo$ESP5400 <= maxNormalPopulationFreq)
  
  mutInfo <- mutInfo[isSomatic, ]
  mutData <- mutData[isSomatic, ]
  #--------------------------------------------------------------------------------------
  
  #----[compute binary, gene-specific mutation profiles]---------------------------------
  mutData <- mutData/100
  genes <- sort(unique(mutInfo$Gene))
  binMutData <- matrix(0, nrow = length(genes), ncol = ncol(mutData))
  rownames(binMutData) <- genes
  colnames(binMutData) <- colnames(mutData)
  
  for (i in (1:length(genes))){
    qw <- which(mutInfo$Gene==genes[i])
    
    if (length(qw) > 1){
      x <- mutData[qw, ]
      binMutData[i, ] <- 1 - apply(1 - x, 2, prod, na.rm=TRUE)
    } else if (length(qw) == 1){
      binMutData[i, ] <- mutData[qw, ]
    }
  }
  
  qw <- which(binMutData > 0)
  binMutData[qw] <- 1
  #--------------------------------------------------------------------------------------
  
  return(binMutData)
}

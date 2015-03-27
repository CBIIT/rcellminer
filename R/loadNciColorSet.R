#' Returns a 60-element color set that matches the color set used on
#' http://discover.nci.nih.gov/
#' 
#' @param returnDf a boolean if a data.frame with tissue names and abbreviations
#'   should be returned (default: FALSE)
#' @return a vector of colors as strings or a data.frame with tissues, tissue
#'   abbreviations, cell line abbreviations and colors
#'   
#' @examples
#' loadNciColorSet()  
#'   
#' @concept rcellminer
#' @export
loadNciColorSet <- function(returnDf=FALSE) {
  tissues <- c(rep("renal", 8), rep("prostrate", 2), rep("ovarian", 7), rep("lung", 9), rep("melanoma", 10), rep("leukemia", 6), rep("colon", 7), rep("central nervous system", 6), rep("breast", 5))
  
  abbrTissues <- c(rep("RE", 8), rep("PR", 2), rep("OV", 7), rep("LC", 9), rep("ME", 10), rep("LE", 6), rep("CO", 7), rep("CNS", 6), rep("BR", 5))
  
#   abbrCellLines <- c("RE:CAKI1", "RE:RXF393", "RE:SN12C", "RE:TK10", "RE:UO31", 
#                      "PR:PC3", "PR:DU145", "RE:7860", "RE:A498", "RE:ACHN", 
#                      "OV:OVCAR4", "OV:OVCAR5", "OV:OVCAR8", "OV:SKOV_3", "OV:NCIADR_RES", 
#                      "LC:NCIH322M", "LC:NCIH460", "LC:NCIH522", "OV:IGROV1", "OV:OVCAR3", 
#                      "LC:EKVX", "LC:HOP62", "LC:HOP92", "LC:NCIH226", "LC:NCIH23", 
#                      "ME:UACC257", "ME:UACC62", "ME:MDAMB_435", "ME:MDA_N", "LC:A549", 
#                      "ME:MALME3M", "ME:M14", "ME:SKMEL_2", "ME:SKMEL_28", "ME:SKMEL_5", 
#                      "LE:K562", "LE:MOLT4", "LE:RPMI8226", "LE:SR", "ME:LOXIMVI", 
#                      "CO:HT29", "CO:KM12", "CO:SW620", "LE:CCRFCEM", "LE:HL60", 
#                      "CNS:U251", "CO:COLO205", "CO:HCC2998", "CO:HCT116", "CO:HCT15", 
#                      "CNS:SF268", "CNS:SF295", "CNS:SF539", "CNS:SNB19", "CNS:SNB75",
#                      "BR:MCF7", "BR:MDAMB_231", "BR:HS578T", "BR:BT549", "BR:T47D")    

	abbrCellLines <- c("RE:UO_31", "RE:TK_10", "RE:SN12C", "RE:RXF_393", "RE:CAKI_1", 
										 "RE:ACHN", "RE:A498", "RE:786_0", "PR:DU_145", "PR:PC_3", 
										 "OV:NCI_ADR_RES", "OV:SK_OV_3", "OV:OVCAR_8", "OV:OVCAR_5", 
										 "OV:OVCAR_4", "OV:OVCAR_3", "OV:IGROV1", "LC:NCI_H522", 
										 "LC:NCI_H460", "LC:NCI_H322M", "LC:NCI_H23", "LC:NCI_H226", 
										 "LC:HOP_92", "LC:HOP_62", "LC:EKVX", "LC:A549", "ME:MDA_N", 
										 "ME:MDA_MB_435", "ME:UACC_62", "ME:UACC_257", "ME:SK_MEL_5",
										 "ME:SK_MEL_28", "ME:SK_MEL_2", "ME:M14", "ME:MALME_3M", 
										 "ME:LOXIMVI", "LE:SR", "LE:RPMI_8226", "LE:MOLT_4", "LE:K_562",  
										 "LE:HL_60", "LE:CCRF_CEM", "CO:SW_620", "CO:KM12", "CO:HT29", 
										 "CO:HCT_15", "CO:HCT_116", "CO:HCC_2998", "CO:COLO205", 
										 "CNS:U251", "CNS:SNB_75", "CNS:SNB_19", "CNS:SF_539", 
										 "CNS:SF_295", "CNS:SF_268", "BR:T47D", "BR:BT_549", "BR:HS578T", 
										 "BR:MDA_MB_231", "BR:MCF7") 
  
  #NCI60 color set
	colors <- c(rep("red", 8), rep("yellow", 2), rep("purple", 7), rep("steelblue", 9), rep("darkolivegreen", 10), rep("yellowgreen", 6), rep("orange", 7), rep("sienna", 6), rep("darkblue", 5))
  #colors <- c(rep("red", 8), rep("yellow", 2), rep("maroon", 7), rep("steelblue", 9), rep("darkolivegreen", 10), rep("yellowgreen", 6), rep("orange", 7), rep("sienna", 6), rep("darkblue", 5))
  #colors <- c(rep("red", 8), rep("yellow", 2), rep("brown4", 7), rep("steelblue3", 9), rep("chartreuse4", 9), rep("olivedrab3", 7), rep("orange1", 7), rep("orange4", 6), rep("dodgerblue4", 5))
	
  df <- data.frame(tissues=tissues, abbrTissues=abbrTissues, abbrCellLines=abbrCellLines, colors=colors, stringsAsFactors=FALSE)
	df <- df[rev(rownames(df)),]

  if(returnDf) {
    return(df)
  } else {
  	colors <- rev(colors)
    return(colors)     
  }
}

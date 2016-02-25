#' CellMiner Drug Response Values 
#' 
#' @details A list containing response values and annotations: 
#' \itemize{
#'  \item{act} {Z-scores of the averaged negative log GI (growth inhibition) 50 values across repeats
#'     for the NCI-60; assay described here: \url{http://dtp.nci.nih.gov/branches/btb/ivclsp.html }}
#'  \item{annot} { 
#'    \itemize{
#'      \item{id} {Dataset identifier; NOTE: DO NOT use this column; the NSC is 
#'        the primary drug identifier}
#'      \item{nsc} {National Service Center identifier; the primary drug identifier}
#'      \item{name} {Compound name}
#'      \item{brand_name} {Brand name for the compound, if sold commericially}
#'      \item{formula} {Compound chemical formula}
#'      \item{testing_status} {Information on whether it is known if the 
#'        compound is FDA approved or undergoing testing in clinical trials}
#'      \item{source} {TODO} 
#'      \item{smiles} {Compound chemical structure as a SMILES string} 
#'      \item{weight} {Compound chemical weight in g/mol} 
#'      \item{mechanism} {Pharmacological mechanism of action} 
#'      \item{confidential_flag} {A flag to indicate if compound information is public} 
#'      \item{total_probes} {TODO} 
#'      \item{total_good_probes} {TODO} 
#'      \item{low_correlations} {TODO} 
#'      \item{failure_reason} {TODO} 
#'      \item{cas} {CAS Registry Number; NOTE: Due to data restrictions PubChem IDs 
#'        are the preferred mapping ID to other datasets} 
#'      \item{pubchem_id} {PubChem ID}
#'  	 }
#'   }
#' }
#'
#' @name drugDB
#' @docType data
#' @author Vinodh Rajapakse \email{vinodh.rajapakse AT nih.gov}
#' @references \url{http://discover.nci.nih.gov/cellminer/loadDownload.do}
#' @keywords data
#' @concept rcellminer 
NULL

#' CellMiner Version
#' 
#' @details The version of CellMiner used 
#'
#' @name cmVersion
#' @docType data
#' @author Vinodh Rajapakse \email{vinodh.rajapakse AT nih.gov}
#' @references \url{http://discover.nci.nih.gov/cellminer}
#' @keywords data
#' @concept rcellminer 
NULL

#' NCI60 Molecular Data
#' 
#' Z-scores of values for a variety of assays conducted on the NCI-60 to facilitate comparison.
#' Z-scores calculated over the 60 cell lines for the given feature. 
#' 
#' @details A list containing various assay values:  
#' \itemize{
#'  \item{cop} { Copy number values; Described in Pubmed ID: 24670534 }
#'  \item{exp} { Expression values; Obtained from "RNA: 5 Platform Gene Transcript" 
#'    \url{http://discover.nci.nih.gov/cellminer/loadDownload.do}; Missing values 
#'    imputed using the R package "impute" }
#'  \item{mut} { Mutation values; Deleterious mutations obtained from TODO }
#'  \item{mir} { MicroRNA values; Obtained from "RNA: Agilent Human microRNA (V2)"
#'    \url{http://discover.nci.nih.gov/cellminer/loadDownload.do} }
#'  \item{pro} { Reverse protein lysate array values; Obtain from "Protein: Lysate Array"
#'    \url{http://discover.nci.nih.gov/cellminer/loadDownload.do} }
#'  \item{mda} { NCI-60 metadata. 
#'    \itemize{ 
#'      \item{CNV_GAIN} {Proportion of genome copy number gains; Described in Pubmed ID: 24670534} 
#'      \item{CNV_LOSS} {Proportion of genome copy number losses; Described in Pubmed ID: 24670534} 
#'      \item{CNV_TOTAL} {Sum of CNV_GAIN and CNV_LOSS} 
#'      \item{P53_BIN} {Binary TP53 profile curated by William Reinhold } 
#'      \item{MSI_OGAN_BIN} {Binary microsatellite instability (MSI) profile curated by Ogan Abaan using COSMIC data; Obtained from Supplementary Table 1 - Ogan Whole Exome Sequencing (WES) paper in Cancer Res. } 
#'      \item{EPITHELIAL} {Epithelial by tissue of origin - pattern extracted from the CellMiner cell line metadata \url{http://discover.nci.nih.gov/cellminer/celllineMetadata.do} } 
#'      \item{EPITHELIAL_KURT} {Kurt Kohn curation for epithelial-like cell lines based on molecular parameters described in Pubmed ID: 24940735} 
#'      \item{DELETERIOUS} {Total deleterious variants from WES dataset; Fabricio Sousa curation} 
#'      \item{MISSENSE} {Total missense variants from WES dataset; Fabricio Sousa curation} 
#'      \item{SILENT} {Total silent variants from WES dataset; Fabricio Sousa curation} 
#'      \item{TOTAL_AA} {Total amino acid changing variants from WES dataset; Fabricio Sousa curation} 
#'      \item{CELL-CELL} {Cell-to-cell adhesion curated by William Reinhold } 
#'      \item{DOUBLINGTIME} {The doubling time pattern was extracted from the CellMiner cell line metadata 
#'        \url{http://discover.nci.nih.gov/cellminer/celllineMetadata.do} } 
#'    }
#'  }
#' }
#' 
#' @name elNetMolDataNCI60
#' @docType data
#' @author Vinodh Rajapakse \email{vinodh.rajapakse AT nih.gov}
#' @references \url{http://discover.nci.nih.gov/cellminer}
#' @keywords data
#' @concept rcellminer 
NULL



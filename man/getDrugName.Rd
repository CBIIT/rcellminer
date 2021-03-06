% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDrugName.R
\name{getDrugName}
\alias{getDrugName}
\title{Get the drug names for a set of NSC identifiers.}
\usage{
getDrugName(nscSet)
}
\arguments{
\item{nscSet}{A character vector of NSC strings}
}
\value{
A named character vector indicating the compound names
for each NSC in nscSet (with an empty string returned if no
such information is available, and an NA returned if the
NSC is not included in the CellMiner database).
}
\description{
Get the drug names for a set of NSC identifiers.
}
\examples{
nscSet <- c("609699", "94600")
getDrugName(nscSet)

}
\concept{rcellminer}

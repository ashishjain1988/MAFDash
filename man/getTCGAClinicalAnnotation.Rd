% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMAFDataTCGA.R
\name{getTCGAClinicalAnnotation}
\alias{getTCGAClinicalAnnotation}
\title{Function to extract the clinical annotations from TCGA}
\usage{
getTCGAClinicalAnnotation(
  cancerCodes = c("ACC"),
  outputFolder = tempdir(),
  plotdata = NULL
)
}
\arguments{
\item{cancerCodes}{A character vector of TCGA cancer codes}

\item{outputFolder}{path of the file containing the clinical
annotations from TCGA}

\item{plotdata}{Flag to plot the annotations}
}
\value{
A list containing the TCGA clinical annotations
}
\description{
This function download and extract the clinical
annotations from TCGA.
}
\examples{
library(MAFDash)
cancerCodes <- "ACC"
#clinicalData <- getTCGAClinicalAnnotation(cancerCodes = cancerCodes)

}
\author{
Mayank Tandon, Ashish Jain
}

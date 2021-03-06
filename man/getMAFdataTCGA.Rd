% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMAFDataTCGA.R
\name{getMAFdataTCGA}
\alias{getMAFdataTCGA}
\title{Function to extract the mutation data in MAF format from TCGA}
\usage{
getMAFdataTCGA(
  cancerCode = "ACC",
  outputFolder = file.path("data"),
  variant_caller = "mutect2"
)
}
\arguments{
\item{cancerCode}{The TCGA cancer code}

\item{outputFolder}{The path of the file containing the mutation
information in the MAF format}

\item{variant_caller}{The type of variant caller in TCGA}
}
\value{
A string containing the path of the downloaded mutation annotation file
}
\description{
This function download and extract the mutation
 data in MAF format from TCGA.
}
\examples{
library(MAFDash)
cancerCode <- "ACC"
outputFolderPath <- tempdir()
\donttest{maf <- getMAFdataTCGA(cancerCode = cancerCode,outputFolder = outputFolderPath)}
}
\author{
Mayank Tandon, Ashish Jain
}

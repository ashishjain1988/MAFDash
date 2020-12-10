#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param filePath The path of the file containing the mutation
#' information in the MAF format
#' @param outputFilePath The path of the output html dashboard
#' file
#' @param outputFileName The name of the output html dashboard
#' file
#' @param outputFileTitle The title of the output html dashboard
#' file
#' @export
#' @return The dashboard html file
#'
#' @examples
#' library(MAFDashRPackage)
#' #MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
getMAFDashboard<-function(filePath,outputFilePath,outputFileName="output",outputFileTitle="MAF Dash"){

  MAFRmdfile <- system.file('extdata', 'MAFDash.Rmd', package = 'MAFDashRPackage')
  html_filename=paste0(outputFileName,".",gsub(".Rmd",".html",basename(MAFRmdfile)))
  rmarkdown::render(MAFRmdfile,
                    knit_root_dir=getwd(),
                    output_format="all", output_file=html_filename,
                    params = list(
                      maffile=filePath,
                      titletext=outputFileTitle
                    ))
  ### rmarkdown::render doesn't let you select output destination (it uses the path of the Rmd file)
  ##  So this bit will move the report to the path in the 'out_dir' variable
  if (!dir.exists(outputFilePath)) { outputFilePath = getwd()}#dir.create(outputFilePath, recursive = T) }
  file.rename(file.path(dirname(MAFRmdfile),html_filename), file.path(outputFilePath,html_filename))
}

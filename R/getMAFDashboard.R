## To Supress Note
utils::globalVariables(c(".","tempdir"))

#' Function to generate a dashboard from a MAF file.
#' @description This function creates an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tandon, Ashish Jain
#' @param maf MAF object or MAF file path
#' information in the MAF format
#' @param plotList A named list denoting which plots to draw.
#' List elements can be:
#' - boolean, if the name matches one of built-in plots: "summary_plot","burden","oncoplot","cooccurence","heatmap"
#' - ggplot, ComplexHeatmap, or plotly object
#' - file path to image (must be absolute path)
#' The default value (if set to NULL) is 'list("summary_plot"=TRUE,"burden"=TRUE,"oncoplot"=TRUE,"cooccurence"=TRUE,"heatmap"=TRUE)'
#' The order of the list determines the order of the tabs, and list names are used as tab titles.
#' @param outputFileName The name of the output html dashboard
#' file
#' @param outputFileTitle The title of the output html dashboard
#' @param outputFilePath The path of the output html dashboard
#' @export
#' @return No return value, the MAF dashboard html file is created in the given output folder
#'
#' @examples
#' library(MAFDash)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' \donttest{getMAFDashboard(maf = maf,outputFilePath=tempdir())}
#' @importFrom rmarkdown render
#' @importFrom knitr knit
#' @importFrom plotly plot_ly ggplotly plotly layout
#' @import maftools
#' @import htmltools
#' @import bsplus
#' @import crosstalk
#' @import circlize
#' @import canvasXpress
#' @importFrom DT datatable formatStyle JS
#' @import flexdashboard

getMAFDashboard<-function(maf=NULL,plotList=NULL,outputFileName="dashboard.html",outputFileTitle="MAF Dash",outputFilePath=tempdir()){

  if (all(is.null(c(maf,plotList)))) {
    stop("Need to define at least a MAF file or a plot list.")
  }
  # if (is.null(masterRmdFile)) {
  #   masterRmdFile <- system.file('extdata', 'MAFDash.Rmd', package = 'MAFDashRPackage')
  # }
  ### Add checks for the conditions
  outputFilePath <- ensurer::ensure_that(outputFilePath,
                              !is.null(.) && (file.exists(outputFilePath)),
                              err_desc = "Need to give the right output file path.")

  html_filename=paste0(basename(outputFileName))
  # html_filename=gsub(".Rmd",".html",basename(masterRmdFile))
  masterRmdFile <- system.file("extdata", "MAFDash.Rmd", package = "MAFDash")
  rmarkdown::render(masterRmdFile,
                    output_format="all", output_file=html_filename,
                    output_dir = outputFilePath,
                    intermediates_dir = tempfile(),
                    params = list(
                      maffile=maf,
                      titletext=outputFileTitle,
                      plot_list=plotList
                      )
                    )
  ##This is fixed now in the 2.5 version
  ### rmarkdown::render doesn't let you select output destination (it uses the path of the Rmd file)
  ##  So this bit will move the report to the path in the 'out_dir' variable
  # outputFilePath=dirname(outputFileName)
  # if (!dir.exists(outputFilePath)) { dir.create(outputFilePath, recursive = T) }
  # file.rename(file.path(dirname(masterRmdFile),html_filename), file.path(outputFilePath,html_filename))
}

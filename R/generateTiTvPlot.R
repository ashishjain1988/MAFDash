#' Function to generate the frequency of Transitions and Transversions of gene mutations
#' @description This function generate the frequency of Transitions and Transversions of gene mutations
#' @author Mayank Tandon, Ashish Jain
#' @param maf A MAF object
#' @param use_silent_mutations Include synonymous variants in analysis. Defaults to FALSE.
#' @param plotType Can be 'bar', 'box' or 'both'. Defaults to 'both'
#' @param color named vector of colors for each coversion class.
#' @param showBarcodes Whether to include sample names for barplot
#' @param sampleOrder Sample names in which the barplot should be ordered. Default NULL
#' @param textSize fontsize if showBarcodes is TRUE. Deafult 2.
#' @param baseFontSize font size. Deafult 1.
#' @param axisTextSize text size x and y tick labels. Default c(1,1).
#' @param plotNotch logical. Include notch in boxplot.
#' @export
#' @return No return value. If 'save_name' is not provided, then the plot is printed to the current graphics device, otherwise a PDF is created at the given path.
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#'
generateTiTvPlot<-function(maf,use_silent_mutations = FALSE,plotType = 'both', sampleOrder = NULL,
                           color = NULL, showBarcodes = FALSE, textSize = 0.8, baseFontSize = 1,
                           axisTextSize = c(1, 1), plotNotch = FALSE,save_name=NULL){

  maf <- ensurer::ensure_that(maf,!is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  use_silent_mutations <- ensurer::ensure_that(use_silent_mutations,
                                               !is.null(.) && (class(.) == "logical"),
                                               err_desc = "Please enter the use_silent_mutations flag in correct format.")
  plotType <- ensurer::ensure_that(plotType,
                                   !is.null(.) && (class(.) == "character") && ((.) %in% c("bar","box","both")),
                                   err_desc = "Please enter plotType in correct format.")

  tiTvResults <- maftools::titv(maf,useSyn=use_silent_mutations,plot = FALSE)

  if(is.null(color)){
    col = c("#F44336", "#3F51B5", "#2196F3", "#4CAF50", "#FFC107", "#FF9800")
    col = grDevices::adjustcolor(col = col, alpha.f = 0.8)
    names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  }else{
    col = color
  }
  titv.frac = tiTvResults$fraction.contribution
  titv.frac.melt = data.table::melt(data = titv.frac, id = 'Tumor_Sample_Barcode')
  conv.class = c('Ti', 'Ti', 'Tv', 'Tv', 'Tv', 'Tv')
  names(conv.class) = c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")
  titv.frac.melt$TiTv = conv.class[as.character(titv.frac.melt$variable)]

  data.table::setDT(x = tiTvResults$TiTv.fractions)
  titv.contrib = suppressMessages(data.table::melt(tiTvResults$TiTv.fractions, id = 'Tumor_Sample_Barcode'))
  titv.frac.melt$variable = factor(x = titv.frac.melt$variable,
                                   levels = c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G"))

  titv.order = titv.frac.melt[,mean(value), by = .(variable)]
  titv.order = titv.order[order(V1, decreasing = TRUE)]
  orderlvl = as.character(titv.order$variable)
  titv.frac.melt$variable = factor(x = titv.frac.melt$variable, levels = rev(orderlvl))

  tf = tiTvResults$TiTv.fractions
  data.table::setDF(x = tf)
  rownames(tf) = tf$Tumor_Sample_Barcode
  tf = tf[,-1]

  # if(plotType == 'bar'){

  ##TODO: Make the plot look better and increase font sizes
  barplot <- ggplot2::ggplot(titv.frac.melt,aes(x=Tumor_Sample_Barcode,y=value, fill=variable)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_fill_manual(values=col) +
    ggplot2::scale_y_continuous(expand = c(0, 0),limits = c(0,100.2)) +
    ggplot2::xlab("") + ggplot2::ylab("% of Mutations") +
    ggplot2::guides(fill = guide_legend(title = "TiTv")) +
    ggplot2::theme_classic(base_size = 10)+
    ggplot2::theme(axis.text.x = element_blank())

  # } else if(plotType == 'box'){

  titv.frac.melt$variable = factor(x = titv.frac.melt$variable, levels = orderlvl)
  ##TODO: Make the plot look better and increase font sizes
  boxplot<-ggplot2::ggplot(titv.frac.melt,aes(x=variable,y=value,fill = variable))+
    ggplot2::theme_classic(base_size = 10) +
    #ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 20,colour = "black"),axis.title = element_text(size = 20,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
    scale_fill_manual(values=col) +
    #ggplot2::guides(fill = guide_legend(title = "TiTv")) +
    #ggplot2::geom_point(size=1.5)+
    ggplot2::guides(fill="none") +
    ggplot2::labs(x="", y = "% of Mutations")
  #ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
  #theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))

  return(list(boxplot=boxplot,barplot=barplot))
}

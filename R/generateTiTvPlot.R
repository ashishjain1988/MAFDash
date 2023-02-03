## To Supress Note
utils::globalVariables(c("value","variable","V1"))

#' Function to plot the frequency of Transitions and Transversions of gene mutations
#' @description This function plot the frequency of Transitions and Transversions of gene mutations
#' @author Ashish Jain, Mayank Tandon
#' @param maf A MAF object
#' @param use_silent_mutations Include synonymous variants in analysis. Defaults to FALSE.
#' @param sampleOrder Sample names in which the barplot should be ordered. Default NULL.
#' @param color named vector of colors for each coversion class.
#' @param showBarcodes Whether to include sample names for barplot.
#' @param textSize fontsize if showBarcodes is TRUE. Default 15 (in pts).
#' @param baseFontSize font size of axis title. Default 15 (in pts).
#' @param axisTextSize text size x and y tick labels. Default 10 (in pts).
# @param save_name The name and path of the output file. Default NULL.
#' @export
#' @return List of objects consisting of different TiTv plots.
# If 'save_name' is provided, a PDF is created at the given path.
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' plots<-generateTiTvPlot(read.maf(maf))
#'
#' @importFrom tidyr pivot_longer
#' @importFrom plotly subplot
#'
generateTiTvPlot<-function(maf,use_silent_mutations = FALSE, sampleOrder = NULL,
                           color = NULL, showBarcodes = FALSE, textSize = 15, baseFontSize = 15,
                           axisTextSize = 10){#,save_name=NULL){

  maf <- ensurer::ensure_that(maf,!is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  use_silent_mutations <- ensurer::ensure_that(use_silent_mutations,
                                               !is.null(.) && (class(.) == "logical"),
                                               err_desc = "Please enter the use_silent_mutations flag in correct format.")
  # plotType <- ensurer::ensure_that(plotType,
  #                                  !is.null(.) && (class(.) == "character") && ((.) %in% c("bar","box","both")),
  #                                  err_desc = "Please enter plotType in correct format.")

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
  titv.melt = tf %>% tidyr::pivot_longer(!Tumor_Sample_Barcode, names_to = "variable", values_to = "value")

  if(!is.null(sampleOrder))
  {
    titv.frac.melt$Tumor_Sample_Barcode = factor(x = titv.frac.melt$Tumor_Sample_Barcode, levels = sampleOrder)
  }

  ##TODO: Make the plot look better and increase font sizes
  barplot <- ggplot2::ggplot(titv.frac.melt,aes(x=Tumor_Sample_Barcode,y=value, fill=variable)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_fill_manual(values=col) +
    ggplot2::scale_y_continuous(expand = c(0, 0),limits = c(0,100.2)) +
    ggplot2::xlab("") + ggplot2::ylab("% of Mutations") +
    ggplot2::guides(fill = guide_legend(title = "TiTv")) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = axisTextSize,colour = "black"),
                   axis.title = element_text(size = baseFontSize,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::theme(axis.title.x = element_text(size = textSize,face = "bold"))

  if(!showBarcodes)
  {
    barplot <- barplot + ggplot2::theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  }

  titv.frac.melt$variable = factor(x = titv.frac.melt$variable, levels = orderlvl)
  ##TODO: Make the plot look better and increase font sizes
  boxplot<-ggplot2::ggplot(titv.frac.melt,aes(x=variable,y=value,fill = variable))+
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
    scale_fill_manual(values=col) +
    ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = axisTextSize,colour = "black"),
                   axis.title = element_text(size = baseFontSize,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::guides(fill="none") +
    ggplot2::labs(x="", y = "% of Mutations")

  ##TODO: Make the plot look better and increase font sizes
  boxplot1<-ggplot2::ggplot(titv.melt,aes(x=variable,y=value,fill = variable))+
    ggplot2::theme_classic(base_size = 10) +
    #ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 20,colour = "black"),axis.title = element_text(size = 20,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::geom_boxplot(outlier.colour="grey", outlier.shape=16,outlier.size=1, notch=FALSE)+
    scale_fill_manual(values=col) +
    ggplot2::guides(fill="none") +
    ggplot2::labs(x="", y = "% of Mutations") +
    ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = axisTextSize,colour = "black"),
                   axis.title = element_text(size = baseFontSize,face = "bold"),legend.background = element_rect(colour = "black"))

  # if(!is.null(save_name))
  # {
  #   #combinePlot <- cowplot::plot_grid(cowplot::plot_grid(boxplot,boxplot1,ncol = 2,rel_widths = c(2,1)),barplot,nrow = 2)
  #   combinePlot <- plotly::subplot(plotly::subplot(boxplot,boxplot1, nrows = 1, widths = c(0.5, 0.25)),barplot,nrows = 2,show.legend=FALSE)
  #   #ggsave(save_name,combinePlot)
  # }
  return(list(tiTvPatterns=boxplot,TiTv = boxplot1,barplot=barplot))
}

#' Function to generate a ribbon plot depicting co-occurence and mutual exclusivity of gene mutations
#' @description This function provides an alternate visualization for maftools::somaticInteractions()
#' @author Mayank Tandon, Ashish Jain
#' @param maf A MAF object
#' @param onco_genes A list of genes to restrict the analysis.  Passed to maftools::somaticInteractions()
#' @param save_name The name and path of the output PDF
#' @param pval_high All interactions with less than this p-value will be shown
#' @param pval_low Links with p-value less than this will be shown in a darker color
#' @param plot_type 'ribbon' returns a customized chord diagram, 'matrix' returns the default somaticInteractions() plot
#' @param shrink_factor Higher values = more shrinkage; use to control whitespace (or lack thereof) around figure.  Mostly useful in 0.5 - 1.5 range.
#' @param rotate_plot_degrees Rotate default layout by this many degrees
#' @param plot_frac_mut_axis Whether or not to draw a numerical axis on the perimeter
#' @param scale_ribbon_to_fracmut Whether or not to scale ribbon widths to their frequency
#' @param sig_colors Vector of 4 colors for coloring significance
#' @param gene_colors Color(s) for gene segments.  By default, they're colored randomly.
#' @export
#' @return No return value. If 'save_name' is not provided, then the plot is printed to the current graphics device, otherwise a PDF is created at the given path.
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' \donttest{generateRibbonPlot(read.maf(maf),save_name=paste0(tempdir(),"/ribbonPlot.pdf"))}
#'
generateRibbonPlot<-function(maf, onco_genes=NULL, save_name=NULL,
                                  pval_high=0.1,  ## All interactions with less than this p-value will be shown
                                  pval_low=0.05,  ## Links with p-value less than this will be highlighted with a dashed border
                                  plot_type="ribbon",  ## 'ribbon' returns a customized chord diagram, 'matrix' returns maftools's somaticInteractions() plot
                                  plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
                                  rotate_plot_degrees=0,   ## For custom rotation
                                  shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure
                                  scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
                                  sig_colors=NULL,   ## Vector of 4 colors for coloring significance
                                  gene_colors=NULL   ## color(s) for gene blocks
){
  # pval_low <- 0.05
  # browser()
  #require(circlize)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  } else {
    if (!plot_type=="matrix") {
      pdf(file = NULL)
    }
  }
  # browser()
  # if (is.null(onco_genes)) {
  #   onco_genes = maf@gene.summary$Hugo_Symbol
  # }

  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high))
  if (!plot_type=="matrix") {
    dev.off()
  } else {
    return(invisible())
  }
  # browser()
  cooccur_data <- som_int
  cooccur_data$pair_string <- apply(cooccur_data[,1:2], 1, function(x) {paste0(sort(x), collapse="_")})
  cooccur_data$popfrac <- NA
  cooccur_idx <- cooccur_data$Event=="Co_Occurence"
  mut_excl_idx = cooccur_data$Event=="Mutually_Exclusive"

  if (scale_ribbon_to_fracmut) {
    cooccur_data$popfrac1[cooccur_idx] <- unlist(cooccur_data[cooccur_idx,"11"]/as.numeric(maf@summary$summary[3]))
    cooccur_data$popfrac2[cooccur_idx] <- cooccur_data$popfrac1[cooccur_idx]
    cooccur_data$popfrac1[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"10"]/as.numeric(maf@summary$summary[3]))
    cooccur_data$popfrac2[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"01"]/as.numeric(maf@summary$summary[3]))
  } else {
    cooccur_data$popfrac1 <- 1
    cooccur_data$popfrac2 <- 1
  }
  chord_data <- cooccur_data[,c("gene1","gene2","popfrac1","popfrac2","pValue","Event")]
  chord_data[which(is.na(chord_data[,3])),3] <- 0

  if (is.null(sig_colors)) {
    sig_colors = RColorBrewer::brewer.pal(5, "BrBG")
    sig_colors <- sig_colors[-3]
  }
  names(sig_colors) <- paste0(rep(c("Co-occurence", "Mutually Exclusive"), each=2), " p-val < ",c(pval_low, pval_high, pval_high, pval_low ))
  color_legend <- ComplexHeatmap::Legend(labels=names(sig_colors),
                         legend_gp = grid::gpar(fill = sig_colors, col=sig_colors),background = sig_colors,size=unit(0.08,"npc"),
                         type="points",direction="vertical")


  chord_data$color_category <- paste0(ifelse(cooccur_data$Event=="Co_Occurence", "Co-occurence", "Mutually Exclusive"),
                                      paste0( " p-val < ", ifelse(cooccur_data$pValue < pval_low, pval_low, pval_high)))
  chord_data$color_val <- sig_colors[chord_data$color_category]

  #require(RColorBrewer)
  # browser()
  interacting_genes <- unique(unlist(chord_data[,1:2]))
  if (is.null(gene_colors)) {
    gene_colors <- colorRampPalette(brewer.pal(8,"Accent"))(length(interacting_genes))
    # gene_colors <- colorRampPalette(brewer.pal(8,"Dark2"))(length(interacting_genes))
    # gene_colors <- colorRampPalette(brewer.pal(8,"Set1"))(length(interacting_genes))
    # gene_colors <- rainbow((length(interacting_genes)))
  }
  if (length(gene_colors) != length(interacting_genes)) {
    # tmpcolors <- rep("grey90", length(interacting_genes))
    tmpcolors <- rep(gene_colors, length(interacting_genes))
    # tmpcolors[1:length(interacting_genes)] <- gene_colors
    gene_colors <- tmpcolors
  }
  names(gene_colors) <- interacting_genes

  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    pdf(file = save_name,height=7,width=7)
  }

  circos.clear()
  circos.par(canvas.xlim=c(-shrink_factor,shrink_factor),
             canvas.ylim=c(-shrink_factor,shrink_factor),
             start.degree = rotate_plot_degrees,
             message=FALSE)
  # circos.par$message = FALSE
  chordDiagram(chord_data[,1:4],grid.col = gene_colors,
               annotationTrack = c("grid",ifelse(plot_frac_mut_axis, "axis", "")),
               col=chord_data$color_val,
               # transparency = link_alpha,
               link.lty = 0,
               link.border = "black",
               link.sort = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  ComplexHeatmap::draw(color_legend, x = unit(0.05, "npc"), y = unit(0.5, "npc"), just = c("left"))
  # draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))

  if (!is.null(save_name)) {
    dev.off()
  }
}

#' Function to generate a dashboard from a MAF file
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param maf The MAF object
#' @param cohort_freq_thresh Cohort Frequency Threshold
#' @param auto_adjust_threshold Auto Adjusted Threshold flag
#' @param oncomat_only Oncomat only flag
#' @param clin_data Clinical data
#' @param clin_data_colors Clinical data colors
#' @export
#' @return The onco plot object.
#'
#' @examples
#' library(MAFDashRPackage)
#' #g <- generateOncoPlot(maf)
#' #g
#'
generateOncoPlot<-function(maf, cohort_freq_thresh = 0.01, auto_adjust_threshold=T,
                           oncomat_only=F,
                           clin_data=NULL, clin_data_colors=NULL){

  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                                !is.null(.) && (class(.) == "MAF"),
                                err_desc = "Please enter correct MAF object")
  cohort_freq_thresh <- ensurer::ensure_that(cohort_freq_thresh,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the cohort_freq_thresh in correct format.")
  auto_adjust_threshold <- ensurer::ensure_that(auto_adjust_threshold,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the auto_adjust_threshold flag in correct format.")
  oncomat_only <- ensurer::ensure_that(oncomat_only,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the oncomat_only flag in correct format.")
  #require(ComplexHeatmap)
  ### Read in MAF file
  # maf <- read.maf(maf_file)

  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf@gene.summary$Hugo_Symbol,
                         frac_mut=(maf@gene.summary$MutatedSamples/as.numeric(maf@summary$summary[3])),
                         stringsAsFactors = F)

  ngene_max=25
  target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
  if (auto_adjust_threshold) {
    cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
  }
  ### Select genes based on the frequency threshold
  freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
  freq_genes <- freq_genes[1:ngene_max]
  if (length(freq_genes) == 0) {
    stop("No genes to plot; change the frequency threshold to include more genes.")
  }
  if (length(freq_genes) > 100) {
    target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))],2)
    # stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    warning(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    # return(NA)
  }
  gene_list <- list(freq_genes)
  reasons <- paste0("Cohort Freq > ",cohort_freq_thresh)

  ### Collect genes to plot
  ##AJ: Optimized the code
  #genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
  genes_for_oncoplot_list<-list()
  for (i in 1:length(gene_list)) {
    if (is.na(gene_list[[i]][1])) {
      next
    }
    genes_for_oncoplot_list[[i]]<-data.frame(Hugo_Symbol=gene_list[[i]],
                                             reason=reasons[i],stringsAsFactors = TRUE)
    #genes_for_oncoplot <- rbind(genes_for_oncoplot,data.frame(Hugo_Symbol=gene_list[[i]],reason=reasons[i]))
  }
  genes_for_oncoplot <-do.call("rbind",genes_for_oncoplot_list)
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])

  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]

  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=genes_for_oncoplot$reason
  ##AJ: Updated the code here
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
  split_colors <- list(Reason=split_colors)

  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  oncomat <- createOncoMatrix(maf, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
  onco_genes <- rownames(oncomat)

  if (oncomat_only) {
    return(oncomat)
  }
  oncomat.plot <- oncomat

  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
  }

  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)

  ### Column labels get cluttered if too many samples
  show_sample_names=T
  if (ncol(oncomat.plot) > 20) {
    show_sample_names=F
  }

  # browser()
  myanno=NULL
  if (!is.null(clin_data)) {

    myanno <- make_column_annotation(clin_data,colnames(oncomat.plot), clin_data_colors)
    # print(myanno)
  }

  ## Show total burden for top annotation
  variant_type_data <- data.frame(maf@variant.classification.summary)
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  top_ha = HeatmapAnnotation("Total\nMutations" = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F),
                             annotation_name_side = "left",annotation_name_rot=90,annotation_name_gp = gpar(cex=0.7))

  # browser()

  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(cex=0.7)), show_annotation_name=F)
  # print(oncomat.plot)
  ### Make the oncoplot
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 show_pct = F,
                                 row_split=split_idx,
                                 row_title = NULL,
                                 bottom_annotation = myanno,
                                 top_annotation = top_ha,
                                 left_annotation = left_ha,
                                 show_column_names = show_sample_names)#,

  ### Return the oncoplot
  return(onco_base_default)
}

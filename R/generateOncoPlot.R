## To Supress Note
utils::globalVariables(c(".", "..anno_columns"))
#' Function to generate a dashboard from a MAF file
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tandon, Ashish Jain
#' @param maf The MAF object
#' @param cohort_freq_thresh Fraction of cohort that a gene must be mutated to select for display
#' @param auto_adjust_cohort_freq Whether or not to automatically adjust the frequen
#' @param genes_to_plot Character vector, data frame, or tab-delimited file name with genes to plot.
#' Data frame or file should contain a column named "Hugo_Symbol" with gene symbols, and optionally,
#' a column named "Reason" for labeling the plot
#' @param include_all Flag to include all the samples including the missing one (Default: FALSE)
#' @param oncomat_only Whether or not to return just the oncoplot matrix
#' @param title_text The title of the plot
#' @param custom_column_order A list containing the order of samples to show in the plot (Optional)
#' @param add_clinical_annotations Whether or not to try to plot column annotations from the 'clinical.data' slot of the MAF object
#' @param clin_data_colors Named list of colors for clinical annoations
#' @export
#' @return A ComplexHeatmap object if 'oncomat_only' is FALSE or a character matrix if 'oncomat_only' is TRUE.
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' generateOncoPlot(read.maf(maf))
#'
generateOncoPlot<-function(maf, cohort_freq_thresh = 0.01, auto_adjust_cohort_freq=TRUE,
                           genes_to_plot=NULL, include_all=FALSE,
                           oncomat_only=FALSE, title_text="",
                           custom_column_order=NULL,
                           add_clinical_annotations=FALSE, clin_data_colors=NULL){

  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                                !is.null(.) && (class(.) == "MAF"),
                                err_desc = "Please enter correct MAF object")
  cohort_freq_thresh <- ensurer::ensure_that(cohort_freq_thresh,
                                       is.null(.) || (class(.) == "numeric"),
                                       err_desc = "Please enter the cohort_freq_thresh in correct format.")
  auto_adjust_cohort_freq <- ensurer::ensure_that(auto_adjust_cohort_freq,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the auto_adjust_cohort_freq flag in correct format.")
  oncomat_only <- ensurer::ensure_that(oncomat_only,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the oncomat_only flag in correct format.")

  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf@gene.summary$Hugo_Symbol,
                         frac_mut=(maf@gene.summary$MutatedSamples/as.numeric(maf@summary$summary[3])),
                         stringsAsFactors = FALSE)

  selected_genes <- geneSelectParser(genes_to_plot)
  selected_genes <- selected_genes[selected_genes$Hugo_Symbol %in% frac_mut$Hugo_Symbol,]

  freq_genes_df <- data.frame(Hugo_Symbol=c(), Reason=c(), stringsAsFactors = FALSE)
  if (!is.null(cohort_freq_thresh)) {
    ngene_max=25
    target_frac = sort(frac_mut$frac_mut, decreasing = TRUE)[min(ngene_max,nrow(frac_mut))]
    if (auto_adjust_cohort_freq) {
      cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
    }
    ### Select genes based on the frequency threshold
    freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
    freq_genes <- freq_genes[1:ngene_max]
    freq_genes <- freq_genes[!is.na(freq_genes)]
    if (length(freq_genes) > 0) {
      if (length(freq_genes) > 100) {
        target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))],2)
        # stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
        warning(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
        # return(NA)
      }
      freq_genes_df <- data.frame(Hugo_Symbol=freq_genes,
                               Reason=paste0("Cohort Freq > ",round(cohort_freq_thresh, digits = 3)),
                               stringsAsFactors = F)
    }

  }
  # browser()
  genes_for_oncoplot <- rbind(freq_genes_df, selected_genes)
  if (nrow(genes_for_oncoplot) < 1) {
    stop("No genes to plot; Try selecting different genes or adjust frequency parameters")
  }
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol, nomatch=0)])

  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$Reason, -genes_for_oncoplot$frac),]

    ### Split the oncoplot based on the reason for picking the gene
  split_idx=genes_for_oncoplot$Reason
  ##AJ: Updated the code here
  split_colors <- rainbow(length(unique(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$Reason[!duplicated(genes_for_oncoplot$Reason)])
  split_colors <- list(Reason=split_colors)

  ### Make matrix to plot, and order it correctly
  oncomat <- createOncoMatrix(maf, g=genes_for_oncoplot$Hugo_Symbol, add_missing = include_all)$oncoMatrix
  ##AJ Update: drop=FALSE for #samples=1
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ,drop=FALSE]
  onco_genes <- rownames(oncomat)

  if (include_all) {
    missing_samples <- setdiff(maf@variants.per.sample$Tumor_Sample_Barcode, colnames(oncomat))
    emptydata <- matrix(data="",nrow=nrow(oncomat), ncol=length(missing_samples))
    colnames(emptydata) <- missing_samples
    rownames(emptydata) <- onco_genes
    oncomat <- cbind(oncomat,emptydata)
  }
  if (oncomat_only) {
    return(oncomat)
  }
  # browser()
  if (is.null(custom_column_order)) {
    custom_column_order <- colnames(oncomat)
  }
  # browser()
  col_order_idx <- match(custom_column_order, colnames(oncomat), nomatch=0)
  col_order_idx <- col_order_idx[col_order_idx>0]
  ##AJ Update: drop=FALSE for #samples=1
  oncomat.plot <- oncomat[,col_order_idx,drop=FALSE]
  # oncomat.plot <- oncomat


  ### Set the height of the plot based on number of genes
  # onco_height=NULL
  # if (is.null(onco_height)) {
  #   onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
  # }

  oncomat.plot <- gsub("_"," ",oncomat.plot)

  ### Column labels get cluttered if too many samples
  show_sample_names=TRUE
  if (ncol(oncomat.plot) > 20) {
    show_sample_names=FALSE
  }

  # browser()
  myanno=NULL
  clin_anno_data <- maf@clinical.data
  anno_columns <- colnames(clin_anno_data)
  if (is.character(add_clinical_annotations)) {
    anno_columns <- intersect(c(add_clinical_annotations,"Tumor_Sample_Barcode"), anno_columns)
    if (length(anno_columns) < 2) {
      warning("Annotation columns not found in clinical data, skipping...")
    }
    add_clinical_annotations=TRUE
  }

  if (is.logical(add_clinical_annotations) && add_clinical_annotations) {
    if (!is.null(clin_data_colors)) {
      continuous_color_columns <- names(clin_data_colors)[unlist(lapply(clin_data_colors, class))=="function"]
      for (col_name in continuous_color_columns) {
        clin_anno_data[[col_name]] <- as.numeric(as.character(clin_anno_data[[col_name]]))
      }
    }
    # browser()
    myanno <- make_column_annotation(my_clin_dat = clin_anno_data[,..anno_columns],names_to_match = colnames(oncomat.plot), my_colors = clin_data_colors)
  }


  ## Show total burden for top annotation
  variant_type_data <- data.frame(maf@variant.classification.summary)
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]

  mutation_colors <- my_mutation_colors()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  top_ha = ComplexHeatmap::HeatmapAnnotation("Total\nMutations" = ComplexHeatmap::anno_barplot(variant_type_data, gp = grid::gpar(fill = var_anno_colors), border = FALSE),
                             annotation_name_side = "left",annotation_name_rot=90,annotation_name_gp = grid::gpar(cex=0.7))

  # browser()

  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = ComplexHeatmap::rowAnnotation("Cohort Pct"=ComplexHeatmap::anno_text(pct_anno,gp = grid::gpar(cex=0.7)), show_annotation_name=FALSE)
  # print(oncomat.plot)
  ### Make the oncoplot
  ##TODO: Where oncoplot_annotation_func is coming or defined?
  onco_plot <- ComplexHeatmap::oncoPrint(oncomat.plot, alter_fun = oncoplot_annotation_func(),
                         col=mutation_colors,
                         row_order=1:nrow(oncomat.plot),
                         column_order=1:ncol(oncomat.plot),
                         # column_order=col_order_idx,
                         remove_empty_columns=!include_all,
                         name="oncoplot",
                         row_title=title_text,
                         show_pct = FALSE,
                         row_split=split_idx,
                         bottom_annotation = myanno,
                         top_annotation = top_ha,
                         left_annotation = left_ha,
                         show_column_names = show_sample_names)

  ### Return the oncoplot
  return(onco_plot)
}

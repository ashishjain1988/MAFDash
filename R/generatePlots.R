#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param file The path of the file containing the mutation
#' information in the MAF format
#' @export
#' @return The output is the html file.
#'
#' @examples
#' library(MAFDashRPackage)
#' MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
generateOncoplot<-function(maf.filtered, cohort_freq_thresh = 0.01, auto_adjust_threshold=T,
                           oncomat_only=F,
                           clin_data=NULL, clin_data_colors=NULL){

  require(ComplexHeatmap)
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)

  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$MutatedSamples/as.numeric(maf.filtered@summary$summary[3])),
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
  genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
  for (i in 1:length(gene_list)) {
    if (is.na(gene_list[[i]][1])) {
      next
    }
    genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                data.frame(Hugo_Symbol=gene_list[[i]],
                                           reason=reasons[i]))
  }
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])

  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]

  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=genes_for_oncoplot$reason
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
  split_colors <- list(Reason=split_colors)

  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
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
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
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



#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param file The path of the file containing the mutation
#' information in the MAF format
#' @export
#' @return The output is the html file.
#'
#' @examples
#' library(MAFDashRPackage)
#' MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
make_burden_plot<-function(maf.filtered, plotType=NULL, mb_covered=NULL, save_data_to_file=NULL){


  num_var_data <- maf.filtered@variants.per.sample
  colnames(num_var_data) <- c("Tumor_Sample_Barcode","Variants_filtered")
  num_var_data$mut_burden_count <- num_var_data$Variants_filtered
  num_var_data$mut_burden <- num_var_data$mut_burden_count
  y_label_text="Mutation Count"
  if (is.numeric(mb_covered)) {
    print("Normalizing mutation count by covered bases...")
    num_var_data$mut_burden <- num_var_data$mut_burden/mb_covered
    y_label_text="Mutation Burden (mutations/Mb)"
  }

  nsamples=nrow(num_var_data)
  if (is.null(plotType)) {
    plotType <- ifelse(nrow(num_var_data) > 15, "Dotplot", "Barplot")
    print(paste0("Using plot type: ", plotType))
  }
  # browser()
  ## Re-jigger the factor levels so they're ordered by decreasing mutation burden (for the plotting)
  num_var_data$Tumor_Sample_Barcode <- factor(num_var_data$Tumor_Sample_Barcode,
                                              levels=num_var_data$Tumor_Sample_Barcode[order(num_var_data$mut_burden, decreasing = T)])

  ########################################################
  #### 5. Generate plots for mutation burden

  ## Pick colors
  median_mut_burdens <- num_var_data %>% summarise(median=median(mut_burden))

  num_var_data$xlabel <- factor(num_var_data$xlabel,
                                levels=num_var_data$xlabel[order(num_var_data$mut_burden, decreasing = T)])
  num_var_data$hoverlabel <- paste0("Sample: ",num_var_data$Tumor_Sample_Barcode,"\nMutations: ", num_var_data$mut_burden)

  ### Mutation burden stacked with variant classification counts
  ### Works better for smaller cohorts, say < 20
  variant_type_per_sample <- as.data.frame(maf.filtered@variant.classification.summary)
  var_type.melt <- reshape2::melt(variant_type_per_sample, id.vars="Tumor_Sample_Barcode",variable.name="classification",value.name="mutation_count")
  var_type.melt$mut_burden <- var_type.melt$mutation_count
  if (is.numeric(mb_covered)) {
    var_type.melt$mut_burden <- var_type.melt$mut_burden/mb_covered
  }
  median_mut_burdens <- data.frame(median=median(var_type.melt[var_type.melt$classification== "total","mut_burden"]))

  plotdata <- var_type.melt[var_type.melt$classification != "total",]
  plotdata$Tumor_Sample_Barcode <- factor(as.character(plotdata$Tumor_Sample_Barcode),
                                          levels=variant_type_per_sample$Tumor_Sample_Barcode[order(variant_type_per_sample$total, decreasing = T)])
  plotdata$classification <- gsub("_"," ",plotdata$classification)

  class_means <- plotdata %>% group_by(classification) %>% summarise(mean=mean(mut_burden))
  plotdata$classification <- factor(as.character(plotdata$classification),
                                    levels=class_means$classification[order(class_means$mean, decreasing = F)])

  my_class_colors <- mutation_colors

  plotdata$hoverlabel <- paste0("Sample: ",plotdata$Tumor_Sample_Barcode,"\nMutations: ", plotdata$mut_burden)

  if (plotType=="Barplot") {
    if (length(unique(plotdata$Tumor_Sample_Barcode)) <= 20) {
      xaxis_text <- element_text(angle=30, hjust=1)
    } else {
      xaxis_text <- element_blank()
    }
    burden_plot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y=mut_burden, text=hoverlabel)) +
      geom_bar(aes(fill=classification), stat="identity",width=1,size=0.3, color="black") +
      scale_fill_manual(values=my_class_colors) +
      theme_linedraw(base_size = 12) +
      xlab("") + ylab(y_label_text) +
      # geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60") +  ### Is screwed up with ggplotly
      theme(
        axis.text.x = xaxis_text,
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.key.height = unit(0.01,"npc"),
        legend.key.width =  unit(0.02,"npc"),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())
    if (nsamples > 1) {
      burden_plot <- burden_plot + geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60")
    }
  } else {

    require(ggbeeswarm)
    ### Mutation Burden - Scatter/Dot plot
    ### Works better for larger cohorts
    alpha_val=1
    point_cex=2
    if (nrow(num_var_data) > 200) {
      alpha_val=0.5
    } else if (nrow(num_var_data) > 20) {
      alpha_val=0.8
    }
    burden_plot <- ggplot(num_var_data, aes(x=1, y=mut_burden, text=hoverlabel)) +
      # geom_beeswarm(color="blue3",cex=2,size=5,dodge.width=0.2,priority='density', alpha=alpha_val) +
      geom_quasirandom(color="blue3",width=0.3,size=5,alpha=alpha_val, method="quasirandom", bandwidth=0.1) +
      # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
      #              geom = "crossbar", width = 0.7, color="gray70", size = 0.2) +
      scale_y_log10()+
      theme_linedraw(base_size = 12) +
      ggtitle("Mutation Burden") +
      ylab(y_label_text) + xlab("") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }



  ## Write data to a file for external plotting if desired
  if (!is.null(save_data_to_file)) {
    if (dir.exists(dirname(save_data_to_file))) {
      outdata <- as.data.frame(maf.filtered@variant.classification.summary)
      outdata$total_per_mb <- outdata$total/mb_covered
      outdata$mb_covered <- mb_covered
      print(paste0("Saving plot data to ", save_data_to_file))
      write.table(outdata, file = save_data_to_file, sep="\t", quote=F,row.names = F)
    } else {
      warning("Path for output data file not found! Skipping...")
    }
  }


  return(burden_plot)
}


#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param file The path of the file containing the mutation
#' information in the MAF format
#' @export
#' @return The output is the html file.
#'
#' @examples
#' library(MAFDashRPackage)
#' MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
plot_silent_nonsilent<-function(mymaf, savename=NULL, returndata=F){

  nonsilent_summary <- mymaf@variant.classification.summary[,c("Tumor_Sample_Barcode","total")]
  nonsilent_summary$type <- "Non-Silent"
  silent_classif_data <- mymaf@maf.silent %>% group_by(Tumor_Sample_Barcode, Variant_Classification) %>% summarise(count=n())
  silent_classif_data <- reshape2::dcast(silent_classif_data,Tumor_Sample_Barcode ~ Variant_Classification, value.var = "count")
  silent_summary <- data.frame(Tumor_Sample_Barcode = silent_classif_data$Tumor_Sample_Barcode,
                               total = rowSums(silent_classif_data[,-1], na.rm=T),
                               type = "Silent"
  )

  # browser()
  plotdata <- rbind(nonsilent_summary, silent_summary)
  tots <- plotdata %>% group_by(Tumor_Sample_Barcode) %>% summarise(tot=sum(total))
  plotdata$Tumor_Sample_Barcode <- factor(plotdata$Tumor_Sample_Barcode,
                                          levels=as.character(tots$Tumor_Sample_Barcode)[order(tots$tot,decreasing = T)])
  myplot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y = total, fill=type)) +
    geom_col() + scale_fill_brewer(palette="Set1") +
    theme_linedraw(base_size = 12) +
    xlab("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle=30, hjust = 1, size=10),
          axis.ticks.x = element_blank())

  if (!is.null(savename)) {
    if (! dir.exists(dirname(savename))) {dir.create(dirname(savename), recursive = T)}
    ggsave(savename,width=6, height=6)
  }
  return_val <- myplot
  if (returndata) {
    return_val <- list(plot=myplot, data=plotdata)
  }
  return(return_val)
}


#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param file The path of the file containing the mutation
#' information in the MAF format
#' @export
#' @return The output is the html file.
#'
#' @examples
#' library(MAFDashRPackage)
#' MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
make_overlap_plot<-function(mymaf, use_silent_mutations=F,
                            summarize_by="gene",
                            plotType=c("ribbon","heatmap"),
                            savename="overlap_plot.pdf",
                            savewidth=8, saveheight=8){

  # browser()
  mafdata <- mymaf@data
  if (use_silent_mutations) {
    mafdata <- rbind(mafdata, mymaf@maf.silent)
  }

  mafdata_by_sample <- split(mafdata, mafdata$Tumor_Sample_Barcode)

  if (summarize_by=="gene") {
    id_cols <- c("Chromosome","Start_Position","End_Position","Hugo_Symbol")
  } else if (summarize_by=="mutation") {
    id_cols <- c("Chromosome","Start_Position","End_Position","Hugo_Symbol","HGVSp_Short")
  } else {
    stop("'summarize_by' must be either 'gene' or 'mutation'")
  }

  mutations_list <- lapply(mafdata_by_sample, function(currmafdata, mycols) {
    id_string <- apply(currmafdata[,..mycols],1,paste0, collapse="_")
    return(id_string)
  }, id_cols)

  pw_combinations <- matrix(0,nrow = length(mutations_list),ncol = length(mutations_list))
  colnames(pw_combinations) <- names(mutations_list)
  rownames(pw_combinations) <- names(mutations_list)
  for ( row_idx in 1:nrow(pw_combinations) ) {
    for (col_idx in 1:ncol(pw_combinations) ) {
      if (!row_idx==col_idx) {
        # overlap_val=sum(mutations_list[[ rownames(pw_combinations)[row_idx] ]]$id_str %in% mutations_list[[colnames(pw_combinations)[col_idx] ]]$id_str)
        overlap_val=length(intersect(mutations_list[[ rownames(pw_combinations)[row_idx] ]],mutations_list[[colnames(pw_combinations)[col_idx] ]]))
        pw_combinations[row_idx,col_idx]=overlap_val
      }
    }
  }


  hm_data <- pw_combinations

  pdf(savename, width=savewidth, height=saveheight)
  if ("heatmap" %in% plotType) {
    library(pheatmap)
    pheatmap(hm_data,cluster_rows = T, cluster_cols = T,
             clustering_method = "complete",
             # annotation_colors = hm_anno_colors,
             # annotation_col = hm_anno_data,annotation_row = hm_anno_data,
             main = "Hierarchically Clustered")
  }

  # pheatmap(log2(pw_combinations+1e-6),cluster_rows = F, cluster_cols = F)
  # dev.off()

  if ("ribbon" %in% plotType) {
    count_clustering <- hclust(dist(pw_combinations), method = "complete")
    cluster_order <- count_clustering$labels[count_clustering$order]

    chordData <- pw_combinations
    grid.col <- rainbow(n=nrow(chordData))
    names(grid.col) <- cluster_order

    chordDiagram(chordData,order=cluster_order,grid.col=grid.col,
                 annotationTrack = c("grid","axis"),
                 preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chordData))))))
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(-0.2, 0))
    }, bg.border = NA) # here set bg.border to NA is important
  }
  dev.off()
}

#' Function to generate a dashboard from a MAF file.
#' @description This function created an HTML file containing the
#' different figures and plots explaining the MAF dataset.
#' @author Mayank Tondon, Ashish Jain
#' @param file The path of the file containing the mutation
#' information in the MAF format
#' @export
#' @return The output is the html file.
#'
#' @examples
#' library(MAFDashRPackage)
#' MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #t <- getMAFDashboard(file = MAFfilePath)
#'
make_single_ribbon_plot<-function(maf, onco_genes=NULL, save_name=NULL, ribbon_color=NULL,
                                  pval_high=0.1,  ## All interactions with less than this p-value will be shown
                                  pval_low=0.05,  ## Links with p-value less than this will be highlighted with a dashed border
                                  plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
                                  rotate_plot_degrees=0,   ## For custom rotation
                                  shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure
                                  scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
                                  sig_colors=NULL,   ## Vector of 4 colors for coloring significance
                                  gene_colors=NULL   ## color(s) for gene blocks
){
  # pval_low <- 0.05
  # browser()
  require(circlize)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  } else {
    pdf(file = NULL)
  }
  # browser()
  # if (is.null(onco_genes)) {
  #   onco_genes = maf@gene.summary$Hugo_Symbol
  # }

  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high))
  dev.off()
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
  color_legend <- Legend(labels=names(sig_colors),
                         legend_gp = gpar(fill = sig_colors, col=sig_colors),background = sig_colors,size=unit(0.08,"npc"),
                         type="points",direction="vertical")


  chord_data$color_category <- paste0(ifelse(cooccur_data$Event=="Co_Occurence", "Co-occurence", "Mutually Exclusive"),
                                      paste0( " p-val < ", ifelse(cooccur_data$pValue < pval_low, pval_low, pval_high)))
  chord_data$color_val <- sig_colors[chord_data$color_category]

  require(RColorBrewer)
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
             message=F)
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
  draw(color_legend, x = unit(0.05, "npc"), y = unit(0.5, "npc"), just = c("left"))
  # draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))

  if (!is.null(save_name)) {
    dev.off()
  }
}

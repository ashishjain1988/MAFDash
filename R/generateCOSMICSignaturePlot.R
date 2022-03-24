## To Supress Note
utils::globalVariables(c("as.dist","colorRamp","rgb","gpar"))

#' Computes and plot the cosine similarity of each individual signature is computed against each COSMIC signature
#' from COMIC V3.2
#'
#' @description This function computes and plot the cosine similarity of each individual signature is computed against each COSMIC signature
#' from COMIC V3.2. The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors
#' are alike.
#'
#' @param mymaf mutation count matrix (dimensions: a mutation features X n samples)
#' @param use_silent_mutations 96 mutation count matrix (dimensions: a mutation features X m samples)
#' @param full_output return full output including the etiology matrix and plot data
#' @param show_broad_categories To show broad etiology categories
#' @param clin_data Clinical data to be plotted in the heatmap
#' @param clin_data_colors Clinical data colors
#' @param add_sample_names Whether or not to add column labels; if set to NULL, will add labels only if # samples less than 10
#' @param savename file name of the plot
#' @param fig_height Output height (inches); set to NULL to size automatically; only used if savename is set.
#' @param fig_width Output width (inches); set to NULL to size automatically; only used if savename is set.
#'
#' @return Complex Heatmap object. If full_output is TRUE it will consist of a list including heatmap object,
#' etiology matrix and plot data.
#'
#' @export
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' library(ComplexHeatmap)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' \donttest{val<-generateCOSMICMutSigSimHeatmap(read.maf(maf));draw(val)}
#'
generateCOSMICMutSigSimHeatmap <- function(mymaf,use_silent_mutations=FALSE,
                                       full_output=FALSE,
                                       show_broad_categories = TRUE,
                                       # genome_build="hg19",
                                       clin_data=NULL, clin_data_colors=NULL,
                                       #signatures_file=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_Mutational_Signatures_v3.1.xlsx"),
                                       #etio_data_xlsx=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_signature_etiology.xlsx"),
                                       add_sample_names = NULL,
                                       savename=NULL,
                                       fig_height=NULL,
                                       fig_width=NULL) {
  mymaf <- ensurer::ensure_that(mymaf,!is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  use_silent_mutations <- ensurer::ensure_that(use_silent_mutations,
                                               !is.null(.) && (class(.) == "logical"),
                                               err_desc = "Please enter the use_silent_mutations flag in correct format.")

  genome_build_res <- detectMAFGenome(mymaf)
  genome_build <- genome_build_res[[1]]
  genome_package <- genome_build_res[[3]]
  prefix_value=ifelse(genome_build_res[[2]],"chr","")
  tnm = maftools::trinucleotideMatrix(maf = mymaf,ref_genome = genome_package,prefix = prefix_value,useSyn = use_silent_mutations)
  mut_mat = t(tnm$nmf_matrix)
  #mut_mat <-  make_tnm(mymaf, use_silent_mutations = use_silent_mutations)

  signatures = "SBS"
  signatures_file=system.file("extdata", paste0("COSMIC_v3.2_",signatures,".txt"), package = "MAFDash")#file.path(paste0("/Users/jaina13/myPART/cosmic_etilogy_scripts/COSMIC_v3.2_",signatures,".txt"))
  etio_data_file=system.file("extdata", paste0("etiology_clustering.raw_",signatures,"_COSMICV3.2.txt"), package = "MAFDash")#file.path(paste0("/Users/jaina13/myPART/cosmic_etilogy_scripts/etiology_clustering.raw_",signatures,"_COSMICV3.2.txt"))
  #signatures_file=file.path(paste0("/Users/jaina13/myPART/cosmic_etilogy_scripts/COSMIC_v3.2_",signatures,".txt"))
  #etio_data_file=file.path(paste0("/Users/jaina13/myPART/cosmic_etilogy_scripts/etiology_clustering.raw_",signatures,"_COSMICV3.2.txt"))

  cosmic_signatures = read.table(signatures_file,header = T,sep = "\t")
  cosmic_signatures$Somatic.Mutation.Type <- cosmic_signatures$Type
  etiology_data <- data.frame(Etiology=NA,Etiology_category = NA,row.names=c(),stringsAsFactors = F)


  etiology_data_raw <- read.table(etio_data_file,header = T,sep = "\t")#read.xlsx(etio_data_xlsx, sheet=1)
  etiology_data <- data.frame(Etiology=etiology_data_raw$etiology,
                              Etiology_category = etiology_data_raw$etiologyBroad,
                              row.names=etiology_data_raw$signature,
                              stringsAsFactors = F)

  #l<-list();for(sig in c("SBS","DBS","ID")){etio_data_file=file.path(paste0("/Users/jaina13/myPART/cosmic_etilogy_scripts/etiology_clustering.raw_",sig,"_COSMICV3.2.txt"));etiology_data_raw <- read.table(etio_data_file,header = T,sep = "\t");l[[sig]]<-etiology_data_raw;}
  #df<-do.call(rbind,l)
  #names(table(df$etiologyBroad))
  categories <- sort(unique(etiology_data$Etiology_category))
  category_colors <- setNames(rainbow(length(categories)), categories)
  etiology_colors <- list(Etiology=unlist(lapply(categories, function(mycat) {
    subcats <- unique(etiology_data$Etiology[etiology_data$Etiology_category==mycat])
    catcolor=category_colors[mycat]
    myrgb <- colorRamp(c(catcolor, "white"))(0:length(subcats)/length(subcats))[1:(length(subcats)), , drop=F]
    rownames(myrgb) <- subcats
    hexcolor <- apply(myrgb, 1, function(x) {rgb(x[1], x[2], x[3], maxColorValue = 255)})
    return(hexcolor)
  })))
  rowOrder=order(etiology_data$Etiology_category)
  etiology_data <- etiology_data[rowOrder, ,drop=F]
  etiology_data$annotation_color <- etiology_colors$Etiology[match(etiology_data$Etiology, names(etiology_colors$Etiology))]
  etiology_data$category_color <- category_colors[match(etiology_data$Etiology_category, names(category_colors))]
  cosmiclist <- list(signatures=cosmic_signatures,etio=etiology_data)

  # sp_url <- file.path(data_dir,"cosmic","sigProfiler_exome_SBS_signatures.csv")
  # cosmic_signatures = read.table(sp_url, sep = ",", header = TRUE, stringsAsFactors = F)
  # cosmic_signatures$Somatic.Mutation.Type <- paste0(substr(cosmic_signatures$SubType, 1, 1),
  #                                                   "[",cosmic_signatures$Type, "]",
  #                                                   substr(cosmic_signatures$SubType, 3, 3))
  #
  cosmic_signatures <- cosmiclist[["signatures"]]

  # Match the order of the mutation types to MutationalPatterns standard
  new_order = match(row.names(mut_mat), cosmic_signatures$Somatic.Mutation.Type)
  # Reorder cancer signatures dataframe
  cosmic_signatures = cosmic_signatures[as.vector(new_order),]
  # Add trinucletiode changes names as row.names
  row.names(cosmic_signatures) = cosmic_signatures$Somatic.Mutation.Type
  # Keep only 96 contributions of the signatures in matrix
  cosmic_signatures = as.matrix(cosmic_signatures[,grep("SBS*", colnames(cosmic_signatures))])

  hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
  # store signatures in new order
  cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
  # plot(hclust_cosmic)

  # if (is.function(progress_func)) {
  #   progress_func(value=50, detail = "Computing similarity to COSMIC...")
  # }
  cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cosmic_signatures)

  # fit_res <- fit_to_signatures(mut_mat, cosmic_signatures)

  # data_dir="data"
  # etio_data_file = file.path(data_dir, "cosmic","COSMIC_signature_etiology.xlsx")
  # etiology_data_raw <- read.xlsx(etio_data_file, sheet="final categories")
  # etiology_data <- as.character(etiology_data_raw$CATEGORY)
  # etiology_data <- data.frame(Etiology=etiology_data_raw$CATEGORY, row.names=etiology_data_raw$signature)
  etio <- cosmiclist[["etio"]]

  if(show_broad_categories)
  {
    etiology_data <- data.frame(Etiology=paste0(etio$Etiology_category),
                                row.names=rownames(etio))
  }else
  {
    etiology_data <- data.frame(Etiology=paste0("[",etio$Etiology_category,"] ", etio$Etiology),
                                row.names=rownames(etio))
  }

  # names(etiology_data) <- etiology_data_raw$signature
  # etiology_colors <-  list(Etiology=c("APOBEC" = "#fce116",
  #                                     "Defective DNA Mismatch Repair" = "#31A354",
  #                                     "Defective DNA Repair" = "#A1D99B",
  #                                     "Exonuclease Domain" = "#E5F5E0",
  #                                     "Exposure to Alfatoxin" = "#DE2D26",
  #                                     "Exposure to Aristolochic Acid" = "#FC9272",
  #                                     "Exposure to Haloalkanes" = "#FEE0D2",
  #                                     "Tobacco - Chewing" = "#6d3617",
  #                                     "Tobacco - Smoking" = "#a85423",
  #                                     "Tobacco - Smoking Associated" = "#d87841",
  #                                     "Prior Therapy - Alkylating Agents" = "#2171B5",
  #                                     "Prior Therapy - Platinum Drugs" = "#6BAED6",
  #                                     "Prior Therapy - Immunosuppression" = "#BDD7E7",
  #                                     "Prior Therapy Associated" = "#EFF3FF",
  #                                     "ROS Damage" = "#BCBDDC",
  #                                     "UV Light" = "#756BB1",
  #                                     "UV Light Associated" = "#a29bca",
  #                                     "Unknown" = "grey70",
  #                                     "Sequencing Artifact" = "grey50")
  # )
  etiology_colors <- list(Etiology=setNames(etio$annotation_color, etiology_data$Etiology))
  plot_matrix <- t(cos_sim_samples_signatures)

  # if (is.function(progress_func)) {
  #   progress_func(value=80, detail = "Making heatmap")
  # }

  # browser()
  rowOrder=order(etiology_data)
  etiology_data <- etiology_data[rowOrder,1,drop=F]
  mycolors <- etiology_colors
  signature_anno <- rowAnnotation(df=etiology_data,
                                  name="Signature Anno", col=etiology_colors, show_annotation_name = FALSE)

  plot_matrix <- plot_matrix[match(rownames(etiology_data),rownames(plot_matrix)),]

  myanno=NULL
  if (!is.null(clin_data)) {
    anno_data <- data.frame(clin_data[match(colnames(plot_matrix), clin_data$Tumor_Sample_Barcode),],stringsAsFactors = F)
    row.names(anno_data) <- anno_data$Tumor_Sample_Barcode
    anno_data <- anno_data[,!colnames(anno_data) %in% "Tumor_Sample_Barcode"]
    myanno <- ComplexHeatmap::HeatmapAnnotation(df=anno_data,col = clin_data_colors)
  }

  if (is.null(add_sample_names)) {
    add_sample_names <- ifelse(ncol(plot_matrix)>10, F, T)
  }
  plot_matrix[is.nan(plot_matrix)] <- 0
  #brewer.pal(9, "BuGn")
  myHM <- ComplexHeatmap::Heatmap(plot_matrix,
                  col=circlize::colorRamp2(seq(min(plot_matrix), max(plot_matrix), length.out = 20),colorRampPalette(c("#F7FCFD","#E5F5F9","#CCECE6","#99D8C9","#66C2A4","#41AE76","#238B45","#006D2C","#00441B"))(20)),
                  left_annotation = signature_anno,
                  bottom_annotation = myanno,
                  cluster_rows = F, row_order = rowOrder,
                  clustering_method_rows = "median",
                  clustering_method_columns = "median",
                  heatmap_height = unit(6, "inches"),
                  heatmap_legend_param = list(
                    # at = c(-2, 0, 2),
                    # labels = c("low", "zero", "high"),
                    title = "Cosine Similarity",
                    # legend_height = unit(4, "cm"),
                    legend_direction = "horizontal"
                  ),
                  show_row_names=T, #row_names_gp = gpar(fontsize = 5),
                  show_column_names = add_sample_names)

  if ( ! is.null(savename) ) {
    if (is.null(fig_height)) {
      anno_height=ifelse(!is.null(clin_data), min(c(4, 0.5*ncol(clin_data))), 0)
      fig_height=max(round(0.1*nrow(plot_matrix),0),4) + anno_height
    }
    if (is.null(fig_width)) {
      fig_width=fig_height*0.75 + anno_height*1.2
    }
    pdf(file = savename,height=fig_height,width=fig_width)
    draw(myHM)
    dev.off()
  }

  ### Return the oncoplot (if function is pointed to a variable)
  returnval <- myHM
  if (full_output) {
    output_list <- list(plot_matrix=plot_matrix,
                        signature_annotations=signature_anno,
                        etiology_data=etiology_data,
                        heatmap_obj=myHM)
    returnval <- output_list
  }

  return(returnval)
  #invisible(myHM)
}

# make_signature_plot <- function(tnm, savepath=NULL, signature_type="SBS") {
#   # Basically stolen from: https://github.com/UMCUGenetics/MutationalPatterns/blob/master/R/plot_96_profile.R
#   require(ggplot2)
#   mydf <- data.frame(substitution = gsub(".*\\[(.*)\\].*","\\1",rownames(tnm)),
#                      context=gsub("\\[.*\\]","\\.",rownames(tnm)),
#                      apply(tnm,2,function(x){round(x/sum(x),digits = 5)}),
#                      stringsAsFactors = F)
#   plotdf <- reshape2::melt(mydf, id.vars=c("substitution","context"))
#   plotdf <- plotdf[!is.nan(plotdf$value),]
#   ymax=max(plotdf$value)
#   COLORS6 <- c(
#     "#2EBAED", "#000000", "#DE1C14",
#     "#D4D2D2", "#ADCC54", "#F0D0CE"
#   )
#   # browser()
#   plot_96 = ggplot(data = plotdf, aes(x = context, y = value,
#                                       fill = substitution, width = 1)) +
#     geom_bar(stat = "identity", colour = "black", size = 0.2) +
#     scale_fill_manual(values = COLORS6) +
#     facet_grid(variable ~ substitution) + ylab("Relative contribution") +
#     coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0,
#                                                                          ymax, 0.1)) + guides(fill = FALSE) + theme_bw() +
#     theme(axis.title.y = element_text(size = 12, vjust = 1),
#           axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12),
#           axis.text.x = element_text(size = 5, angle = 90,
#                                      vjust = 0.4), strip.text.x = element_text(size = 9),
#           strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(),
#           panel.spacing.x = unit(0, "lines"))
#   # plot
#   if (!is.null(savepath)) {
#     myheight=ncol(mydf)*0.8
#     ggsave(savepath,plot=plot_96,height=myheight, width=8, limitsize=F)
#     return(savepath)
#   } else {
#     return(plot_96)
#   }
# }
#
#
# signature_plot_colorpalette <- function(signature_type="SBS") {
#
#   COLORS6 <- c(
#     "#2EBAED", "#000000", "#DE1C14",
#     "#D4D2D2", "#ADCC54", "#F0D0CE"
#   )
#
# }

#' Code copied from MutationalPatterns package
#' Compute all pairwise cosine similarities between mutational profiles/signatures
#'
#' Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices.
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#'
#' @param mut_matrix1 mutation count matrix (dimensions: a mutation features X n samples)
#' @param mut_matrix2 96 mutation count matrix (dimensions: a mutation features X m samples)
#' @return Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
#'
#' @noRd
cos_sim_matrix <- function(mut_matrix1, mut_matrix2) {

  # Check that both inputs are numeric.
  if (!all(apply(mut_matrix1, 2, is.numeric))){
    stop("The first input contains non-numeric columns, while all columns should be numeric.")
  }
  if (!all(apply(mut_matrix2, 2, is.numeric))){
    stop("The second input contains non-numeric columns, while all columns should be numeric.")
  }

  # Determine number of samples
  n_samples1 <- ncol(mut_matrix1)
  n_samples2 <- ncol(mut_matrix2)
  res_matrix <- matrix(nrow = n_samples1, ncol = n_samples2)

  # Loop over the columns of both input matrices,
  # to determine the cosine similarities.
  for (s in seq_len(n_samples1))
  {
    signal1 <- mut_matrix1[, s, drop = TRUE]
    cos_sim_vector <- c()
    for (i in seq_len(n_samples2))
    {
      signal2 <- mut_matrix2[, i, drop = TRUE]
      cos_sim_vector[i] <- as.numeric(signal1 %*% signal2 / (sqrt(signal1 %*% signal1) * sqrt(signal2 %*% signal2)))#cos_sim(signal1, signal2)
    }
    res_matrix[s, ] <- cos_sim_vector
  }
  rownames(res_matrix) <- colnames(mut_matrix1)
  colnames(res_matrix) <- colnames(mut_matrix2)

  return(res_matrix)
}

#' Code copied from MutationalPatterns package
#' Signature clustering function
#'
#' Hierarchical clustering of signatures based on cosine similarity
#'
#' @param signatures Matrix with 96 trinucleotides (rows) and any number of
#' signatures (columns)
#' @param method     The agglomeration method to be used for hierarchical
#' clustering. This should be one of "ward.D", "ward.D2", "single", "complete",
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC). Default = "complete".
#' @return hclust object
#'
#' @noRd
cluster_signatures <- function(signatures, method = "complete") {
  # construct cosine similarity matrix
  sim <- cos_sim_matrix(signatures, signatures)
  # transform to distance
  dist <- as.dist(1 - sim)
  # perform hierarchical clustering
  hc_sig_cos <- hclust(dist, method = method)
  return(hc_sig_cos)
}

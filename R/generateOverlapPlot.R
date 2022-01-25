#' Function to generate Overlap plot
#' @description This function generates an overlap plot
#' using the MAF data.
#' @author Mayank Tandon, Ashish Jain
#' @param mymaf The MAF object
#' @param use_silent_mutations Flag to use the silent mutations
#' in the plot
#' @param summarize_by Paramter to summarize the data
#' (Either 'gene' or 'mutation'")
#' @param plotType The type of plot generated ("ribbon"
#' or "heatmap" or both)
#' @param savename The name and path of the output file
#' @param savewidth Width of plot
#' @param saveheight Height of plot
#' @export
#' @return No return value, the plot is saved as a pdf
#'
#' @examples
#'
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' \donttest{generateOverlapPlot(read.maf(maf))}
generateOverlapPlot<-function(mymaf, use_silent_mutations=FALSE,
                            summarize_by="gene",
                            plotType=c("ribbon","heatmap"),
                            savename=NULL,
                            savewidth=8, saveheight=8){

  ### Add checks for the conditions
  mymaf <- ensurer::ensure_that(mymaf,
                              !is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  use_silent_mutations <- ensurer::ensure_that(use_silent_mutations,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the use_silent_mutations flag in correct format.")
  summarize_by <- ensurer::ensure_that(summarize_by,
                                     !is.null(.) && ((. == "gene") || (. == "mutation")),
                                     err_desc = "Please enter the summarize_by in correct format.")
  plotType <- ensurer::ensure_that(plotType,
                                       !is.null(.) && (class(.) == "character"),
                                       err_desc = "Please enter the plotType in correct format.")
  savewidth <- ensurer::ensure_that(savewidth,
                                    !is.null(.) && (class(.) == "numeric"),
                                    err_desc = "Please enter the savewidth in correct format.")
  saveheight <- ensurer::ensure_that(saveheight,
                                    !is.null(.) && (class(.) == "numeric"),
                                    err_desc = "Please enter the saveheight in correct format.")
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

  if (!is.null(savename)) {pdf(savename, width=savewidth, height=saveheight)}
  if ("heatmap" %in% plotType) {
    #library(pheatmap)
    pheatmap::pheatmap(hm_data,cluster_rows = TRUE, cluster_cols = TRUE,
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
  if (!is.null(savename)) {dev.off()}
}

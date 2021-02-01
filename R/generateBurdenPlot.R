## To Supress Note
utils::globalVariables(c(".", "mut_burden", "classification",
                         "Tumor_Sample_Barcode", "hoverlabel"))
#' Function to generate Burden plot
#' @description This function generates burden plot
#' using the MAF data.
#' @author Mayank Tondon, Ashish Jain
#' @param mymaf The MAF object
#' @param plotType Type of plot (Either "Dotplot" or "Barplot")
#' @param mb_covered Total covered bases for mutation count
#' normalization
#' @param save_data_to_file The name and path of the output file
#'
#' @export
#' @return A ggplot object containing the burden plot
#'
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' generateBurdenPlot(read.maf(maf), plotType="Dotplot")
#' generateBurdenPlot(read.maf(maf), plotType="Barplot")
#'
#' @importFrom dplyr summarise group_by %>%
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#'
generateBurdenPlot<-function(mymaf, plotType=NULL, mb_covered=NULL, save_data_to_file=NULL){

  num_var_data <- mymaf@variants.per.sample
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
                                              levels=num_var_data$Tumor_Sample_Barcode[order(num_var_data$mut_burden, decreasing = TRUE)])

  ########################################################
  #### 5. Generate plots for mutation burden

  ## Pick colors
  median_mut_burdens <- num_var_data %>% dplyr::summarise(median=median(mut_burden))

  num_var_data$xlabel <- factor(num_var_data$xlabel,
                                levels=num_var_data$xlabel[order(num_var_data$mut_burden, decreasing = TRUE)])
  num_var_data$hoverlabel <- paste0("Sample: ",num_var_data$Tumor_Sample_Barcode,"\nMutations: ", num_var_data$mut_burden)

  ### Mutation burden stacked with variant classification counts
  ### Works better for smaller cohorts, say < 20
  variant_type_per_sample <- as.data.frame(mymaf@variant.classification.summary)
  var_type.melt <- reshape2::melt(variant_type_per_sample, id.vars="Tumor_Sample_Barcode",variable.name="classification",value.name="mutation_count")
  var_type.melt$mut_burden <- var_type.melt$mutation_count
  if (is.numeric(mb_covered)) {
    var_type.melt$mut_burden <- var_type.melt$mut_burden/mb_covered
  }
  median_mut_burdens <- data.frame(median=median(var_type.melt[var_type.melt$classification== "total","mut_burden"]))

  plotdata <- var_type.melt[var_type.melt$classification != "total",]
  plotdata$Tumor_Sample_Barcode <- factor(as.character(plotdata$Tumor_Sample_Barcode),
                                          levels=variant_type_per_sample$Tumor_Sample_Barcode[order(variant_type_per_sample$total, decreasing = TRUE)])
  plotdata$classification <- gsub("_"," ",plotdata$classification)

  class_means <- plotdata %>% dplyr::group_by(classification) %>% dplyr::summarise(mean=mean(mut_burden))
  plotdata$classification <- factor(as.character(plotdata$classification),
                                    levels=class_means$classification[order(class_means$mean, decreasing = FALSE)])

  my_class_colors <- my_mutation_colors()

  plotdata$hoverlabel <- paste0("Sample: ",plotdata$Tumor_Sample_Barcode,"\nMutations: ", plotdata$mut_burden)

  if (plotType=="Barplot") {
    if (length(unique(plotdata$Tumor_Sample_Barcode)) <= 20) {
      xaxis_text <- element_text(angle=30, hjust=1)
    } else {
      xaxis_text <- element_blank()
    }
    burden_plot <- ggplot2::ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y=mut_burden, text=hoverlabel)) +
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

    # require(ggbeeswarm)
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
      outdata <- as.data.frame(mymaf@variant.classification.summary)
      outdata$total_per_mb <- outdata$total/mb_covered
      outdata$mb_covered <- mb_covered
      print(paste0("Saving plot data to ", save_data_to_file))
      write.table(outdata, file = save_data_to_file, sep="\t", quote=FALSE,row.names = FALSE)
    } else {
      warning("Path for output data file not found! Skipping...")
    }
  }
  return(burden_plot)
}

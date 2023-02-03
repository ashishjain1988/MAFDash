## To Supress Note
utils::globalVariables(c(".", "boxplot.stats", ".N",
                         "V2", "TCGA","Cohort","Pval",
                         "Median_Mutations","Sample",
                         "total_perMB","plot_total","site",
                         "cohort"))
#' Plot the comparison of the mutation load against TCGA cohorts
#' @description This function computes and plot the mutation load of the input MAF against all of 33 TCGA cohorts derived from MC3 project.
#' @author Ashish Jain, Mayank Tandon
#' @param maf A MAF object
#' @param capture_size capture size for input MAF in MBs. Default NULL. If provided plot will be scaled to mutations per MB
#' @param tcga_capture_size capture size for TCGA cohort in MB. Default 35.8.
#' @param cohortName name for the input MAF cohort. Default "Input"
#' @param tcga_cohorts restrict tcga data to these cohorts.
#' @param primarySite If TRUE uses primary site of cancer as labels instead of TCGA project IDs. Default FALSE.
#' @param col color vector for length 2 TCGA cohorts and input MAF cohort. Default gray70 and black.
#' @param medianCol color for median line. Default red.
#' @param logscale Default TRUE
#' @param decreasing Default FALSE. Cohorts are arranged in increasing mutation burden.
#' @param rm_hyper Remove hyper mutated samples (outliers)? Default FALSE
#' @param rm_zero Remove samples with zero mutations? Default TRUE
#' @return A list consisting of median mutations per cohort, mutation burden per cohort,
#' significant differences between cohorts, and ggplot object to show mutation burden plot
#' @examples
#' library(maftools)
#' library(MAFDash)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' maf <- read.maf(maf = maf)
#' l<-generateTCGAComparePlot(maf = maf, cohortName = "test")
#' l$tcga_compare_plot
#' @export
#' @details Tumor mutation burden for TCGA cohorts is obtained from TCGA MC3 study. For consistency TMB is estimated by restricting variants within Agilent Sureselect capture kit of size 35.8 MB.
#' @source TCGA MC3 file was obtained from  https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc. See TCGAmutations R package for more details. Further downstream script to estimate TMB for each sample can be found in ‘inst/scripts/estimate_tcga_tmb.R’
#' @references Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines Kyle Ellrott, Matthew H. Bailey, Gordon Saksena, et. al. Cell Syst. 2018 Mar 28; 6(3): 271–281.e7. https://doi.org/10.1016/j.cels.2018.03.002
#' @importFrom grDevices boxplot.stats
#' @importFrom stats pairwise.t.test
#' @importFrom data.table fread setDT rbindlist data.table

generateTCGAComparePlot <- function(maf, capture_size = NULL, tcga_capture_size = 35.8, cohortName = NULL, tcga_cohorts = NULL, primarySite = FALSE, col = c('gray70', 'black'), medianCol = 'red', decreasing = FALSE, logscale = TRUE, rm_hyper = FALSE, rm_zero = TRUE){

  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                              !is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  capture_size <- ensurer::ensure_that(capture_size,
                                             is.null(.) || (class(.) == "numeric"),
                                             err_desc = "Please enter the input sample capture size in correct format.")
  tcga_capture_size <- ensurer::ensure_that(tcga_capture_size,
                                       is.null(.) || (class(.) == "numeric"),
                                       err_desc = "Please enter the tcga_capture_size in correct format.")
  cohortName <- ensurer::ensure_that(cohortName,
                                            is.null(.) || (class(.) == "character"),
                                            err_desc = "Please enter correct cohort name")


  tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
  tcga.cohort = data.table::fread(file = tcga.cohort, sep = '\t', stringsAsFactors = FALSE)

  if(primarySite){
    tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, site)]
    colnames(tcga.cohort)[3] = 'cohort'
  }else{
    tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, cohort)]
  }

  if(!is.null(tcga_cohorts)){
    tcga.cohort = tcga.cohort[cohort %in% tcga_cohorts]
    if(nrow(tcga.cohort) == 0){
      stop("Something went wrong. Provide correct names for 'tcga_cohorts' arguments")
    }
  }

  if(length(maf) == 1){
    maf = list(maf)
  }
  maf.mutload = lapply(maf, function(m){
    x = getSampleSummary(m)[,.(Tumor_Sample_Barcode, total)]
    if(rm_zero){
      if(nrow(x[x$total == 0]) > 0)
      {
        warning(paste0("Removed ", nrow(x[x$total == 0]), " samples with zero mutations."))
      }
      x = x[!total == 0]
    }
    x
  })


  if(is.null(cohortName)){
    cohortName = paste0('Input', seq_len(length(maf)))
  }else if(length(cohortName) != length(maf)){
    stop("Please provide names for all input cohorts")
  }

  names(maf.mutload) = cohortName
  maf.mutload = data.table::rbindlist(l = maf.mutload, idcol = "cohort")

  #maf.mutload[,cohort := cohortName]
  tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
  maf.mutload$total = as.numeric(as.character(maf.mutload$total))

  if(rm_hyper){
    #Remove outliers from tcga cohort
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    tcga.cohort = lapply(tcga.cohort, function(x){
      xout = boxplot.stats(x = x$total)$out
      if(length(xout) > 0){
        message(paste0("Removed ", length(xout), " outliers from ", x[1,"cohort"]))
        x = x[!total %in% xout]
      }
      x
    })
    tcga.cohort = data.table::rbindlist(l = tcga.cohort)
    #Remove outliers from input cohort
    xout = boxplot.stats(x = maf.mutload$total)$out
    if(length(xout) > 0){
      message(paste0("Removed ", length(xout), " outliers from Input MAF"))
      maf.mutload = maf.mutload[!total %in% xout]
    }
  }

  if(is.null(capture_size)){
    tcga.cohort = rbind(tcga.cohort, maf.mutload)
    pt.test = stats::pairwise.t.test(x = tcga.cohort$total, g = tcga.cohort$cohort, p.adjust.method = "fdr")
    pt.test.pval = as.data.frame(pt.test$p.value)
    data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
    colnames(pt.test.pval)[1] = 'Cohort'
    #message("Performing pairwise t-test for differences in mutation burden..")
    pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
    colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
    tcga.cohort$plot_total = tcga.cohort$total
    plotYAxisText <- "TMB"
  }else{
    message(paste0("Capture size [TCGA]:  ", tcga_capture_size))
    message(paste0("Capture size [Input]: ", capture_size))
    maf.mutload[,total_perMB := total/capture_size]
    tcga.cohort[,total_perMB := total/tcga_capture_size]
    tcga.cohort = rbind(tcga.cohort, maf.mutload)
    #message("Performing pairwise t-test for differences in mutation burden (per MB)..")
    pt.test = pairwise.t.test(x = tcga.cohort$total_perMB, g = tcga.cohort$cohort, p.adjust.method = "fdr")
    pt.test.pval = as.data.frame(pt.test$p.value)
    data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
    colnames(pt.test.pval)[1] = 'Cohort'
    pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
    colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
    #tcga.cohort[, plot_total := total_perMB]
    tcga.cohort$plot_total = tcga.cohort$total_perMB
    plotYAxisText <- "TMB (per MB)"
  }

  #Median mutations
  tcga.cohort.med = tcga.cohort[,.(.N, median(plot_total)),cohort][order(V2, decreasing = decreasing)]
  tcga.cohort$cohort = factor(x = tcga.cohort$cohort,levels = tcga.cohort.med$cohort)
  colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
  tcga.cohort$TCGA = ifelse(test = tcga.cohort$cohort %in% cohortName,
                            yes = 'Input', no = 'TCGA')
  tcga.cohort.med$Median_Mutations_log10 = log10(tcga.cohort.med$Median_Mutations)

  #Plotting data
  tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
  plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
    x = tcga.cohort[[i]]
    x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                               x[order(plot_total, decreasing = TRUE), plot_total],
                               x[,TCGA])
    x
  })
  names(plot.dat) = names(tcga.cohort)

  # if(logscale){
  #   y_lims = range(log10(data.table::rbindlist(l = plot.dat)[V2 != 0][,V2]))
  # }else{
  #   y_lims = range(data.table::rbindlist(l = plot.dat)[,V2])
  # }
  #
  # #y_lims = range(log10(unlist(lapply(plot.dat, function(x) max(x[,V2], na.rm = TRUE)))))
  # y_min = floor(min(y_lims))
  # y_max = ceiling(max(y_lims))
  # y_lims = c(y_min, y_max)
  # y_at = pretty(y_lims)

  ###New code for ggplot ones
  temp<-lapply(names(plot.dat),function(x){d=plot.dat[[x]];d$Cohort=x;d$Sample=c(nrow(d):1);return(d)})
  t<-do.call(rbind,temp)
  t$color <- col[1]
  t$color[t$Cohort==cohortName]<-col[2]
  t$Cohort<-factor(t$Cohort,levels = tcga.cohort.med$Cohort)
  tcga.cohort.med$Cohort <- factor(tcga.cohort.med$Cohort,levels = tcga.cohort.med$Cohort)
  tcga.cohort.med$color<-col[1]
  tcga.cohort.med$color[tcga.cohort.med$Cohort==cohortName]<-col[2]
  tcgaComparePlot<-ggplot(t,aes(x=Sample, y=V2,colour = Cohort)) +
    geom_point(size=1) +
    scale_colour_manual(values=setNames(tcga.cohort.med$color,as.character(tcga.cohort.med$Cohort))) +
    guides(color="none") +
    labs(x="",y="Tumor Mutation Burden") +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
    geom_hline(data = tcga.cohort.med, aes(yintercept = Median_Mutations),color=medianCol) +
    #geom_rect(data = subset(tp,day == 'Fri'),aes(fill = day),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
    facet_wrap(~Cohort,nrow=1,scales="free_x", strip.position = "bottom") +
    theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank())
  if(logscale)
  {
    tcgaComparePlot<-tcgaComparePlot+scale_y_continuous(trans='log10')
  }

  tcga.cohort = data.table::rbindlist(l = tcga.cohort)
  tcga.cohort[, plot_total := NULL]
  tcga.cohort[, TCGA := NULL]
  pt.test.pval = pt.test.pval[!is.na(Pval)][order(Pval, decreasing = FALSE)]
  return(list(median_mutation_burden = tcga.cohort.med, mutation_burden_perSample = tcga.cohort, pairwise_t_test = pt.test.pval,tcga_compare_plot=tcgaComparePlot))
}

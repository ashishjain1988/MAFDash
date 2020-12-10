#' Function to filter the mutations
#' @description This function filter the mutations in the MAF format using thresholds on various features
#' @author Mayank Tondon, Ashish Jain
#' @param mafFilePath The path of the file containing the mutation
#' information in the MAF format
#' @param flag_genes
#' @param save_name
#' @param no_filter
#' @param norm_alt_max
#' @param t_alt_min
#' @param t_depth_min
#' @param tumor_freq_min
#' @param norm_freq_max
#' @param gnomAD_AF_max
#' @param AF_max
#' @param ExAC_AF_max
#' @param n_callers
#' @param variant_caller
#'
#' @export
#' @return The filtered MAF object
#'
#' @examples
#' library(MAFDashRPackage)
#' #MAFfilePath <- system.file('extdata', 'test.maf', package = 'MAFDashRPackage')
#' #filteredMAF <- filterMAF(mafFilePath = MAFfilePath)

filterMAF<-function(mafFilePath, flag_genes="default",save_name=NULL,no_filter=F,
                    norm_alt_max=1,t_alt_min=1,t_depth_min=20,
                    tumor_freq_min=0.05, norm_freq_max=0.02,
                    gnomAD_AF_max=0.001, AF_max=0.001, ExAC_AF_max=0.001,
                    n_callers=2, variant_caller=NULL){

  # browser()
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  maf_df.raw <- read.table(mafFilePath, sep="\t", header=T, fill = T, quote="\"", stringsAsFactors = F)
  maf_df.raw <- maf_df.raw[maf_df.raw$Hugo_Symbol != "Hugo_Symbol",]
  filter_genes=!maf_df.raw$Hugo_Symbol %in% flag_genes
  maf_df.raw <- maf_df.raw[filter_genes,]


  if (!"tumor_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  }
  if (!"norm_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  }

  filter_tumor_depth=rep(TRUE,nrow(maf_df.raw))
  filter_tumor_alt=rep(TRUE,nrow(maf_df.raw))
  filter_norm_alt=rep(TRUE,nrow(maf_df.raw))
  filter_pop_freq=rep(TRUE,nrow(maf_df.raw))

  if (!no_filter) {
    options(warn=-1)
    filter_tumor_depth=as.numeric(maf_df.raw$t_depth) > t_depth_min
    if (!sum(is.na(maf_df.raw$norm_freq)) == nrow(maf_df.raw)){
      filter_norm_alt=maf_df.raw$norm_freq < norm_freq_max
    }
    filter_tumor_alt=maf_df.raw$tumor_freq > tumor_freq_min
    if (! is.null(t_alt_min)){
      filter_tumor_alt <- filter_tumor_alt & maf_df.raw$t_alt_count > t_alt_min
    }

    filter_gnomad=rep(TRUE,nrow(maf_df.raw))
    filter_1000G=rep(TRUE,nrow(maf_df.raw))
    filter_exac=rep(TRUE,nrow(maf_df.raw))
    if (!is.null(maf_df.raw$gnomAD_AF)) {
      filter_gnomad=maf_df.raw$gnomAD_AF %in% c("-","") | is.na(maf_df.raw$gnomAD_AF) | as.numeric(maf_df.raw$gnomAD_AF) < min(gnomAD_AF_max,1)
    }
    if (!is.null(maf_df.raw$AF)) {
      filter_1000G=maf_df.raw$AF %in% c("-","") | is.na(maf_df.raw$AF)  | as.numeric(maf_df.raw$AF) < min(AF_max,1)
    }
    if (!is.null(maf_df.raw$ExAC_AF)) {
      filter_exac=maf_df.raw$ExAC_AF %in% c("-","") | is.na(maf_df.raw$ExAC_AF) | as.numeric(maf_df.raw$ExAC_AF) < min(ExAC_AF_max,1)
    }
    filter_pop_freq=filter_gnomad & filter_1000G & filter_exac
    options(warn=0)
  }
  filter_caller=rep(TRUE,nrow(maf_df.raw))
  if (! is.null(variant_caller)) {       ### Set 'variant_caller' to NULL to skip any filtering based on caller
    maf_df.raw$set[maf_df.raw$set=="" & maf_df.raw$Hugo_Symbol=="Hugo_Symbol"] <- "set"
    maf_df.raw$set[maf_df.raw$set==""] <- "N.A."
    if (variant_caller == "consensus") {   ### Set 'variant_caller' to 'consensus' to keep variants by two or more callers
      # filter_caller <- grepl("-|Intersection", maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {length(x)>=n_callers | "Intersection" %in% x}))
    } else {                             ### Set 'variant_caller' to one of the callers (mutect, mutect2, vardict, or strelka) to get only that caller
      # filter_caller <- grepl(paste0(variant_caller,"[|-]|Intersection"), maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {any(c(variant_caller,"Intersection") %in% x)}))
    }
  }

  maf_df.rawest <- maf_df.raw
  maf_df.raw <- maf_df.raw[filter_tumor_depth & filter_norm_alt & filter_tumor_alt & filter_pop_freq & filter_caller,]
  # browser()
  maf_df.raw <- maf_df.raw[rowSums(is.na(maf_df.raw))!=ncol(maf_df.raw),]

  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    write.table(maf_df.raw, sep="\t", quote=F, file = save_name, row.names = F)
    print(paste0("Saving filtered maf to ",save_name))
    return(save_name)
  } else {
    return(maf_df.raw)
  }
}

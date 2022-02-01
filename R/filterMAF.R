## To Supress Note
utils::globalVariables(c("%>%",".","suppressWarnings"))

#' Function to filter the mutations
#' @description This function filter the mutations in the MAF format using thresholds on various features
#' @author Mayank Tandon, Ashish Jain
#' @param mafFilePath The path of the file containing the mutation
#' information in the MAF format
#' @param flag_genes The list of genes used as flag genes
#' @param save_name The name and path of the output file to save filtered MAFs
#' @param no_filter Flag to filter the MAF (Default no_filter=FALSE)
#' @param norm_alt_max Alt norm max (Default norm_alt_max=1)
#' @param t_alt_min Alt t min (Default t_alt_min=1)
#' @param t_depth_min Depth t min (Default t_depth_min=20)
#' @param tumor_freq_min Tumor Frequency Minimum (Default tumor_freq_min=0.05)
#' @param norm_freq_max norm_freq_max (Default norm_freq_max=0.02)
#' @param gnomAD_AF_max Maximum allele frequency in gnomAD database  (Default gnomAD_AF_max=0.001)
#' @param AF_max Maximum allele frequency in 1000 genome database (Default AF_max=0.001)
#' @param ExAC_AF_max Maximum allele frequency in ExAC database (Default ExAC_AF_max=0.01)
#' @param n_callers Minimum number of callers identified mutation. (Default n_callers=2)
#' @param variant_caller Name of variant caller to be used or "consensus"
#' to apply filter based on `n_callers` (Default variant_caller=NULL)
#'
#' @export
#' @return An object of class MAF with the filtered mutations
#'
#' @examples
#' library(MAFDash)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' filteredMAF <- filterMAF(mafFilePath = maf)
#'
#' @importFrom ensurer ensure_that
#' @importFrom dplyr pull

filterMAF<-function(mafFilePath, flag_genes="default",save_name=NULL,no_filter=FALSE,
                    norm_alt_max=1,t_alt_min=1,t_depth_min=20,
                    tumor_freq_min=0.05, norm_freq_max=0.02,
                    gnomAD_AF_max=0.001, AF_max=0.001, ExAC_AF_max=0.001,
                    n_callers=2, variant_caller=NULL){

  # grep_vcf_filter_col="PASS",
  # non_silent_only=F,
  # t_alt_max=1e12,
  # t_depth_max=1e12,
  # tumor_freq_max=1,
  # n_alt_min=0,
  # n_depth_min=0,
  # n_depth_max=1e12,
  # norm_freq_min=0,
  # gnomAD_AF_min=0,
  # AF_min=0,
  # ExAC_AF_min=0,

  ### Add checks for the conditions
  mafFilePath <- ensurer::ensure_that(mafFilePath,
                                     !is.null(.) && file.exists(.),
                                     err_desc = "Please enter correct file path.")
  flag_genes <- ensurer::ensure_that(flag_genes,
                                      is.null(.) || (class(.) == "character"),
                                      err_desc = "Please enter the gene list in correct format.")
  no_filter <- ensurer::ensure_that(no_filter,
                                     !is.null(.) && (class(.) == "logical"),
                                     err_desc = "Please enter the no_filter flag in correct format.")
  norm_alt_max <- ensurer::ensure_that(norm_alt_max,
                                    !is.null(.) && (class(.) == "numeric"),
                                    err_desc = "Please enter the norm_alt_max in correct format.")
  t_alt_min <- ensurer::ensure_that(t_alt_min,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the t_alt_min in correct format.")
  t_depth_min <- ensurer::ensure_that(t_depth_min,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the t_depth_min in correct format.")
  tumor_freq_min <- ensurer::ensure_that(tumor_freq_min,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the tumor_freq_min in correct format.")
  norm_freq_max <- ensurer::ensure_that(norm_freq_max,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the norm_freq_max in correct format.")
  gnomAD_AF_max <- ensurer::ensure_that(gnomAD_AF_max,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the gnomAD_AF_max in correct format.")
  AF_max <- ensurer::ensure_that(AF_max,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the AF_max in correct format.")
  ExAC_AF_max <- ensurer::ensure_that(ExAC_AF_max,
                                 !is.null(.) && (class(.) == "numeric"),
                                 err_desc = "Please enter the ExAC_AF_max in correct format.")
  n_callers <- ensurer::ensure_that(n_callers,
                                 !is.null(.) && (class(.) == "numeric"),
                                 err_desc = "Please enter the n_callers in correct format.")


  if (is.null(flag_genes) || length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  maf_df.raw <- read.table(mafFilePath, sep="\t", header=TRUE, fill = TRUE, quote="\"", stringsAsFactors = FALSE)
  #filtered_df <- read_tsv_chunked(maf_file,chunk_size = chunk_lines, col_types = cols(CLIN_SIG = col_character()), callback = DataFrameCallback$new(readr_filterfunc), comment="#")
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
    #options(warn=-1)
    #suppressWarnings()
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
    write.table(maf_df.raw, sep="\t", quote=FALSE, file = save_name, row.names = FALSE)
    print(paste0("Saving filtered maf to ",save_name))
  }
  # } else {
  #   return(maf_df.raw)
  # }
  return(maf_df.raw)
}


#' Function to filter the mutations
#' @description This function reads the the big MAF files in chunks and filter
#' the mutations using various features
#' @author Mayank Tondon, Ashish Jain
#' @param maf_file The path of the file containing the mutation
#' information in the MAF format
#' @param flag_genes The list of genes used as flag genes
#' @param chunk_lines The number of lines to be read at once
#' @param save_name The name and path of the output file to save filtered MAFs
#' @param no_filter Flag to filter the MAF (Default no_filter=FALSE)
#' @param grep_vcf_filter_col FILTER column (Default grep_vcf_filter_col="PASS")
#' @param n_alt_max Alt norm max (Default n_alt_max=1)
#' @param non_silent_only Flag to filter non slient SNVs only
#' (Default non_silent_only=FALSE)
#' @param t_alt_min Alt t min (Default t_alt_min=1)
#' @param t_depth_min Depth t min (Default t_depth_min=20)
#' @param tumor_freq_min Tumor Frequency Minimum (Default tumor_freq_min=0.05)
#' @param norm_freq_max norm_freq_max (Default norm_freq_max=0.02)
#' @param gnomAD_AF_max gnomAD_AF_max (Default gnomAD_AF_max=0.001)
#' @param AF_max 1000 genome data AF_max (Default AF_max=0.001)
#' @param ExAC_AF_max ExAC_AF_max (Default ExAC_AF_max=0.01)
#' @param n_callers n_callers (Default n_callers=2)
#' @param variant_caller variant_caller
#'
#' @export
#' @return An object of class MAF with the filtered mutations
#'
#' @examples
#' library(MAFDash)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' filteredMAF <- filter_maf_chunked(maf_file = maf)
#'
#' @importFrom ensurer ensure_that
#' @importFrom readr read_tsv_chunked
#' @importFrom readr cols
#' @importFrom readr col_character
#' @importFrom readr DataFrameCallback
#' @importFrom tibble as_tibble
filter_maf_chunked <- function(maf_file, chunk_lines=10000, flag_genes="default",save_name=NULL,no_filter=FALSE,grep_vcf_filter_col="PASS",
                               non_silent_only=FALSE,t_alt_min=2,t_depth_min=5,tumor_freq_min=0.01,n_alt_max=1,
                               norm_freq_max=0.01,gnomAD_AF_max=0.001,AF_max=0.001,ExAC_AF_max=0.001,n_callers=2,
                               variant_caller=NULL) {
  ### Add checks for the conditions
  maf_file <- ensurer::ensure_that(maf_file,
                                      !is.null(.) && file.exists(.),
                                      err_desc = "Please enter correct file path.")
  flag_genes <- ensurer::ensure_that(flag_genes,
                                     is.null(.) || (class(.) == "character"),
                                     err_desc = "Please enter the gene list in correct format.")
  save_name <- ensurer::ensure_that(save_name,
                                     is.null(.) || (class(.) == "character"),
                                     err_desc = "Please enter correct filename.")
  no_filter <- ensurer::ensure_that(no_filter,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the no_filter flag in correct format.")
  grep_vcf_filter_col <- ensurer::ensure_that(grep_vcf_filter_col,
                                    (class(.) == "character"),
                                    err_desc = "Please enter correct FILTER column pattern.")
  non_silent_only <- ensurer::ensure_that(non_silent_only,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the non_silent_only flag in correct format.")
  n_alt_max <- ensurer::ensure_that(n_alt_max,
                                       !is.null(.) && (class(.) == "numeric"),
                                       err_desc = "Please enter the norm_alt_max in correct format.")
  t_alt_min <- ensurer::ensure_that(t_alt_min,
                                    !is.null(.) && (class(.) == "numeric"),
                                    err_desc = "Please enter the t_alt_min in correct format.")
  t_depth_min <- ensurer::ensure_that(t_depth_min,
                                      !is.null(.) && (class(.) == "numeric"),
                                      err_desc = "Please enter the t_depth_min in correct format.")
  tumor_freq_min <- ensurer::ensure_that(tumor_freq_min,
                                         !is.null(.) && (class(.) == "numeric"),
                                         err_desc = "Please enter the tumor_freq_min in correct format.")
  norm_freq_max <- ensurer::ensure_that(norm_freq_max,
                                        !is.null(.) && (class(.) == "numeric"),
                                        err_desc = "Please enter the norm_freq_max in correct format.")
  gnomAD_AF_max <- ensurer::ensure_that(gnomAD_AF_max,
                                        !is.null(.) && (class(.) == "numeric"),
                                        err_desc = "Please enter the gnomAD_AF_max in correct format.")
  AF_max <- ensurer::ensure_that(AF_max,
                                 !is.null(.) && (class(.) == "numeric"),
                                 err_desc = "Please enter the AF_max in correct format.")
  ExAC_AF_max <- ensurer::ensure_that(ExAC_AF_max,
                                      !is.null(.) && (class(.) == "numeric"),
                                      err_desc = "Please enter the ExAC_AF_max in correct format.")
  n_callers <- ensurer::ensure_that(n_callers,
                                    !is.null(.) && (class(.) == "numeric"),
                                    err_desc = "Please enter the n_callers in correct format.")
  variant_caller <- ensurer::ensure_that(variant_caller,
                                         is.null(.) || (class(.) == "character"),
                                              err_desc = "Please enter list of variants callers in correct format.")


  filtered_df <- readr::read_tsv_chunked(maf_file,chunk_size = chunk_lines, col_types = readr::cols(CLIN_SIG = readr::col_character()),
                                         callback = readr::DataFrameCallback$new(filter_maf_tbl(flag_genes=flag_genes,no_filter=no_filter,
                                                                                                grep_vcf_filter_col=grep_vcf_filter_col,non_silent_only=non_silent_only,
                                                                                                t_alt_min=t_alt_min,t_depth_min=t_depth_min,tumor_freq_min=tumor_freq_min,
                                                                                                n_alt_max=n_alt_max,norm_freq_max=norm_freq_max,gnomAD_AF_max=gnomAD_AF_max,
                                                                                                AF_max=AF_max,ExAC_AF_max=ExAC_AF_max,n_callers=n_callers,
                                                                                                variant_caller=variant_caller)), comment="#")
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    write.table(filtered_df, sep="\t", quote=F, file = save_name, row.names = F)
    print(paste0("Saving filtered maf to ",save_name))
    return(save_name)
  } else {
    return(filtered_df)
  }
}

filter_maf_tbl <- function(flag_genes="default",
                           no_filter=FALSE,
                           grep_vcf_filter_col="PASS",
                           non_silent_only=FALSE,
                           t_alt_min=2,
                           t_alt_max=1e12,
                           t_depth_min=5,
                           t_depth_max=1e12,
                           tumor_freq_min=0.01,
                           tumor_freq_max=1,
                           n_alt_min=0,
                           n_alt_max=1,
                           n_depth_min=0,
                           n_depth_max=1e12,
                           norm_freq_min=0,
                           norm_freq_max=0.02,
                           gnomAD_AF_min=0,
                           gnomAD_AF_max=0.001,
                           AF_min=0,
                           AF_max=0.001,
                           ExAC_AF_min=0,
                           ExAC_AF_max=0.001,
                           n_callers=2,
                           variant_caller=NULL) {
  function(maftbl,pos){

    if (length(flag_genes)==0) {
      flag_genes <- c()
    } else if (flag_genes[1]=="default") {
      flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
    }
    # browser()
    #require(tibble)
    #require(dplyr)
    df <- tibble::as_tibble(maftbl)
    maf_df.raw <- df[df$Hugo_Symbol != "Hugo_Symbol",]

    if ("FILTER" %in% colnames(maf_df.raw)) {
      if (!is.null(grep_vcf_filter_col)) {
        # maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,pull(maf_df.raw,FILTER)),]
        maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,maf_df.raw$FILTER),]
      }
    } else {
      warning("FILTER column not found; skipping...")
    }

    #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
    vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                     "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                     "In_Frame_Ins", "Missense_Mutation")
    if ("Variant_Classification" %in% colnames(maf_df.raw)) {
      if (non_silent_only) {
        # maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,pull(maf_df.raw,FILTER)),]
        maf_df.raw <- maf_df.raw[maf_df.raw$Variant_Classification %in% vc.nonSilent,]
      }
    } else {
      warning("Variant_Classification column not found; skipping...")
    }

    if (!"tumor_freq" %in% colnames(maf_df.raw)) {
      if (! all(c("t_alt_count","t_depth") %in% colnames(maf_df.raw))) {
        stop("Can't find t_alt_count or t_depth columns")
      }
      maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
    }
    if (!"norm_freq" %in% colnames(maf_df.raw)) {
      if (! all(c("n_alt_count","n_depth") %in% colnames(maf_df.raw))) {
        maf_df.raw$norm_freq <- rep(0, nrow(maf_df.raw))
      } else {
        maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
      }
    }

    maf_num_filter_columns <- list("t_alt_count"=c(min=t_alt_min, max=t_alt_max),
                                   "t_depth"=c(min=t_depth_min, max=t_depth_max),
                                   "tumor_freq"=c(min=tumor_freq_min, max=tumor_freq_max),
                                   "n_alt_count"=c(min=n_alt_min, max=n_alt_max),
                                   "n_depth"=c(min=n_depth_min, max=n_depth_max),
                                   "norm_freq"=c(min=norm_freq_min, max=norm_freq_max),
                                   "gnomAD_AF"=c(min=gnomAD_AF_min, max=gnomAD_AF_max),
                                   "AF"=c(min=AF_min, max=AF_max),
                                   "ExAC_AF"=c(min=ExAC_AF_min, max=ExAC_AF_max)
    )

    numfilter_columns <- names(maf_num_filter_columns)[names(maf_num_filter_columns) %in% colnames(maf_df.raw)]
    notfound <- setdiff(numfilter_columns, names(maf_num_filter_columns))
    if (length(notfound) > 0 ) {
      warning(paste0("Couldn't find these columns; skipping filtering for these: ", paste0(notfound, collapse=",")))
    }

    all_num_filters <- lapply(numfilter_columns, function(col_name) {
      currdata <- as.numeric(dplyr::pull(maf_df.raw,col_name))
      currdata[is.na(currdata)] <- 0
      filter_vec <- currdata >= maf_num_filter_columns[[col_name]]["min"] & currdata <= maf_num_filter_columns[[col_name]]["max"]
      return(filter_vec)
    })
    final_num_filters <- Reduce("&", all_num_filters)

    return_df <- maf_df.raw[final_num_filters,]
    return(return_df)
  }
}

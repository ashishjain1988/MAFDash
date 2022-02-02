## To Supress Note
utils::globalVariables(c(".", ":=", "Tumor_Sample_Barcode",
                         "Hugo_Symbol", "Variant_Classification"))
#' A function to detect MAF genome
#' @description A function to detect MAF genome
#' @author Mayank Tandon, Ashish Jain
#' @param maf The MAF object
#' @export
#' @return A list containing the genome information
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' detectMAFGenome(read.maf(maf))
#'
#' @importFrom grDevices colorRampPalette dev.off pdf rainbow
#' @importFrom graphics strwidth
#' @importFrom stats dist hclust median setNames
#' @importFrom utils count.fields read.table write.table globalVariables
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#' @importFrom IRanges width
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation oncoPrint Legend draw
#' @importFrom pheatmap pheatmap
#' @importFrom plotly plot_ly ggplotly
#' @importFrom dplyr mutate select all_of
#' @importFrom reshape2 dcast

detectMAFGenome<-function(maf){
  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                              !is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  if (! "NCBI_Build" %in% colnames(maf@data)) {
    warning("No genome information in MAF obj.")
    return(NA)
  }

  my_genome = unique(maf@data$NCBI_Build)
  if (length(my_genome) > 1) {
    warning("Multiple genomes listed in MAF obj. Trying the first one")
    my_genome <- my_genome[1]
  }

  return_genome <- switch(my_genome,GRCh38="hg38",GRCh37="hg19",GRCm38="mm10", NA)

  my_chrs <- unique(maf@data$Chromosome)
  add_chr = sum(grepl("^chr", my_chrs)) < length(my_chrs)

  pkg_prefix=ifelse(return_genome=="mm10","BSgenome.Mmusculus.UCSC.","BSgenome.Hsapiens.UCSC.")
  genome_package=paste0(pkg_prefix,return_genome)

  return(list(genome=return_genome, add_chr=add_chr, bsgenome_pkg=genome_package))
}


### Creates matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
#' Creates matrix for oncoplot
#' @description Creates matrix for oncoplot from maf file
#' @author Mayank Tandon, Ashish Jain
#' @param maf The MAF object
#' @param g g
#' @param add_missing add_missing
#' @export
#' @return A list of objects to be used as an input for \link[MAFDash]{generateOncoPlot} function
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' oncoMatrix<-createOncoMatrix(read.maf(maf),g=c("GNA11","MACF1"))
createOncoMatrix = function(maf, g = NULL, add_missing = FALSE){

  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                              !is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  ### Add checks for the conditions
  g <- ensurer::ensure_that(g,
                              !is.null(.) && (class(.) == "character") && (length(.) >= 2),
                              err_desc = "Please provde at least two genes")
  add_missing <- ensurer::ensure_that(add_missing,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the add_missing flag in correct format.")

  subMaf = subsetMaf(maf = maf, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }

  oncomat = reshape2::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                  # xvc = paste0(xvc, collapse="|")
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

#' Make variant table from maf file
#' @description Make variant table from maf file
#' @author Mayank Tandon, Ashish Jain
#' @param maf The MAF object
#' @param use_syn Whether or not to include synonymous variants (default is FALSE, i.e. returns non-synonymous mutations only)
#' @param extra_cols Vector of column names to include from the MAF file.  If it's a named vector, the names will be used in the output table.
#' @export
#' @return A data frame containing the variant information
#' @examples
#' library(MAFDash)
#' library(maftools)
#' maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#' variantTable<-generateVariantTable(read.maf(maf))
generateVariantTable <- function(maf, use_syn=FALSE, extra_cols=c()) {
  ### Add checks for the conditions
  maf <- ensurer::ensure_that(maf,
                              !is.null(.) && (class(.) == "MAF"),
                              err_desc = "Please enter correct MAF object")
  use_syn <- ensurer::ensure_that(use_syn,
                                    !is.null(.) && (class(.) == "logical"),
                                    err_desc = "Please enter the use_syn flag in correct format.")
  # extra_cols <- ensurer::ensure_that(extra_cols,
  #                                 !is.null(.) && (class(.) == "character"),
  #                                 err_desc = "Please provide a character vector to include extra MAF columns.")


  output_data <- maf@data
  if (use_syn) {
    output_data <- rbind(output_data, maf@maf.silent)
  }

  if (all(c("Tumor_Seq_Allele1","Tumor_Seq_Allele2") %in% colnames(output_data))) {
    output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
  }
  if (all(c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2") %in% colnames(output_data))) {
    output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")
  }

  if (! "tumor_freq" %in% colnames(output_data)) {
    # browser()
    if (all(c("t_depth","t_alt_count")%in% colnames(output_data))) {
      output_data$tumor_freq <- as.numeric(as.character(output_data$t_alt_count))/as.numeric(as.character(output_data$t_depth))
    }
  }
  cols_for_table <- c("Hugo Symbol" = "Hugo_Symbol",
                      "Sample ID" = "Tumor_Sample_Barcode",
                      "Variant Classification"="Variant_Classification",
                      "Variant Type"="Variant_Type",
                      "Consequence"="Consequence",
                      "Chromosome"="Chromosome","Start Position" ="Start_Position","End Position"="End_Position","Strand"="Strand",
                      "Reference Allele"="Reference_Allele",
                      "Tumor Genotype"="tumor_genotype",
                      "Normal Genotype"="normal_genotype",
                      "Known Effects ClinVar"="CLIN_SIG",
                      "Transcript Change"="HGVSc",
                      "Protein Change"="HGVSp_Short",
                      "Normal Depth"="n_depth",
                      "Normal Ref Depth"="n_ref_count",
                      "Normal Alt Depth"="n_alt_count",
                      "Tumor Depth"="t_depth",
                      "Tumor Ref Depth"="t_ref_count",
                      "Tumor Alt Depth"="t_alt_count",
                      "Tumor Alt Frequency"="tumor_freq",
                      "Existing Annotation"="Existing_variation",
                      "gnomAD Frequency"="gnomAD_AF",
                      "ExAC Frequency"="ExAC_AF",
                      "1000Genomes Frequency"="AF",
                      "Effect Prediction - SIFT"="SIFT",
                      "Effect Prediction - PolyPhen"="PolyPhen"
  )
  # browser()
  cols_for_table <- c(cols_for_table, extra_cols)
  # variant_info <- as.data.frame(output_data)[,cols_for_table]
  # norm_info_cols <- grep("^n_",cols_for_table, value=T)
  # # mydat <- apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)})
  # if (sum(rowSums(apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)}), na.rm=T), na.rm = T)==0) {
  #   cols_for_table <- cols_for_table[!cols_for_table %in% norm_info_cols]
  # }
  output_cols <- colnames(output_data)[match(cols_for_table, colnames(output_data), nomatch=0)]
  not_output <- cols_for_table[!cols_for_table %in% output_cols]
  if (length(not_output) > 0) {
    message(paste0("Not outputting these columns: ", paste(not_output, collapse=", ")))
  }
  variant_info <- as.data.frame(output_data)[,output_cols]
  colnames(variant_info) <- names(cols_for_table)[match(colnames(variant_info),cols_for_table)]
  return(variant_info)
}

#' Compute exome coverage from a region file
#' @description This function will take a bed file, and return the sum of the lengths of unique regions
#' @author Mayank Tandon, Ashish Jain
#' @param targets_bed_file Path to a bed file with exome target regions
#' @param out_file A file name to which the number of covered bases will be written, instead of returning the value
#' @export
#' @return An integer value of the sum of the length of the covered regions
#' @examples
#' library(MAFDash)
#' bedFile <- system.file("extdata", "test.bed", package = "MAFDash")
#' coverage<-compute_exome_coverage(bedFile)
compute_exome_coverage <- function(targets_bed_file, out_file=NULL) {
  ##### This function will read the target regions BED file and
  #####  compute the sum of the lengths of the regions
  # require(GenomicRanges)

  ## This bit will only read in the first three columns
  num_fields <- max(count.fields(targets_bed_file, sep = "\t"))
  my_classes <- c("character","integer","integer", rep("NULL", num_fields-3))

  ## Read the BED file as a table
  bed_data <- read.table(targets_bed_file, sep="\t",colClasses = my_classes,
                         stringsAsFactors = FALSE)
  colnames(bed_data) <- c("chr","start","end")

  ## Convert to a GenomicRanges object
  bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed_data)

  ## Collapse any overlapping features
  bed.gr.collapse <- GenomicRanges::reduce(bed.gr)

  ## Sum up the widths of each range
  total_exome_coverage = sum(IRanges::width(bed.gr.collapse))

  if (! is.null(out_file)) {
    ## Write to a file
    write.table(total_exome_coverage,file = out_file, col.names =FALSE, row.names = FALSE)
    invisible()
  } else {
    return(total_exome_coverage)
  }
}


#' Function to the extact the Gene symbol from the input genes
#' @description Function to the extact the Gene symbol from the input genes
#' @author Mayank Tandon, Ashish Jain
#' @param genes_arg genes_arg
#' @return A character vector containing the list of genes with gene symbols
#' @noRd
#' @examples
#' library(MAFDash)
#' geneSelectParser()
geneSelectParser <- function(genes_arg=NULL) {

  genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), Reason=c(), stringsAsFactors = FALSE)

  if (! is.null(genes_arg)) {
    if (class(genes_arg)=="character") {
      if (length(genes_arg)==1) {
        ## Then it's either a file name or a single gene
        if (file.exists(genes_arg)) {
          ### Need to parse file type and read accordingly; assuming tsv for now
          gene_data <- read.table(genes_arg,sep="\t", header=TRUE)
          if (sum(c("Hugo_Symbol","Reason") %in% colnames(gene_data)) != 2) {
            stop("Can't find Hugo Symbol or Reason in custom gene input.")
          }
          genes_for_oncoplot <- gene_data
        } else {
          stop(paste0("Can't find file: ",genes_arg))
        }
      } else {
        genes_for_oncoplot <- data.frame(Hugo_Symbol=genes_arg,Reason="Selected Genes", stringsAsFactors = FALSE)
      }
    } else if (class(genes_arg)=="data.frame") {
      genes_for_oncoplot <- genes_arg
      if (! "Reason" %in% colnames(genes_for_oncoplot)) {
        genes_for_oncoplot$Reason <- "Selected Genes"
      }
    } else {
      stop(paste0("Don't know what to do with 'genes_arg' of class: ",class(genes_arg)))
    }

    genes_for_oncoplot <- genes_for_oncoplot[,c("Hugo_Symbol","Reason")]
    genes_for_oncoplot <- data.frame(apply(genes_for_oncoplot,2,as.character), stringsAsFactors = FALSE)
    # genes_for_oncoplot <- genes_for_oncoplot[genes_for_oncoplot$Hugo_Symbol %in% maf.filtered@gene.summary$Hugo_Symbol, ]
  }

  return(genes_for_oncoplot)
}

#' Make the annotation data frame from the TCGA clinical dataset
#' @description This function creates a annotation data frame from the TCGA clinical dataset
#' @author Mayank Tandon, Ashish Jain
#' @param my_clin_dat Clinical dataset in a data frame
#' @param names_to_match The list containing the matched patient's name
#' @param my_colors my_colors
#' @noRd
make_column_annotation <- function(my_clin_dat, names_to_match, my_colors=NULL) {

  myanno <- NULL
  # browser()
  # my_clin_dat <- as.data.frame(my_clin_dat)
  if ("Tumor_Sample_Barcode" %in% colnames(my_clin_dat)) {
    tsb_idx=which(colnames(my_clin_dat)=="Tumor_Sample_Barcode")
  } else {
    ## Try to guess column containing sample IDs
    tsb_idx <- which(apply(my_clin_dat,2, function(x) { sum(names_to_match %in% x)/length(names_to_match)})>0.9)[1]
  }

  if ( is.na(tsb_idx) ) {
    warning("No Tumor Sample Barcode match found")
  } else {
    anno_data <- as.data.frame(my_clin_dat, stringsAsFactors = FALSE)
    # anno_data <- my_clin_dat
    # browser()
    colnames(my_clin_dat)[tsb_idx] <- "Tumor_Sample_Barcode"
    rownames(anno_data) <- anno_data$Tumor_Sample_Barcode
    anno_data <- anno_data[,colnames(anno_data)!="Tumor_Sample_Barcode", drop=FALSE]

    anno_data <- anno_data[match(names_to_match,rownames(anno_data)),,drop=FALSE]
    anno_data <- anno_data[,unlist(lapply(anno_data,function(x){!all(is.na(x))}))]

    if (ncol(anno_data) > 0) {
      if (ncol(anno_data) > 10) {
        warning("Too many columns for annotation, using only the first 10...")
        anno_data <- anno_data[,1:10]
      }
      make_legends <- unlist(lapply(anno_data,function(x) {
        ret_val=TRUE
        if (is.factor(x) && length(levels(x))>12) {
          ret_val=FALSE
        }
        return(ret_val)
      }))
      if (sum(!make_legends) > 0) {
        warning(paste0("Suppressing legends for: ", paste0(names(make_legends)[!make_legends], collapse=", ")))
      }

      # browser()
      if (!is.null(my_colors)) {
        # testcolors <- my_colors[1]
        myanno <- ComplexHeatmap::HeatmapAnnotation(df=anno_data,
                                    which="column",
                                    col = my_colors,
                                    # col = testcolors,
                                    show_legend = make_legends,
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "right")
        # draw(myanno)

      } else {
        myanno <- ComplexHeatmap::HeatmapAnnotation(df=anno_data,
                                    which="column",
                                    show_legend = make_legends,
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "right")
      }
    }

  }

  return(myanno)
}

#' Returns the colors for each mutation
#' @description This function returns the colors for each mutation
#' @author Mayank Tandon, Ashish Jain
#' @noRd
my_mutation_colors <- function() {
  mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
                       In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
                       In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
                       no_variants="#d6d6d6")
  names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  return(mutation_colors)
}

#' Returns the mutation colors for oncoplot function
#' @description This function returns the mutation colors for oncoplot function
#' @author Mayank Tandon, Ashish Jain
#' @noRd
oncoplot_annotation_func <- function() {

  mutation_colors <- my_mutation_colors()
  alter_fun = list(
    background = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = "#CCCCCC", col = NA))
    },
    # "0" = function(x, y, w, h) {
    #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
    #             gp = gpar(fill = "#CCCCCC", col = NA))
    # },
    "Nonsense Mutation" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
    },
    "Missense Mutation" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Missense Mutation"], col = NA))
    },
    "Frame Shift Del" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
    },
    "In Frame Ins" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["In Frame Ins"], col = NA))
    },
    "Splice Site" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Splice Site"], col = NA))
    },
    "Multi Hit" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Multi Hit"], col = NA))
    },
    "Frame Shift Ins" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
    },
    "In Frame Del" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["In Frame Del"], col = NA))
    },
    "Nonstop Mutation" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
    },
    "Translation Start Site" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = grid::gpar(fill = mutation_colors["Translation Start Site"], col = NA))
    },
    "no variants" = function(x, y, w, h) {
      grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                # gp = gpar(fill = "#e0e0e0", col = NA))
                gp = grid::gpar(fill = "#CCCCCC", col = NA))
    }
  )
  return(alter_fun)
}

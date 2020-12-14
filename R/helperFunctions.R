#' A function to detect MAF genome
#' @description A function to detect MAF genome
#' @param maf The MAF object
#' @export
#' @return The list of object containing the genome
#' information
#' @examples
#' library(MAFDashRPackage)
#' #g<-detectMAFGenome(maf)
detectMAFGenome<-function(maf){
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


### Cretaes matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
#' Creates matrix for oncoplot
#' @description Creates matrix for oncoplot from maf file
#' @param maf The MAF object
#' @param g g
#' @param chatty chatty
#' @param add_missing add_missing
#' @export
#' @return The list of objects required for oncoplot function
#' @examples
#' library(MAFDashRPackage)
#' #g<-createOncoMatrix(maf)
createOncoMatrix = function(maf, g = NULL, chatty = TRUE, add_missing = FALSE){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

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

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
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
#' @param maf.filter The MAF object
#' @param use_syn use_syn
#' @param extra_cols extra_cols
#' @export
#' @return Data table containing the variant information
#' @examples
#' library(MAFDashRPackage)
#' #g<-generateVariantTable(maf)
generateVariantTable <- function(maf.filter, use_syn=F, extra_cols=c()) {

  output_data <- maf.filter@data
  if (use_syn) {
    output_data <- rbind(output_data, maf.filter@maf.silent)
  }

  output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
  output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")

  # pheno_info <- sample_info.exome[match(output_data$Tumor_Sample_Barcode, sample_info.exome$Tumor_Sample_Barcode),]
  # pheno_info <- cbind(pheno_info[,"Tumor_Sample_Barcode"],pheno_info[,-c("Tumor_Sample_Barcode")])
  # pheno_columns <- colnames(pheno_info)
  # names(pheno_columns) <- make.names(pheno_columns, unique = T)
  if (! "tumor_freq" %in% colnames(output_data)) {
    # browser()
    if (all(c("t_depth","t_alt_count")%in% colnames(output_data))) {
      output_data$tumor_freq <- as.numeric(as.character(output_data$t_alt_count))/as.numeric(as.character(output_data$t_depth))
    }
  }
  # output_data <- cbind(output_data,pheno_info)
  cols_for_table <- c("Hugo Symbol" = "Hugo_Symbol",
                      "Sample ID" = "Tumor_Sample_Barcode",
                      "Variant Classification"="Variant_Classification",
                      "Variant Type"="Variant_Type",
                      "Consequence"="Consequence",
                      # pheno_columns,
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
  norm_info_cols <- grep("^n_",cols_for_table, value=T)
  # mydat <- apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)})
  if (sum(rowSums(apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)}), na.rm=T), na.rm = T)==0) {
    cols_for_table <- cols_for_table[!cols_for_table %in% norm_info_cols]
  }
  output_cols <- colnames(output_data)[match(cols_for_table, colnames(output_data), nomatch=0)]
  not_output <- cols_for_table[!cols_for_table %in% output_cols]
  if (length(not_output) > 0) {
    print("Not outputting these columsn: ")
    print(not_output)
  }
  variant_info <- as.data.frame(output_data)[,output_cols]
  colnames(variant_info) <- names(cols_for_table)[match(colnames(variant_info),cols_for_table)]
  return(variant_info)
}


### Define colors for mutation types
mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
                     In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
                     In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
                     no_variants="#d6d6d6")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  }
)

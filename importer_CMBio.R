#' This function imports taxonomy and relative abundance data from a cmBio Excel workbook. 
#' Assays are extracted from Excel file sheets and taxonomy information is added to rowData. 
#' 
#' @param file Character. Excel file path. 
#' @param assay.file Named list. Excel file sheet name that has the assay data. 
#' @param row.file Character. Name of the taxonomy sheet in Excel. 
#' @param col.Data NULL
#' @param ID.prefix 
#' @param filter_out Character vector. IDs that need to be filtered out. 
#' 
#' @return A TreeSummarizedExperiment object
#' 
#' 
import_cmBio_taxonomy <- function(file, assay.file=list(), row.file, col.Data=NULL, type=NULL ,ID.prefix=NULL, filter_out=NULL){

  # Check inputs.
  if(!is.character(file) || length(file) != 1){
    stop("'file' must be a single character value.")
  }
  if(!is.character(row.file) || length(row.file)!=1){
    stop("'row.file' must be a single character value.")
  }
  if(!is.null(col.Data) && !(is.data.frame(col.Data) || is.matrix(col.Data) || methods::is(col.Data, "DataFrame"))){
    stop("'col.data' must be data.frame or matrix", call.=FALSE)
  }
  if (!is.list(assay.file) || is.null(names(assay.file))) {
    stop("'assay.file' must be a named list of sheet names.", call. = FALSE)
  }

  #taxonomy table
  taxonomy<- .read_cmbio_taxonomy(file, row.file) # creates taxonomy file
  
  #Assay table
  assay.information<- .read_cmbio_assays(file, assay.file, ID.prefix, filter_out) # Assays
  
  assays <- assay.information$assays 
  taxa_level <- assay.information$taxa_levels_detected
  
  #Check that everything is same size. 
  metadata<- .check_compability(assays, taxonomy, col.Data, taxa_level)
  
  taxonomy <- metadata$row.Data
  col.Data <- metadata$col.Data
  assay_matrices <- metadata$assays
  
  names(assay_matrices) <- names(assay.file)
  #Form TreeSE
  if(is.null(col.Data)){
  tse <- TreeSummarizedExperiment(
      assays = assay_matrices,
      rowData = taxonomy
  )
  }else{
  tse <- TreeSummarizedExperiment(
    assays = assay_matrices,
    rowData = taxonomy,
    colData = col.Data
  )
}
return(tse)
} 
#' Import vOTU profiles into a TreeSE object. 
#' 
import_cmBio_votu <- function(file, assay.file, col.Data=NULL, ID.prefix=NULL, filter_out=NULL){
  
  #Check 
  if(!is.character(file) || length(file) != 1){
    stop("'file' must be a single character value.")
  }
  
  if(!is.null(col.Data) && !(is.data.frame(col.Data) || is.matrix(col.Data) || methods::is(col.Data, "DataFrame"))){
    stop("'col.data' must be data.frame or matrix", call.=FALSE)
  }
  if (!is.list(assay.file) || is.null(names(assay.file))) {
    stop("'assay.file' must be a named list of sheet names.", call. = FALSE)
  }
  
  
  votu_read<- .read_cmbio_vOTU(file, assay.file,ID.prefix, filter_out) 
  
  assays_list<-votu_read$assay
  taxa_profiles <- votu_read$taxa_profiles
  taxa_level <- votu_read$taxa_level
  
  votu_filtered<- .check_compability(assays_list, taxa_profiles, col.Data, taxa_level)
  
  assays<-votu_filtered$assay
  row.Data <- votu_filtered$row.Data
  col.Data <- votu_filtered$col.Data
  
  names(assays) <- names(assay.file)
  
  #Form TreeSE
  if(is.null(col.Data)){
    tse <- TreeSummarizedExperiment(
      assays = assays,
      rowData = row.Data
    )
  }else{
    tse <- TreeSummarizedExperiment(
      assays = assays,
      rowData = row.Data,
      colData = col.Data
    )
  }
  return(tse)
} 
import_cmBio_arg <- function(file, assay.file, col.Data=NULL, ID.prefix=NULL, filter_out=NULL){
  #Check 
  if(!is.character(file) || length(file) != 1){
    stop("'file' must be a single character value.")
  }
  
  if(!is.null(col.Data) && !(is.data.frame(col.Data) || is.matrix(col.Data) || methods::is(col.Data, "DataFrame"))){
    stop("'col.data' must be data.frame or matrix", call.=FALSE)
  }
  if (!is.list(assay.file) || is.null(names(assay.file))) {
    stop("'assay.file' must be a named list of sheet names.", call. = FALSE)
  }
  
  read_arg<- .read_cmbio_ARG(file, assay.file, ID.prefix, filter_out)
  
  assay_list<-read_arg$assays
  row.Data <- read_arg$row.Data
  level<- read_arg$level
  
  
  arg_filtered<- .check_compability(assay_list, row.Data, col.Data, level)
  
  assays<-arg_filtered$assays
  row.Data <- arg_filtered$row.Data
  col.Data <- arg_filtered$col.Data
  
  names(assays) <- names(assay.file)
  
  #Form TreeSE
  if(is.null(col.Data)){
    tse <- TreeSummarizedExperiment(
      assays = assays,
      rowData = row.Data
    )
  }else{
    tse <- TreeSummarizedExperiment(
      assays = assays,
      rowData = row.Data,
      colData = col.Data
    )
  }
  return(tse)
}
#'
#'
#' Helper functions. 
.file_and_sheet_exists <- function(file, sheet){
  
  if (!file.exists(file)) {
    stop("The input file does not exist: ", file, call. = FALSE)
  }
  all_sheets <- readxl::excel_sheets(file)
  if (!sheet %in% all_sheets) {
    stop("Sheet '", sheet, "' not found in Excel file. Available sheets are: ",
         paste(all_sheets, collapse = ", "), call. = FALSE)
  }
  
} # Check that the file and sheet exist. 
#'
.read_cmbio_taxonomy <- function(file, row.file){
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The 'readxl' package is required to read Excel files.", call. = FALSE)
  }
  
  #Check
  .file_and_sheet_exists(file, row.file)
  taxonomies<- readxl::read_excel(file, sheet=row.file)

  # Possible taxonomies
  taxa_columns <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  matched_columns <- sapply(taxa_columns, function(x) {
    match <- grep(paste0("^", x, "$"), names(taxonomies), ignore.case = TRUE, value = TRUE)
    if (length(match) > 0) return(match[1]) else return(NA)
  }) 
  

  if (any(is.na(matched_columns))) {
    
    warning("Some taxa levels are missing: ",
            paste(taxa_columns[is.na(matched_columns)], collapse = ", ")
    )
    
  }
  matched_columns <- matched_columns[!is.na(matched_columns)] 
  
  taxa_profiles <- taxonomies[, matched_columns, drop = FALSE] 
  
  return(taxa_profiles)
} 
.read_cmbio_assays <- function(file, assay.file, ID.prefix = NULL, filter_out = NULL) {
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The 'readxl' package is required to read Excel files.", call. = FALSE)
  }
  

  results <- lapply(assay.file, function(sheet) {
    .file_and_sheet_exists(file, sheet)
    
    
    dat <- readxl::read_excel(file, sheet = sheet)
    dat <- as.data.frame(dat)
    
    taxa <- colnames(dat)[1]

    
    rownames(dat) <- dat[[1]]
    dat <- dat[, -1, drop = FALSE]
    dat <- as.matrix(dat)
    
    if (!is.null(ID.prefix)) {
      colnames(dat) <- sub(paste0("^", ID.prefix), "", colnames(dat))
    }
    
    if (!is.null(filter_out)) {
      dat <- dat[, !colnames(dat) %in% filter_out, drop = FALSE]
    }
    
    return(list(assay = dat, taxa = taxa))
  })
  
  assay_matrices <- lapply(results, `[[`, "assay")
  taxa_levels_detected <- unique(unlist(lapply(results, `[[`, "taxa")))
                                 
  return(list(
    assays = assay_matrices,
    taxa_levels_detected = unique(taxa_levels_detected)
  ))
} 
.check_compability <- function(assays, rowdata, col.Data, taxa_level){
  
  reference_taxa <- rownames(assays[[1]])
  reference_samples <- colnames(assays[[1]])
  
  #Check that all assays have the same taxa and samples. 
  inconsistent <- !sapply(assays, function(mat) identical(rownames(mat), reference_taxa))
  if (any(inconsistent)) {
    stop("Not all assays contain the same taxa: ",
         paste(names(assays)[inconsistent], collapse = ", "))
  }
  
  inconsistent <- !sapply(assays, function(mat) identical(colnames(mat), reference_samples))
  
  if (any(inconsistent)) {
    stop("Not all assays contain the same samples: ",
         paste(names(assays)[inconsistent], collapse = ", "))
  }
  
  #Is assay taxa level in RowData?
  if (!(taxa_level %in% colnames(rowdata))) {
    stop("Taxa level '", taxa_level, "' not found in taxonomy.")
  }
  
  # Keep only taxa present in the assays. 
  taxa_profiles_filtered <- rowdata[
    rowdata[[taxa_level]] %in% reference_taxa, , drop = FALSE
  ]
  
  taxa_profiles_filtered <- taxa_profiles_filtered[
    !duplicated(taxa_profiles_filtered[[taxa_level]]),, drop = FALSE ]#filter down 
  
  rownames(taxa_profiles_filtered) <- taxa_profiles_filtered[[taxa_level]]  #rownames()
  
  taxa_profiles_filtered <- taxa_profiles_filtered[reference_taxa, , drop = FALSE]
  
  if(!is.null(col.Data)){
    
    reference_samples <- colnames(assays[[1]])
    sample_meta       <- rownames(col.Data)
    
    common_samples <- reference_samples[reference_samples %in% sample_meta]
    
    if (length(common_samples) == 0) {
      stop("Yhtään näytettä ei löydy sekä assayssa että col.Data:ssa.")
    }
    
    if (length(common_samples) < length(reference_samples)) {
      warning(
        sprintf(
          "%d/%d näytettä assays-listasta puuttuu col.Data:ssa ja ne pudotetaan.",
          length(reference_samples) - length(common_samples),
          length(reference_samples)
        )
      )
    }
    assays <- lapply(assays, function(mat) {
      mat[, common_samples, drop = FALSE]
    })
    col.Data <- col.Data[common_samples, , drop = FALSE]
  }
  
  
  return(list(
    row.Data=taxa_profiles_filtered ,
    col.Data=col.Data,
    assays=assays
  ))
}
#'
.read_cmbio_ARG <- function(file, assay.file, ID.prefix = NULL, filter_out = NULL){
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The 'readxl' package is required to read Excel files.", call. = FALSE)
  }
  ARG_columns <- c("Measure","ARG_class","ARG","ARG_gene_ID")
  
  results <- lapply(assay.file, function(sheet) {
    
    .file_and_sheet_exists(file, sheet)
    
    dat <- readxl::read_excel(file, sheet = sheet)
    dat <- as.data.frame(dat)  
    
    
    present_taxa <- names(dat)[tolower(names(dat)) %in% tolower(ARG_columns)]
    if (length(present_taxa) == 0) {
      stop("No taxonomy columns found in sheet '", sheet, "'.")
    }
    level <- tail(intersect(ARG_columns, present_taxa), 1) #taxa_level. 

    row.Data <- dat[, present_taxa, drop = FALSE] #taxonomic table
    
    if (is.null(level)) {
      stop("Taxa level column '", level, "' not found in data.")
    }
    rownames(dat)<- dat[[level]]
    
    #Drop taxa columns from assay tables. 
    dat <- dat[, !names(dat) %in% present_taxa, drop = FALSE]
    
    if (!is.null(ID.prefix)) {
      colnames(dat) <- sub(paste0("^", ID.prefix), "", colnames(dat))
    }
    if (!is.null(filter_out)) {
      dat <- dat[, !colnames(dat) %in% filter_out, drop = FALSE]
    }
    dat <- as.matrix(dat)
    list(
      assay = dat,
      row.Data = row.Data,
      level =level
    )
  })
  
  assay_list <- lapply(results, "[[", "assay")
  
  row.Data <- results[[1]]$row.Data
  level <- results[[1]]$level
  
  return(list(
    assays        = assay_list,
    row.Data = row.Data,
    level= level
  ))
}
#'
#'
.read_cmbio_vOTU <- function(file, assay.file, ID.prefix = NULL, filter_out = NULL){
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The 'readxl' package is required to read Excel files.", call. = FALSE)
  }
  
  all_taxa_cols <- c(
    "virus_type", "realm", "superkingdom", "domain", "kingdom", "phylum",
    "class", "order", "family", "genus", "species",
    "lifestyle", "vOTU"
  )
  
  results <- lapply(assay.file, function(sheet) {
    .file_and_sheet_exists(file, sheet)
    dat <- readxl::read_excel(file, sheet = sheet)
    dat <- as.data.frame(dat)
    
    present_taxa <- names(dat)[tolower(names(dat)) %in% tolower(all_taxa_cols)]
    if (length(present_taxa) == 0) {
      stop("No taxonomy columns found in sheet '", sheet, "'.")
    }
    taxa_level <- tail(intersect(all_taxa_cols, present_taxa), 1) #taxa_level. 
    
    taxa_profiles <- dat[, present_taxa, drop = FALSE] #taxonomic table
    
    if (is.null(taxa_level)) {
      stop("Taxa level column '", taxa_level, "' not found in data.")
    }
    rownames(dat)<- dat[[taxa_level]]
    
    #Drop taxa columns from assay tables. 
    dat <- dat[, !names(dat) %in% present_taxa, drop = FALSE]
    
    if (!is.null(ID.prefix)) {
      colnames(dat) <- sub(paste0("^", ID.prefix), "", colnames(dat))
    }
    if (!is.null(filter_out)) {
      dat <- dat[, !colnames(dat) %in% filter_out, drop = FALSE]
    }
    
    dat <- as.matrix(dat)
    
    list(
      assay         = dat,
      taxa_profiles = taxa_profiles,
      taxa_level =taxa_level
    )
  })
  
  assay_list <- lapply(results, "[[", "assay")

  taxa_profile <- results[[1]]$taxa_profiles
  taxa_level <- results[[1]]$taxa_level
  
  return(list(
    assays        = assay_list,
    taxa_profiles = taxa_profile,
    taxa_level= taxa_level
  ))
}


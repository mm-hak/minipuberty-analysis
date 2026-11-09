# Data cleaning functions
clean_continuous_variables <- function(df, variables) {
  improper_values <- c("###", ".", "?", "#", "tarkasta", "-", "ei")
  
  for (i in variables) {
    if (is.character(df[[i]]) || is.factor(df[[i]])) {
      df[[i]] <- as.character(df[[i]])
      
      df[[i]] <- gsub(",", ".", df[[i]])  
      df[[i]] <- gsub("\u2212", "-", df[[i]])  
      df[[i]] <- gsub("\u00A0", " ", df[[i]])  
      df[[i]] <- gsub("[[:space:]]+", " ", df[[i]])
      df[[i]] <- trimws(df[[i]])
      
      
      df[[i]] <- gsub("(?i)(-?\\d+(\\.\\d+)?)\\s*(sd|%)", "\\1", df[[i]], perl = TRUE)
      
      
      df[[i]][df[[i]] %in% improper_values] <- NA
      
      df[[i]] <- as.numeric(df[[i]])
    }
  }
  
  return(df)
}

clean_binary_variables <- function(df, variables) {
  for (var in variables) {
    
    df[[var]] <- as.character(df[[var]])
    
    df[[var]] <- gsub("[[:space:]]+", " ", df[[var]])
    df[[var]] <- trimws(df[[var]])
    
    df[[var]][!df[[var]] %in% c("1", "2")] <- NA
    
    df[[var]] <- factor(df[[var]], levels = c("1", "2"))
  }
  return(df)
}

clean_categoric_variables <- function(df, variables,allowed){
  for(var in variables){
    df[[var]] <- as.character(df[[var]])
    
    df[[var]] <- trimws(df[[var]])
    df[[var]] <- gsub("[[:space:]]+", " ", df[[var]])
    df[[var]] <- tolower(df[[var]]) 
    
    if( var %in% names(allowed)){
      df[[var]][!df[[var]] %in% allowed[[var]]] <- NA
    }
    df[[var]] <- as.factor(df[[var]])
  }
  return(df)
}

# Adds alpha diversity 
process_experiment <- function(x, label,
                               detection = 1e-4,
                               pseudocount = 1e-6) {
  
  if (!"relabundance" %in% assayNames(x)) {
    stop("Assay 'relabundance' puuttuu tasolta: ", label)
  }
  
  # Transformations from relabundance
  x <- transformAssay(
    x,
    method      = "clr",
    assay.type  = "relabundance",
    pseudocount = pseudocount,
    name        = "clr"
  )
  
  x <- transformAssay(
    x,
    method      = "log",
    assay.type  = "relabundance",
    pseudocount = pseudocount,
    name        = "log"
  )
  
  # Alpha relabundancesta (Shannon/Simpson + richness detection-limit)
  x <- addAlpha(
    x,
    assay.type = "relabundance",
    detection  = detection
  )
  
  # Betadiversity (db-RDA)
  ## Bray–Curtis relabundancesta
  x <- runRDA(
    x,
    FUN        = getDissimilarity,
    method     = "bray",
    assay.type = "relabundance",
    name       = paste0("dbrda_bray_", label)
  )
  
  
  ## Aitchison (euclidean in CLR-space)
  x <- runRDA(
    x,
    FUN        = getDissimilarity,
    method     = "euclidean",
    assay.type = "clr",
    name       = paste0("dbrda_aitchison_", label)
  )
  
  x
}

clean_tax_label <- function(x){
  x <- sub("^species:", "", x)     
  x <- sub("^s__","", x)           
  x <- sub("^g__","", x)           
  x <- gsub(";", "_", x)           
  x <- gsub("\\s+", " ", x)        
  trimws(x)
  
  return(x)
}

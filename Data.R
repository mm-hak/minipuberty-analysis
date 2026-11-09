#' This Data.R script imports the data and makes the necessary mutation and additions to it. 

# Run this script with 01_data_processing as the working directory.
# From the project root, use source("01_data_processing/Data.R", chdir = TRUE).

# Set this path to your participant metadata CSV file.
metadata_path <- "/path/to/Minipuberty_colData3.csv"

# Set this path to your cmBio taxonomic profiles Excel workbook containing non-rarefied taxonomic profiles. 
# Required sheets: 001_Taxonomies and 002_species-relative-abundance. 
cmbio_taxonomic_profiles_path <- "/path/to/cmbio_taxonomic_profiles.xlsx"

# Keep this relative path: processed data are saved inside the project.
output_dir <- file.path("..", "data")

required_input_paths <- c(metadata_path, cmbio_taxonomic_profiles_path)
missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]
if (length(missing_input_paths) > 0) {
  stop(
    "Required input file(s) not found: ",
    paste(missing_input_paths, collapse = ", "),
    call. = FALSE
  )
}

# Libraries.
pkgs <- c(
  "tidyverse",        # includes dplyr, readr, tibble, etc.
  "readxl",
  "mia",              # Data wrangling & analysis
  "miaTime",          # longitudinal tools
  "miaViz",           # visualizations
  "S4Vectors",        # Bioconductor infrastructure
  "ape",        
  "TreeSummarizedExperiment",
  "MultiAssayExperiment",
  "SummarizedExperiment"
)

missing_packages <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Required package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running Data.R.",
    call. = FALSE
  )
}

invisible(lapply(pkgs, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE)
  )
}))

# Helper and import functions---------------------------------------------------
source("helper_functions.R", local = FALSE)
source("importer_CMBio.R", local = FALSE) # Import functions
# ------------------------------------------------------------------------------
# Rawdata from csv file that is imported from excel
colnames <- make.names(names(read_csv(metadata_path, n_max = 0, show_col_types = FALSE))) #column names


actual_data <- read_csv(metadata_path,skip=2, col_names = colnames,show_col_types = FALSE, locale = locale(decimal_mark = ","))
# ------------------------------------------------------------------------------

#' The data was missing the growth velocity information on the 0 to 6 month period
#' Adding missing time period value. 
actual_data<- actual_data %>% 
  mutate(
    Growth.velocity.0.6M = (X6M.height-Birth.height)/0.5
  )
#' Creating data that contains information on baby boys and girls. Filtering out the 
#' rows that contain missing values, selecting necessary columns and making the 
#' transformations (log).
#' Also adding early life antibiotic exposure var that tells if the infant has been exposed to antibiotics 
#' during pregnancy, birth or before 7 days of life.

#For the gender comparison in stage 1.
g1_gender <- actual_data %>%
  dplyr::filter(!if_any(c("Birth.height","X3M.height", "X6M.height", "Birth.weight","X3M.weight","X6M.weight"), is.na)) %>%
  column_to_rownames("Original.study") %>% 
  dplyr::mutate(
    Gender=Child.Sex
  ) %>% 
  dplyr::select(Gender,Intervention,Pregnancy.antibiotic,
         Labour.antibiotic,Child.antibiotic.under.7.days,
         Child.antibiotic.7.days.to.3M, Birth.mode, Mothers.age,
         Mother.BMI.before.pregnancy, Weight.gain.during.pregnancy,BMI.after.pregnancy,
         Pregnancy.weeks, Birth.height,Birth.weight, Growth.velocity.0.3M,Growth.velocity.3.6M,
         Growth.velocity.0.6M, Weight.gain.0.3M,Weight.gain.3.6M,Weight.gain.0.6M) %>%
  dplyr::mutate(
    Growth.velocity.0.3M_log=log(Growth.velocity.0.3M),
    Growth.velocity.3.6M_log=log(Growth.velocity.3.6M),
    Growth.velocity.0.6M_log=log(Growth.velocity.0.6M),
    Weight.gain.0.3M_log=log(Weight.gain.0.3M),
    Weight.gain.3.6M_log=log(Weight.gain.3.6M),
    Weight.gain.0.6M_log=log(Weight.gain.0.6M),
    early_life_ab =as.factor(ifelse(Pregnancy.antibiotic==2| Labour.antibiotic==2| Child.antibiotic.under.7.days==2,2,1))
  )%>% 
  as.data.frame()

# Making sure all values are proper in the data. 
Categoric_variables<- c('Intervention','Gender','Birth.mode')
Binary_variables <- c('Pregnancy.antibiotic','Labour.antibiotic','Child.antibiotic.under.7.days','Child.antibiotic.7.days.to.3M')
Continous_variables <- c('Mothers.age','Pregnancy.weeks','Mother.BMI.before.pregnancy','BMI.after.pregnancy','Weight.gain.during.pregnancy', 'Birth.height','Birth.weight','Growth.velocity.0.3M','Growth.velocity.3.6M','Growth.velocity.0.6M','Weight.gain.0.3M','Weight.gain.3.6M', 'Weight.gain.0.6M')
Outcome_variables <- c("Growth.velocity.0.3M","Growth.velocity.3.6M","Growth.velocity.0.6M",
                       "Weight.gain.0.3M","Weight.gain.3.6M","Weight.gain.0.6M")
allowed <- list(
  "Intervention"= c("probiotic","placebo"),
  "Gender"= c("male","female"),
  "Birth.mode"=c("1","2")
)  
g1_gender<- clean_continuous_variables(g1_gender, Continous_variables)
g1_gender <- clean_binary_variables(g1_gender,Binary_variables)
g1_gender <- clean_categoric_variables(g1_gender, Categoric_variables, allowed) 

# Check
outcome_zero <- g1_gender %>%
  filter(if_any(any_of(Outcome_variables), ~ .x == 0)) 

g1_gender <- na.omit(g1_gender) # fullrows. 
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(g1_gender, file = file.path(output_dir, "g1_gender.Rds"))

# Create the blood-sample dataset used in Stage 2. 
# In the source data, blood samples are available only for male infants.
data_reduced <- actual_data %>%
  column_to_rownames("Original.study") %>%
  filter(!if_any(c("Birth.height", "X3M.height", "X6M.height"), is.na)) %>%
  filter(!if_any(c("Birth.weight", "X3M.weight", "X6M.weight"), is.na)) %>%
  filter(!is.na(Blood.Sample)) %>%
  dplyr::select(
    Intervention, Pregnancy.antibiotic, Labour.antibiotic,
    Child.antibiotic.under.7.days, Child.antibiotic.7.days.to.3M,
    Birth.mode, Mothers.age, Mother.BMI.before.pregnancy,
    Weight.gain.during.pregnancy, BMI.after.pregnancy, Pregnancy.weeks,
    Birth.height, Birth.weight, S.TestoMS, P.LH, P.FSH, P.SHBG,
    Growth.velocity.0.3M, Growth.velocity.3.6M, Growth.velocity.0.6M,
    Weight.gain.0.3M, Weight.gain.3.6M, Weight.gain.0.6M
  ) %>%
  as.data.frame()

# Clean the blood-sample dataset.
data_reduced_categorical_variables <- c("Intervention", "Birth.mode")
data_reduced_binary_variables <- c(
  "Pregnancy.antibiotic", "Labour.antibiotic",
  "Child.antibiotic.under.7.days", "Child.antibiotic.7.days.to.3M"
)
data_reduced_continuous_variables <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy",
  "Birth.height", "Birth.weight", "Growth.velocity.0.3M",
  "Growth.velocity.3.6M", "Growth.velocity.0.6M", "Weight.gain.0.3M",
  "Weight.gain.3.6M", "Weight.gain.0.6M", "S.TestoMS", "P.LH",
  "P.FSH", "P.SHBG"
)
data_reduced_allowed <- list(
  "Intervention" = c("probiotic", "placebo"),
  "Birth.mode" = c("1", "2")
)

data_reduced <- clean_continuous_variables(
  data_reduced, data_reduced_continuous_variables
)
data_reduced <- clean_binary_variables(
  data_reduced, data_reduced_binary_variables
)
data_reduced <- clean_categoric_variables(
  data_reduced,
  data_reduced_categorical_variables,
  data_reduced_allowed
)
data_reduced <- na.omit(data_reduced)
data_reduced <- data_reduced %>%
  dplyr::filter(
    !if_any(
      all_of(c("S.TestoMS", "P.LH", "P.FSH", "P.SHBG")),
      ~ .x == 0
    )
  ) %>%
  dplyr::mutate(
    S.TestoMS_log = log(S.TestoMS),
    P.LH_log = log(P.LH),
    P.FSH_log = log(P.FSH),
    P.SHBG_log = log(P.SHBG),
    Growth.velocity.0.3M_log = log(Growth.velocity.0.3M),
    Growth.velocity.3.6M_log = log(Growth.velocity.3.6M),
    Growth.velocity.0.6M_log = log(Growth.velocity.0.6M),
    Weight.gain.0.3M_log = log(Weight.gain.0.3M),
    Weight.gain.3.6M_log = log(Weight.gain.3.6M),
    Weight.gain.0.6M_log = log(Weight.gain.0.6M),
    early_life_ab = as.factor(
      ifelse(
        Pregnancy.antibiotic == 2 |
          Labour.antibiotic == 2 |
          Child.antibiotic.under.7.days == 2,
        2,
        1
      )
    )
  )
saveRDS(data_reduced, file = file.path(output_dir, "data_reduced.Rds"))

# Prepare metadata for the microbiome analysis using infants with both
# blood and stool samples. Blood samples were available only for male infants.
sample_metadata<- actual_data %>% 
  column_to_rownames("Original.study") %>%
  filter(!if_any(c("Birth.height","X3M.height", "X6M.height"), is.na)) %>% 
  filter(!if_any(c("Birth.weight","X3M.weight","X6M.weight"), is.na)) %>% 
  filter(!is.na(Blood.Sample)) %>% 
  filter(!is.na(Fecal.Sample) & !grepl("low", Faecal.coding, ignore.case = TRUE)) %>% 
  dplyr::select(Intervention,Pregnancy.antibiotic,
         Labour.antibiotic,Child.antibiotic.under.7.days, 
         Child.antibiotic.7.days.to.3M, Birth.mode, Mothers.age, 
         Mother.BMI.before.pregnancy, Weight.gain.during.pregnancy,BMI.after.pregnancy,
         Pregnancy.weeks, Birth.height,Birth.weight,S.TestoMS,P.LH,P.FSH, P.SHBG,
         Growth.velocity.0.3M,Growth.velocity.3.6M,Growth.velocity.0.6M,
         Weight.gain.0.3M,Weight.gain.3.6M, Weight.gain.0.6M) %>% 
  as.data.frame()



# Clean the sample metadata. 
Categoric_variables<- c('Intervention','Birth.mode')
Binary_variables <- c('Pregnancy.antibiotic','Labour.antibiotic','Child.antibiotic.under.7.days','Child.antibiotic.7.days.to.3M')
Continous_variables <- c('Mothers.age','Pregnancy.weeks','Mother.BMI.before.pregnancy','BMI.after.pregnancy','Weight.gain.during.pregnancy', 'Birth.height','Birth.weight','Growth.velocity.0.3M','Growth.velocity.3.6M','Growth.velocity.0.6M','Weight.gain.0.3M','Weight.gain.3.6M', 'Weight.gain.0.6M','S.TestoMS','P.LH','P.FSH','P.SHBG')
allowed <- list(
  "Intervention"= c("probiotic","placebo"),
  "Birth.mode"=c("1","2") #1=Vb, 2=Cs
)

# Cleaning the data types correct and checking to see that . and , are marked correct. 
sample_metadata<- clean_continuous_variables(sample_metadata, Continous_variables)
sample_metadata <- clean_binary_variables(sample_metadata,Binary_variables)
sample_metadata <- clean_categoric_variables(sample_metadata, Categoric_variables, allowed)
sample_metadata <- na.omit(sample_metadata)
sample_metadata <- sample_metadata %>%
  dplyr::filter(!if_any(all_of(c("S.TestoMS", "P.LH", "P.FSH", "P.SHBG")), ~ .x == 0)) %>% 
  dplyr::mutate(
    S.TestoMS_log=log(S.TestoMS),
    P.LH_log  = log(P.LH),
    P.FSH_log = log(P.FSH),
    P.SHBG_log= log(P.SHBG),
    Growth.velocity.0.3M_log=log(Growth.velocity.0.3M),
    Growth.velocity.3.6M_log=log(Growth.velocity.3.6M),
    Growth.velocity.0.6M_log=log(Growth.velocity.0.6M),
    Weight.gain.0.3M_log=log(Weight.gain.0.3M),
    Weight.gain.3.6M_log=log(Weight.gain.3.6M),
    Weight.gain.0.6M_log=log(Weight.gain.0.6M), 
    early_life_ab =as.factor(ifelse(Pregnancy.antibiotic==2| Labour.antibiotic==2| Child.antibiotic.under.7.days==2,2,1))
  )
sample_metadata<- sample_metadata %>%S4Vectors::DataFrame()

# Collect the IDs of stool samples marked as "low" for exclusion during import.

low_faecal_coding_samples <- actual_data %>%
  filter(grepl("low", Faecal.coding, ignore.case = TRUE)) %>%
   pull(Original.study)


# Create treeSE. 
path <- cmbio_taxonomic_profiles_path

assay_sheets <- list(
  "relabundance" = "002_species-relative-abundance"
) 

# Set this to the sample ID prefix used in your cmBio workbook.
cmbio_id_prefix <- "your_project_prefix_"

tse <- import_cmBio_taxonomy(
  file = path,
  assay.file = assay_sheets,
  row.file = "001_Taxonomies",
  col.Data = sample_metadata,
  ID.prefix = cmbio_id_prefix, # Remove the prefix from assay sample IDs.
  filter_out = low_faecal_coding_samples
)

# Scale from 0-1. 
assay(tse, "relabundance") <- assay(tse, "relabundance") / 100


taxonomyRanks(tse)
tse<- addHierarchyTree(tse)


phy  <- rowTree(tse, "phylo")
phy$tip.label <- clean_tax_label(phy$tip.label)
rowTree(tse, "phylo") <- phy

tse<- agglomerateByRanks(tse)

names<- altExpNames(tse)

alpha_detection <- 1e-4
pseudocount     <- 1e-6

tse <- process_experiment(
  tse,
  label       = "rela",
  detection   = alpha_detection,
  pseudocount = pseudocount
)

for (rk in altExpNames(tse)) {
  ax <- altExp(tse, rk)
  
  ax <- process_experiment(
    ax,
    label       = rk,
    detection   = alpha_detection,
    pseudocount = pseudocount
  )
  
  altExp(tse, rk) <- ax
}

# save
saveRDS(tse, file = file.path(output_dir, "tse.Rds"))

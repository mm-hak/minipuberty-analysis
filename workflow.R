# Control document

# Set this path to the downloaded project root (the folder containing workflow.R).
setwd("/path/to/folder")

# Data processing.

# Creates the datasets needed for the analyses. 
source("01_data_processing/Data.R", encoding = "UTF-8", chdir = TRUE)

# Data description
source("02_analysis/00_data_description.R", encoding = "UTF-8", chdir = TRUE)

# Growth and hormone analyses 
source("02_analysis/01_stage1_gender.R", encoding = "UTF-8", chdir = TRUE)
source("02_analysis/02_stage1.R", encoding = "UTF-8", chdir = TRUE)
source("02_analysis/03_stage2.R", encoding = "UTF-8", chdir = TRUE)

# Differential-abundance analyses 

# Might take some time to run...
source(
  "02_analysis/04_differential_abundance_relative.R",
  encoding = "UTF-8",
  chdir = TRUE
)
source(
  "02_analysis/05_differential_abundance_clr.R",
  encoding = "UTF-8",
  chdir = TRUE
)

# Stage 3 microbiome analyses

# These scripts use the differential-abundance results created above!
source("02_analysis/06_stage3_species.R", encoding = "UTF-8", chdir = TRUE)
source("02_analysis/07_stage3_genus.R", encoding = "UTF-8", chdir = TRUE)
source("02_analysis/08_stage3_phylum.R", encoding = "UTF-8", chdir = TRUE)

# Create the final article figures and multi-panel layouts.
source(
  "03_reports/02_create_article_figures.R",
  encoding = "UTF-8",
  chdir = TRUE
)

# Render the reports

setwd("03_reports")

quarto::quarto_render("reports/00_data_description_report.qmd")
quarto::quarto_render("reports/01_stage1_gender_report.qmd")
quarto::quarto_render("reports/02_stage1_report.qmd")
quarto::quarto_render("reports/03_stage2_report.qmd")
quarto::quarto_render("reports/04_stage3_species_report.qmd")
quarto::quarto_render("reports/05_stage3_genus_report.qmd")
quarto::quarto_render("reports/06_stage3_phylum_report.qmd")
quarto::quarto_render("reports/07_summary_report.qmd")

# Render all report pages as one multi-page website. (will take some time )
quarto::quarto_render(input = ".", profile = "combined")

# Open the rendered reports in the default web browser.
browseURL(normalizePath("output/reports/00_data_description_report.html"))
browseURL(normalizePath("output/reports/01_stage1_gender_report.html"))
browseURL(normalizePath("output/reports/02_stage1_report.html"))
browseURL(normalizePath("output/reports/03_stage2_report.html"))
browseURL(normalizePath("output/reports/04_stage3_species_report.html"))
browseURL(normalizePath("output/reports/05_stage3_genus_report.html"))
browseURL(normalizePath("output/reports/06_stage3_phylum_report.html"))
browseURL(normalizePath("output/reports/07_summary_report.html"))

browseURL(
  normalizePath(
    "output/combined_report/index.html",
    mustWork = TRUE
  )
)
# Open the multi-page report in the default web browser.
browseURL(normalizePath("output/combined_report/index.html"))


# Return to the folder.
setwd("..")

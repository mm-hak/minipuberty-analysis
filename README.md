# Minipuberty analysis

This repository accompanies the manuscript *Perinatal antibiotic exposure may have an impact on male infant testosterone levels during the critical period of intense growth* (short title: *Perinatal antibiotics and infant testosterone*).

The repository contains R scripts and Quarto documents used to process the study data, run the analyses, and generate manuscript figures. The study examines associations between perinatal antibiotic exposure and testosterone levels in male infants during minipuberty. Related analyses examine early infant growth, other hormone measures, and the gut microbiota.

## Analysis overview

1. **Infant growth.** Examines growth during the first six months of life, using two intervals (0–3 and 3–6 months). The workflow includes analyses in the larger growth cohort and in a subset of male infants.
2. **Hormones and growth.** Examines associations of early-life exposures and maternal and infant characteristics with hormone levels in male infants. It also examines associations between hormone levels and growth. The hormone measures are testosterone, luteinizing hormone (LH), follicle-stimulating hormone (FSH), and sex hormone-binding globulin (SHBG).
3. **Gut microbiota.** Examines associations between early-life factors and gut microbiota composition, with a focus on bifidobacteria, and explores how microbiota measures relate to growth and hormone levels.

Blood and stool samples used for the hormone and microbiota analyses were collected at the three-month visit.

## Repository structure

| Path | Contents |
| --- | --- |
| `workflow.R` | Main script that runs the analysis steps in order. |
| `01_data_processing/` | Scripts that import and clean the source data and prepare the analysis datasets. |
| `02_analysis/` | Scripts for descriptive summaries and the three analysis stages, along with helper functions. |
| `03_reports/` | Quarto documents and scripts used to generate reports and manuscript figures. |

## Running the workflow

The workflow uses R. Quarto is required to render the reports.

1. Obtain authorized access to the required study data.
2. Open `01_data_processing/Data.R` and set the input paths to the locations of the source files on your machine.
3. Set the working directory in `workflow.R` to the project folder containing `workflow.R` and the numbered subfolders. Keep the folder structure shown above.
4. Run `workflow.R` in R.

## Data availability

This repository provides analysis code. Participant-level study data are required to reproduce the analyses and are kept separately. See the manuscript’s Data Availability Statement for information on access to the study data.

## License

The analysis code in this repository is licensed under the Artistic License 2.0. See [LICENSE](LICENSE).

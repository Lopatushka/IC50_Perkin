# Enter data
name_of_dir <- "10.12.21_MTS_MISIS/HEK293_MS"
cell_line_names <- "HEK293"
setwd("C:/Users/User/Google Диск/R_scipts/IC50_Perkin")

# Download functions
source("functions.R")

# Paths to files
working_dir <- getwd()
path_data <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = cell_line_names,
                        recursive=TRUE)
path_names <- list.files(path = name_of_dir, ignore.case=TRUE,
                         full.names = TRUE,pattern = "names", recursive=TRUE)

path_conc <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = "concentrations", recursive=TRUE)

path_export <- paste(name_of_dir, paste(cell_line_names, "results", sep="_"), sep="/")
path_CC50_lm <- paste(name_of_dir,
                          paste(paste("CC50", "HEK", sep="_"), "xlsx", sep="."),
                          sep='/')

# Download and process row data from one cell line
data_1 <- ImportDataFile_MISIS(path_data[1])
data_2 <- ImportDataFile_MISIS(path_data[2])
data_2[1] <- "Планшет 2"

data <- CombineTwoDataFiles(data_1, data_2)

data <- SubstractBackground_MISIS(df=data, wlength=490, backwlength=700)
data <- AddDrugNamesManual_MISIS(df=data, path_names=path_names, n_replicates=1)
data <- AddConcentrations_MISIS(df=data, path_conc=path_conc)
data <- DropNull(df=data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)
drug_names

# Subset and plot control (DMSO)
sb_drugs <- SubsetManual_MISIS(df=data, name="DMSO_drugs")
Plot(sb_drugs)
sb_MMAE <- SubsetManual_MISIS(df=data, name="DMSO_MMAE")
Plot(sb_MMAE)

# Find control medians and replace outliers
control_medians_drugs <- RmOutliersFromControl(sb_drugs)
control_medians_drugs

control_medians_MMAE <- RmOutliersFromControl(sb_MMAE)
control_medians_MMAE

# Fit curves: bunch processing
curves <- DRC_bunch_MISIS_new(df=data, drug_names=drug_names,
                                  controls=control_medians_MMAE,
                                  normilized=TRUE, start_dose=100,
                                  step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                                  path_export=path_export, export=FALSE,
                                  need_CCX=TRUE,
                                  manual_drugs_add=TRUE)

CC50s <- CC50_slope_bunch_MISIS(df=data, controls=control_medians_MMAE,
                                    path_to_table=path_CC50_lm,
                                    merge=TRUE, merge_with=curves)

Export_xlsx(df=CC50s,
            path=paste(paste(path_export, cell_line_names, sep="/"),
                       ".xlsx", sep=""))

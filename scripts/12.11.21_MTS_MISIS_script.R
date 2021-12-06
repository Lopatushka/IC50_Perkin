# Enter data
name_of_dir <- "12.11.21_MTS_MISIS"
cell_line_names <- "HEK"
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

path_export_HEK <- paste(name_of_dir, paste("HEK", "results", sep="_"), sep="/")
path_export_PC3 <- paste(name_of_dir, paste("PC3", "results", sep="_"), sep="/")

path_CC50_lm_HEK <- paste(name_of_dir,
                      paste(paste("CC50", "HEK", sep="_"), "xlsx", sep="."),
                      sep='/')

path_CC50_lm_PC3 <- paste(name_of_dir,
                          paste(paste("CC50", "PC3", sep="_"), "xlsx", sep="."),
                          sep='/')
# Download and process row data from one cell line
data <- ImportDataFile_MISIS(path_data[[1]])
data <- SubstractBackground_MISIS(df=data, wlength=490, backwlength=700)
data <- AddDrugNamesManual_MISIS(df=data, path_names=path_names, n_replicates=1)
data <- AddConcentrations_MISIS(df=data, path_conc=path_conc)
data <- DropNull(df=data)

# Change plate names
unique(data$Планшет)

data[data$Планшет == "Plate", 1] <-  "HEK293_1"
data[data$Планшет == "Plate (2)", 1] <-  "HEK293_2"
data[data$Планшет == "Plate (3)", 1] <-  "PC3_1"
data[data$Планшет == "Plate (4)", 1] <-  "PC3_2"

unique(data$Планшет)

# Split single file by cell lines
HEK <- data[grep("HEK293", data$Планшет),]
PC3 <- data[grep("PC3", data$Планшет),]

# List of all drugs in experiment
drug_names_HEK <- unique(HEK$Drug)
drug_names_PC3 <- unique(PC3$Drug)
cat("HEK:", drug_names_HEK)
cat("PC3:", drug_names_PC3)

# Subset and plot control (DMSO)
sb_HEK <- SubsetManual_MISIS(df=HEK, name="DMSO")
sb_PC3 <- SubsetManual_MISIS(df=PC3, name="DMSO")

Plot(sb_HEK)
Plot(sb_PC3)

# Find control medians and replace outliers
control_medians_HEK <- RmOutliersFromControl(sb_HEK)
control_medians_HEK

control_medians_PC3 <- RmOutliersFromControl(sb_PC3)
control_medians_PC3

# Fit curves: bunch processing
curves_HEK <- DRC_bunch_MISIS_new(df=HEK, drug_names=drug_names_HEK,
                    controls=control_medians_HEK,
                    normilized=TRUE, start_dose=100,
                    step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                    path_export=path_export_HEK, export=TRUE, need_CCX=TRUE,
                    manual_drugs_add=TRUE)

curves_PC3 <- DRC_bunch_MISIS_new(df=PC3, drug_names=drug_names_PC3,
                                  controls=control_medians_PC3,
                                  normilized=TRUE, start_dose=100,
                                  step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                                  path_export=path_export_PC3, export=TRUE,
                                  need_CCX=TRUE, manual_drugs_add=TRUE)

# CC50_lm for several drugs + merge to final table
CC50s_HEK <- CC50_slope_bunch_MISIS(df=HEK, controls=control_medians_HEK,
                          path_to_table=path_CC50_lm_HEK,
                          merge=TRUE, merge_with=curves_HEK)

CC50s_PC3 <- CC50_slope_bunch_MISIS(df=HEK, controls=control_medians_PC3,
                              path_to_table=path_CC50_lm_PC3,
                              merge=TRUE, merge_with=curves_PC3)
# Export
Export_xlsx(df=CC50s_HEK, path=paste(name_of_dir, "final_tables", "HEK_results.xlsx",
                                    sep='/'))

Export_xlsx(df=CC50s_PC3, path=paste(name_of_dir, "final_tables", "PC3_results.xlsx",
                                     sep='/'))

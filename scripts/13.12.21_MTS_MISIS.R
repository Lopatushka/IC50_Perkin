# Enter data
name_of_dir <- "13.12.21_MTS"
cell_line_names <- "PC3"
setwd("C:/Users/User/Google Диск/R_scipts/IC50_Perkin")

# Paths to files
working_dir <- getwd()
path_data <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = cell_line_names,
                        recursive=TRUE)[2]
path_names <- list.files(path = name_of_dir, ignore.case=TRUE,
                         full.names = TRUE,pattern = "names", recursive=TRUE)[3]

path_conc <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = "concentrations", recursive=TRUE)[3]

path_export <- paste(name_of_dir, paste("PC3", "results", sep="_"), sep="/")

path_CC50_lm <- list.files(path = name_of_dir, ignore.case=TRUE,
           full.names = TRUE,pattern = "CC50", recursive=TRUE)[3]


# Download and process row data from one cell line
data <- ImportDataFile_MISIS(path_data)
data <- SubstractBackground_MISIS(df=data, wlength=490, backwlength=700)
data <- AddDrugNamesManual_MISIS(df=data, path_names=path_names, n_replicates=1)
data <- AddConcentrations_MISIS(df=data, path_conc=path_conc)
data <- DropNull(df=data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)
drug_names

# Subset and plot control (DMSO)
sb_drugs <- SubsetManual_MISIS(df=data, name="DMSO")
Plot(sb_drugs)

# Find control medians and replace outliers
control_medians_drugs <- RmOutliersFromControl(sb_drugs)
control_medians_drugs

# Fit curves: bunch processing
curves <- DRC_bunch_MISIS_new(df=data, drug_names=drug_names,
                              controls=control_medians_drugs,
                              normilized=TRUE, start_dose=100,
                              step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                              path_export=path_export, export=FALSE,
                              need_CCX=TRUE,
                              manual_drugs_add=TRUE)

CC50s <- CC50_slope_bunch_MISIS(df=data, controls=control_medians_drugs,
                                path_to_table=path_CC50_lm,
                                merge=TRUE, merge_with=curves)

Export_xlsx(df=CC50s,
            path=paste(paste(path_export, cell_line_names, sep="/"),
                       ".xlsx", sep=""))

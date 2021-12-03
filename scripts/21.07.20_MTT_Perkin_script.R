# Enter data
name_of_dir <- "21.07.20_MTT_Perkin_data"
cell_line_names <- "Huh7"
setwd("C:/Users/User/Google Диск/R_scipts/IC50_Perkin/row_data")

# Download functions
source("functions.R")

# Paths to files
working_dir <- getwd()
path_data <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = cell_line_names, recursive=TRUE)
path_names <- list.files(path = name_of_dir, ignore.case=TRUE,
                         full.names = TRUE,pattern = "names", recursive=TRUE)
path_conc <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = "concentrations", recursive=TRUE)
path_CC50_lm <- paste(name_of_dir,
                      paste(paste("CC50", cell_line_names, sep="_"), "xlsx", sep="."),
                      sep='/')

# TODO Create a dir
path_export <- paste(name_of_dir, paste(cell_line_names, "results", sep="_"), sep="/")

# Download and process row data from one cell line
data <- ImportDataFile(path_data)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(df = data, path_conc = path_conc)
data <- DropNull(data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)
drug_names

# Subset and plot control (DMSO)
sb <- Subset(df = data, name = "DMSO")
Plot(df = sb)

# Find control medians and replace outliers
control_medians <- RmOutliersFromControl(sb)
control_medians

# Fit curves: bunch processing
curves <- DRC_bunch(df=data, drug_names=drug_names,
                    controls=control_medians,
                    normilized=TRUE, start_dose=100,
                    step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                    path_export=path_export, export=TRUE, need_CCX=TRUE)

drugs_of_interest <- c("TT19", "TT20", "TT21", "TT29", "TT30", "TT31",
                       "TT35", "TT36", "TT37", "TT41", "TT42", "TT43")
length(drugs_of_interest)

# CC50_lm for several drugs + merge to final table
CC50s <- CC50_slope_bunch(df=data, controls=control_medians,
                          path_to_table=path_CC50_lm,
                          merge=TRUE, merge_with=curves)

Export_xlsx(CC50s, paste(name_of_dir,
                        paste(cell_line_names, "results_merged.xlsx"), sep='/'))

# Merge final tables into a single one and export to .xlsx file
combined <- MergeFilesFun(path=paste(name_of_dir, "final_tables", sep='/'),
                          split_by=" ")

Export_xlsx(df=combined, path=paste(name_of_dir, "final_tables", "results_combined.xlsx",
                                    sep='/'))

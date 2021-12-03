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

# Download and process row data from one cell line
data <- ImportDataFile_MISIS(path_data)
data <- SubstractBackground_MISIS(df=data, wlength=490, backwlength=700)
data <- AddDrugNamesManual_MISIS(df=data, path_names=path_names, n_replicates=1)

data_c <- AddConcentrationsManual_MISIS(df=data)

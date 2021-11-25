# Enter data
name_of_dir <- "01.04.19_MTT_Perkin_data"
cell_line_names <- "HEK293"

# Download functions
source("functions.R")

# Paths to files
working_dir <- getwd()
path_data <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = cell_line_names)
path_names <- list.files(path = name_of_dir, ignore.case=TRUE,
                         full.names = TRUE,pattern = "names")
path_conc <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = "concentrations")

# TODO Create a dir
path_export <- paste(name_of_dir, "results", sep="/")

# Download and process row data from one cell line
data <- ImportDataFile(path_data)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(df = data, path_conc = path_conc)
data <- DropNull(data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)

# Subset and plot control (DMSO)
sb <- Subset(df = data, name = "DMSO")
Plot(df = sb)

# Find control medians and replace outliers
control_medians <- RmOutliersFromControl(sb)



todor::todor()
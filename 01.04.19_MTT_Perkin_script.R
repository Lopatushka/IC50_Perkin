# Enter data
name_of_dir <- "01.04.19_MTT_Perkin_data"
cell_line_names <- "HEK293"

# Download functions
source("functions.R")
todor::todor()

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

# Download row data
data <- ImportDataFile(path_data)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(df = data, path_conc = path_conc)
data <- DropNull(data)
source("functions.R")

# Paths to files
path_data1 <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/28.05.21_PC3.xls"
path_names <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/names.xlsx"
path_conc <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/concentrations.xlsx"
path_export <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/PC3"

# Download and process row data
data <- ImportDataFile(path_data1)
#data2 <- ImportDataFile(path_data2)
#data <- CombineTwoDataFiles(data1, data2)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(data)
data <- DropNull(data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)

# Subset and plot control (DMSO)
sb1 <- Subset(data, "DMSO")
Plot(sb1)

# Find control medians and replace outliers
control_medians_1 <- RmOutliersFromControl(sb1)

# Fit curve: bunch processing
curves <- DRC_bunch(df=data, drug_names=drug_names,
                    controls=control_medians_1,
                    normilized=TRUE, start_dose=100,
                    step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                    path_export=path_export, export=TRUE, CCX=TRUE)

# Fit CC50: bunch procassing
PC3_from <- c(13, 1, 6, 1, 2, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 9, 6, 2, 1, 5, 2, 13, 2, 2, 1, 4, 1, 1, 1, 3, 1)
PC3_to <- c(18, 24, 14, 24, 5, 9, 6, 6, 6, 24, 24, 6, 8, 24, 4, 24, 11, 12, 11, 7, 7, 12, 10, 18, 10, 10, 5, 9, 24, 24, 9, 7, 7)

CC50s <- CC50_slope_bunch(df=data, drug_names=drug_names,
                          controls=control_medians_1,
                          from=PC3_from,
                          to=PC3_to,
                          response=c(50))



# Construct final table and export it
final_table <- cbind(curves, CC50s)
write_xlsx(curves,
           paste(path_export, "/", "PC3_final_table.xlsx", sep=""))

slope <- CC50_slope(df=data, name = "EC52",
                    controls = control_medians_1,
                    from=1, to=4, response=c(50))
slope

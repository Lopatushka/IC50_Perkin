source("functions.R")

# Paths to files
path_data1 <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/HEK293 lena.xls"
path_data2 <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/HEK293 lena+skv.xls"
path_names <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/names.xlsx"
path_conc <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/concentrations.xlsx"

# Download and process row data
data1 <- ImportDataFile(path_data1)
data2 <- ImportDataFile(path_data2)
data <- CombineTwoDataFiles(data1, data2)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(data)
data <- DropNull(data)

# List of all drugs in experiment
drug_names <- unique(data$Drug)

# Subset and plot control (DMSO)
sb1 <- Subset(data, "DMSO_1")
sb2 <- Subset(data, "DMSO_2")
Plot(sb1)
Plot(sb2)

# Find control medians and replace outliers
control_medians_1 <- RmOutliersFromControl(sb1)
#control_medians_2 <- RmOutliersFromControl(sb2)

# Subset particular drug, normilize data and plot results
drug <- Subset(data, "GK140p")
Plot(drug)
drug <- Normalization(drug, control_medians_1)
Plot(drug, y=drug$D555_N)
#drug[-c(11:19), ]

# Fit the model and plot it
statistics <- DRC(df=drug, normilized=TRUE,
                  start_dose=100, step_dose=0.02,
                  X=50, plot=TRUE)

# Create an empty data frame for bind resuls
GKs <- data.frame(matrix(NA, ncol=20, nrow=0))
colnames(GKs) <- c('Drug', 'F val', 'p-val',
                   'Slope', 'LL','UL', 'ED50',
                   'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                   'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                   'Slope p-val', 'LL p-val','UL p-val','ED50 p-val',
                   'CC50')

# Add statistics to summary dataframe
GKs <- rbind(GKs, statistics)

# Export summary data frame
write_xlsx(GKs,"C:/Users/acer/Desktop/Work/Data/MTT/SKV/plots_new/HEK293/HEK293_results.xlsx")

drug_names

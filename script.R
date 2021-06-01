source("functions.R")

# Paths to files
path_data1 <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/28.05.21_HEK293.xls"
path_names <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/names.xlsx"
path_conc <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/row_data/concentrations.xlsx"
path_export <- "C:/Users/User/Documents/Work/Data/MTT/28.05.2021_MTT/HEK293"

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
                    path_export=path_export, export=TRUE, need_CCX=TRUE)

exclude=c("PAV7", )

boundaries <- list(name = drug_names[!is.element(drug_names, exclude)],
                   from=c(12, 2, 3, 3, 1, 1, 3, 4, 6, 3, 1, 5, 5, 6, 1, 1, 2, 1),
                   to=c(17, 8, 6, 8, 5, 6, 8, 9, 11, 7, 8, 10, 10, 11, 14, 14, 5, 6),
                   response=rep(50, length(drug_names)))

for(i in 1:length(boundaries$from))
{
  cat(sprintf("%s %s %s\n", boundaries$name[i], boundaries$from[i], boundaries$to[i]))
}


CC50s <- CC50_slope_bunch(df=data, controls=control_medians_1,
                          boundaries=boundaries)


# Construct final table and export it
final_table <- merge(curves, CC50s, by="Drug", all = T)

write_xlsx(curves,
           paste(path_export, "/", "HEK293_final_table.xlsx", sep=""))



slope <- CC50_slope(df=data, name = "EC51",
                    controls = control_medians_1,
                    from=2, to=4, response=c(50))

slope

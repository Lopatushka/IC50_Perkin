source("functions.R")

# Paths to files
path_data1 <- "C:/Users/User/Documents/Work/Data/MTT/mtt30_11_20/row_data/30.11.20_PC3.xls"
path_names <- "C:/Users/User/Documents/Work/Data/MTT/mtt30_11_20/row_data/names.xlsx"
path_conc <- "C:/Users/User/Documents/Work/Data/MTT/mtt30_11_20/row_data/concentrations.xlsx"
path_export <- "C:/Users/User/Documents/Work/Data/MTT/mtt30_11_20/PC3"

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

exclude=c("PAV15","PAV14-2","PAV14-1","PAV13","EC22", "PAV19-1","OA","DXK","XK", "TT59","DOC", "DMSO","LXK" )
drugs <- drug_names[!is.element(drug_names, exclude)]

boundaries <- list(name = drugs,
                   from=c(7, 3, 7, 4, 2, 1, 1, 4, 4, 4, 7, 4, 3, 6, 2, 1, 2, 1, 1, 1, 7),
                   to=c(12, 11, 12, 9, 7, 6, 12, 9, 9, 9, 14, 12, 10, 11, 7, 6, 7, 9, 6, 6, 11),
                   response=rep(50, length(drugs)))

for(i in 1:length(boundaries$from))
{
  cat(sprintf("%s %s %s\n", boundaries$name[i], boundaries$from[i], boundaries$to[i]))
}


CC50s <- CC50_slope_bunch(df=data, controls=control_medians_1,
                          boundaries=boundaries)


# Construct final table and export it
final_table <- merge(curves, CC50s, by="Drug", all = T)

write_xlsx(final_table,
           paste(path_export, "/", "PC3_merged_final.xlsx", sep=""))



slope <- CC50_slope(df=data, name = "BA",
                    controls = control_medians_1,
                    from=1, to=9, response=c(50))

slope

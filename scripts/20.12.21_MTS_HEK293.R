# Enter data
name_of_dir <- "20.12.21_MTS_HEK293"
cell_line_names <- "HEK293"
setwd("C:/Users/User/Google Диск/R_scipts/IC50_Perkin")

# Paths to files
working_dir <- getwd()
path_data <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = cell_line_names,
                        recursive=TRUE)
path_names <- list.files(path = name_of_dir, ignore.case=TRUE,
                         full.names = TRUE,pattern = "names", recursive=TRUE)

path_conc <- list.files(path = name_of_dir, ignore.case=TRUE,
                        full.names = TRUE,pattern = "concentrations", recursive=TRUE)

path_export <- paste(name_of_dir, paste(cell_line_names, "results", sep="_"), sep="/")

path_CC50_lm <- list.files(path = name_of_dir, ignore.case=TRUE,
                           full.names = TRUE,pattern = "CC50", recursive=TRUE)

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
sb_drugs <- SubsetManual_MISIS(df=data, name="DMSO_drugs")
Plot(sb_drugs)

sb_MMP86 <- SubsetManual_MISIS(df=data, name="DMSO_MMP86")
Plot(sb_MMP86)

# Find control medians and replace outliers
control_medians_drugs <- RmOutliersFromControl(sb_drugs)
control_medians_drugs

control_medians_MMP86 <- RmOutliersFromControl(sb_MMP86)
control_medians_MMP86

# Fit curves: bunch processing
curves <- DRC_bunch_MISIS_new(df=data, drug_names=drug_names,
                              exclude=c("MMP86", "DMSO_MMP86"),
                              controls=control_medians_drugs,
                              normilized=TRUE, start_dose=100,
                              step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                              path_export=path_export, export=FALSE,
                              need_CCX=TRUE,
                              manual_drugs_add=TRUE)

  
CC50s <- CC50_slope_bunch_MISIS(df=data, controls=control_medians_drugs,
                                path_to_table=path_CC50_lm,
                                merge=TRUE, merge_with=curves)

MMP86 <- data[data$Drug=="DMSO_MMP86" | data$Drug=="MMP86",]

curves_MMP86 <- DRC_bunch_MISIS_new(df=MMP86, drug_names=c("MMP86","DMSO_MMP86"),
                              controls=control_medians_MMP86,
                              normilized=TRUE, start_dose=100,
                              step_dose=0.02, X=50, plot=TRUE, save_plot=TRUE,
                              path_export=path_export, export=FALSE,
                              need_CCX=TRUE,
                              manual_drugs_add=TRUE)

CC50s_MMP86 <- CC50_slope_MISIS(df=MMP86, name="MMP86", from=15, to=20,
                                controls=control_medians_MMP86, normalized=TRUE,
                                response=c(50), manual_drugs_add=TRUE)


Export_xlsx(df=curves_MMP86,
            path=paste(paste(path_export, "MMP86_HEK293", sep="/"),
                       ".xlsx", sep=""))

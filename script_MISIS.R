source("functions.R")

# Paths
path_data <- "C:/Users/User/Documents/Work/Data/MTS/09.04.21_MTS/row_data/09.04.MTS_MCF7.xlsx"
export <- "C:/Users/User/Documents/Work/Data/MTS/09.04.21_MTS/MCF7"


# Function callings
data <- ImportDataFile_MISIS(path_data=path_data)
data <- SubstractBackground_MISIS(data, 490, 700)
data <- DropBlank_MISIS(data)
drug_names <- unique(data$Образец)
drug_names
conc_info_HEK <- list(drug=drug_names,
                  stock_conc=c(20.07, 20.08, 20.08, 20.05, 19.99, 19.97, 20, 20),
                  first_dilution=rep(200, 8),
                  step_dilution=c(3, 3, 2, 2, 3, 2, 2, 3),
                  n_dilutions=rep(8, 8),
                  n_replicates=rep(3, 8))

conc_info_HEK_2 <- list(drug=drug_names,
                      stock_conc=c(20.08, 20.07, 19.97, 20, 20),
                      first_dilution=rep(200, 5),
                      step_dilution=c(3, 3, 2, 2, 3),
                      n_dilutions=rep(8, 5),
                      n_replicates=c(3, 3, 2, 2, 3))

conc_info_MCF7_09.04 <- list(drug=drug_names,
                             stock_conc=rep(20, 10),
                             first_dilution=rep(200, 10),
                             step_dilution=rep(2, 10),
                             n_dilutions=rep(8, 10),
                             n_replicates=rep(3, 10))


data <- AddConcentrations_MISIS(data, conc_info_MCF7_09.04)

controls_dil2 <- RmOutliersFromControl(Subset_MISIS(data, "DMSO-dil2"), n_replicates=3)
controls_dil3 <- RmOutliersFromControl(Subset_MISIS(data, "DMSO-dil3"))

summary <- DRC_bunch_MISIS(df=data, drug_names=c("DG4CkSek","DG4ClSe","DG605k","DG606k","DG618k","DGAllC2"),
                           controls=controls_dil2, conc=conc_info_HEK,
                           normilized=TRUE, step_dose=0.02, X=50,
                           path_export=export, export=FALSE,
                           plot=TRUE, save_plot=FALSE, need_CCX=TRUE)



boundaries_MCF7 <- list(name=c("DG4ClSe","DG603k","DG605k","DG618k"),
                   from=c(1, 6, 7, 4),
                   to=c(6, 9, 11, 8),
                   response=rep(50, 4))

boundaries_HEK293 <- list(name=c("DG4CkSek","DG4ClSe", "DG605k", "DG606k", "DG618k", "DGAllC2"),
                        from=c(10, 1, 5, 1, 4, 13),
                        to=c(14, 4, 9, 5, 9, 18),
                        response=rep(50, 6))

boundaries_HEK293_2 <- list(name=c("4ClSe", "4ClSek", "DGAlC2"),
                            from=c(1, 13, 9),
                            to=c(4, 17, 12),
                            response=rep(50, 3))

CC50s <- CC50_slope_bunch_MISIS(df=data, controls=controls_dil3, conc=conc_info_HEK_2,
                                boundaries=boundaries_HEK293_2, normalized=TRUE,
                                exclude=c())

# Construct final table and export it
write_xlsx(summary,
           paste(export, "/", "HEK293_2_final_table.xlsx", sep=""))

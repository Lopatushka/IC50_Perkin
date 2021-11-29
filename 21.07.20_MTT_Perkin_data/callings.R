source('C:/Users/acer/Desktop/work/14.07.20_MTT_functions.R')

getwd()
setwd("C:/Users/acer/Desktop/work/21.07.20_mtt")

names <- ImportNames("./names.xlsx")

df <- ImportData("./21.07.20_pc3.xls", nplates_1file=9)

conc <- ImportConcentrations("./concentrations.xlsx")

predrc <- DataProcessing(df = df, nblocks = 3, control_name="DMSO",
                         df_c=conc, export=TRUE, path_export="./pc3_row_data")

results <- DrcOutliers(df_drc = predrc, exclude="DMSO", 
                       plot=TRUE,export_plot=TRUE, path_export="./pc3_plots",
                       short=TRUE, max_value=20)

ExportResults(df = results)

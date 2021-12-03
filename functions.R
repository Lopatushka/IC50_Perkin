# Import packages
library(readxl)
library(drc)
library(outliers)
library(writexl)
library(purrr)
library(dplyr)
library(writexl)
library (berryFunctions)

# Import .xls file with D555 data
ImportDataFile <- function(path_data)
{
  df <- read_excel(path_data)
  colnames(df)[6] <- "D555"
  df$D555 <- sapply(df$D555, as.double)
  return(df)
}

# Export .xlsx files
Export_xlsx <- function(df, path="./results.xlsx")
{
  write_xlsx(df, path)
}

# Combine 2 dataframe with row data into single one
CombineTwoDataFiles <- function(df1, df2)
{
  for (i in 1:nrow(df2))
  {
    df2$Plate[i] <- df2$Plate[i]+ df1$Plate[as.numeric(dim(df1)[1])]
  }
  df <- rbind(df1,df2)
  return(df)
}

# Add Drug name col to dataframe
AddDrugNames <- function(df, path_names, plate_type=96, n_replicates=3)
{
  names <- read_excel(path_names)
  df$Block <- NA
  block_names <- colnames(names)
  df$Block <- rep(block_names, each=plate_type*n_replicates)
  df$Drug <- NA
  for (row in 1:dim(df)[1])
  {
    name_well <- as.numeric(substr(df$Well[row], start=2, stop=3))
    name_block <- df$Block[row]
    df$Drug[row] <- names[colnames(names) == name_block][[1]][name_well]
  }
  return(df)
}


# Add Concentration col to dataframe
AddConcentrations <- function(df, path_conc, first_dilution=200,step_dilution=3,
                              n_dilutions=8,  n_replicates=3)
{
  out <- tryCatch(
    expr = {
      # Create new col in dataframe and add null to Drug null
      df$C_mkM <- NA
      df$C_mkM[df$Drug == 'null'] <- 'null'
      
      # Make a list of drugs
      list_of_drugs <- unique(df$Drug)
      
      # Read concentration file
      concentrations <- read_excel(path_conc)
      
      # Big for cycle
      for (drug in list_of_drugs)
      {
        if (drug!='null')
        {
          #print(drug)
          # Create a temp dataframe
          temp <- data.frame(matrix(NA, ncol=1, nrow=n_dilutions))
          rownames(temp) <- paste("C_mkM", 1:n_dilutions,  sep="_")
          colnames(temp) <- drug
          
          # Fill the temp dataframe: First dilution, turn to mkM
          temp[1,] <- 1000*concentrations[concentrations$Drug == drug, ][2] /first_dilution
          # Other dilutions
          for (row in 2:nrow(temp)) temp[row,] <- temp[row-1, ]/step_dilution
          
          # Fill C_mkM col in df
          df$C_mkM[df$Drug==drug] <- rep(temp[, 1], n_replicates)
        }
      }
      
      df$C_mkM <- lapply(df$C_mkM, as.numeric)
      return(df)
    },
    error = function(e){
      message(paste("There is no drug in the concentration file:", drug))
      message(paste('The original message from R is below:'))
      message(e)
    }
  )
  return(out)
}

# Drop rows with null drug name
DropNull <- function(df)
{
  return(as.data.frame(subset(df, df$Drug!='null')))
}

# Subset rows with particular drug name
Subset <- function(df, name)
{
  out <- tryCatch(
    expr = {
      drug <- subset(df, df$Drug == name)
      drug <- drug[, c(6, 9, 8)]
      drug <- drug[order(unlist(drug$C_mkM), decreasing = TRUE), ]
    },
    error = function(e){
      message(paste("There is no drug in the data:", name))
      message(paste('The original message from R is below:'))
      message(e)
    }
  )
  return(out)
}

# Find plateau in control (DMSO) subset using linear regression
FindPlateuForControl <- function(df, alpha=0.05)
{
  notPlateu <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(notPlateu) <- colnames(df)
  
  lm <- lm(df$D555 ~ unlist(df$C_mkM))
  p_val <- summary(lm)$coefficients[2,4]
  
  while (p_val < alpha)
  {
    notPlateu <- rbind(notPlateu, df[(1:3), ])
    df <- df[-(1:3), ]
    lm <- lm(df$D555 ~ unlist(df$C_mkM))
    p_val <- summary(lm)$coefficients[2,4]
  }
  
  return(list(df, notPlateu))
}

# Plot subset
Plot <- function(df, x=df$C_mkM, y=df$D555, log=TRUE)
{
  x = unlist(x)
  if(log==TRUE) x=log(x, base=10)
  
  plot(x=x,
       y=y,
       xlab='Log10[C], mkM', ylab='D555', main=df[1, 3])
}

# Find medians of D555 data. Return vector
FindMedians <- function(df, n_dilutions=8,  n_replicates=3)
{
  medians <- c()
  i <- 1
  while (i < n_dilutions*n_replicates + 1)
  {
    median <- median(df$D555[i:(i+2)])
    medians <- c(medians, median)
    i <- i + 3
  }
  medians <- medians[!is.na(medians)]
  return(medians)
}

# Remove outliers in control medians and replace them with median values
# Return vector
RmOutliersFromControl <- function(df, alpha=0.05, n_dilutions=8,  n_replicates=3)
{
  pl <- FindPlateuForControl(df, alpha=alpha)[[1]]
  non_pl <- FindPlateuForControl(df, alpha=alpha)[[2]]
  
  pl_medians <- FindMedians(pl, n_dilutions=n_dilutions,  n_replicates=n_replicates)
  pl_medians <- rm.outlier(x=pl_medians, fill=TRUE, median=TRUE)
  
  non_pl_medians <- FindMedians(non_pl, n_dilutions=n_dilutions,  n_replicates=n_replicates)
  if(length(non_pl_medians)==0){
    non_pl_medians_rm <- NULL
  }
  if(length(non_pl_medians)==1){
    non_pl_medians_rm <- non_pl_medians
  }
  if(length(non_pl_medians)>1){
    non_pl_medians_rm <- rm.outlier(x=non_pl_medians, fill=TRUE, median=TRUE)
  }
  
  return(c(non_pl_medians_rm, pl_medians))
  
  
}

RmOutliers <- function(df, var="D555")
{
  return(df[!outlier(df[[var]], logical = TRUE), ])
}

# Create a new col D555_N in data subset by normalizing D555 by control vector
# Return subset
Normalization <- function(df, controls, n_replicates=3)
{
  df$D555_N <- 100*df$D555/rep(controls, each=n_replicates)
  return(df)
}

# Find CC50 (or others) using predict fun. Return CC50 concentration without SE
CCX <- function(model, start_dose=100, step_dose=0.01, X=50)
{
  dose <- start_dose
  result <- predict(model, data.frame(dose), se.fit=TRUE)
  while(result[[1]]<X & dose > step_dose)
  {
    dose <- dose-step_dose
    result <- predict(model, data.frame(dose), se.fit=TRUE)
  }
  return(dose)
}

# Fit curve for one drug
DRC <- function(df, normilized=TRUE,
                start_dose=100, step_dose=0.02, X=50,
                plot=TRUE, save_plot=TRUE, path_export=".", need_CCX=TRUE)
{
  if(need_CCX==TRUE)
  {
    results <- data.frame(matrix(NA, ncol=20, nrow=1))
    colnames(results) <- c('Drug', 'F val', 'p-val',
                           'Slope', 'LL','UL', 'ED50',
                           'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                           'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                           'Slope p-val', 'LL p-val','UL p-val','ED50 p-val', 'CC50')
  }
  
  if(need_CCX==FALSE)
  {
    results <- data.frame(matrix(NA, ncol=19, nrow=1))
    colnames(results) <- c('Drug', 'F val', 'p-val',
                           'Slope', 'LL','UL', 'ED50',
                           'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                           'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                           'Slope p-val', 'LL p-val','UL p-val','ED50 p-val')
  }
  
  drug_name <- df[1, 3]
  
  if(normilized==TRUE)
  {
    df=df[, c(2,4)]
  }
  if(normilized==FALSE)
  {
    df=df[, c(2,1)]
  }
  
  colnames(df) <- c("C_mkM", "D555")
  
  model <- try_model_fun(code=(drm(D555 ~ unlist(C_mkM),
               data = df,
               fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))))
  #Logic variable: TRUE - converge error, FALSE - normals
  converge_error <- is.logical(model)
  
  if (converge_error==FALSE) # No converge error
  {
    model_fit <- data.frame(modelFit(model))
    coeffs <- summary(model)$coefficients
    
    results$Drug <- drug_name
    results$`F val` <- round(model_fit$"F.value"[2], digits=3)
    results$`p-val` <- round(model_fit$"p.value"[2], digits=3)
    
    i <- 4
    for(coef in coeffs)
    {
      results[1, i] <- round(coef, digits=3)
      i <- i + 1
    }
  }
  
  if (converge_error==TRUE) # Converge error
  {
    results$Drug <- drug_name
    sprintf('Converged error for %s', drug_name)
  }
  
  
  if(need_CCX==TRUE & converge_error==FALSE)
  {
    results$CC50 <- CCX(model=model, start_dose=start_dose,
                        step_dose=step_dose, X=X)
  }
  
  
  if(plot==TRUE & converge_error==FALSE)
  {
    plot(model, broken=TRUE, bty="l",
         xlab="Log(drug concentration)", ylab="Response",
         main=drug_name, type="all")
  }
  
  if(save_plot==TRUE)
  {
    path_to_image <- paste(path_export, '/', drug_name,'.jpeg', sep="")
    dev.copy(jpeg, filename=path_to_image)
    dev.off() # Close plots
  }
  
  return(results)
}

# Fit curves for drugs in list drug_names
DRC_bunch <- function(df, drug_names, controls,
                      normilized=TRUE, start_dose=100,
                      step_dose=0.02, X=50,
                      path_export=".", export=TRUE,
                      plot=TRUE, save_plot=TRUE, need_CCX=TRUE)
{
  message("Waiting...")
  # Create an empty data frame for bind resuls
  if(need_CCX==TRUE)
  {
    GKs <- data.frame(matrix(NA, ncol=20, nrow=0))
    colnames(GKs) <- c('Drug', 'F val', 'p-val',
                       'Slope', 'LL','UL', 'ED50',
                       'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                       'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                       'Slope p-val', 'LL p-val','UL p-val','ED50 p-val',
                       'CC50')
  }
  
  if(need_CCX==FALSE)
  {
    GKs <- data.frame(matrix(NA, ncol=19, nrow=0))
    colnames(GKs) <- c('Drug', 'F val', 'p-val',
                       'Slope', 'LL','UL', 'ED50',
                       'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                       'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                       'Slope p-val', 'LL p-val','UL p-val','ED50 p-val')
  }
  
  
  for (name in drug_names)
  {
    drug <- Subset(df, name)
    drug <- Normalization(drug, controls)
    drug <- RmOutliers(drug)
    statistics <- DRC(df=drug, normilized=normilized,
                      start_dose=start_dose, step_dose=step_dose,
                      path_export=path_export,
                      X=X, plot=plot, save_plot=save_plot, need_CCX=need_CCX)
    GKs <- rbind(GKs, statistics)
    if(export==TRUE)
    {
      path <- paste(path_export, '/', name,'.xlsx', sep="")
      write_xlsx(drug, path)
    }
  }
  message("Done!")
  return(GKs)
}

# Fit linear model and find CC50, SE, CIs for one drug
# name, from, to are str and int
CC50_slope <- function(df, name, from, to, controls, normalized=TRUE,
                       response=c(50))
{
  drug <- Subset(df, name)
  if(normalized==TRUE)
  {
    drug <- Normalization(drug, controls=controls)
  }
  
  drug <- RmOutliers(drug)
  drug <- drug[from:to, ]
  
  
  if(normalized==TRUE)
  {
    model <- lm(unlist(C_mkM) ~ D555_N, data=drug)
    predition <- predict(model, data.frame(D555_N=response),
                         se.fit=TRUE, interval = "confidence")
  }
  
  if(normalized==FALSE)
  {
    model <- lm(unlist(C_mkM) ~ D555, data=drug)
    predition <- predict(model, data.frame(D555=response),
                         se.fit=TRUE, interval = "confidence")
  }
  
  results <- data.frame(matrix(NA, ncol=5, nrow=1))
  colnames(results) <- c('Drug', 'CC50', 'Lower', 'Upper', 'SE')
  results[1] <- drug[1,3]
  results[2:4] <- predition$fit
  results[5] <- predition$se.fit
  
  return(results)
}

# Fit linear model and find CC50, SE, CIs for several drugs
# drug_names, from, to are vectors
CC50_slope_bunch <- function(df, path_to_table, controls,
                             normalized=TRUE, merge=FALSE, merge_with=curves)
{
  # Create new dataframe
  results <- data.frame(matrix(NA, ncol=5, nrow=0))
  colnames(results) <- c('Drug', 'CC50', 'Lower', 'Upper', 'SE')
  
  # Download a table
  boundaries <- read_excel(path_to_table)
  
  # Number of drugs
  n_drugs <- nrow(boundaries)
  for (i in 1:n_drugs)
  {
    drug_name <- boundaries$Drug[i]
    print(drug_name)
    temp <- CC50_slope(df=df, name=drug_name,
                         controls=controls,
                         from=boundaries[boundaries$Drug==drug_name, ][[2]],
                         to=boundaries[boundaries$Drug==drug_name, ][[3]],
                         response=boundaries[boundaries$Drug==drug_name, ][[4]])
    results <- rbind(results, temp)
    
  }
  
  if (merge==TRUE) return(merge(merge_with, results, by="Drug", all = T))
  
  return(results)
}

try_model_fun <- function(code)
{
  temp <- try(code, silent = TRUE)
  if (class(temp)=="try-error") temp <- NA
  return(temp)
}

MergeFilesFun <- function(path, split_by=" ", files_type=".xlsx")
{
  # List of files
  paths_files <- list.files(path, full.names = TRUE,  pattern = files_type)
  list_of_files <- list.files(path, full.names = FALSE, pattern = files_type)
  n_files <- length(list_of_files)
  vector_names <- c()
  for (name in list_of_files)
  {
    vector_names <- c(vector_names, strsplit(name, split_by)[[1]][1])
  }
  
  temp <- list()
  for (file in paths_files)
  {
    df <- read_xlsx(file)
    df <- df[c(1, 7, 11, 19, 20:24)]
    df <- df[!is.na(df$CC50.y),]
    df$CC50.x[df$CC50.x==100] <- ">100"
    temp <- append(temp, df)
  }
  # Combine files into a single dataframe
  combined_data <- as.data.frame(do.call(cbind, temp))
  # Add emptu row above
  combined_data <- insertRows(combined_data, 1 , new = NA)
  # Add cell line names
  j <- 1
  for (i in 1:(n_files*9))
  {
    if(i%%9==1) 
    {
      combined_data[1, i] <- vector_names[j]
      j <- j + 1
    }
  }
  
  return(combined_data)
}

### MISIS
DropBlank_MISIS_2 <- function(df) 
{
  df_drop <- df[df[4]!="Бланк",]
  return(df_drop)
}

AddDrugNamesManual_MISIS <- function(df, path_names, plate_type=96, n_replicates=3)
{
  names <- read_excel(path_names)
  df$Block <- NA
  block_names <- colnames(names)
  df$Block <- rep(block_names, each=plate_type*n_replicates)
  df$Drug <- NA
  for (row in 1:dim(df)[1])
  {
    name_well <- as.numeric(substr(df$Лунка[row], start=2, stop=3))
    name_block <- df$Block[row]
    df$Drug[row] <- names[colnames(names)==name_block][[1]][name_well]
  }
  return(df)
}

split_string <- function(s) return(strsplit(s, "_")[[1]][1])
convert_to_numeric <- function(s) return(as.double(s))

ImportDataFile_MISIS <- function(path_data)
{
  df <- read_excel(path_data)
  colnames(df) <- df[2,]
  df <- df[-c(1, 2), ]
  df <- df[, 1:8]
  df[5] <- unlist(map(strsplit(dplyr::pull(df[5]), split="_"), 1))
  df[7] <- sapply(df[7], convert_to_numeric)
  return(df)
}

SubstractBackground_MISIS <- function(df, wlength=490, backwlength=700)
{
  temp <- df[df[6]==wlength,]
  temp[7] <- subset(df[7], df[6]==wlength)-subset(df[7], df[6]==backwlength)
  return(temp)
}

SubstractBackground_2files_MISIS <- function(df1, df2)
{
  df1[7] <- df1[7]-df2[7]
  return(df1)
}

DrugList_MISIS <- function(df) return(as.vector(unlist(unique(df[5]))))

AddConcentrations_MISIS <- function(df, conc)
{
  temp <- data.frame(matrix(NA, ncol=length(conc$drug), nrow=max(conc$n_dilutions)))
  rownames(temp) <- paste("C_mkM", 1:max(conc$n_dilutions),  sep="_")
  colnames(temp) <- conc$drug
  temp[1, ] <- 1000*conc$stock_conc/conc$first_dilution
  for(i in 1:(max(conc$n_dilutions)-1)) temp[i+1, ] <- temp[i, ]/conc$step_dilution
  
  df$C_mkM <- NA
  drug_list <- DrugList_MISIS(df)
  for (i in 1:length(drug_list))
  {
    drug <- drug_list[i]
    drug_conc <- temp[colnames(temp) == drug][[1]]
    df$C_mkM[df[5]==drug] <- rep(drug_conc, conc$n_replicates[i])
  }
  
  # Convert concentrations to numbers
  df$C_mkM <- sapply(df$C_mkM, as.double)
  
  return(df)
}

# Subset rows with particular drug name
Subset_MISIS <- function(df, name)
{
  drug <- subset(df, df[5]==name)
  drug <- drug[, c(9, 7, 5)]
  colnames(drug) <- c("C_mkM", "D555", "Drug")
  drug <- drug[order(drug$C_mkM, decreasing = TRUE), ]
  drug <- as.data.frame(drug)
  return(drug)
}

# Fit curve for one drug
DRC_MISIS <- function(df, normilized=TRUE,
                      step_dose=0.02, X=50,
                      plot=TRUE, save_plot=TRUE, path_export=".", need_CCX=TRUE)
{
  if(need_CCX==TRUE)
  {
    results <- data.frame(matrix(NA, ncol=20, nrow=1))
    colnames(results) <- c('Drug', 'F val', 'p-val',
                           'Slope', 'LL','UL', 'ED50',
                           'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                           'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                           'Slope p-val', 'LL p-val','UL p-val','ED50 p-val', 'CC50')
  }
  
  if(need_CCX==FALSE)
  {
    results <- data.frame(matrix(NA, ncol=19, nrow=1))
    colnames(results) <- c('Drug', 'F val', 'p-val',
                           'Slope', 'LL','UL', 'ED50',
                           'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                           'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                           'Slope p-val', 'LL p-val','UL p-val','ED50 p-val')
  }
  
  drug_name <- df[1, 3]
  
  if(normilized==TRUE)
  {
    df=df[, c(1,4)]
  }
  if(normilized==FALSE)
  {
    df=df[, c(1,2)]
  }
  
  colnames(df) <- c("C_mkM", "D555")
  
  model <- try_model_fun(code=(drm(D555 ~ C_mkM,
                                   data = df,
                                   fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))))
  #Logic variable: TRUE - converge error, FALSE - normals
  converge_error <- is.logical(model)
  
  if (converge_error==FALSE) # No converge error
  {
    model_fit <- data.frame(modelFit(model))
    coeffs <- summary(model)$coefficients
    
    results$Drug <- drug_name
    results$`F val` <- round(model_fit$"F.value"[2], digits=3)
    results$`p-val` <- round(model_fit$"p.value"[2], digits=3)
    
    i <- 4
    for(coef in coeffs)
    {
      results[1, i] <- round(coef, digits=3)
      i <- i + 1
    }
  }
  
  if (converge_error==TRUE) # Converge error
  {
    results$Drug <- drug_name
    sprintf('Converged error for %s', drug_name)
  }
  
  
  if(need_CCX==TRUE & converge_error==FALSE)
  {
    results$CC50 <- CCX(model=model, start_dose=df$C_mkM[1],
                        step_dose=step_dose, X=X)
  }
  
  
  if(plot==TRUE & converge_error==FALSE)
  {
    plot(model, broken=TRUE, bty="l",
         xlab="Log(drug concentration)", ylab="Response",
         main=drug_name, type="all")
  }
  
  if(save_plot==TRUE)
  {
    path_to_image <- paste(path_export, '/', drug_name,'.jpeg', sep="")
    dev.copy(jpeg, filename=path_to_image)
    dev.off() # Close plots
  }
  
  return(results)
}

# Fit curves for drugs in list drug_names
DRC_bunch_MISIS <- function(df, drug_names, controls, conc,
                            normilized=TRUE, step_dose=0.02, X=50,
                            path_export=".", export=TRUE,
                            plot=TRUE, save_plot=TRUE, need_CCX=TRUE)
{
  # Create an empty data frame for bind resuls
  if(need_CCX==TRUE)
  {
    GKs <- data.frame(matrix(NA, ncol=20, nrow=0))
    colnames(GKs) <- c('Drug', 'F val', 'p-val',
                       'Slope', 'LL','UL', 'ED50',
                       'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                       'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                       'Slope p-val', 'LL p-val','UL p-val','ED50 p-val',
                       'CC50')
  }
  
  if(need_CCX==FALSE)
  {
    GKs <- data.frame(matrix(NA, ncol=19, nrow=0))
    colnames(GKs) <- c('Drug', 'F val', 'p-val',
                       'Slope', 'LL','UL', 'ED50',
                       'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                       'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                       'Slope p-val', 'LL p-val','UL p-val','ED50 p-val')
  }
  
  for (name in drug_names)
  {
    drug <- Subset_MISIS(df, name)
    
    replics <- conc$n_replicates[match(name, conc$drug)]
    drug <- Normalization(drug, controls, n_replicates = replics)
    
    drug <- RmOutliers(drug)
    statistics <- DRC_MISIS(df=drug, normilized=normilized,
                            step_dose=step_dose, path_export=path_export,
                            X=X, plot=plot, save_plot=save_plot, need_CCX=need_CCX)
    GKs <- rbind(GKs, statistics)
    if(export==TRUE)
    {
      path <- paste(path_export, '/', name,'.xlsx', sep="")
      write_xlsx(drug, path)
    }
  }
  return(GKs)
}

# Fit linear model and find CC50, SE, CIs for one drug
# name, from, to are str and int
CC50_slope_MISIS <- function(df, name, controls, normalized=TRUE,
                             from, to, response=c(50), conc)
{
  drug <- Subset_MISIS(df, name)
  if(normalized==TRUE)
  {
    replics <- conc$n_replicates[match(name, conc$drug)]
    drug <- Normalization(drug, controls=controls, n_replicates = replics)
  }
  
  drug <- RmOutliers(drug)
  drug <- drug[from:to, ]
  
  if(normalized==TRUE)
  {
    model <- lm(C_mkM ~ D555_N, data=drug)
  }
  
  if(normalized==FALSE)
  {
    model <- lm(C_mkM ~ D555, data=drug)
  }
  
  
  predition <- predict(model, data.frame(D555_N=response),
                       se.fit=TRUE, interval = "confidence")
  
  results <- data.frame(matrix(NA, ncol=5, nrow=1))
  colnames(results) <- c('Drug', 'CC50', 'Lower', 'Upper', 'SE')
  results[1] <- drug[1,3]
  results[2:4] <- predition$fit
  results[5] <- predition$se.fit
  
  return(results)
}


# Fit linear model and find CC50, SE, CIs for several drugs
# drug_names, from, to are vectors
CC50_slope_bunch_MISIS <- function(df, controls, conc,
                                   boundaries, normalized=TRUE,
                                   exclude=c())
{
  results <- data.frame(matrix(NA, ncol=5, nrow=0))
  colnames(results) <- c('Drug', 'CC50', 'Lower', 'Upper', 'SE')
  
  n_drugs <- length(boundaries$name) - length(exclude)
  for (i in 1:n_drugs)
  {
    #print(c(i, boundaries$name[i], boundaries$from[i], boundaries$to[i]))
    if(!is.element(boundaries$name[i], exclude))
    {
      temp <- CC50_slope_MISIS(df=df, name=boundaries$name[i],
                               controls=controls, conc=conc,
                               from=boundaries$from[i], to=boundaries$to[i],
                               response=boundaries$response[i], normalized=normalized)
      results <- rbind(results, temp)
    }
    
  }
  return(results)
}


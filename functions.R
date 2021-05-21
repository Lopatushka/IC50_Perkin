# Import packages
library(drc)
library(outliers)
library(writexl)

# Import .xls file with D555 data
ImportDataFile <- function(path_data)
{
  library(readxl)
  df <- read_excel(path_data)
  colnames(df)[6] <- "D555"
  df$D555 <- sapply(df$D555, as.double)
  return(df)
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
AddConcentrations <- function(df, first_dilution=200,step_dilution=3,
                              n_dilutions=8,  n_replicates=3)
{
  concentrations <- read_excel(path_conc)
  temp <- data.frame(matrix(NA, ncol=ncol(concentrations), nrow=n_dilutions))
  rownames(temp) <- paste("C_mkM", 1:n_dilutions,  sep="_")
  colnames(temp) <- colnames(concentrations)
  concentrations <- as.data.frame(rbind(concentrations, temp))
  rownames(concentrations)[1] <- 'C_stock_mM'
  
  concentrations[2,] <- 1000*as.numeric(concentrations[1, ])/first_dilution
  for (row in 3:nrow(concentrations))
  {
    concentrations[row,] <- concentrations[row-1, ]/step_dilution
  }
  
  #log10 transformation
  #concentrations[2:nrow(concentrations), ] <- log10(concentrations[2:nrow(concentrations),])
  
  # Add null to Drug null
  df$C_mkM <- NA
  df$C_mkM[df$Drug == 'null'] <- 'null'
  
  for (i in 1:length(unique(df$Drug)))
  {
    drug <- unique(df$Drug)[i]
    if (drug!='null')
    {
      conc <- concentrations[colnames(concentrations) == drug][[1]][-1]
      df$C_mkM[df$Drug==drug] <- rep(conc, n_replicates)
    }
  }
  
  
  df$C_mkM <- sapply(df$C_mkM, as.double)
  return(df)
}

# Drop rows with null drug name
DropNull <- function(df)
{
  return(as.data.frame(subset(df, df$Drug!='null')))
}

# Subset rows with particular drug name
Subset <- function(df, name)
{
  drug <- subset(df, df$Drug == name)
  drug <- drug[, c(6, 9, 8)]
  drug <- drug[order(drug$C_mkM, decreasing = TRUE), ]
  return(drug)
}

# Find plateau in control (DMSO) subset using linear regression
FindPlateuForControl <- function(df, alpha=0.05)
{
  notPlateu <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(notPlateu) <- colnames(df)
  
  lm <- lm(df$D555 ~ df$C_mkM)
  p_val <- summary(lm)$coefficients[2,4]
  
  while (p_val < alpha)
  {
    notPlateu <- rbind(notPlateu, df[(1:3), ])
    df <- df[-(1:3), ]
    lm <- lm(df$D555 ~ df$C_mkM)
    p_val <- summary(lm)$coefficients[2,4]
  }
  
  return(list(df, notPlateu))
}

# Plot subset
Plot <- function(df, x=df$C_mkM, y=df$D555, log=TRUE)
{
  if(log==TRUE)
  {
    x=log10(x)
  }
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

# Remove outliers in control medians and replce them with median values
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

# Find CC50 (or others) using predict fun. Return CC50 concentration
CCX <- function(model, start_dose=100, step_dose=0.01, X=50)
{
  dose <- start_dose
  result <- predict(model, data.frame(dose), se.fit=TRUE)
  while(result[[1]]<X)
  {
    dose <- dose-step_dose
    result <- predict(model, data.frame(dose), se.fit=TRUE)
  }
  return(dose)
}

DRC <- function(df, normilized=TRUE,
                start_dose=100, step_dose=0.02, X=50, plot=TRUE)
{
  results <- data.frame(matrix(NA, ncol=20, nrow=1))
  colnames(results) <- c('Drug', 'F val', 'p-val',
                         'Slope', 'LL','UL', 'ED50',
                         'Slope SE', 'LL SE', 'UL SE', 'ED50 SE',
                         'Slope t-val', 'LL t-val','UL t-val','ED50 t-val',
                         'Slope p-val', 'LL p-val','UL p-val','ED50 p-val',
                         'CC50')
  
  
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
  
  model <- drm(D555 ~ C_mkM,
               data = df,
               fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))
  
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
  
  results$CC50 <- CCX(model=model, start_dose=start_dose,
                      step_dose=step_dose, X=X)
  
  if(plot==TRUE)
  {
    plot(model, broken=TRUE, bty="l",
         xlab="Log(drug concentration)", ylab="Response",
         main=drug_name, type="all")
  }
  return(results)
}

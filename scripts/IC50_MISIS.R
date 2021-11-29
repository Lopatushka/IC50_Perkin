# Import packages
library(readxl)
library(drc)
library(outliers)
library(writexl)

# Paths
path_data <- "C:/Users/User/Documents/Work/Data/MTS/12.04.21_MTS/row data/12.04.21_MTS_MCF7.xlsx"

# Functions
split_string <- function(s) return(strsplit(s, "_")[[1]][1])
convert_to_numeric <- function(s) return(as.double(s))

ImportDataFile_MISIS <- function(path_data)
{
  df <- read_excel(path_data)
  colnames(df) <- df[2,]
  df <- df[-c(1, 2), 1:8]
  df$Образец <- sapply(df$Образец, split_string)
  df$Погл. <- sapply(df$Погл., convert_to_numeric)
  return(df)
}

SubstractBackground_MISIS <- function(df, wlength=490, backwlength=700)
{
  temp <- df[df$`Длина волны`==wlength,]
  temp$Погл. <- df$Погл.[df$`Длина волны`==wlength]-df$Погл.[df$`Длина волны`==backwlength]
  return(temp)
}

DropBlank_MISIS <- function(df) return(df[!df$Тип=="Бланк",])

AddConcentrations_MISIS <- function(df, conc)
{
  temp <- data.frame(matrix(NA, ncol=length(conc$drug), nrow=max(conc$n_dilutions)))
  rownames(temp) <- paste("C_mkM", 1:max(conc$n_dilutions),  sep="_")
  colnames(temp) <- conc$drug
  temp[1, ] <- 1000*conc$stock_conc/conc$first_dilution
  for(i in 1:(max(conc$n_dilutions)-1)) temp[i+1, ] <- temp[i, ]/conc$step_dilution
  
  df$C_mkM <- NA
  for (i in 1:length(unique(df$Образец)))
  {
    drug <- unique(df$Образец)[i]
    drug_conc <- temp[colnames(temp) == drug][[1]]
    df$C_mkM[df$Образец==drug] <- rep(drug_conc, conc$n_replicates[i])
  }
  
  # Convert concentrations to numbers
  df$C_mkM <- sapply(df$C_mkM, as.double)
  
  return(df)
}

# Subset rows with particular drug name
Subset_MISIS <- function(df, name)
{
  drug <- subset(df, df$Образец==name)
  drug <- drug[, c(9, 7, 5)]
  colnames(drug) <- c("C_mkM", "D555", "Drug")
  drug <- drug[order(drug$C_mkM, decreasing = TRUE), ]
  return(drug)
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

try_model_fun <- function(code)
{
  temp <- try(code, silent = TRUE)
  if (class(temp)=="try-error") temp <- NA
  return(temp)
}

# Fit curves for drugs in list drug_names
DRC_bunch_MISIS <- function(df, drug_names, controls,
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
    drug <- Normalization(drug, controls)
    drug <- RmOutliers(drug)
    statistics <- DRC_MISIS(df=drug, normilized=normilized,
                      step_dose=step_dose,
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
                       from, to, response=c(50))
{
  drug <- Subset_MISIS(df, name)
  if(normalized==TRUE)
  {
    drug <- Normalization(drug, controls=controls)
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
CC50_slope_bunch_MISIS <- function(df, controls,
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
                         controls=controls,
                         from=boundaries$from[i], to=boundaries$to[i],
                         response=boundaries$response[i])
      results <- rbind(results, temp)
    }
    
  }
  return(results)
}

# Function callings
data1 <- ImportDataFile_MISIS(path_data=path_data)
data1 <- SubstractBackground_MISIS(data1, 490, 700)
data1 <- DropBlank_MISIS(data1)
drug_names <- unique(data1$Образец)
drug_names
conc_info <- list(drug=drug_names,
             stock_conc=c(20.08, 20.02, 20.08, 19.99, 20, 20),
             first_dilution=rep(200, 6),
             step_dilution=c(3, 3, 2, 2, 2, 3),
             n_dilutions=rep(8, 6),
             n_replicates=rep(3, 6))

data1 <- AddConcentrations_MISIS(data1, conc_info)

controls_dil2 <- RmOutliersFromControl(Subset_MISIS(data1, "DMSO-dil2"))
controls_dil3 <- RmOutliersFromControl(Subset_MISIS(data1, "DMSO-dil3"))


# Draft
df <- drug



drug <- Subset_MISIS(data1, "DG605k")
Plot(drug)

drug <- Normalization(drug, controls_dil2)
Plot(drug, x=drug$C_mkM, y=drug$D555_N)

result <- DRC_MISIS(df=drug, normilized=TRUE,step_dose=0.02, X=50,
    plot=TRUE, save_plot=FALSE, path_export=".", need_CCX=TRUE)

summary <- DRC_bunch_MISIS(df=data1, drug_names=c("DG4ClSe","DG603k","DG605k","DG618k"),
                           controls=controls_dil2,
                           normilized=TRUE, step_dose=0.02, X=50,
                           path_export=".", export=FALSE,
                           plot=TRUE, save_plot=FALSE, need_CCX=TRUE)

CC50_slope_MISIS(data1, name="DG605k",
                 controls=controls_dil2, normalized=TRUE, from=6, to=10)

boundaries <- list(name=c("DG4ClSe","DG603k","DG605k","DG618k"),
                   from=c(1, 6, 7, 4),
                   to=c(6, 9, 11, 8),
                   response=rep(50, 4))

CC50s <- CC50_slope_bunch_MISIS(df=data1, controls=controls_dil2,
                                boundaries=boundaries, normalized=TRUE,
                                exclude=c())
# Construct final table and export it
final_table <- cbind(summary, CC50s)

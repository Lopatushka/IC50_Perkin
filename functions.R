library(drc)
library(outliers)
library(compare)

# Paths
path_data1 <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/mcf7 lena.xls"
path_data2 <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/mcf7 lena+skv.xls"
path_names <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/names.xlsx"
path_conc <- "C:/Users/acer/Desktop/Work/Data/MTT/SKV/200519/concentrations.xlsx"

ImportDataFile <- function(path_data)
{
  library(readxl)
  df <- read_excel(path_data)
  colnames(df)[6] <- "D555"
  df$D555 <- sapply(df$D555, as.double)
  return(df)
}

CombineTwoDataFiles <- function(df1, df2)
{
  for (i in 1:nrow(df2))
  {
    df2$Plate[i] <- df2$Plate[i]+ df1$Plate[as.numeric(dim(df1)[1])]
  }
  df <- rbind(df1,df2)
  return(df)
}

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

DropNull <- function(df)
{
  return(as.data.frame(subset(df, df$Drug!='null')))
}

Subset <- function(df, name)
{
  drug <- subset(df, df$Drug == name)
  drug <- drug[, c(6, 9, 8)]
  drug <- drug[order(drug$C_mkM, decreasing = TRUE), ]
  return(drug)
}

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

Normalization <- function(df, controls, n_replicates=3)
{
  df$D555_N <- 100*df$D555/rep(controls, each=n_replicates)
  return(df)
}


# Function callings
data1 <- ImportDataFile(path_data1)
data2 <- ImportDataFile(path_data2)
data <- CombineTwoDataFiles(data1, data2)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(data)
data <- DropNull(data)

sb1 <- Subset(data, "DMSO_1")
sb2 <- Subset(data, "DMSO_2")

drug <- Subset(data, "GK149p")

Plot(sb1)
Plot(sb2)
Plot(drug)

control_medians <- RmOutliersFromControl(sb1)


drug_names <- unique(data$Drug)
#names(summary(lm))
drug$D555_N <- 100*drug$D555/rep(control_medians, each=3)


model <- drm(D555_N ~ C_mkM,
             data = drug[, c(2,4)],
             fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))

summary(model)
plot(model)

median(drug$D555_N)

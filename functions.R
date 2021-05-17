# Paths
path_data1 <- "C:/Users/User/Documents/Work/data/MTT/SKV/200519/mcf7 lena.xls"
path_data2 <- "C:/Users/User/Documents/Work/data/MTT/SKV/200519/mcf7 lena+skv.xls"
path_names <- "C:/Users/User/Documents/Work/data/MTT/SKV/200519/names.xlsx"
path_conc <- "C:/Users/User/Documents/Work/data/MTT/SKV/200519/concentrations.xlsx"

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
  concentrations[2:nrow(concentrations), ] <- log10(concentrations[2:nrow(concentrations),])
  
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
  
  
  df$C_mkM <- sapply(df$C_mkM[df$Drug == 'null'], as.double)
  return(df)
}

DropNull <- function(df)
{
  return(as.data.frame(subset(df, df$Drug!='null')))
}

# Function callings
data1 <- ImportDataFile(path_data1)
data2 <- ImportDataFile(path_data2)
data <- CombineTwoDataFiles(data1, data2)
data <- AddDrugNames(data, path_names)
data <- AddConcentrations(data)
data <- DropNull(data)


# DRC
library(drc)
library(boot)
library(lmtest)
library(metafor)
library(sandwich)

drug_names <- unique(data$Drug)

drug <- subset(data, data$Drug == "GK149p")
drug <- drug[, c(6, 9)]
drug <- drug[order(drug$C_mkM, decreasing = TRUE), ]
plot(x=drug$C_mkM, y=drug$D555, xlab='Log10[C], mkM', ylab='D555')


model <- drm(D555 ~ C_mkM,
             data = drug,
             fct = LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "ED50")))

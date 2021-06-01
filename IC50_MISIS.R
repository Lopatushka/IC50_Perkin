# Libraries
library(readxl)

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
Plot_MISIS <- function(df, x=df$C_mkM, y=df$D555, log=TRUE)
{
  if(log==TRUE)
  {
    x=log10(x)
  }
  plot(x=x,
       y=y,
       xlab='Log10[C], mkM', ylab='D555', main=df[1, 3])
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

drug <- Subset_MISIS(data1, "DG605k")

Plot(drug)



# Draft
df <- drug
name <- "DG4ClSe"

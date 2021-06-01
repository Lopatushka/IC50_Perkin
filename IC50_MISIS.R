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

# Function callings
data1 <- ImportDataFile_MISIS(path_data=path_data)
data1 <- SubstractBackground_MISIS(data1, 490, 700)
data1 <- DropBlank_MISIS(data1)



AddConcentrations <- function(df, path_conc, first_dilution=200,step_dilution=3,
                              n_dilutions=8,  n_replicates=3)
{
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
      print(drug)
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
  
  df$C_mkM <- lapply(df$C_mkM, as.double)
  return(df) 
}


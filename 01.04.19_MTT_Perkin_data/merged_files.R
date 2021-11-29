combined <- MergeFilesFun(path="C:/Users/User/Google Диск/R_scipts/IC50_Perkin/01.04.19_MTT_Perkin_data/final_data",
              split_by=" ")


Export_xlsx <- function(df, path="./results.xlsx", append = FALSE)
{
  write_xlsx(df, path, append=append)
}









#' import data from files
#'
#' this function will return one data frame from aa-tRNA-seq-pipeline
#' @return A data frame containing sequencing metrics
#' @export
import_df <- function(sample_ids, replicates, filter_files_by){
  library(fs)
  library(tidyverse)
  library(here)
  library(janitor)

  ### this takes sample ids, replicate number (vector), and a string to filter the files by ###
  
  sample_ids <- sample_ids
  replicates <- replicates 

  # recursively add paths for each file from each directory into tibble 
  files <- crossing(id = sample_ids, repl = sprintf("%02d", replicates)) |>
      mutate(path = map2(id, repl, ~ list.files(fs::path(here(), 'cleaned-data-folders/', glue("{.x}-{.y}")), 
                                              pattern = "charging.cpm.tsv.gz|bcerror.tsv.gz", 
                                              full.names = TRUE))) |>
      unnest(cols = c(path))

  # determine data type by file name
  files <- files |> filter(str_detect(path, filter_files_by))
  
  # recursively go through tibble containing file paths and load each one into the data frame.
  df <- map_dfr(files$path, ~ read_tsv(gzfile(.x)) |> 
      mutate(sample = basename(.x)))
  
  # clean sample name
  if (filter_files_by == "charging.cpm") {
      df <- df |> mutate(sample = str_remove(sample, "\\.charging.cpm\\.tsv\\.gz$")) # Strip file suffix
  } 
  
  if (filter_files_by == "bcerror") {
     df <- df |> mutate(sample = str_remove(sample, "\\.bcerror\\.tsv\\.gz$"))
  }
  
  
  # clean column names
  df <- clean_names(df)

  return(df)
}


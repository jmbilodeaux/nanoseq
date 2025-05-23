#' this file will summarize the cpm data frame by separating out the tRNA reference name and summing the reads together where 'uncharged' was included in the refereence file name.
#' @param cpm_df Counts per million dataframe from import_df()
#' @return clean tRNA reference names and summarize redundant names
#' @export
clean_cpm_df <- function(cpm_df){
  library(tidyverse)

  cleaned_cpm_df <- cpm_df |> 
      # separate reference and drop "uncharged" label
      separate_wider_delim(
          t_rna,
          delim = "-",
          names = c("origin", "trna", "aa", "anticodon", "family", "species"),
          too_few = "align_start",
          too_many = "drop"
      ) |> 
      # combine the references that were labeled uncharged
      group_by(origin, aa, anticodon, family, species, sample) |> 
      summarize(
      charged = sum(charged),
      uncharged = sum(uncharged),
      charged_cpm = sum(charged_cpm),
      uncharged_cpm = sum(uncharged_cpm),
      .groups = 'drop'
      )
  
  return(cleaned_cpm_df)
}
#' @param cleaned_cpm_df counts per million df with cleaned reference names
#' @return a dataframe with isodecoder count per million means (average tRNA species)
#' @export
isodecoder_means <- function(cleaned_cpm_df){
  # get mean cpm and chrg for each isodecoder
  iso_means <- cleaned_cpm_df |> 
    group_by(origin, aa, anticodon, sample) |> 
    summarise(
      charged_cpm = mean(charged_cpm),
      uncharged_cpm = mean(uncharged_cpm),
      total_cpm = charged_cpm + uncharged_cpm,
      pct_charged = ifelse(total_cpm > 0, 100 * (charged_cpm / total_cpm), NA),
      .groups = 'drop'
    ) |> 
    # separate out sample into relevant parts
    separate_wider_delim(
        sample,
        delim = "-",
        names = c("genotype", "timepoint", "infection_status", "replicate")
    ) |> 
      mutate(isodecoder = paste(aa, anticodon, sep = "-"))

  return(iso_means)
}

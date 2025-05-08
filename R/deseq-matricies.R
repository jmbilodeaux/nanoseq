#'this file runs DESEQ2 analysis
#' @param clean_cpm_df count per million dataframe with cleaned reference names
#' @param organism either host or phage. filters relevant colummns to generate single organism dataframe
#' @param metric either total_counts (abunance) or pct_charged (charging) (or just charged reads?)
#' @return generate count matrix from clean_cpm_df
#' @export
count_matrix <- function(clean_cpm_df, organism, metric) {
  
  # Ensure valid organism
  if (!organism %in% c("host", "phage")) {
    stop("Organism must be 'host' or 'phage'")
  }

  # filter and summarize cpm dataframe (handle replicates)
  df <- clean_cpm_df |> 
    dplyr::rename(trna_id = t_rna) |> 
    mutate(
      trna_id = str_remove(trna_id, "-uncharged")
    ) |> 
    # If phage, only keep infected samples
    filter(
      str_starts(trna_id, organism),
      if (organism == "phage") grepl("inf", sample) else TRUE
    ) |> 
    group_by(trna_id, sample) |> 
    summarize(
      charged = sum(charged),
      uncharged = sum(uncharged),
      total_counts = charged + uncharged,
      pct_charged = (charged / total_counts) * 100,
      .groups = 'drop'
    )

  # Pivot to make count matrix
  count_matrix <- df |> 
    select(trna_id, sample, metric) |> 
    pivot_wider(
      names_from = sample, 
      values_from = metric
    ) |> 
    column_to_rownames("trna_id") |> 
    as.matrix()

  return(count_matrix)
}
#' @param count_matrix use count matrix to generate metadata
#' @return metadata dataframe
#' @export
coldata <- function(count_matrix){
  # Extract metadata from the column names
  sample_names <- colnames(count_matrix)
  metadata <- data.frame(
    sample = gsub("_charged|_uncharged", "", sample_names)) |> 
    mutate(
      condition = ifelse(grepl("inf", sample), "infected", "control"),
      time = ifelse(grepl("15", sample), "15", "30"),
      genotype = sub("-.*", "", sample_names),
      replicate = gsub(".*-(\\d+)$", "r\\1", sample)
    )

  # Convert to factors
  metadata$time <- factor(metadata$time, levels = c("15", "30"))
  metadata$condition <- factor(metadata$condition, levels = c("control", "infected"))
  metadata$genotype <- factor(metadata$genotype, levels = c("wt", "ta", "tb", "tbta", "thi"))

  rownames(metadata) <- sample_names

  return(metadata)

}
#' @param count_matrix count matrix
#' @param coldata metadata
#' @param design_formula statistical formula 
#' @return dds dataframe containing differencial expression results
#' @export
run_deseq <- function(count_matrix, coldata, design_formula, min_count = 10) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = as.formula(design_formula)
  )
  dds <- dds[rowSums(counts(dds)) > min_count, ]
  dds <- DESeq(dds)
  return(dds)
}

#' @param dds dds dataframe
#' @param coef_name name of coefficient 
#' @param padj_cutoff 0.05 by default
#' @return results dataframe
#' @export
get_deseq_results <- function(dds, coef_name, padj_cutoff = 0.05) {
  res <- results(dds, name = coef_name)
  sig_res <- subset(res, padj < padj_cutoff)
  return(list(
    res = res,
    sig = sig_res,
    summary = summary(res)
  ))
}

#' @param res_object results from get_res_results
#' @param title title of volcano plot
#' @return volcano plot object
#' @export
vol_plot <- function (res_object, title, pCutoff = 0.05, FCcutoff = 0.5){
  # Convert results to data frames for visualization
  # Create a data frame with rownames converted to a column
  res_df <- as.data.frame(res_object) %>% 
    rownames_to_column("tRNA")

  # Extract just the isodecoder name for cleaner labels
  res_df <- res_df |> mutate(short_name = str_remove(tRNA, "tRNA-"))

  # Create enhanced volcano plot with customizations
  EnhancedVolcano(res_df,
                  lab = res_df$short_name,  # Use short name for cleaner labels
                  x = 'log2FoldChange',
                  y = 'padj',            
                  title = title,
                  subtitle = 'DESeq2 results',
                  caption = paste0('Total tRNAs: ', nrow(res_df)),
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,       
                  pointSize = 3,
                  labSize = 2.5,
                  drawConnectors = T,   # Draw lines to connect points to labels
                  widthConnectors = 0.5,
                  colConnectors = 'grey30',
                  maxoverlapsConnectors = 50,
                  legendPosition = 'bottom',
                  legendLabSize = 10,
                  legendIconSize = 4.0,
                  axisLabSize = 12,
                  # Show more labels, not just significant ones
                  selectLab = res_df$short_name[order(res_df$pvalue)[1:15]]
                  )

}

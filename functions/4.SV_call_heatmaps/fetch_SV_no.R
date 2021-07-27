fetch_SV_no <- function(
  samplenames, 
  SV_dir, 
  patient_df, 
  dilution_df
) {
  
  library(GenomicRanges)
  simplify_id <- dget(paste0(func_dir, "simplify_id.R"))
  
  # load SV results and fetch SV numbers:
  SV_no <- lapply(samplenames, function(x) {
    
    print(x)
    
    SV_results <- readRDS(paste0(SV_dir, x, "/detected_SVs.Rdata"))
    
    res <- data.frame(
      stringent_true_positives = SV_results$high_conf_bp$true_positives$SV_nos,
      stringent_false_positives = SV_results$high_conf_bp$false_positives$SV_nos,
      less_stringent_true_positives = SV_results$low_conf_bp$true_positives$SV_nos,
      less_stringent_false_positives = SV_results$low_conf_bp$false_positives$SV_nos
    )
    colnames(res) <- c(
      "Stringent_true_positives", "Stringent_false_positives",
      "Less_stringent_true_positives", "Less_stringent_false_positives"
    )
    return(res)
    
  })
  names(SV_no) <- samplenames
  
  SV_no <- do.call("rbind", SV_no)
  
  # simplify sample ids:
  SV_no$Sample <- simplify_id(rownames(SV_no), divider = "_")
  
  # separate into patient and dilution SVs:
  return(
    list(
      patient = SV_no[SV_no$Sample %in% patient_df$Sample,],
      dilution = SV_no[SV_no$Sample %in% dilution_df$Sample,]
    )
  )
  
}
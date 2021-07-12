fetch_fusion_no <- function(
  samplenames, 
  fusion_dir, 
  patient_df, 
  dilution_df
) {
  
  library(GenomicRanges)
  
  # load fusion results and fetch fusion numbers:
  fusion_no <- lapply(samplenames, function(x) {
    
    print(x)
    
    fresults <- readRDS(paste0(fusion_dir, x, "/EWSR1_GOI_fusions.Rdata"))
    
    # format list structures:
    if (is.na(fresults$high_conf_bp$true_positives$fusion_nos)) {
      fresults$high_conf_bp$true_positives$fusion_nos$FLI1 <- 
        fresults$high_conf_bp$true_positives$fusion_nos <- list(FLI1 = NA)
    }
    if (is.na(fresults$high_conf_bp$true_positives$fusions)) {
      fresults$high_conf_bp$true_positives$fusions$FLI1 <- 
        fresults$high_conf_bp$true_positives$fusions <- list(FLI1 = NA)
    }
    if (is.na(fresults$high_conf_bp$false_positives$fusion_nos)) {
      fresults$high_conf_bp$false_positives$fusion_nos <- list(FLI1 = NA)
    }
    if (is.na(fresults$high_conf_bp$false_positives$fusion_nos)) {
      fresults$high_conf_bp$false_positives$fusions <- list(FLI1 = NA)
    }
    
    if (is.na(fresults$low_conf_bp$true_positives$fusion_nos)) {
      fresults$low_conf_bp$true_positives$fusion_nos$FLI1 <- 
        fresults$low_conf_bp$true_positives$fusion_nos <- list(FLI1 = NA)
    }
    if (is.na(fresults$low_conf_bp$true_positives$fusions)) {
      fresults$low_conf_bp$true_positives$fusions$FLI1 <- 
        fresults$low_conf_bp$true_positives$fusions <- list(FLI1 = NA)
    }
    if (is.na(fresults$low_conf_bp$false_positives$fusion_nos)) {
      fresults$low_conf_bp$false_positives$fusion_nos <- list(FLI1 = NA)
    }
    if (is.na(fresults$low_conf_bp$false_positives$fusions)) {
      fresults$low_conf_bp$false_positives$fusions <- list(FLI1 = NA)
    }
    
    res <- data.frame(
      stringent_true_positives = fresults$high_conf_bp$true_positives$fusion_nos$FLI1,
      stringent_false_positives = fresults$high_conf_bp$false_positives$fusion_nos$FLI1,
      less_stringent_true_positives = fresults$low_conf_bp$true_positives$fusion_nos$FLI1,
      less_stringent_false_positives = fresults$low_conf_bp$false_positives$fusion_nos$FLI1
    )
    colnames(res) <- c(
      "Stringent_true_positives", "Stringent_false_positives",
      "Less_stringent_true_positives", "Less_stringent_false_positives"
    )
    return(res)
    
  })
  names(fusion_no) <- samplenames
  
  fusion_no <- do.call("rbind", fusion_no)
  
  # simplify sample ids:
  fusion_no$Sample <- gsub(
    "_.*$", "", 
    sub(".*?_(.+)", "\\1", rownames(fusion_no))
  )
    
  
  # separate into patient and dilution fusions:
  return(
    list(
      patient = fusion_no[
        gsub(
          "_[A-Z].*$", "",
          gsub("409_", "", rownames(fusion_no))
        ) %in% patient_df$Sample,
      ],
      dilution = fusion_no[
        gsub(
          "_[A-Z].*$", "",
          gsub("409_", "", rownames(fusion_no))
        ) %in% dilution_df$Sample,
      ]
    )
  )
  
}
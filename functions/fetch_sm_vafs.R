fetch_sm_vafs <- function(sample_df, roi) {
  
  library(plyr)
  print(sample_df$Library_id)
  if (exists("vcf")) rm(vcf)
  
  # read in variants:
  try(
    vcf <- read.table(
      paste0(
        variant_path, "/", sample_df$Library_id, "/", 
        sample_df$Library_id, ".smCounter.anno.vcf" ),
      ), silent = TRUE)
  
  if (exists("vcf")) {
    
    # change TRUE to T where should be nucleotide:
    if (class(vcf$V4) == "logical") {
      vcf$V4 <- revalue(as.character(vcf$V4), c("TRUE" = "T")) 
    }
    
    # keep only those which passed filtering as gr, and add effect size and 
    # vaf columns:
    pass_var <- vcf[vcf$V7 == "PASS",]
    pass_gr <- GRanges(
      seqnames = pass_var$V1,
      ranges = IRanges(start = pass_var$V2, end = pass_var$V2),
      ref = pass_var$V4,
      alt = pass_var$V5,
      qual = pass_var$V6,
      effect_size = sapply(strsplit(pass_var$V8, "\\|"), function(x) x[3]),
      VAF = round(as.numeric(
        gsub("^.*=", "", gsub(",.*$", "", sapply(strsplit(pass_var$V8, ";"), function(x) x[6])))
      )*100, 1),
      info = pass_var$V8
    )
    
    if (sample_df$Treatment == "tumour" | sample_df$Treatment == "100") {
      # remove tumour gDNA variants with VAF == 100 as likely germline:
      pass_gr <- pass_gr[pass_gr$VAF != 100]
    } else {
      # remove ctDNA variants with VAF > 95 as likely germline:
      pass_gr <- pass_gr[pass_gr$VAF <= 95]
    }
    
    # convert effect size to numeric score:
    pass_gr$effect_size <- revalue(as.character(pass_gr$effect_size), 
      c("HIGH" = "3", "MODERATE" = "2", "LOW" = "1","MODIFIER" = "0") )
    
    # keep variant in roi with highest effect and top quality score, 
    # but remove those with 'MODIFIER' effect:
      
    # keep variants within roi:
    pint <- pintersect(pass_gr, roi)
    keep_var <- pint[pint$hit]
    
    # keep variants with highest effect size:
    keep_var <- keep_var[which.max(keep_var$effect_size)]
    top_var <- keep_var[which.max(keep_var$qual)]
    
    # change effect size back:
    top_var$effect_size <- revalue(as.character(keep_var$effect_size), 
      c("3" = "HIGH", "2" = "MODERATE", "1" = "LOW", "0" = "MODIFIER" ) )
    
    # keep best quality variant call:
    if (length(top_var) > 0) {
      if (top_var$effect_size == "MODIFIER" | top_var$effect_size == "LOW") {
        return(GRanges(NULL))
      } else {
        return(top_var)
      }
    } else {
      return(GRanges(NULL))
    }
    
  } else {
    return(GRanges(NULL))
  }
  
}
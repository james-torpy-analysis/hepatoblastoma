fetch_reads <- function(reads) {
  
  # remove supplementary alignments:
  filt_reads <- reads[
    as.numeric(reads$flag) < 2000
  ]
    
  # fetch paired (duplicated) reads:
  res <- list(
    all = reads,
    pairs = filt_reads[
      filt_reads$qname %in% filt_reads$qname[
        duplicated(filt_reads$qname)
      ]
    ]
  )
  
  # store supplementary reads:
  supp_reads <- reads[
    as.numeric(reads$flag) >= 2000
  ]
  
  # fetch unpaired (non-duplicated) reads:
  filt_reads <- filt_reads[
    !(filt_reads$qname %in% filt_reads$qname[
      duplicated(filt_reads$qname)
    ])
  ]
  # fetch unpaired (non-duplicated) supplementary reads:
  supp_reads <- supp_reads[
    !(supp_reads$qname %in% supp_reads$qname[
      duplicated(supp_reads$qname)
    ])
  ]
  
  # add to list and return:
  return(
    c(
      res,
      singles = c(filt_reads, supp_reads)
    )
  )
  
}
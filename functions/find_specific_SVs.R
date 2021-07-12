find_specific_SVs <- function(
  breakpoints, 
  ROI
) {

  # load function and define ROI regions:
  find_SVs <- dget(paste0(func_dir, "find_SVs.R"))

  # find breakpoint overlaps with ROI:
  specific_SVs <- find_SVs(
    query_coord = breakpoints, 
    subject_coord = ROI
  )
  
  # count ranges:
  print("Number of SVs:")
  if (length(specific_SVs) < 1) {
    print(0)
  } else {
    print(length(specific_SVs))
  }
  
  # determine whether joining ranges overlap ROI:
  if (length(specific_SVs) > 0) {

    seqnams <- gsub(
      ":.*$", "", 
      gsub("^.*chr", "chr", specific_SVs$join)
    )
    coord <- as.numeric(
      gsub(
        "[^0-9.-]", "", 
        gsub("^.*chr.*:", "", specific_SVs$join)
      )
    )
    
    join_ranges <- GRanges(
      seqnames = specific_SVs$join_chr,
      ranges = IRanges(
        start = specific_SVs$join_coord, 
        end = specific_SVs$join_coord
      ),
      strand = "*",
      join_chr = seqnames(specific_SVs),
      join_coord = start(specific_SVs)
    )

  } else {
    join_ranges <- GRanges(NULL)
  }
  
  # identify ROI SVs: 
  all_SVs <- list(
    true_positives = list(
      SVs = find_SVs(join_ranges, ROI)
    )
  )
    
  # count ROI SVs:
  if (all(length(all_SVs$true_positives$SVs) < 1)) {
    all_SVs$true_positives$SV_nos <- GRanges(NULL)
  } else {
    all_SVs$true_positives$SV_nos <- length(all_SVs$true_positives$SVs)
  }

  # identify false SVs:
  all_SVs$false_positives$SVs <- find_SVs(
    query_coord = join_ranges, 
    subject_coord = ROI, 
    invert = T
  )
  
  # count ROI SVs:
  if (all(length(all_SVs$false_positives$SVs) < 1)) {
    all_SVs$false_positives$SV_nos <- GRanges(NULL)
  } else {
    all_SVs$false_positives$SV_nos <- length(all_SVs$false_positives$SVs)
  }
  
  return(all_SVs)

}

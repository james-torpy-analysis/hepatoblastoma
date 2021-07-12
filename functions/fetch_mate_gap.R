fetch_mate_gap <- function(split_read_pair) {
  
  library(rtracklayer)
  library(GenomicRanges)
  
  print(
    paste0("Finding contig for read id ", unique(split_read_pair$qname), "...")
  )
  
  if (length(unique(seqnames(split_read_pair))) == 1) {
    
    if (length(split_read_pair) == 2) {
      # find gap between the reads:
      mgap <- gaps(split_read_pair)
      return(mgap[start(mgap) != 1])
    } else if (length(split_read_pair) == 1) {
      print(
        paste0(
          split_read_pair$qname, " does not have a corresponding mate, removing..."
        )
      )
      return(NA)
    } else if (length(split_read_pair) > 2) {
      print(
        paste0(split_read_pair$qname, " has more than one corresponding mate, removing...")
      )
      return(NA)
    }
    
  } else {
    print(
      paste0(split_read_pair$qname, " has mates in different chromosomes, removing...")
    )
    return(NA)
  }
  
}
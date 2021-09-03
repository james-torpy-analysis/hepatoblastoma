fetch_mate_contig <- function(split_read) {
  
  library(rtracklayer)
  library(GenomicRanges)
  
  print(paste0("Finding contig for read id ", split_read$qname, "..."))
  
  if (length(unique(seqnames(split_read))) == 1) {
    
    # ensure granges is in the correct order:
    split_read <- sortSeqlevels(split_read)
    split_read <- sort(split_read)
    
    if (length(split_read) == 2) {
      contig <- GRanges(
        seqnames = unique(seqnames(split_read)),
        ranges = IRanges(start = start(split_read[1]), end = end(split_read[2])),
        strand = "*",
        qname = unique(split_read$qname)
      )
      return(contig)
    } else if (length(split_read) == 1) {
      print(
        paste0(
          split_read$qname, " does not have a corresponding mate, removing..."
        )
      )
      return(NA)
    } else if (length(split_read) > 2) {
      print(
        paste0(split_read$qname, " has more than one corresponding mate, removing...")
      )
    }
    
  } else {
    print(
      paste0(split_read$qname, " has mates in different chromosomes, removing...")
    )
  }
  
}
filter_overlaps <- function(reads, min_overlap, chromosomes) {
  
  if (any(chromosomes %in% seqnames(reads))) {
    # split by chromosome:
    split_reads <- split(reads, seqnames(reads))
    split_reads <- split_reads[chromosomes]
    
    # calculate lengths from start of read to SV, and from SV to end:
    split_reads <- lapply(split_reads, function(y) {
      if (length(y) > 0) {
        if (length(unique(seqnames(y))) == 1) {
          y$start_to_SV <- y$SV_coord - start(y)
          y$SV_to_end <- end(y) - y$SV_coord
        }
      }
      return(y)
    })
    
    if (length(split_reads) == 2) {
      # concatentate reads back together:
      prefilt_reads <- c(split_reads[[1]], split_reads[[2]])
    } else {
      prefilt_reads <- split_reads[[1]]
    }
    
    # split by qname:
    prefilt_split <- split(prefilt_reads, prefilt_reads$qname)
    
    # filter out reads without overlaps of at least min_overlap on both sides of
    # SV (need at least 1/2 read to satisfy this condition):
    filt_split <- lapply(prefilt_split, function(x) {
      remove_vec <- !(
        x$start_to_SV >= min_overlap & x$SV_to_end >= min_overlap
      )
      remove_vec[is.na(remove_vec)] <- FALSE
      if (any(remove_vec)) {
        return(NULL)
      } else {
        return(x)
      }
    })
    
    res <- unlist(
      as(
        filt_split[sapply(filt_split, function(x) !is.null(x))],
        "GRangesList"
      )
    )
    names(res) <- NULL
    
    return(res)
    
  } else {
    return(reads)
  }
 
  
}
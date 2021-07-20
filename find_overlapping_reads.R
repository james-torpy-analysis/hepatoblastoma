find_overlapping_reads <- function(
  read_pair,
  SV,
  alt_coords
) {
  
  # if finding overlaps with alternate SV coord, make this coord main one:
  if (alt_coords) {
    
    SV <- GRanges(
      seqnames = SV$join_chr,
      ranges = IRanges(start = SV$join_coord, end = SV$join_coord),
      strand = "*",
      join_chr = seqnames(SV),
      join_coord = start(SV)
    )
    
  }
  
  # find overlaps with SV:
  olaps <- findOverlaps(read_pair, SV)
  
  if (length(olaps) > 0) {
    
    # add overlapping breakpoint information:
    read_pair$SV_chr <- NA
    read_pair$SV_coord <- NA
    read_pair$SV_mate_chr <- NA
    read_pair$SV_mate_coord <- NA
    
    read_pair[queryHits(olaps)]$SV_chr <- seqnames(SV)
    read_pair[queryHits(olaps)]$SV_coord <- start(SV)[
      subjectHits(olaps)
    ]
    
    read_pair[queryHits(olaps)]$SV_mate_chr <- as.character(
      unique(SV$join_chr)
    )
    read_pair[queryHits(olaps)]$SV_mate_coord <- SV$join_coord[
      subjectHits(olaps)
    ]
    
  } else {
    read_pair <- NULL
  }
  
  return(read_pair)
  
}
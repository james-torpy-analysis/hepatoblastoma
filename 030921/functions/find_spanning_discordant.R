find_spanning_discordant <- function(
  split_read_pair,
  fusions,
  exp_window = 200
) {
  
  # separate chr11 and chr2 fusion bps:
  chr22_fusions <- fusions
  chr11_fusions <- GRanges(
    seqnames = fusions$join_chr,
    ranges = IRanges(start = fusions$join_coord, end = fusions$join_coord),
    strand = "*",
    join_chr = seqnames(fusions),
    join_coord = start(fusions)
  )
  
  # make breakpoint windows:
  chr11_windows <- resize(chr11_fusions, exp_window+1)
  chr22_windows <- resize(chr22_fusions, exp_window+1)
  
  # remove strand info from read pairs:
  temp_pair <- split_read_pair
  strand(temp_pair) <- "*"
  
  # find overlaps of read mates:
  olaps <- list(
    chr11 = findOverlaps(temp_pair, chr11_windows),
    chr22 = findOverlaps(temp_pair, chr22_windows)
  )
  
  if (all(sapply(olaps, length) > 0)) {
    
    # annotate the chr11 mate with breakpoints:
    split_read_pair$fusion_chr <- NA
    split_read_pair$fusion_coord <- NA
    split_read_pair$fusion_mate_chr <- NA
    split_read_pair$fusion_mate_coord <- NA
    
    split_read_pair$fusion_chr[queryHits(olaps$chr11)] <- seqnames(chr11_fusions)[
      subjectHits(olaps$chr11)
    ]
    split_read_pair$fusion_coord[queryHits(olaps$chr11)] <- start(chr11_fusions)[
      subjectHits(olaps$chr11)
    ]
    
    # annotate the chr11 mate with mate breakpoints:
    split_read_pair$fusion_mate_chr[queryHits(olaps$chr11)] <- as.character(
      chr11_fusions$join_chr[subjectHits(olaps$chr11)]
    )
    split_read_pair$fusion_mate_coord[queryHits(olaps$chr11)] <- as.character(
      chr11_fusions$join_coord[
        subjectHits(olaps$chr11)
      ]
    )
    
    # annotate the chr22 mate with breakpoints:
    split_read_pair$fusion_chr[queryHits(olaps$chr22)] <- seqnames(chr22_fusions)[
      subjectHits(olaps$chr22)
    ]
    split_read_pair$fusion_coord[queryHits(olaps$chr22)] <- start(chr22_fusions)[
      subjectHits(olaps$chr22)
    ]
    
    # annotate the chr22 mate with mate breakpoints:
    split_read_pair$fusion_mate_chr[queryHits(olaps$chr22)] <- as.character(
      chr22_fusions$join_chr[subjectHits(olaps$chr22)]
    )
    split_read_pair$fusion_mate_coord[queryHits(olaps$chr22)] <- as.character(
      chr22_fusions$join_coord[
        subjectHits(olaps$chr22)
      ]
    )
    
    # check if reads on positive strand are to the left of the
    # corresponding breakpoint, allowing for non-split breakpoint overlaps:
    split_read_pair$correct_orientation <- FALSE
    split_read_pair$correct_orientation[
      as.logical(strand(split_read_pair) == "+" & 
                   (split_read_pair$fusion_coord - start(split_read_pair)) > 0)
    ] <- TRUE
    
    # check if reads on positive strand are to the left of the
    # corresponding breakpoint, allowing for non-split breakpoint overlaps:
    split_read_pair$correct_orientation[
      as.logical(strand(split_read_pair) == "-" & 
                   (end(split_read_pair) - split_read_pair$fusion_coord) > 0)
    ] <- TRUE
    
    if (all(split_read_pair$correct_orientation)) {
      return(split_read_pair)
    } else {
      return(NULL)
    }
    
  } else {
    return(NULL)
  }
  
}
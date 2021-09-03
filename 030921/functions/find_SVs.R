### find overlaps between gene SV breakpoints and other genes ###
find_SVs <- function(query_coord, subject_coord, invert = FALSE) {

  if (length(query_coord) > 0) {

    olaps <- findOverlaps(query_coord, subject_coord) 

    if (invert) {

      SV <- query_coord[-queryHits(olaps)]
      if (length(SV) < 1) {
        SV <-GRanges(NULL)
      }

    } else {

      if (length(olaps) > 0) {
        SV <- query_coord[queryHits(olaps)]
      } else {
        SV <-GRanges(NULL)
      }

    }
    
  } else {
    SV <-GRanges(NULL)
  }

  if (length(SV) < 1) {
    # remove duplicates:
    SV <- SV[!(duplicated(start(SV)) & duplicated(SV$join_coord))]
  }

  return(SV)
  
}
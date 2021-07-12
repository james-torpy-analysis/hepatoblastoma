create_gviz <- function(
  bps,
  gene_df,
  start_coord,
  end_coord,
  genome_build = "hg19"
) {

  # create genomic coordinate track:
  gtrack <- GenomeAxisTrack()

  # create chromosome ideogram track:
  itrack <- IdeogramTrack(
    genome = genome_build, 
    chromosome = unique(as.character(gene_df$chromosome))
  )

  # split breakpoint tracks by tag:
  split_bps <- split(bps, bps$filter)

  # create breakpoint tracks:
  bp_tracks = lapply(split_bps, function(x) {

  	if (unique(x$filter) == "PASS") {
  	  bp_col <- "#7CBA61"
  	} else if (unique(x$filter) == "LOWMAPQDISC") {
  	  bp_col <- "#4A70D1"
  	} else {
  	  bp_col <- "#E7298A"
  	}
    return(
      AnnotationTrack(
        x, 
        name = "SV bp",
        col = bp_col,
        fill = bp_col,
        background.title = bp_col,
        genome = genome_build, 
        chromosome = unique(as.character(gene_df$chromosome))
      )
    )
  })

   # create protein coding gene model:
  gr_track <- GeneRegionTrack(
    gene_df, 
    genome = genome_build, 
    chromosome = unique(as.character(gene_df$chromosome)),
    name = unique(as.character(gene_df$symbol)),
    transcriptAnnotation = "symbol",
    background.title = "#F9D480"
  )

  # return plot paramaters:
  return(
    list(
      gviz = append(
        list(
          itrack = itrack,
          gtrack = gtrack,
          gr_track = gr_track
        ),
        bp_tracks
      ),
      start_coord = start_coord,
      end_coord = end_coord
    )
  )

}
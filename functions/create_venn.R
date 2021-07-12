create_venn <- function(read_list, venn_cols) {

  library(ggvenn)
  
  # create venn diagrams of read breakdowns:
  read_vector <- lapply(read_list, function(x) x$qname)
  venn_vector <- read_vector[order(names(read_vector))]
  
  return(
    ggvenn(
      venn_vector, 
      fill_color = venn_cols[1:length(venn_vector)],
      stroke_size = 0.5, set_name_size = 4
    )
  )
  
}
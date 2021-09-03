simplify_id <- function(sample_id, divider) {
  return(
    gsub(
      paste0(divider, ".*$"), "", 
      sub(paste0(".*?", divider, "(.+)"), "\\1", sample_id)
    )
  )
}
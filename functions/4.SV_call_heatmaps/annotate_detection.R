annotate_detection <- function(column, label) {
  column[suppressWarnings(as.numeric(column)) > 0] <- label
  column[suppressWarnings(as.numeric(column)) == 0] <- paste0("no_", label)
  return(column)
}
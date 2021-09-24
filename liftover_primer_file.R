projectname <- "hepatoblastoma_240921"

home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects", projectname)
ref_dir <- file.path(project_dir, "refs")

library(GenomicRanges)
library(rtracklayer)

# read in original primer file:
orig <- read.table(
  file.path(ref_dir, "CDHS-33412Z-324.primers_orig_hg19.txt"),
  header = FALSE,
  sep = "\t" )
orig$strand <- "+"
orig$strand[orig$V3 == 0] <- "-"

# convert to GRanges:
gr <- GRanges(
  seqnames = orig$V1,
  ranges = IRanges(start = orig$V2, end = orig$V2),
  strand = orig$strand )

# retrieve hg19 to 38 chain:
ch = import.chain(file.path(ref_dir, "hg19ToHg38.over.chain"))
gr19 = unlist(liftOver(gr, ch))

# convert to bedfile format:
hg38 <- as.data.frame(gr19)
hg38$strand_no <- 1
hg38$strand_no[hg38$strand == "-"] <- 0
hg38 <- subset(hg38, select = -c(end, width, strand))
hg38$seq <- orig$V4

# save new primer file:
write.table(
  hg38,
  file.path(ref_dir, "CDHS-33412Z-324.primers.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE )

args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
min_overlap_R1 <- args[3]
min_overlap_R2 <- args[4]
#projectname <- "ewing_ctDNA"
#samplename <- "409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001"
#min_overlap_R1 <- 19
#min_overlap_R2 <- 19

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
VAF_dir <- paste0(result_dir, "VAF_calculation/", samplename, "/")
in_dir <- paste0(VAF_dir, "Rdata/")

library(GenomicAlignments)

# load filtered reads:
freads <- readRDS(paste0(in_dir, "filtered_reads.Rdata"))
R1 <- freads$R1
R2 <- freads$R2

# define ROI:
ROI <- list(
  first = GRanges(
    seqnames = "chr22",
    ranges = IRanges(start = 29663998, end = 29696515),
    strand = "*"
  ),
  second = GRanges(
    seqnames = "chr11",
    ranges = IRanges(start = 128554430, end = 128685162),
    strand = "*"
  )
)

# find fusion supporting reads with read 1 in first ROI:
supp1 <- intersect(
  names(which(sum(width(pintersect(R1, ROI$first))) >= min_overlap_R1)),
  names(which(sum(width(pintersect(R2, ROI$second))) >= min_overlap_R2)))

# find fusion supporting reads with read 1 in second ROI:
supp2 <- intersect(
  names(which(sum(width(pintersect(R1, ROI$second))) >= min_overlap_R1)),
  names(which(sum(width(pintersect(R2, ROI$first))) >= min_overlap_R2)))

# merge and remove duplicates:
supp <- c(supp1, supp2)
supp <- supp[!duplicated(supp)]

# save supporting read no:
write.table(
  length(supp),
  paste0(VAF_dir, "non_specific_fusion_supporting_reads.tsv"),
  row.names = F,
  col.names = T,
  quote = F
)



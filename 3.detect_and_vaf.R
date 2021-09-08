
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
min_bp_coverage <- args[3]
#projectname <- "hepatoblastoma"
#samplename <- "324_056_D9YWF_TAAGGCGA-CTCTCTAT_L001"
#min_bp_coverage <- 20


home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
bam_dir <- paste0(result_dir, "picard/")
out_path <- paste0(result_dir, "detection_and_VAF/", samplename, "/")


Robject_dir <- paste0(out_path, "Rdata/")
table_dir <- paste0(out_path, "tables/")
out_bam_dir <- paste0(out_path, "bams/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", out_bam_dir))

func_dir <- paste0(project_dir, "scripts/functions/")

library(GenomicAlignments)
library(Rsamtools)
library(scales)

source(paste0(func_dir, "vaf_functions.R"))

min_overlap_R1 <- 19
min_overlap_R2 <- 19

CTNNB1_gr <- GRanges(
  seqname = "chr3",
  ranges = IRanges(start = 41199422, end = 41240445),
  strand = "*"
)


## 1) read bam file

file_bam <- paste0(
    bam_dir,
    samplename, "/", samplename, ".dedup.sorted.by.coord.bam" )

what <- c(
    "qname",
    "flag",
    "rname",
    "strand",
    "pos",
    "qwidth",
    "mapq",
    "cigar",
    "mrnm",
    "mpos",
    "isize",
    "seq",
    "qual" )

param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = F),
    what = what )

ga <- readGAlignments(file_bam, use.names = TRUE, param = param)
gr <- granges(ga, use.mcols = TRUE)

# only keepreads within CTNNB1 for efficiency:
gr <- gr[queryHits(findOverlaps(gr, CTNNB1_gr, type="within"))]
mcols(gr) <- subset(mcols(gr), select = -strand)

## check all reads paired
stopifnot(all(mcols(gr)$flag %% 2 >= 1))

## check no unmapped reads
stopifnot(all(mcols(gr)$flag %% 8 < 4))

## check no multi-mappers
stopifnot(all(mcols(gr)$flag %% 512 < 256))

## add whether alignment is for first or second read in template
mcols(gr)$R1 <- mcols(gr)$flag %% 128 >= 64
mcols(gr)$R2 <- mcols(gr)$flag %% 256 >= 128
stopifnot(all(mcols(gr)$R1 + mcols(gr)$R2 == 1L))

## fetch read pairs with supp alignments:
mcols(gr)$supp <- mcols(gr)$flag >= 2048
supp_gr <- gr[mcols(gr)$supp]
filt_gr <- gr[gr$qname %in% supp_gr$qname]

# keep read pairs with only 1-2 supplementary alignments, all mapped to same chromosome:
spl_gr <- split(filt_gr, filt_gr$qname)
split_gr <- unlist(as(spl_gr[
  lapply(spl_gr, function(x) {
    (length(x) == 3 | length(x) == 4) & length(unique(x$rname)) == 1
  }) 
], "GRangesList"))
names(split_gr) <- NULL

# save as RDS:
saveRDS(split_gr, paste0(Robject_dir, "CTNNB1_filtered_split_reads.rds"))

# isolate and keep non-split reads:
non_split_gr <- gr[!(gr$qname %in% split_gr$qname)]
names(non_split_gr) <- NULL
length(split_gr) + length(non_split_gr) == length(gr)

saveRDS(non_split_gr, paste0(Robject_dir, "CTNNB1_filtered_non_split_reads.rds"))


## 2) identify deletion breakpoints from split reads:

if (length(split_gr) > 0) {

  # split by qname:
  spl_gr <- split(split_gr, split_gr$qname)
  
  # determine deletion ranges from primary and supp alignments
  deletions <- lapply(spl_gr, function(x) {
    
    # split R1 and 2:
    read_spl <- split(x, x$R1)
  
    deletion <- lapply(read_spl, function(y) {
      # if supp present, determine breakpoints depending on orientation:
      if (any(y$supp)) {
        if (end(y[y$supp]) > end(y[!y$supp])) {
          bps <- c(end(y[!y$supp])-width(y[y$supp]), start(y[y$supp]))
        } else {
          bps <- c(end(y[y$supp]), start(y[!y$supp])+width(y[y$supp]))
        }
        return(GRanges(
          seqnames = seqnames(y)[1],
          ranges = IRanges(start = bps[1], end = bps[2]),
          qname = x$qname[1]
        ))
      }
    })
    return(deletion[sapply(deletion, function(y) length(y) > 0)][[1]])
  
  })
  initial_deletions <- unlist(as(deletions, "GRangesList"))
  
  # remove duplicate deletions:
  deletions <- initial_deletions[!duplicated(initial_deletions)]
  
  # separate starts and ends:
  deletion_starts <- GRanges(
    seqnames = seqnames(deletions),
    ranges = IRanges(start = start(deletions), width = 1),
    strand = "*"
  )
  deletion_ends <- GRanges(
    seqnames = seqnames(deletions),
    ranges = IRanges(start = end(deletions), width = 1),
    strand = "*"
  )
  
  
  ## 3) calculate VAFs for all deletions
  
  # isolate read 2 in forward orientation for upstream breakpoint coords, 
  # as it contains gene-speific primer:
  fwd_split <- split_gr[split_gr$R2 & strand(split_gr) == "+"]
  fwd_non_split <- non_split_gr[non_split_gr$R2 & strand(non_split_gr) == "+"]
  
  # isolate read 2 in reverse orientation for downstream breakpoint coords, 
  # as it contains gene-speific primer:
  rev_split <- split_gr[split_gr$R2 & strand(split_gr) == "-"]
  rev_non_split <- non_split_gr[non_split_gr$R2 & strand(non_split_gr) == "-"]
  
  deletions$upstream_supporting <- NA
  deletions$upstream_non_supporting <- NA
  deletions$downstream_supporting <- NA
  deletions$downstream_non_supporting <- NA
  deletions$upstream_VAF <- NA
  deletions$downstream_VAF <- NA
  deletions$combined_VAF <- NA
  
  for (i in seq_along(deletions)) {
  
    # fetch qnames of split reads which overlap either breakpoint coord of 
    # deletion:
    supporting_reads <- list(
      upstream = pintersect(fwd_split, deletion_starts[i]), 
      downstream = pintersect(rev_split, deletion_ends[i]) )
    supporting_reads <- lapply(supporting_reads, function(x) {
      unique(x$qname[x$hit])
    })
  
    # fetch qnames of non-split reads which overlap either breakpoint coord of 
    # deletion:
    non_supporting_reads <- list(
      upstream = pintersect(fwd_non_split, deletion_starts[i]), 
      downstream = pintersect(rev_non_split, deletion_ends[i]) )
    non_supporting_reads <- lapply(non_supporting_reads, function(x) {
      unique(x$qname[x$hit])
    })
  
    # fill in numbers:
    deletions$upstream_supporting[i] <- length(supporting_reads$upstream)
    deletions$upstream_non_supporting[i] <- length(non_supporting_reads$upstream)
    deletions$downstream_supporting[i] <- length(supporting_reads$downstream)
    deletions$downstream_non_supporting[i] <- length(non_supporting_reads$downstream)
    
    # calculate deletion VAF:
    deletions$upstream_VAF[i] <- deletions$upstream_supporting[i] /
      (deletions$upstream_supporting[i] + deletions$upstream_non_supporting[i])
    deletions$downstream_VAF[i] <- deletions$downstream_supporting[i]/
      (deletions$downstream_supporting[i] + deletions$downstream_non_supporting[i])
    deletions$combined_VAF[i] <- (deletions$upstream_supporting[i]+deletions$downstream_supporting[i])/
      (deletions$upstream_supporting[i]+deletions$downstream_supporting[i] + 
      (deletions$upstream_non_supporting[i]+deletions$downstream_non_supporting[i])/2)
  
  }
  names(deletions) <- NULL

  # format:
  upstream_VAFs <- subset(
    as.data.frame(deletions), 
    select = -c(strand, qname, downstream_supporting, downstream_non_supporting, 
      downstream_VAF, combined_VAF ) )
  colnames(upstream_VAFs) <- gsub("upstream_", "", colnames(upstream_VAFs))
  upstream_VAFs$type <- "upstream"

  downstream_VAFs <- subset(
    as.data.frame(deletions), 
    select = -c(strand, qname, upstream_supporting, upstream_non_supporting, 
      upstream_VAF, combined_VAF ) )
  colnames(downstream_VAFs) <- gsub("downstream_", "", colnames(downstream_VAFs))
  downstream_VAFs$type <- "downstream"

  both_VAFs <- rbind(upstream_VAFs, downstream_VAFs)

  # remove entires with less than min_bp_coverage total reads:
  filtered_VAFs <- both_VAFs[
    both_VAFs$supporting + both_VAFs$non_supporting >= min_bp_coverage, ]
  # keep VAF with max supporting reads:
  selected_VAF <- filtered_VAFs[
    filtered_VAFs$supporting == max(filtered_VAFs$supporting), ]

  colnames(selected_VAF) <- gsub("upstream_", "", colnames(selected_VAF))
  colnames(selected_VAF) <- gsub("downstream_", "", colnames(selected_VAF))
  selected_VAF <- selected_VAF[which.max(selected_VAF$VAF),]
  
  # keep and save deletions and VAFs:
  saveRDS(selected_VAF, paste0(Robject_dir, "selected_VAF.rds"))
  write.table(
    as.data.frame(selected_VAF),
    paste0(table_dir, "selected_VAF.txt"),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t" )
  
  saveRDS(deletions, paste0(Robject_dir, "deletion_VAFs.rds"))
  write.table(
    as.data.frame(deletions),
    paste0(table_dir, "deletion_VAFs.txt"),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t" )

} else {
  # create output file for snakemake:
  selected_VAF <- NULL
  saveRDS(selected_VAF, paste0(Robject_dir, "selected_VAF.rds"))
}



args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
# projectname <- "hepatoblastoma"
# samplename <- "324_057_DB674_CGTACTAG-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
result_dir <- paste0(project_dir, "results/")
bam_dir <- paste0(result_dir, "BWA_and_picard/bams/")
SV_dir <- paste0(result_dir, "SVs/", samplename, "/")
out_path <- paste0(result_dir, "VAF_calculation/", samplename, "/")

Robject_dir <- paste0(out_path, "Rdata/")
hist_dir <- paste0(out_path, "histograms/")
out_bam_dir <- paste0(out_path, "bams/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", hist_dir))
system(paste0("mkdir -p ", out_bam_dir))

func_dir <- paste0(project_dir, "scripts/functions/3.vaf/")

library(GenomicAlignments)
library(Rsamtools)
library(scales)

source(paste0(func_dir, "vaf_functions.R"))

min_overlap_R1 <- 19
min_overlap_R2 <- 19

## 1) read bam file

file_bam <- paste0(
    bam_dir,
    samplename, "/", samplename, ".consensus.bam"
)

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
    "qual")

param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = F),
    what = what)

ga <- readGAlignments(file_bam, use.names = TRUE, param = param)
gr <- granges(ga, use.mcols = TRUE)

## check all reads paired
stopifnot(all(mcols(gr)$flag %% 2 >= 1))

# ## remove reads with segments not properly aligned
# remove <- unique(names(gr)[which(mcols(gr)$flag %% 4 < 2)])
# gr <- gr[which(!names(gr) %in% remove), ]

## check no unmapped reads
stopifnot(all(mcols(gr)$flag %% 8 < 4))

## check no multi-mappers
stopifnot(all(mcols(gr)$flag %% 512 < 256))

## add whether alignment is for first or second read in template
mcols(gr)$R1 <- mcols(gr)$flag %% 128 >= 64
mcols(gr)$R2 <- mcols(gr)$flag %% 256 >= 128
stopifnot(all(mcols(gr)$R1 + mcols(gr)$R2 == 1L))

## extract R1 + R2
## R1 ~ gene-specific primer
## R2 ~ universal primer
tmp  <- gr[mcols(gr)$R1]
mcols(tmp) <- NULL
R1 <- split(tmp, names(tmp))
tmp  <- gr[mcols(gr)$R2]
mcols(tmp) <- NULL
R2 <- split(tmp, names(tmp))

# save as RDS:
saveRDS(list(R1=R1, R2=R2), paste0(Robject_dir, "filtered_reads.Rdata"))


## 2) read SV information

# load in SVs:
both_SVs <- readRDS(paste0(SV_dir, "detected_SVs.Rdata"))

# annotate high and low confidence breakpoints:
try(both_SVs$high_conf_bp$true_positives$SVs$conf <- "high")
try(both_SVs$low_conf_bp$true_positives$SVs$conf <- "low")

# merge high confidence and low confidence breakpoints:
SVs <- c(
  both_SVs$high_conf_bp$true_positives$SVs, 
  both_SVs$low_conf_bp$true_positives$SVs)

# remove SVs with join_coord before start coord:
SVs <- SVs[start(SVs) < SVs$join_coord]

if (length(SVs) >= 1) {
  # remove duplicates:
  SVs <- SVs[
    !duplicated(start(SVs)) & !duplicated(SVs$join_coord)
  ]
  
  ## need to add orientation for identified translocations
  ## if not output by SVABA, can infer from gene annotation
  ## i.e. driver translocation in sense orientation for EWSR1 and SV partner
  strand(SVs) <- "+"
  mcols(SVs)$join_strand <- "+"
  
  for (i in seq_along(SVs)) {
    
    print(i)
    
    SV <- SVs[i]
    
    # calculate length of deletion, to use as read windows:
    read_window_length <- abs(SV$join_coord - start(SV))
    
    gene_b_breakpoint <- GRanges(
      SV$join_chr,
      IRanges(
        SV$join_coord,
        SV$join_coord),
      SV$join_strand)
    gene_a_breakpoint <- SV
    mcols(gene_a_breakpoint) <- NULL
    
    breakpoint <- gsub(":", "_", as.character(gene_a_breakpoint), fixed = TRUE)
    print(breakpoint)
    
    ## consider read_window_length upstream or downstream of EWSR1 breakpoint
    gene_a_upstream <- flank(gene_a_breakpoint, read_window_length, start = TRUE)
    gene_a_upstream <- resize(gene_a_upstream, read_window_length + 1) ## add breakpoint position
    gene_a_dnstream <- flank(gene_a_breakpoint, read_window_length, start = FALSE)
    gene_a_dnstream_rev <- reverseStrand(gene_a_dnstream)
    
    ## consider read_window_length upstream or downstream of SV partner breakpoint
    gene_b_upstream <- flank(gene_b_breakpoint, read_window_length, start = TRUE)
    gene_b_upstream <- resize(gene_b_upstream, read_window_length + 1) ## add breakpoint position
    gene_b_dnstream <- flank(gene_b_breakpoint, read_window_length, start = FALSE)
    gene_b_dnstream_rev <- reverseStrand(gene_b_dnstream)
    
    ## non-supporting reads that satisfy overlap criteria for driver SV
    nonsupp_fwd <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_a_dnstream_rev))) >= min_overlap_R2)))
    ## supporting reads that satisfy overlap criteria for driver SV
    supp_fwd <- intersect(
      names(which(sum(width(pintersect(R1, gene_a_upstream))) >= min_overlap_R1)),
      names(which(sum(width(pintersect(R2, gene_b_dnstream_rev))) >= min_overlap_R2)))
    
    ## diagnost plots for overlap
    R1_nonsupp_fwd <- R1[which(names(R1) %in% nonsupp_fwd)]
    R2_nonsupp_fwd <- R2[which(names(R2) %in% nonsupp_fwd)]
    R1_supp_fwd <- R1[which(names(R1) %in% supp_fwd)]
    R2_supp_fwd <- R2[which(names(R2) %in% supp_fwd)]
    
    pdf(paste0(hist_dir, "hist_overlap_", breakpoint, "_fwd.pdf"))
      par(mfrow = c(2, 2))
      try(
        hist(-sum(width(pintersect(R1_nonsupp_fwd, gene_a_upstream))),
             xlim = c(-180, 0), xlab = "Overlap [bp]",
             main = "SV non-supporting EWSR1 upstream")
      )
      try(
        hist(sum(width(pintersect(R2_nonsupp_fwd, gene_a_dnstream_rev))),
           xlim = c(0, 180), xlab = "Overlap [bp]",
           main = "SV non-supporting EWSR1 dnstream")
      )
      try(
        hist(-sum(width(pintersect(R1_supp_fwd, gene_a_upstream))),
             xlim = c(-180, 0), xlab = "Overlap [bp]",
             main = "SV supporting EWSR1 upstream"))
      try(
        hist(sum(width(pintersect(R2_supp_fwd, gene_b_dnstream_rev))),
             xlim = c(0, 180), xlab = "Overlap [bp]",
             main = "SV supporting FLI1 dnstream"))
    dev.off()
    
    ## calculate VAF
    print(length(nonsupp_fwd))
    print(length(supp_fwd))
    VAF <- length(supp_fwd) / (length(nonsupp_fwd) + length(supp_fwd))
    print(VAF)
    
    
    if (i==1) {
      VAFs <- data.frame(VAF = VAF)
    } else {
      VAFs <- rbind(VAFs, data.frame(VAF = VAF))
    }
    rownames(VAFs)[i] <- paste0("SV_", SVs[i]$join_coord)
    
    ## write reads to SAM for inspection
    writeSam(file_bam, nonsupp_fwd, paste0(out_bam_dir, "/reads_", breakpoint, "_nonsupp_fwd.sam"))
    writeSam(file_bam, supp_fwd, paste0(out_bam_dir, "/reads_", breakpoint, "_supp_fwd.sam"))
   
  }
  
  save.image(paste0(Robject_dir, "VAFs_calculated.Rdata"))
  
  # add SV confidences to and write VAFs
  VAF_df <- do.call("rbind", VAFs)
  VAFs$conf <- SVs$conf
  
  write.table(
    VAFs,
    paste0(out_path, "VAFs.tsv"),
    sep = "\t",
    quote = F,
    row.names = T,
    col.names = T)
  
  saveRDS(VAFs, paste0(Robject_dir, "VAFs.Rdata"))
} else {
  # create dummy file for snakemake:
  system(paste0("touch ", Robject_dir, "VAFs.Rdata"))
}



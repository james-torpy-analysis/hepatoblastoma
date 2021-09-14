
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
min_bp_coverage <- as.numeric(args[3])
min_overlap <- as.numeric(args[4])
min_supporting <- as.numeric(args[5])
#projectname <- "hepatoblastoma"
#samplename <- "324_022_D9HGF_CGAGGCTG-CTCTCTAT_L001"
#min_bp_coverage <- 10
#min_overlap <- 19
#min_supporting <- 2

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

CTNNB1_gr <- GRanges(
  seqname = "chr3",
  ranges = IRanges(start = 41199422, end = 41240445),
  strand = "*" )


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

# only keep reads within CTNNB1 for efficiency:
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

# write samfile:
writeSam(file_bam, filt_gr$qname, 
  paste0(out_bam_dir, "unfiltered_split_read_pairs.sam") )

## fetch read pairs without supp alignments:
non_supp_gr <- gr[!(gr$qname %in% supp_gr$qname)]

# write samfile:
writeSam(file_bam, non_supp_gr$qname, 
  paste0(out_bam_dir, "unfiltered_non_split_read_pairs.sam") )


## 2) identify deletion breakpoints from split reads:

# following doesn't work - get error:
# unable to find an inherited method for function ‘psetdiff’ for signature ‘"CompressedGRangesList", "CompressedGRangesList"’
# spit into GrangesList:
#grl <- split(gr, paste(mcols(gr)$qname, mcols(gr)$R))
# find gaps between primary and supplementary alignments:
#del <- psetdiff(range(grl), grl)
# also tried on each element with lapply, but psetdiff only works on GRanges objects of same length:
#del <- lapply(split_gr, function(x) {
#  return(psetdiff(range(x), x)
#})

# instead, for R1 or R2 with 2 ranges of same strand, did psetdiff for 
# first range of GRanges  elements then pintersected the result with the 
# other range.
# made R1 or R2 GRanges with length != 2 return empty GRanges

# annotate reads:
gr$R <- "R1"
gr$R[gr$R2] <- "R2"
# spit into GrangesList:
grl <- split(gr, paste(gr$qname, gr$R))
# remove entries with supp alignments on different strand to corresponding
# primary alignment:
grl <- lapply(grl, function(x) x[length(unique(strand(x))) == 1])
filt_reads <- unlist(as(grl, "GRangesList"))

# find gaps between primary and supplementary alignments:
del <- lapply(grl, function(x) {
  if (length(x) == 2) {
    diffrange <- psetdiff(range(x), x[1])
    return(psetdiff(diffrange, x[2]))
  } else {
    return(GRanges(NULL))
  }
})

# keep only unique deletion ranges:
del <- unlist(as(del[sapply(del, function(x) length(x) > 0)], "GRangesList"))
strand(del) <- "*"
del <- unique(del)
names(del) <- NULL
# expand del range by 1 to include reads split at breakpoints:
del <- del+1

table(as.character(del))


## 3) calculate VAFs for all deletions

if (length(del) > 0) {

  spl_gr <- split(filt_reads, filt_reads$qname)

  # only keep read pairs containing both primary alignments:
  filt_gr <- spl_gr[sapply(spl_gr, function(x) {
    any(x$R == "R1" & x$flag < 2000) & any(x$R == "R2" & x$flag < 2000)
  })]

  # separate reads with R2 upstream, and those with R2 downstream:
  upstream_R2 <- filt_gr[sapply(filt_gr, function(x) {
    strand(x[x$R2 & x$flag < 2000]) == "+" & 
      end(x[x$R2 & x$flag < 2000]) <= end(x[x$R1 & x$flag < 2000])
  })]

  dnstream_R2 <- filt_gr[sapply(filt_gr, function(x) {
    strand(x[x$R2 & x$flag < 2000]) == "-" & 
      end(x[x$R2 & x$flag < 2000]) >= end(x[x$R1 & x$flag < 2000])
  })]

  # separate split and non-split reads:
  upstream_split <- unlist(as(
    upstream_R2[sapply(upstream_R2, function(x) any(x$flag > 2000))], "GRangesList" ))
  upstream_non_split <- unlist(as(
    upstream_R2[sapply(upstream_R2, function(x) all(x$flag < 2000))], "GRangesList" ))
  dnstream_split <- unlist(as(
    dnstream_R2[sapply(dnstream_R2, function(x) any(x$flag > 2000))], "GRangesList" ))
  dnstream_non_split <- unlist(as(
    dnstream_R2[sapply(dnstream_R2, function(x) all(x$flag < 2000))], "GRangesList" ))

  del_grl <- split(del, paste(start(del), end(del)))
  VAF <- lapply(del_grl, function(x) {

    # separate deletion breakpoints:
    up_del <- x
    end(up_del) <- start(x)
    down_del <- x
    start(down_del) <- end(x)
  
    # fetch split reads which contain deletion breakpoint:
    upstream_supp <- pintersect(upstream_split, up_del)
    upstream_supp <- unique(upstream_supp$qname[upstream_supp$hit])
    dnstream_supp <- pintersect(dnstream_split, down_del)
    dnstream_supp <- unique(dnstream_supp$qname[dnstream_supp$hit])
  
    # define up and downstream regions:
    del_upstream <- flank(x, 1e6, start = TRUE)
    del_dnstream  <- flank(x, 1e6, start = FALSE)
  
    # fetch non-split reads flanking upstream breakpoint:
    upstream_non_supp <- unique(upstream_non_split$qname[intersect(
      which(width(pintersect(upstream_non_split, del_upstream)) >= min_overlap),
      which(width(pintersect(upstream_non_split, x)) >= min_overlap) )])
    dnstream_non_supp <- unique(dnstream_non_split$qname[intersect(
      which(width(pintersect(dnstream_non_split, del_dnstream)) >= min_overlap),
      which(width(pintersect(dnstream_non_split, x)) >= min_overlap) )])
  
    # write to bams:
    writeSam(file_bam, upstream_supp, 
      paste0(out_bam_dir, "deletion_", up_del, "_", down_del, 
        "_upstream_supporting.sam" ))
    writeSam(file_bam, dnstream_supp, 
      paste0(out_bam_dir, "deletion_", up_del, "_", down_del, 
        "_dnstream_supporting.sam" ))
    writeSam(file_bam, upstream_non_supp, 
      paste0(out_bam_dir, "deletion_", up_del, "_", down_del, 
        "_upstream_non_supporting.sam" ))
    writeSam(file_bam, dnstream_non_supp, 
      paste0(out_bam_dir, "deletion_", up_del, "_", down_del, 
        "_dnstream_non_supporting.sam" ))

    # add to deletion record:
    x$upstream_supp <- length(upstream_supp)
    x$upstream_non_supp <- length(upstream_non_supp)
    x$dnstream_supp <- length(dnstream_supp)
    x$dnstream_non_supp <- length(dnstream_non_supp)
    x$total_supp <- x$upstream_supp + x$dnstream_supp

    # calculate VAFs:
    x$up_VAF <- round(x$upstream_supp/(x$upstream_supp + x$upstream_non_supp), 3)
    x$dn_VAF <- round(x$dnstream_supp/(x$dnstream_supp + x$dnstream_non_supp), 3)

    # calculate weighted average VAF:
    up_weight <- round(x$upstream_supp/(x$upstream_supp + x$dnstream_supp), 3)
    dn_weight <- round(x$dnstream_supp/(x$upstream_supp + x$dnstream_supp), 3)
    x$avg_VAF <- round((x$up_VAF*up_weight) + (x$dn_VAF*dn_weight), 3)

    return(x)

  })
  VAF <- unlist(as(VAF, "GRangesList"))

  # make NA values = 0:
  VAF$up_VAF[is.na(VAF$up_VAF)] <- 0
  VAF$dn_VAF[is.na(VAF$dn_VAF)] <- 0
  VAF$avg_VAF[is.na(VAF$avg_VAF)] <- 0


  ## 4) filter and save VAFs:

  # select VAF with max supporting reads:
  selected_VAF <- VAF[VAF$total_supp == max(VAF$total_supp)]

  if (any(!is.na(VAF$avg_VAF))) {
    selected_VAF <- VAF[which.max(VAF$avg_VAF)]
  } else {
    selected_VAF <- VAF[which.max(VAF$VAF)]
  }

  # keep and save deletions and VAFs:
  saveRDS(selected_VAF, paste0(Robject_dir, "selected_VAF.rds"))
  write.table(
    as.data.frame(selected_VAF),
    paste0(table_dir, "selected_VAF.txt"),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t" )

} else {
  # create output file for snakemake:
  selected_VAF <- NULL
  saveRDS(selected_VAF, paste0(Robject_dir, "selected_VAF.rds"))
}


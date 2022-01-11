args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
min_bp_coverage <- as.numeric(args[3])
min_overlap <- as.numeric(args[4])
max_overlap_del <- as.numeric(args[5])

#projectname <- "hepatoblastoma"
#samplename <- "324_054_combined"
#min_bp_coverage <- 10
#min_overlap <- 19
#max_overlap_del <- 10

home_dir <- "/share/ScratchGeneral/jamtor"
project_dir <- file.path(home_dir, "projects", projectname)
func_dir <- file.path(project_dir, "scripts", "functions")
result_dir <- file.path(project_dir, "results")
bam_dir <- file.path(result_dir, "picard")
out_path <- file.path(result_dir, "detection_and_VAF", samplename)

Robject_dir <- file.path(out_path, "Rdata")
table_dir <- file.path(out_path, "tables")
out_bam_dir <- file.path(out_path, "bams")

dir.create(Robject_dir, FALSE, TRUE)
dir.create(table_dir, FALSE, TRUE)
dir.create(out_bam_dir, FALSE, TRUE)

library(GenomicAlignments)
library(Rsamtools)
library(scales)

source(file.path(func_dir, "3.detect_and_vaf_functions.R"))

CTNNB1_gr <- GRanges(
    seqname = "chr3",
    ranges = IRanges(start = 41199422, end = 41240445),
    strand = "*")

## 1) read bam file

file_bam <- file.path(
    bam_dir,
    samplename, paste0(samplename, ".dedup.sorted.by.coord.bam"))

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
mcols(gr)$R <- "R1"
mcols(gr)$R[mcols(gr)$R2] <- "R2"

## remove reads with missing R1 or R2
keep <- intersect(names(gr)[gr$R1], names(gr)[gr$R2])
gr <- gr[names(gr) %in% keep]

## check all reads have primary alignment
stopifnot(all(names(gr) %in% names(gr)[gr$flag < 2048 & gr$R1]))
stopifnot(all(names(gr) %in% names(gr)[gr$flag < 2048 & gr$R2]))

## remove reads with split alignments mapping to different chromosomes or strands
tmp <- split(gr, paste(gr$qname, gr$R))
excl <- unlist(tmp[which(lengths(range(tmp)) > 1)])
excl <- unique(mcols(excl)$qname)
gr <- gr[which(!mcols(gr)$qname %in% excl)]

## only keep reads within CTNNB1 for efficiency
keep <- unique(names(gr)[gr %within% CTNNB1_gr])
gr <- gr[names(gr) %in% keep]

# write to bam:
writeSam(file_bam, gr$qname, file.path(out_bam_dir, "CTNNB1_reads.sam") )

# write split reads to bam for visualisation:
print(length(gr[gr$flag > 2048]))
writeSam(file_bam, gr[gr$flag > 2048]$qname, file.path(out_bam_dir, "filtered_split_CTNNB1_reads.sam") )

## 2) identify deletion breakpoints from split reads

## find gaps between primary and supplementary alignments
grl <- split(gr, paste(gr$qname, gr$R))

## only keep reads within CTNNB1
if (!file.exists(file.path(Robject_dir, "CTNNB1_reads.rds"))) {
  grl <- grl[sapply(grl, function(x) all(x %within% CTNNB1_gr))]
  saveRDS(grl, file.path(Robject_dir, "CTNNB1_reads.rds"))
} else {
  grl <- readRDS(file.path(Robject_dir, "CTNNB1_reads.rds"))
}

# take the range of all elements, and find the difference between
# this and the primary and supp alignments in the elements which have them:
del <- psetdiff(unlist(range(grl)), grl)

# remove empty ranges and convert to gr:
del <- unlist(del)
strand(del) <- "*"

print(max_overlap_del) 

## check deletion supporting reads are well behaved
## i.e. mates not overlapping deletion
del$qname <- sapply(strsplit(names(del), " R(1|2)"), "[", 1)
exclude <- sapply(split(del, names(del)), function (x) {
    read_name <- x$qname
    support_gr <- gr[names(gr) == read_name]
    print(support_gr)
    w <- width(pintersect(support_gr, x))
    print(w)
    if (any(w >= max_overlap_del)) return(read_name)
    else return()
})
del <- del[which(!(del$qname %in% unlist(exclude)))]

## keep only unique deletion ranges
del <- unique(del)
names(del) <- NULL
table(as.character(del))

## 3) calculate VAFs for all deletions

if (length(del) > 0) {

    grl <- split(gr, gr$qname)

    ## upstream breakpoint (A)

    bp_A_grl <- grl[sapply(grl, function(x) {
        r1 <- range(x[x$R1])
        r2 <- range(x[x$R2])
        start(r1) < end(r2) &&
            as.character(strand(r1)) == "+" &&
            as.character(strand(r2)) == "-"
    })]
    bp_A_R1 <- endoapply(bp_A_grl, function (x) x[x$R1])
    bp_A_R1_del <- psetdiff(unlist(range(bp_A_R1)), bp_A_R1)
    bp_A_R2 <- endoapply(bp_A_grl, function (x) x[x$R2])
    bp_A_R2_del <- psetdiff(unlist(range(bp_A_R2)), bp_A_R2)

    ## downstream breakpoint (B)

    bp_B_grl <- grl[sapply(grl, function(x) {
        r1 <- range(x[x$R1])
        r2 <- range(x[x$R2])
        start(r2) < end(r1) &&
            as.character(strand(r1)) == "-" &&
            as.character(strand(r2)) == "+"
    })]
    bp_B_R1 <- endoapply(bp_B_grl, function (x) x[x$R1])
    bp_B_R1_del <- psetdiff(unlist(range(bp_B_R1)), bp_B_R1)
    bp_B_R2 <- endoapply(bp_B_grl, function (x) x[x$R2])
    bp_B_R2_del <- psetdiff(unlist(range(bp_B_R2)), bp_B_R2)

    del_grl <- as(del, "GRangesList")

    VAF <- lapply(del_grl, function(x) {

        ## supporting reads must have R1 or R2 contain deletion
        bp_A_ol_R1 <- findOverlaps(x, bp_A_R1_del,
            type = "equal", ignore.strand = TRUE)
        bp_A_ol_R2 <- findOverlaps(x, bp_A_R2_del,
            type = "equal", ignore.strand = TRUE)
        bp_A_supp <- union(
            names(bp_A_R1_del)[subjectHits(bp_A_ol_R1)],
            names(bp_A_R2_del)[subjectHits(bp_A_ol_R2)])
        bp_B_ol_R1 <- findOverlaps(x, bp_B_R1_del,
            type = "equal", ignore.strand = TRUE)
        bp_B_ol_R2 <- findOverlaps(x, bp_B_R2_del,
            type = "equal", ignore.strand = TRUE)
        bp_B_supp <- union(
            names(bp_B_R1_del)[subjectHits(bp_B_ol_R1)],
            names(bp_B_R2_del)[subjectHits(bp_B_ol_R2)])

        # write bams:
        writeSam(
          file_bam, bp_A_supp, paste0(
            out_bam_dir, "/deletion_", start(x), "_to_",
            end(x),"_upstream_bp_supporting_reads.sam" ) )
        writeSam(
          file_bam, bp_B_supp, paste0(
            out_bam_dir, "/deletion_", start(x), "_to_",
            end(x),"_downstream_bp_supporting_reads.sam" ) )

        ## define up and downstream regions
        bp_A_up <- flank(x, min_overlap, start = TRUE)
        bp_A_dn <- flank(x, -min_overlap, start = TRUE)
        bp_B_up <- flank(x, -min_overlap, start = FALSE)
        bp_B_dn <- flank(x, min_overlap, start = FALSE)

        ## non-supporting reads must have R1 or R2 overlapping breakpoint
        ## with min_overlap upstream and downstream

        ## upstream breakpoint (A)
        bp_A_non_supp_R1 <- intersect(
            names(which(sum(width(pintersect(
                bp_A_R1, bp_A_up, ignore.strand = TRUE))) >= min_overlap)),
            names(which(sum(width(pintersect(
                bp_A_R1, bp_A_dn, ignore.strand = TRUE))) >= min_overlap)))
        bp_A_non_supp_R2 <- intersect(
            names(which(sum(width(pintersect(
                bp_A_R2, bp_A_up, ignore.strand = TRUE))) >= min_overlap)),
            names(which(sum(width(pintersect(
                bp_A_R2, bp_A_dn, ignore.strand = TRUE))) >= min_overlap)))
        bp_A_non_supp <- union(bp_A_non_supp_R1, bp_A_non_supp_R2)

        ## downstream breakpoint (B)
        bp_B_non_supp_R1 <- intersect(
            names(which(sum(width(pintersect(
                bp_B_R1, bp_B_up, ignore.strand = TRUE))) >= min_overlap)),
            names(which(sum(width(pintersect(
                bp_B_R1, bp_B_dn, ignore.strand = TRUE))) >= min_overlap)))
        bp_B_non_supp_R2 <- intersect(
            names(which(sum(width(pintersect(
                bp_B_R2, bp_B_up, ignore.strand = TRUE))) >= min_overlap)),
            names(which(sum(width(pintersect(
                bp_B_R2, bp_B_dn, ignore.strand = TRUE))) >= min_overlap)))
        bp_B_non_supp <- union(bp_B_non_supp_R1, bp_B_non_supp_R2)

        # write bams:
        writeSam(
          file_bam, bp_A_non_supp, paste0(
            out_bam_dir, "/deletion_", start(x), "_to_",
            end(x),"_upstream_bp_non_supporting_reads.sam" ) )
        writeSam(
          file_bam, bp_B_non_supp, paste0(
            out_bam_dir, "/deletion_", start(x), "_to_",
            end(x),"_downstream_bp_non_supporting_reads.sam" ) )

        ## add to deletion record
        x$bp_A_supp <- length(bp_A_supp)
        x$bp_A_non_supp <- length(bp_A_non_supp)
        x$bp_A_total <- length(bp_A_supp) + length(bp_A_non_supp)
        x$bp_B_supp <- length(bp_B_supp)
        x$bp_B_non_supp <- length(bp_B_non_supp)
        x$bp_B_total <- length(bp_B_supp) + length(bp_B_non_supp)
        x$total_supp <- x$bp_A_supp + x$bp_B_supp

        ## calculate VAFs, removing those if total read no or non-supporting
        # read number == 0:
        x$up_VAF <- round(x$bp_A_supp/x$bp_A_total, 4)
        if (x$bp_A_total == 0 | x$bp_A_non_supp == 0) {
          x$up_VAF <- 0
        }
        x$dn_VAF <- round(x$bp_B_supp/x$bp_B_total, 4)
        if (x$bp_B_total == 0 | x$bp_B_non_supp == 0) {
          x$dn_VAF <- 0
        }

        ## calculate weighted average VAF
        up_weight <- x$bp_A_total/(x$bp_A_total + x$bp_B_total)
        dn_weight <- x$bp_B_total/(x$bp_A_total + x$bp_B_total)
        x$avg_VAF <- round((x$up_VAF*up_weight) + (x$dn_VAF*dn_weight), 4)
        if (x$bp_A_total == 0 && x$bp_B_total == 0) x$avg_VAF <- 0

        return(x)

    })

    VAF <- unlist(as(VAF, "GRangesList"))

    # order VAFs according to no supporting reads:
    VAF <- VAF[order(VAF$total_supp)]

    ## keep and save deletions and VAFs
    saveRDS(VAF, file.path(Robject_dir, "/VAF.rds"))
    write.table(
        as.data.frame(VAF),
        file.path(table_dir, "/VAF.txt"),
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t" )

} else {
    ## create output file for snakemake
    VAF <- NULL
    saveRDS(VAF, paste0(Robject_dir, "/VAF.rds"))
}

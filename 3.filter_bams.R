
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
chromosomes <- args[3]
ROI <- args[4]
supps_allowed <- as.numeric(args[5])
SV_type <- args[6]
#projectname <- "hepatoblastoma"
#samplename <- "324_003_DB674_AGGCAGAA-CTCTCTAT_L001"
#chromosomes <- c("chr3")
#ROI <- "41265552_41266717" # region covered by primers
#supps_allowed <- as.numeric("2")
#SV_type <- "deletion"

# split variables into vectors:
chromosomes <- unlist(strsplit(chromosomes, "_")[[1]])
ROI <- as.numeric(
  unlist(
    strsplit(ROI, "_")[[1]]
  )
)

venn_cols <- c("#7C1BE2", "#1B9E77", "#EFC000FF", "blue")

#home_dir <- "/Users/torpor/clusterHome/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")

out_path <- paste0(result_dir, "VAF_calculation/", samplename, 
  "/max_supps_allowed_", supps_allowed, "/")
Robject_dir <- paste0(out_path, "Rdata/")
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))

fusion_dir <- paste0(result_dir, "fusions/", samplename, "/")
bam_path <- paste0(result_dir, "BWA_and_picard/bams/")
col_dir <- paste0(home_dir, "R/colour_palettes/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(cowplot)
library(ggvenn)

create_venn <- dget(paste0(func_dir, "create_venn.R"))
fetch_reads <- dget(paste0(func_dir, "fetch_reads.R"))

plot_cols <- read.table(
  paste0(col_dir, "labelled_colour_palette.txt"),
  sep = "\t",
  header = F,
  comment.char = "",
  fill = TRUE
)$V1
plot_cols <- plot_cols[c(1:3, 5, 4, 6:length(plot_cols))]


####################################################################################
### 1. Load bams ###
####################################################################################

# define file types to import:
bamtypes <- list(
  concordant_pairs = "consensus.concordant.pairs.bam", 
  discordant_pairs = "consensus.discordant.bam",
  split_supp = "consensus.split.bam"
)

# set up scanBamParam to filter out unmapped reads:
param <- ScanBamParam(
  flag = scanBamFlag(isUnmappedQuery = F),
  what = c(
    "rname", "pos", "qwidth", "strand", 
    "qname", "flag", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
  )
)

if (!file.exists(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))) {
  
  # load in bam:
  unfilt_bam <- lapply(bamtypes, function(x) {
    
    bam_obj <- scanBam(
      paste0(bam_path, samplename, "/", samplename, ".", x),
      param = param
    )
    
    # convert to GRanges:
    gr <- GRanges(
      seqnames = bam_obj[[1]]$rname,
      ranges = IRanges(
        start = bam_obj[[1]]$pos, 
        width = bam_obj[[1]]$qwidth
      ),
      strand = bam_obj[[1]]$strand,
      qname = bam_obj[[1]]$qname,
      flag = bam_obj[[1]]$flag,
      mapq = bam_obj[[1]]$mapq,
      cigar = bam_obj[[1]]$cigar,
      rnext = bam_obj[[1]]$mrnm,
      pnext = bam_obj[[1]]$mpos,
      tlen = bam_obj[[1]]$isize,
      seq = bam_obj[[1]]$seq,
      qual= bam_obj[[1]]$qual
    )
    
    return(gr)
    
  })
  
  # calculate read numbers:
  read_numbers <- data.frame(
    unfiltered = sapply(unfilt_bam, length)
  )
  read_numbers["split_pairs",] = NA
  read_numbers["non_split_concordant_pairs",] = NA
  
  
  ####################################################################################
  ### 2a. Filter bams ###
  ####################################################################################

  # remove pairs not mapped to chromosomes of interest:
  chr_filt_bam <- list(
    concordant_pairs = unfilt_bam$concordant_pairs,
    discordant_pairs = unfilt_bam$discordant_pairs,
    split_supp = unfilt_bam$split_supp
  )
  spl <- lapply(chr_filt_bam, function(x) split(x, x$qname))
  
  # initiate cluster:
  cl <- makeCluster(7)
  clusterExport(
    cl, varlist = c("spl", "chromosomes")
  )
  
  # remove concordant pairs not mapped to chromosomes of interest:
  spl$concordant_pairs <- spl$concordant_pairs[
    unlist(
      parLapply(cl, spl$concordant_pairs, function(x) {
        all(seqnames(x) %in% chromosomes)
      })
    )
  ]
  
  # remove discordant pairs not mapped to chromosomes of interest:
  spl$discordant_pairs <- spl$discordant_pairs[
    unlist(
      parLapply(cl, spl$discordant_pairs, function(x) {
        "chr11" %in% seqnames(x) & "chr22" %in% seqnames(x)
      })
    )
  ]
  
  stopCluster(cl)
  
  # remove one mate of double-split read pairs:
  spl$split_supp <- lapply(spl$split_supp, function(x) {
    return(x[1])
  })
  
  # merge each set of reads and separate split supp and primary:
  chr_filt_bam <- lapply(spl, function(x) {
    x <- unlist(
      as(x, "GRangesList")
    )
    names(x) <- NULL
    return(x)
  })
  
  # update read number record:
  read_numbers$non_fusion_chr_removed <- c(
    sapply(chr_filt_bam, length),
    split_pairs = NA,
    non_split_concordant_pairs = NA
  )
  
  temp_bam <- lapply(chr_filt_bam, function(x) {

  	# remove multimapping reads (mapq score = 0):
    mmappers <- unique(x$qname[x$mapq == 0])
    res <- x[!(x$qname %in% mmappers)]
    
    # calculate read numbers:
    mmappers_removed = length(res)
  
    return(
      list(
        bam = res,
        read_no = mmappers_removed
      )
    )

  })


  ####################################################################################
  ### 2b. Filter bams ###
  ####################################################################################

  # update read number record:
  read_numbers$mmappers_removed <- c(
    sapply(temp_bam, function(x) {
      return(x$read_no)
    }),
    split_pairs = NA,
    non_split_concordant_pairs = NA
  )
  
  temp_bam <- lapply(temp_bam, function(x) return(x$bam))

  # filter concordant and discordant pairs:
  paired_bam <- temp_bam[names(temp_bam) %in% c("concordant_pairs", "discordant_pairs")]
  
  filt_bam <- lapply(paired_bam, function(x) {
    
    # remove reads with > supps_allowed supplementary alignments:
    spl <- split(x, x$qname)
    filt_spl <- spl[
      sapply(spl, function(y) {
        return(length(which(y$flag >= 2000)) <= supps_allowed)
      })
    ]
    x <- unlist(filt_spl)
    
    # record read numbers:
    read_no <- list(too_many_supp_removed = length(x))

  	# remove reads with only supplementary alignments:
    spl <- split(x, x$qname)

    filt_spl <- spl[
      !(
        sapply(spl, function(y) {
          return(all(y$flag >= 2000))
        })
      )
    ]
    res <- unlist(filt_spl)
    
    # record read numbers:
    read_no$only_supp_removed <- length(res)

    # remove unpaired reads:
    fetch_reads <- dget(paste0(func_dir, "fetch_reads.R"))
  	singles_vs_pairs <- fetch_reads(res)
  	
  	# record read numbers:
  	read_no$unpaired_removed = length(singles_vs_pairs$pairs)

    return(
      list(
        bam = singles_vs_pairs$pairs,
        read_no = unlist(read_no)
      )
    )

  })
  
  # add read numbers to record:
  temp_read_no <- rbind(
    as.data.frame(
      t(
        sapply(filt_bam, function(x) {
          return(x$read_no)
        })
      )
    ),
    data.frame(
      too_many_supp_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      ), 
      only_supp_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      ),
      unpaired_removed = c(
        read_numbers$mmappers_removed[rownames(read_numbers) == "split_supp"], 
        NA, NA
      )
    )
  )
  rownames(temp_read_no)[3:5] <- c(
    "split_supp", "split_pairs", "non_split_concordant_pairs"
  )
  
  # combine with other read counts:
  read_numbers <- cbind(read_numbers, temp_read_no)

  # combine read objects into list:
  bam <- list(
  	concordant_pairs = filt_bam$concordant_pairs$bam,
  	discordant_pairs = filt_bam$discordant_pairs$bam,
  	split_supp = temp_bam$split_supp
  )
  
  # remove split supps not in either concordant or discordant pair elements:
  bam$split_supp <- bam$split_supp[
    bam$split_supp$qname %in% c(
      bam$concordant_pairs$qname,
      bam$discordant_pairs$qname
    )
  ]
  
  # fetch split primary pairs:
  bam$split_pairs <- c(
  	bam$concordant_pairs[bam$concordant_pairs$qname %in% bam$split_supp$qname],
  	bam$discordant_pairs[bam$discordant_pairs$qname %in% bam$split_supp$qname]
  )
  
  # record read numbers:
  read_numbers$unmatched_split_removed <- c(
    sapply(bam, length),
    non_split_concordant_pairs = NA
  )

  # remove split pairs from discordant pairs:
  bam$discordant_pairs <- bam$discordant_pairs[
  	!(bam$discordant_pairs$qname %in% bam$split_pairs$qname)
  ]
  
  # record read numbers:
  read_numbers$split_removed_from_discordant <- c(
    sapply(bam, length),
    non_split_concordant_pairs = NA
  )
  
  # fetch non_split_concordant_pairs:
  bam$non_split_concordant_pairs <- bam$concordant_pairs[
    !(bam$concordant_pairs$qname %in% bam$split_pairs$qname)
  ]
  
  # remove rownames:
  final_bam <- lapply(bam, function(x) {
    names(x) <- NULL
    return(x)
  })
  
  # record read numbers:
  read_numbers$split_removed_from_discordant <- sapply(final_bam, length)
  
  saveRDS(unfilt_bam, paste0(Robject_dir, "unfiltered_VAF_calculation_reads.Rdata"))
  saveRDS(final_bam, paste0(Robject_dir, "VAF_calculation_reads.Rdata"))
  saveRDS(read_numbers, paste0(Robject_dir, "VAF_calculation_read_nos.Rdata"))
  
} else {

  unfilt_bam <- readRDS(paste0(Robject_dir, "unfiltered_VAF_calculation_reads.Rdata"))
  final_bam <- readRDS(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))
  read_numbers <- readRDS(paste0(Robject_dir, "VAF_calculation_read_nos.Rdata"))

}


####################################################################################
### 3. Fetch and save non-specific SV-supporting reads ###
####################################################################################

# isolate reads within roi for supporting reads on those samples without 
# SV called:
roi <- GRanges(
  seqnames = chromosomes,
  range = IRanges(
    start = ROI[1], 
    end = 
  ),
  strand = "*"
)

if (length(final_bam$discordant_pairs) > 0) {
  # keep roi discordant reads:
  # split by read name to keep pairs together:
  spl <- split(final_bam$discordant_pairs, final_bam$discordant_pairs$qname)
  
  roi_discordant <- lapply(spl, function(x) {

    # return if chr11 read overlaps chr11 region and same for ch22:
    olaps <- findOverlaps(x, roi)
    olap_roi <- roi[subjectHits(olaps)]
    
    if (length(olap_roi) > 1) {
      if (
        as.logical(
          seqnames(olap_roi)[1] == "chr11" & seqnames(olap_roi[2]) == "chr22"
        )
      ) {
        return(x)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
    
  })
  
  # remove NULLs:
  roi_discordant <- roi_discordant[
    sapply(roi_discordant, function(y) !is.null(y))
  ]

  # bind as list:
  roi_discordant <- unlist(
    as(roi_discordant, "GRangesList")
  )
  names(roi_discordant) <- NULL
  
} else {
  roi_discordant <- GRanges(NULL)
}

if (length(final_bam$split_pairs) > 0) {
  
  # fetch roi split reads:
  # split by read name to keep pairs together:
  spl <- split(final_bam$split_pairs, final_bam$split_pairs$qname)
  
  for (i in 1:length(spl)) {
    
    # keep if reads overlap roi:
    olaps <- findOverlaps(spl[[i]], roi)
    olap_roi <- roi[subjectHits(olaps)]
    
    if (length(olap_roi) > 1) {
      
      if (SV_type == "fusion") {
        # add to keep list if chr11 read overlaps chr11 region and same for ch22
        if (
          as.logical(
            seqnames(olap_roi)[1] == "chr11" & seqnames(olap_roi[2]) == "chr22"
          )
        ) {
          
          keep <- TRUE
          
        } else {
          
          supp <- final_bam$split_supp[
            final_bam$split_supp$qname == unique(spl[[i]]$qname)
          ]
          
          if (as.logical(unique(seqnames(olap_roi)) != seqnames(supp))) {
            keep <- TRUE
          } else {
            keep <- FALSE
          }
          
        }
      } else {
        # keep regardless as overlaps with roi:
        keep <- TRUE
      }
      
    } else {
      keep <- FALSE
    }
    
    if (keep) {
      if (i==1) {
        roi_split <- list(spl[[i]])
      } else {
        roi_split <- c(roi_split, spl[[i]])
      }
    } else {
      if (i==1) {
        roi_split <- list(GRanges(NULL))
      } else {
        roi_split <- c(roi_split, GRanges(NULL))
      }
    }

  }

  # bind as list:
  roi_split <- unlist(
    as(roi_split, "GRangesList")
  )
  names(roi_split) <- NULL
  
} else {
  roi_split <- GRanges(NULL)
}

# calculate and save total number of fusion supporting reads:
supporting_read_length <- length(roi_discordant)/2 + length(roi_split)/2
write.table(
  supporting_read_length,
  paste0(table_dir, "non_specific_fusion_supporting_reads.txt"),
  quote = F,
  row.names = F,
  col.names = F
)


####################################################################################
### 4. Plot read numbers ###
####################################################################################

# create venn diagram of filtered vs unfiltered reads:
filt_vs_unfilt <- list(
  unfiltered = c(unfilt_bam$concordant_pairs, unfilt_bam$discordant_pairs), 
  filtered = c(final_bam$concordant_pairs, final_bam$discordant_pairs)
)
filt_vs_unfilt_venn <- create_venn(filt_vs_unfilt, venn_cols)

# create venn diagram of split concordant pairs:
concordant_split <- final_bam[
  names(final_bam) %in% c(
    "concordant_pairs", "split_pairs", "split_supp"
  )
]
concordant_split_venn <- create_venn(concordant_split, venn_cols)

# create venn diagram of non-split concordant pairs:
concordant_non_split <- final_bam[
  names(final_bam) %in% c(
    "concordant_pairs", "non_split_concordant_pairs"
  )
]
concordant_non_split_venn <- create_venn(concordant_non_split, venn_cols)

# create venn diagram of discordant pairs:
discordant <- final_bam[
  names(final_bam) %in% c(
    "discordant_pairs", "split_pairs", "split_supp"
  )
]
discordant_venn <- create_venn(discordant, venn_cols)

# create barplot of read filtering:
read_numbers$type <- gsub("_", " ", rownames(read_numbers))
colnames(read_numbers) <- gsub("_", " ", colnames(read_numbers))
plot_df <- melt(read_numbers)
plot_df$type <- factor(
  plot_df$type, 
  levels = c(
    "concordant pairs", "discordant pairs", "split supp", 
    "split pairs", "non split concordant pairs"
  )
)

# plot:
p <- ggplot(plot_df, aes(x=variable, y=value, fill=type))
p <- p + geom_bar(stat="identity", position = "dodge")
p <- p + scale_fill_manual(values=plot_cols)
p <- p + ylab("No. reads")
p <- p + theme_cowplot(12)
read_summary_plot <- p + theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x=element_blank()
)

# plot log10:
plot_df$value[plot_df$value == 0] <- NA
p <- ggplot(plot_df, aes(x=variable, y=value, fill=type))
p <- p + geom_bar(stat="identity", position = "dodge")
p <- p + scale_y_continuous(trans='log10')
p <- p + scale_fill_manual(values=plot_cols)
p <- p + ylab("No. reads (log10)")
p <- p + theme_cowplot(12)
log10_read_summary_plot <- p + theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x=element_blank()
)

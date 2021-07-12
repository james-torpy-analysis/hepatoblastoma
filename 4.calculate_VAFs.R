args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
chromosomes <- args[3]
supps_allowed <- as.numeric(args[4])
SV_type <- args[5]
min_overlap <- as.numeric(args[6])
disc_read_window <- as.numeric(args[7])
same_SV_window <- as.numeric(args[8])
#projectname <- "hepatoblastoma"
#samplename <- "324_001_DB674_TAAGGCGA-CTCTCTAT_L001"
#chromosomes <- c("chr3")
#supps_allowed <- as.numeric("2")
#SV_type <- "deletion"
#min_overlap <- as.numeric("19")
#disc_read_window <- as.numeric("200")
#same_SV_window <- as.numeric("10")

venn_cols <- c("#7C1BE2", "#1B9E77", "#EFC000FF", "blue")

#home_dir <- "/Users/torpor/clusterHome/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
in_path <- paste0(
  result_dir, "VAF_calculation/", samplename, 
  "/max_supps_allowed_", supps_allowed, "/"
)
Robject_dir <- paste0(in_path, "Rdata/")

SV_dir <- paste0(result_dir, "SVs/", samplename, "/")
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
find_overlapping_reads <- dget(paste0(func_dir, "find_overlapping_reads.R"))

filter_overlaps <- function(reads, min_overlap, chromosomes) {
  
  if (any(chromosomes %in% seqnames(reads))) {
    # split by chromosome:
    split_reads <- split(reads, seqnames(reads))
    split_reads <- split_reads[chromosomes]
    
    # calculate lengths from start of read to SV, and from SV to end:
    split_reads <- lapply(split_reads, function(y) {
      if (length(y) > 0) {
        if (length(unique(seqnames(y))) == 1) {
          y$start_to_SV <- y$SV_coord - start(y)
          y$SV_to_end <- end(y) - y$SV_coord
        }
      }
      return(y)
    })
    
    if (length(split_reads) == 2) {
      # concatentate reads back together:
      prefilt_reads <- c(split_reads[[1]], split_reads[[2]])
    } else {
      prefilt_reads <- split_reads[[1]]
    }
    
    # split by qname:
    prefilt_split <- split(prefilt_reads, prefilt_reads$qname)
    
    # filter out reads without overlaps of at least min_overlap on both sides of
    # SV (need at least 1/2 read to satisfy this condition):
    filt_split <- lapply(prefilt_split, function(x) {
      remove_vec <- !(
        x$start_to_SV >= min_overlap & x$SV_to_end >= min_overlap
      )
      remove_vec[is.na(remove_vec)] <- FALSE
      if (any(remove_vec)) {
        return(NULL)
      } else {
        return(x)
      }
    })
    
    res <- unlist(
      as(
        filt_split[sapply(filt_split, function(x) !is.null(x))],
        "GRangesList"
      )
    )
    names(res) <- NULL
    
    return(res)
    
  } else {
    return(reads)
  }
  
}

fetch_mate_gap <- dget(paste0(func_dir, "fetch_mate_gap.R"))
find_spanning_discordant <- dget(paste0(func_dir, "find_spanning_discordant.R"))

plot_cols <- read.table(
  paste0(col_dir, "labelled_colour_palette.txt"),
  sep = "\t",
  header = F,
  comment.char = "",
  fill = TRUE
)$V1
plot_cols <- plot_cols[c(1:3, 5, 4, 6:length(plot_cols))]


####################################################################################
### 1. Load data ###
####################################################################################

# load in SVs:
both_SVs <- readRDS(paste0(SV_dir, "detected_SVs.Rdata"))

# fetch high confidence and low confidence breakpoints:
hc_SVs <- both_SVs$high_conf_bp$true_positives$SVs
lc_SVs <- both_SVs$low_conf_bp$true_positives$SVs

# label confidence of breakpoints:
try(hc_SVs$svaba_conf <- "high")
try(lc_SVs$svaba_conf <- "low")
  
# merge and remove duplicates:
SVs <- c(hc_SVs, lc_SVs)
SVs <- SVs[
  !duplicated(start(SVs)) | !duplicated(SVs$join_coord)
]

if (length(SVs) > 0) {

  # combine SVs with <= same_SV_window bp difference in coord:
  SVs$remove <- FALSE
  
  for (i in 1:length(SVs)) {
    
    if (!(SVs$remove[i])) {
      
      window1 <- (start(SVs[i])-same_SV_window):
        (start(SVs[i])+same_SV_window)
      window2 <- (SVs[i]$join_coord-same_SV_window):
        (SVs[i]$join_coord+same_SV_window)
      
      # check each other SV to see whether both coords is close
      # enough to those of another SV to justify removing one of them:
      for (j in 1:length(SVs)) {
        
        if (j!=i & start(SVs)[j] %in% window1 & 
            SVs$join_coord[j] %in% window2) {
          SVs$remove[j] <- TRUE
        }
        
        # check whether SVs are duplicated, but start(SV1) = join_coord(SV2)
        # and same for chromosomes:
        if (j!=i & start(SVs)[j] %in% window2 & 
            SVs$join_coord[j] %in% window1) {
          SVs$remove[j] <- TRUE
        }
        
      }
      
    }
    
  }
  
  # remove those marked as the same as another SV:
  SVs <- SVs[!(SVs$remove)]
  mcols(SVs) <- subset(SVs, select = -remove)
  
  # load filtered reads:
  filtered_reads <- readRDS(paste0(Robject_dir, "VAF_calculation_reads.Rdata"))
  
}

save.image(paste0(Robject_dir, "data_loaded_img.Rdata"))


####################################################################################
### 2. find breakpoint-overlapping reads: ###
####################################################################################

if (exists("filtered_reads")) {
    
  if (!file.exists(paste0(Robject_dir, "SV_overlapping_reads.Rdata"))) {
    
    # create empty list to be filled:
    temp_list <- list(
      spanning = NULL,
      overlapping = NULL
    )
    
    for (i in 1:length(SVs)) {
      if (i==1) {
        olap_reads <- list(
          list(
            non_supporting = temp_list,
            supporting = temp_list
          )
        )
      } else {
        olap_reads[[i]] <- list(
          non_supporting = temp_list,
          supporting = temp_list
        )
      }
    }
    names(olap_reads) <- start(SVs)
    
    # find overlaps of non-split concordant pairs with main breakpoint coords, 
    # add to non_supporting$overlapping:
    
    # split concordant reads by qname:
    spl <- split(
      filtered_reads$non_split_concordant_pairs, 
      filtered_reads$non_split_concordant_pairs$qname
    )
    
    # initiate cluster:
    cl <- makeCluster(7)
    clusterExport(
      cl, varlist = c("spl")
    )
    
    for (i in 1:length(SVs)) {

      # find overlaps with SV breakpoints:
      overlapping_reads <- parLapply(
        cl,
        spl,
        find_overlapping_reads,
        SV = SVs[i],
        alt_coords = FALSE
      )
      
      overlapping_reads <- unlist(
        as(
          overlapping_reads[
            sapply(overlapping_reads, function(x) !is.null(x))
          ],
          "GRangesList"
        )
      )
      names(overlapping_reads) <- NULL
      
      # make NULL entries empty GRanges object, to keep consistent with other
      # olap_reads elements:
      if (!is.null(overlapping_reads)) {
        olap_reads[[i]]$non_supporting$overlapping <- overlapping_reads
      } else {
        olap_reads[[i]]$non_supporting$overlapping <- GRanges(NULL)
      }
      
    }
    
    stopCluster(cl)
    
    # find overlaps of split primary reads with major breakpoint coords, 
    # add to supporting_$overlapping:
    spl <- split(filtered_reads$split_pairs, filtered_reads$split_pairs$qname)
    
    # initiate cluster:
    cl <- makeCluster(7)
    clusterExport(
      cl, varlist = c("spl")
    )
    
    for (i in 1:length(SVs)) {
      
      # find overlaps with SV breakpoints:
      overlapping_reads <- parLapply(
        cl,
        spl,
        find_overlapping_reads,
        SV = SVs[i],
        alt_coords = FALSE
      )
      
      overlapping_reads <- unlist(
        as(
          overlapping_reads[
            sapply(overlapping_reads, function(x) !is.null(x))
          ],
          "GRangesList"
        )
      )
      names(overlapping_reads) <- NULL
      
      # make NULL entries empty GRanges object:
      if (!is.null(overlapping_reads)) {
        olap_reads[[i]]$supporting$overlapping <- overlapping_reads
      } else {
        olap_reads[[i]]$supporting$overlapping <- GRanges(NULL)
      }
      
    }
    
    if (SV_type == "fusion") {
      
      # find overlaps of split primary reads with minor breakpoint coords, 
      # add to supporting_reads$overlapping:
      
      # find overlaps with SV breakpoint:
      overlapping_reads <- parLapply(
        cl,
        spl,
        find_overlapping_reads,
        SV = SVs[i],
        chromosome = "chr11"
      )
      
      stopCluster(cl)
      
      olap_reads[[i]]$supporting$overlapping <- c(
        olap_reads[[i]]$supporting$overlapping,
        unlist(
          as(
            overlapping_reads[
              sapply(overlapping_reads, function(x) !is.null(x))
            ],
            "GRangesList"
          )
        )
      )
      names(olap_reads[[i]]$supporting$overlapping) <- NULL
      
    }
    
    # filter reads:
    olap_reads <- lapply(olap_reads, function(x) {
      
      return(
        lapply(x, function(y) {
          
          if (length(y$overlapping) > 0) {
            
            # deduplicate reads:
            spl <- split(y$overlapping, y$overlapping$qname)
            spl <- spl[!duplicated(spl)]
            
            y$overlapping <- unlist(
              as(
                spl[sapply(spl, function(z) !is.null(z))],
                "GRangesList"
              )
            )
            names(y$overlapping) <- NULL
            
            # filter overlapping reads without at least min_overlap on both 
            # sides of SV:
            y$overlapping <- filter_overlaps(
              reads = y$overlapping, 
              min_overlap = min_overlap,
              chromosomes = chromosomes
            )
          }
          
          return(y)
          
        })
      )
      
    })
    
    saveRDS(olap_reads, paste0(Robject_dir, "SV_overlapping_reads.Rdata"))
    
  } else {
    olap_reads <- readRDS(paste0(Robject_dir, "SV_overlapping_reads.Rdata"))
  }

}

save.image(paste0(Robject_dir, "overlaps_found_img.Rdata"))


####################################################################################
### 3. Find breakpoint-spanning reads ###
####################################################################################

if (SV_type == "fusion") {
  if (exists("olap_reads")) {
    
    # redefine read object as 'all_reads' to add spanning reads to:
    all_reads <- olap_reads
    
    if (!file.exists(paste0(Robject_dir, "SV_reads.Rdata"))) {
      
      # find breakpoint-spanning reads from non-split non-discordant reads, 
      # add to spanning$non_supporting:
      
      # fetch read gaps (read + gap coords):
      # split by qname:
      spl <- split(
        filtered_reads$non_split_concordant_pairs, 
        filtered_reads$non_split_concordant_pairs$qname
      )
      
      # initiate cluster:
      cl <- makeCluster(7)
      clusterExport(
        cl, varlist = c("spl")
      )
      
      # fetch read gaps:
      non_split_concordant_gaps <- parLapply(
        cl, spl, fetch_mate_gap
      )
      
      stopCluster(cl)
      
      # merge into granges:
      non_split_concordant_gaps <- unlist(
        as(non_split_concordant_gaps, "GRangesList")
      )
      names(non_split_concordant_gaps) <- NULL
      
      # split by qname:
      spl <- split(non_split_concordant_gaps, non_split_concordant_gaps$qname)
      
      if (length(spl) > 0) {
        
        # initiate cluster:
        cl <- makeCluster(7)
        clusterExport(
          cl, varlist = c("spl")
        )
        
        for (i in 1:length(SVs)) {
          
          # find overlaps of non-supporting gaps with SVs:
          spanning_gaps <- parLapply(
            cl, 
            spl, 
            find_overlapping_reads, 
            SVs = SVs[i],
            chromosome = "chr22"
          )
          
          # merge ranges:
          spanning_gaps <- unlist(
            as(
              spanning_gaps[
                sapply(spanning_gaps, function(x) !is.null(x))
              ],
              "GRangesList"
            )
          )
          names(spanning_gaps) <- NULL
          
          stopCluster(cl)
          
          if (!is.null(spanning_gaps)) {
            
            # fetch corresponding read pairs:
            all_reads[[i]]$non_supporting$spanning <- filtered_reads$non_split_concordant_pairs[
              filtered_reads$non_split_concordant_pairs$qname %in% 
                spanning_gaps$qname
            ]
            m <- match(
              all_reads[[i]]$non_supporting$spanning$qname, spanning_gaps$qname
            )
            all_reads[[i]]$non_supporting$spanning$SV_coord <- 
              spanning_gaps$SV_coord[m]
            
          } else {
            all_reads[[i]]$non_supporting$spanning <- GRanges(NULL)
          }
          
        }
        
      } else {
        for (i in 1:length(SVs)) {
          all_reads[[i]]$non_supporting$spanning <- GRanges(NULL)
        }
      }
      
      if (length(filtered_reads$discordant_pairs) > 0) {
        
        # find breakpoint-spanning reads from discordant reads, 
        # add to spanning$supporting:
        spl <- split(
          filtered_reads$discordant_pairs, 
          filtered_reads$discordant_pairs$qname
        )
        
        for (i in 1:length(SVs)) {
          
          spanning_reads <- lapply(
            spl, 
            find_spanning_discordant, 
            SVs = SVs[i], 
            exp_window = disc_read_window
          )
          
          # merge to granges object:
          spanning_reads <- unlist(
            as(
              spanning_reads[
                sapply(spanning_reads, function(x) !is.null(x))
              ], "GRangesList"
            )
          )
          names(spanning_reads) <- NULL
          
          if (!is.null(overlapping_reads)) {
            all_reads[[i]]$supporting$spanning <- spanning_reads
          } else {
            all_reads[[i]]$supporting$spanning <- GRanges(NULL)
          }
          
        }
        
      } else {
        for (i in 1:length(all_reads)) {
          all_reads[[i]]$supporting$spanning <- GRanges(NULL)
        }
      }
      
      saveRDS(all_reads, paste0(Robject_dir, "SV_reads.Rdata"))
      
    } else {
      
      all_reads <- readRDS(paste0(Robject_dir, "SV_reads.Rdata"))
      
    }
    
  }
} else {
  if (exists("olap_reads")) {
    
    # make spanning read elements empty Granges:
    all_reads <- olap_reads
    all_reads <- lapply(all_reads, function(x) {
      x$non_supporting$spanning <- GRanges(NULL)
      x$supporting$spanning <- GRanges(NULL)
      return(x)
    })
    
  }
}

save.image(paste0(Robject_dir, "spanning_found_img.Rdata"))


####################################################################################
### 4. Calculate proportions of non-supporting vs supporting read pairs for each
# sample ###
####################################################################################

if (exists("all_reads")) {
  
  # combine all supporting reads, and all non-supporting reads, for each SV:
  combined_reads <- lapply(all_reads, function(x) {
    return(
      lapply(x, function(y) {
        return(c(y$spanning, y$overlapping))
      })
    )
  })
  
  # check all reads are in pairs:
  pair_check <- lapply(combined_reads, function(x) {
    sapply(x, function(y) {
      spl <- split(y, y$qname)
      return(
        all(
          sapply(spl, function(z) {
            length(z) == 2
          })
        )
      )
    })
  })
  
  print(
    paste0(
      "Are all supporting reads in pairs? ", 
      all(sapply(pair_check, function(x) x[names(x) == "supporting"]))
    )
  )
  
  print(
    paste0(
      "Are all non-supporting reads in pairs? ", 
      all(sapply(pair_check, function(x) x[names(x) == "non_supporting"]))
    )
  )
  
  # save each group as sam file:
  for (i in 1:length(all_reads)) {
    for (j in 1:length(all_reads[[i]])) {
      for (k in 1:length(all_reads[[i]][[j]])) {
        
        if (length(all_reads[[i]][[j]][[k]]) > 0) {
          
          # define sam cols:
          sam <- data.frame(
            qname = all_reads[[i]][[j]][[k]]$qname,
            flag = all_reads[[i]][[j]][[k]]$flag,
            rname = seqnames(all_reads[[i]][[j]][[k]]),
            pos = start(all_reads[[i]][[j]][[k]]),
            mapq = all_reads[[i]][[j]][[k]]$mapq,
            cigar = all_reads[[i]][[j]][[k]]$cigar,
            rnext = all_reads[[i]][[j]][[k]]$rnext,
            pnext = all_reads[[i]][[j]][[k]]$pnext,
            tlen = all_reads[[i]][[j]][[k]]$tlen,
            seq = all_reads[[i]][[j]][[k]]$seq,
            qual = all_reads[[i]][[j]][[k]]$qual
          )
          
          # write sam to tab-separated file:
          write.table(
            sam,
            paste0(
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam"
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
          )
          
          # add header:
          system(
            paste0(
              "samtools view -H ", 
              bam_path, samplename, "/", samplename, ".consensus.concordant.pairs.bam",
              " > ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam"
            )
          )
          
          # add rest of sam:
          system(
            paste0(
              "cat ", 
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam", 
              " >> ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam"
            )
          )
          
          # convert to bam:
          system(
            paste0(
              "samtools view -bh ", 
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam",
              " > ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
          # sort:
          system(
            paste0(
              "samtools sort -o ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sorted.bam ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
          # index:
          system(
            paste0(
              "samtools index ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sorted.bam"
            )
          )
          
          # clean:
          system(
            paste0(
              "rm ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads_temp.sam ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.sam ",
              SV_dir, "SV_coord_", names(all_reads)[i], "_", 
              names(all_reads[[i]])[j],"_", names(all_reads[[i]][[j]])[k], 
              "_reads.bam"
            )
          )
          
        }
        
      }
    }
  }
  
  # calculate VAFs for all SVs:
  read_counts <- sapply(combined_reads, function(x) {
    return(
      sapply(x, function(y) {
        length(y)
      })
    )
  })
  
  final_reads <- all_reads[[which.max(read_counts["supporting",])]]
  final_combined <- combined_reads[[which.max(read_counts["supporting",])]]
  
  final_combined_no <- data.frame(
    supporting = length(final_combined$supporting),
    non_supporting = length(final_combined$non_supporting)
  )
  
  write.table(
    final_combined_no, 
    paste0(in_path, "final_SV_read_nos.txt"),
    quote = F,
    row.names = F,
    col.names = T
  )
  
  # final_venns <- lapply(final_reads, function(x) {
  #   return(lapply(x, create_venn, venn_cols))
  # })
  
  # calculate VAFs for all SVs:
  VAFs <- round(
    (read_counts["supporting",]/
       (read_counts["supporting",] + read_counts["non_supporting",]))*100
    , 1
  )
  
  # save VAFs:
  write.table(
    as.data.frame(VAFs),
    paste0(in_path, "VAFs.txt"),
    quote = F,
    row.names = T,
    col.names = F
  )
  
}

# create dummy file for Snakemake if no VAF was calculated:
if (!file.exists(paste0(in_path, "VAFs.txt"))) {
  system(paste0("touch ", in_path, "VAFs.txt"))
  system(paste0("touch ", in_path, "final_SV_read_nos.txt"))
} 
  
  
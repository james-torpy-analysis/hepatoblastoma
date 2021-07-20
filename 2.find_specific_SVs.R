
args = commandArgs(trailingOnly=TRUE)

projectname <- args[1]
samplename <- args[2]
# projectname <- "hepatoblastoma"
# samplename <- "324_003_DB674_AGGCAGAA-CTCTCTAT_L001"

home_dir <- "/share/ScratchGeneral/jamtor/"
# home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
in_dir <- paste0(project_dir, "results/svaba/BWA_and_picard/", samplename, "/")
non_collapsed_dir <- paste0(project_dir, "results/svaba/non_collapsed/", samplename, "/")

func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")
ref_dir <- paste0(project_dir, "refs/")
Robject_dir <- paste0(project_dir, "results/SVs/", samplename, "/")

system(paste0("mkdir -p ", Robject_dir))


####################################################################################
### 0. Load packages, SVs and GOI coordinates ###
####################################################################################

library(rtracklayer)
library(dplyr)
library(tibble)
library(ComplexHeatmap)

load_breakpoints <- dget(paste0(func_dir, "load_breakpoints.R"))
find_specific_SVs <- dget(paste0(func_dir, "find_specific_SVs.R"))

if (!file.exists(paste0(Robject_dir, "detected_SVs.Rdata"))) {
  
  # define ROI co-ordinates:
  ROI <- GRanges(
    seqnames = "chr3",
    ranges = IRanges(start = 41265552, end = 41266752), # region covered by primers
    strand = "*"
  )

  ####################################################################################
  ### 1. Load VCFs and define EWSR1 co-ordinates ###
  ####################################################################################
  
  bp <- list(
    high_conf_bp = load_breakpoints(samplename, in_dir),
    low_conf_bp = load_breakpoints(samplename, in_dir, high_conf = FALSE)
  )
  
  
  ####################################################################################
  ### 2. Find overlaps with breakpoints and EWSR1 ###
  ####################################################################################
  
  SV_results <- lapply(
    bp,
    find_specific_SVs,
    ROI = ROI
  )
  
  # save SVs:
  saveRDS(
    SV_results,
    paste0(Robject_dir, "detected_SVs.Rdata")
  )
  
}


#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- file.path(home_dir, "projects/hepatoblastoma")
ref_dir <- file.path(project_dir, "refs")
func_dir <- file.path(project_dir, "scripts/functions")
variant_path <- file.path(project_dir, "results/smcounter2")

plot_dir <- file.path(variant_path, "plots")
table_dir <- file.path(variant_path, "tables")
Robject_dir <- file.path(variant_path, "Rdata")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", Robject_dir))

library(ggplot2)
library(GenomicAlignments)
library(rtracklayer)
library(tibble)
library(cowplot)
library(ggrepel)

fetch_sm_vafs <- dget(file.path(func_dir, "fetch_sm_vafs.R"))


####################################################################################
### 0. Load mutation data ###
####################################################################################

# load Sanger/GeneGlobe SNV info:
meta <- read.table(paste0(ref_dir, "/metadata.tsv"), header = T )

# fetch sample and treatment info:
meta_list <- split(meta, meta$Sample_id)

### check what happened to 021-1 record in metadata file! ###

# define CTNNB1 ROI:
CTNNB1 <- import(file.path(ref_dir, "CDHS-33412Z-324.roi_hg38.bed"))

# fetch and annotate VAFs:
sm_vars <- unlist(as(lapply(meta_list, fetch_sm_vafs, CTNNB1), "GRangesList"))

# add to SNV df:
var_df <- as.data.frame(sm_vars)
var_df$new_SNV <- paste0(var_df$seqnames, ":", var_df$start, "_", 
  var_df$ref, ">", var_df$alt )
var_df <- subset(var_df, select = c(new_SNV, effect_size, VAF) )
colnames(var_df) <- c("smCounter2_SNV", "smCounter2_effect_size", "smCounter2_VAF")
var_df$id <- rownames(var_df)

# merge with prior info:
final_vars <- merge(SNV, var_df, by = "id")
extra_df <- SNV[!(SNV$id %in% var_df$id),]
extra_df <- cbind(extra_df, 
  as.data.frame(matrix("not_identified", nrow = nrow(extra_df), ncol = 3)) )
colnames(extra_df) <- colnames(final_vars)
extra_df$smCounter2_effect_size <- NA

# bind NA values, format:
final_df <- subset(rbind(final_vars, extra_df), 
  select = c(id, Treatment, GeneGlobe_SNV, smCounter2_SNV, GeneGlobe_VAF, 
    smCounter2_VAF, smCounter2_effect_size ))
final_df$id <- gsub("_D.*$", "", final_df$id)
  
write.table(
  final_df,
  file.path(table_dir, "GeneGlobe_vs_smCounter2_SNV_VAFs.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE )

no_samples <- 51

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
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

source(file.path(func_dir, "4.compare_snv_functions.R"))


####################################################################################
### 0. Load mutation data ###
####################################################################################

# load Sanger/GeneGlobe point_mut info:
meta <- read.table(paste0(ref_dir, "/metadata.tsv"), header = T)

stopifnot(nrow(meta) == no_samples)

# fetch sample and treatment info:
meta_list <- split(meta, meta$Library_id)

# define CTNNB1 ROI:
CTNNB1 <- import(file.path(ref_dir, "CDHS-33412Z-324.roi_hg38.bed"))

# fetch and annotate VAFs:
sm_vars <- unlist(as(lapply(meta_list, fetch_sm_vafs, CTNNB1), "GRangesList"))

# format:
var_df <- as.data.frame(sm_vars)
var_df$new_point_mut <- paste0(var_df$seqnames, ":", var_df$start, "_", 
  var_df$ref, ">", var_df$alt )
var_df <- subset(var_df, select = c(new_point_mut, effect_size, VAF) )
colnames(var_df) <- c("smCounter2_point_mut", "smCounter2_effect_size", "smCounter2_VAF")
var_df$Library_id <- rownames(var_df)

# merge with prior info:
final_vars <- merge(meta, var_df, by = "Library_id", all=T)

# bind NA values, format:
final_df <- subset(final_vars, 
  select = c(Patient_id, Sample_id, Library_id, Treatment.dilution, Sanger_point_mut, smCounter2_point_mut, 
    ddPCR_point_mut_VAF, GeneGlobe_point_mut_VAF, smCounter2_VAF, smCounter2_effect_size,
    Sanger_deletion, Andre_deletion, Andre_deletion_confidence, ddPCR_deletion_VAF, 
    Andre_deletion_VAF, Resequenced, 
    Reads, UMIs, Reads_per_UMI ) )

# sort by original metadata order:
final_df <- final_df[match(meta$Sample_id, final_df$Sample_id),]

write.table(
  final_df,
  file.path(table_dir, "sample_summary_with_smcounter_point_mut.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE )

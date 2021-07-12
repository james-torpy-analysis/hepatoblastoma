
#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
ref_dir <- paste0(project_dir, "refs/")

func_dir <- paste0(project_dir, "scripts/functions/")
fusion_dir <- paste0(project_dir, "results/fusions/")
VAF_dir <- paste0(project_dir, "results/VAF_calculation/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
Robject_dir <- paste0(detection_dir, "/Rdata/")
plot_dir <- paste0(VAF_dir, "plots/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

library(dplyr)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_fusion_no <- dget(paste0(func_dir, "fetch_fusion_no.R"))
longitudinal_heatmap <- dget(paste0(func_dir, "longitudinal_heatmap.R"))

hm_cols <- c(
  FISH_detection = "#75EA3D",
  no_FISH_detection = "#D68EB7",
  stringent_detection = "#F7A006",
  less_stringent_detection = "#F4D30B",
  supporting_reads = "white",
  no_supporting_reads = "black",
  unknown = "grey"
)


####################################################################################
### 1. Load patient data, fusion calls and VAFs ###
####################################################################################

# load metadata:
patient_df <- read.table(
  paste0(ref_dir, "ES_samples_by_patient.tsv"),
  sep = "\t",
  header = T
)
# format colnames and Sample column:
colnames(patient_df) <- gsub("\\.", "_", colnames(patient_df))
samplenames <- patient_df$Sample
patient_df$Sample <- gsub(
  "_.*$", "", 
  gsub("409_", "", patient_df$Sample)
)

# order rows:
split_df <- split(patient_df, patient_df$Patient)

split_df <- lapply(split_df, function(x) {
  
  # tumour samples first:
  df <- x[x$Treatment == "tumour",]
  # add resected:
  df <- rbind(df, x[x$Treatment == "resection",])
  # add treatment naive:
  df <- rbind(df,x[x$Treatment == "naive",])
  # add NACT:
  df <- rbind(df, x[grep("NACT", x$Treatment),])
  # add pre-relapse ACT:
  df <- rbind(df,x[x$Treatment == "ACT1",])
  # add targeted:
  df <- rbind(df,x[x$Treatment == "targeted",])
  # add relapses:
  df <- rbind(df, x[x$Treatment == "relapse",])
  # add post-relapse ACT:
  return(rbind(df,x[x$Treatment == "ACT2",]))
  
})
# merge split elements:
patient_df <- do.call("rbind", split_df)

# load dilution samples metadata:
dilution_df <- read.table(
  paste0(ref_dir, "ES_samples_by_dilution.tsv"),
  sep = "\t",
  header = T
)
# order rows and format Sample column:
colnames(dilution_df) <- gsub("\\.", "_", colnames(dilution_df))
dilution_df <- arrange(dilution_df, desc(as.numeric(Dilution)))
samplenames <- c(samplenames, dilution_df$Sample)
samplenames <- sort(samplenames)
dilution_df$Sample <- gsub(
  "_.*$", "", 
  gsub("409_", "", dilution_df$Sample)
)

# define sample names:
samplenames <- as.list(list.files(VAF_dir, pattern = "409"))

# fetch fusion numbers:
fusion_nos <- fetch_fusion_no(
  samplenames, 
  fusion_dir,
  patient_df,
  dilution_df
)

# add fusion numbers to metadata df:
for (i in 1:length(fusion_nos)) {
  
  if (i==1) {
    fusion_dfs <- list(
      patient = merge(patient_df, fusion_nos[[i]], by = "Sample")
    )
  } else {
    fusion_dfs$dilution = merge(dilution_df, fusion_nos[[i]], by = "Sample")
  }
  
}

VAFs <- unlist(
  lapply(samplenames, function(x) {
    
    print(x)
    
    all_VAFs <- tryCatch(
      read.table(
        paste0(VAF_dir, x, "/VAFs.txt"),
        header = F
      ),
      error=function(err) NA
    )
    
    if (suppressWarnings(!is.na(all_VAFs))) {
      return(max(all_VAFs$V2))
    } else {
      return(all_VAFs)
    }
    
  })
)

# turn NAs to 0:
VAFs[is.na(VAFs)] <- 0

VAF_df <- data.frame(
  VAF = VAFs
)
VAF_df$Sample <- gsub(
  "_.*$", "", 
  sub(".*?_(.+)", "\\1", unlist(samplenames))
)

fusion_dfs <- lapply(fusion_dfs, function(x) {
  x <- merge(
    x,
    VAF_df,
    by = "Sample"
  )
  return(x)
})


# fetch fusion-supporting reads:
fusion_read_nos <- sapply(samplenames, function(x) {
  
  print(x)
  
  return(
    read.table(
      paste0(VAF_dir, x, "/tables/non_specific_fusion_supporting_reads.txt"),
      header = F
    )
  )
})
fusion_read_nos <- data.frame(
  Sample = gsub(
    "_.*$", "", 
    sub(".*?_(.+)", "\\1", unlist(samplenames))
  ),
  Supporting_read_pairs = unlist(fusion_read_nos)
)

# add to fusion_dfs:
fusion_dfs <- lapply(fusion_dfs, function(x) {
  return(
    merge(
      x,
      fusion_read_nos,
      by = "Sample"
    )
  )
})

#save.image(paste0(Robject_dir, "temp_img.Rdata"))
#load(paste0(Robject_dir, "temp_img.Rdata"))


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmap_FP <- longitudinal_heatmap(
  fusion_df = fusion_dfs$patient,
  hm_title = "Patient EWSR1/FLI1 fusion detections",
  type = "patient",
  annotation = "false positives",
  hm_cols = hm_cols
)

png(
  paste0(plot_dir, "patient_fusion_detection_heatmap_FP_annotated.png"),
  width = 10,
  height = 6,
  unit = "in",
  res = 300
)
  print(patient_heatmap_FP)
dev.off()

patient_heatmaps_VAF <- longitudinal_heatmap(
  fusion_df = fusion_dfs$patient,
  hm_title = "Patient EWSR1/FLI1 fusion detections",
  type = "patient",
  annotation = "VAF",
  hm_cols = hm_cols
)

png(paste0(plot_dir, "patient_fusion_detection_heatmap_VAFs_annotated.png"),
    width = 10,
    height = 6,
    unit = "in",
    res = 300
)
  print(patient_heatmaps_VAF$VAF_annot)
dev.off()

png(
  paste0(plot_dir, "patient_fusion_detection_heatmap_supporting_reads_annotated.png"),
  width = 10,
  height = 6,
  unit = "in",
  res = 300
)
  print(patient_heatmaps_VAF$sread_annot)
dev.off()

dilution_heatmap_FP <- longitudinal_heatmap(
  fusion_df = fusion_dfs$dilution,
  hm_title = "Cell line EWSR1/FLI1 fusion detections",
  type = "dilution",
  annotation = "false positives",
  hm_cols = hm_cols
)

png(
  paste0(plot_dir, "cell_line_fusion_detection_heatmap_FP_annotated.png"),
  width = 10,
  height = 3,
  unit = "in",
  res = 300
)
  print(dilution_heatmap_FP)
dev.off()

dilution_heatmaps_VAF <- longitudinal_heatmap(
  fusion_df = fusion_dfs$dilution,
  hm_title = "Cell line EWSR1/FLI1 fusion detections",
  type = "dilution",
  annotation = "VAF",
  hm_cols = hm_cols
)

png(
  paste0(plot_dir, "cell_line_fusion_detection_heatmap_VAFs_annotated.png"),
  width = 10,
  height = 3,
  unit = "in",
  res = 300
)
  print(dilution_heatmaps_VAF$VAF_annot)
dev.off()

png(
  paste0(plot_dir, "cell_line_fusion_detection_heatmap_supporting_reads_annotated.png"),
  width = 10,
  height = 3,
  unit = "in",
  res = 300
)
  print(dilution_heatmaps_VAF$sread_annot)
dev.off()



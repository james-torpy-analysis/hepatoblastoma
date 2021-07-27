
project_id <- "324"
SV_type <- "deletion"
condition_order <- c(
  paste0("naive", 1:3), paste0("NACT", 1:3), 
  "resection", paste0("relapse", 1:2))
dilution_order <- c(
  "100", "50", "10", "2", 
  "1", "0.4", "0.2", "0.1", 
  "0.04", "0"
)


#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/hepatoblastoma/")
ref_dir <- paste0(project_dir, "refs/")

func_dir <- paste0(project_dir, "scripts/functions/4.SV_call_heatmaps/")
SV_dir <- paste0(project_dir, "results/SVs/")
VAF_dir <- paste0(project_dir, "results/VAF_calculation/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
Robject_dir <- paste0(detection_dir, "/Rdata/")
plot_dir <- paste0(detection_dir, "plots/")

system(paste0("mkdir -p ", func_dir))
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

library(dplyr)
library(tibble)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

fetch_SV_no <- dget(paste0(func_dir, "fetch_SV_no.R"))
longitudinal_patient_heatmap <- dget(
  paste0(func_dir, "longitudinal_patient_heatmap.R")
)
longitudinal_dilution_heatmap <- dget(
  paste0(func_dir, "longitudinal_dilution_heatmap.R")
)
simplify_id <- dget(paste0(func_dir, "simplify_id.R"))

hm_cols <- c(
  pathology_detection = "#75EA3D",
  no_pathology_detection = "#58B9DB",
  stringent_detection = "#F7A006",
  less_stringent_detection = "#F4D30B",
  no_detection = "black",
  unknown = "grey"
)

annot_cols <- c(
  FP = "#991425",
  VAF = "#221699",
  sread = "#115E0F"
)


####################################################################################
### 1. Load patient data, SV calls and VAFs ###
####################################################################################

# load metadata:
patient_df <- read.table(
  paste0(ref_dir, "hepatoblastoma_samples_by_patient.tsv"),
  sep = "\t",
  header = T
)
# format colnames and Sample column:
colnames(patient_df) <- gsub("\\.", "_", colnames(patient_df))
samplenames <- patient_df$Sample
patient_df$Sample <- simplify_id(patient_df$Sample, divider = "-")

# load dilution samples metadata:
dilution_df <- read.table(
  paste0(ref_dir, "hepatoblastoma_samples_by_dilution.tsv"),
  sep = "\t",
  header = T
)
# # order rows and format Sample column:
# colnames(dilution_df) <- gsub("\\.", "_", colnames(dilution_df))
# dilution_df <- arrange(dilution_df, desc(as.numeric(Condition)))
dilution_df$Sample <- simplify_id(dilution_df$Sample, divider = "-")

# define sample names:
samplenames <- as.list(list.files(VAF_dir, pattern = project_id))

# fetch SV numbers:
SV_nos <- fetch_SV_no(
  samplenames, 
  SV_dir,
  patient_df,
  dilution_df
)

# add SV numbers to metadata df:
for (i in 1:length(SV_nos)) {
  
  if (i==1) {
    SV_dfs <- list(
      patient = merge(patient_df, SV_nos[[i]], by = "Sample")
    )
  } else {
    SV_dfs$dilution = merge(dilution_df, SV_nos[[i]], by = "Sample")
  }
  
}

# fetch VAFs with most supporting reads for each sample, 
# prioritising high conf:
VAFs <- unlist(
  lapply(samplenames, function(x) {
    
    print(x)
    
    all_VAFs <- tryCatch(
      all_VAFs <- readRDS(paste0(VAF_dir, x, "/Rdata/VAFs.Rdata")),
      error=function(err) NA
    )
    
    if (suppressWarnings(!is.na(all_VAFs))) {
      
      # prioritise high confidence VAFs:
      if ("high" %in% all_VAFs$conf) {
        temp_VAFs <- all_VAFs[all_VAFs$conf == "high",]
      } else {
        temp_VAFs <- all_VAFs
      }
      return(round(max(temp_VAFs$VAF)*100, 1))
    
    } else {
      return(0)
    }
    
  })
)

VAF_df <- data.frame(
  VAF = VAFs
)
VAF_df$Sample <- simplify_id(unlist(samplenames), divider = "_")

# merge VAFs with SV df:
SV_dfs <- lapply(SV_dfs, function(x) {
  x <- merge(
    x,
    VAF_df,
    by = "Sample"
  )
  return(x)
})

# save VAFs:
VAF_df <- VAF_df %>% column_to_rownames("Sample")
saveRDS(VAF_df, paste0(Robject_dir, "all_VAFs.Rdata"))

save.image(paste0(Robject_dir, "temp_img.Rdata"))
#load(paste0(Robject_dir, "temp_img.Rdata"))


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmaps_VAF <- longitudinal_patient_heatmap(
  SV_df = SV_dfs$patient,
  hm_title = paste0("Patient ", SV_type, " detections"),
  annot = "VAF",
  condition_order,
  hm_cols = hm_cols,
  func_dir
)

dilution_heatmaps_VAF <- longitudinal_dilution_heatmap(
  SV_df = SV_dfs$dilution,
  hm_title = paste0("Cell line ", SV_type, " detections"),
  annot = "VAF",
  dilution_order,
  hm_cols = hm_cols,
  func_dir
)

all_heatmaps <- list(
  patient_heatmap_FP = longitudinal_patient_heatmap(
    SV_df = SV_dfs$patient,
    hm_title = paste0("Patient ", SV_type, " detections"),
    annot = "false positives",
    condition_order,
    hm_cols = hm_cols,
    func_dir
  ),
  patient_heatmap_VAF = patient_heatmaps_VAF,
  dilution_heatmap_FP = longitudinal_dilution_heatmap(
    SV_df = SV_dfs$dilution,
    hm_title = paste0("Cell line ", SV_type, " detections"),
    annot = "false positives",
    dilution_order,
    hm_cols = hm_cols,
    func_dir
  ),
  dilution_heatmap_VAF = dilution_heatmaps_VAF
)

for (i in seq_along(all_heatmaps)) {
  
  # convert to grob:
  hm_grob <- grid.grabExpr(
    draw(all_heatmaps[[i]], gap = unit(6, "mm"))
  )
  dev.off()
  
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 10,
      height = 6,
      unit = "in",
      res = 300
    )
    annot_coord <-  0.41
    
  } else {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 10,
      height = 3,
      unit = "in",
      res = 300
    )
    annot_coord <-  0.37
    
  }
  
    grid.newpage()
    
    pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                          just = c("left", "bottom")))
      grid.draw(hm_grob)
    popViewport()
    
    if (length(grep("FP", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = 0.99, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
      
        grid.text(
          "x = No. false positives", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "FP"])
        )
      
      popViewport()
      
    } else if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = 0.95, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
        grid.text(
          "x% = VAF", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "VAF"])
        )
      popViewport()
      
    } else {
      
      pushViewport(viewport(x = 0.995, y = annot_coord, width = 0.2, height = 0.1, 
        just = c("right", "top")))
        grid.text(
          "x = No. supporting reads", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "sread"])
        )
      popViewport()
    }
  
  dev.off()
  
}


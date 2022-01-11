
home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/hepatoblastoma/")
ref_dir <- paste0(project_dir, "refs/")

func_dir <- paste0(project_dir, "scripts/functions/")
fusion_dir <- paste0(project_dir, "results/fusions/")
VAF_dir <- paste0(project_dir, "results/detection_and_VAF/")

Robject_dir <- paste0(VAF_dir, "/Rdata/")
plot_dir <- paste0(VAF_dir, "plots/")

dir.create(Robject_dir)
dir.create(plot_dir)

library(dplyr)
library(tibble)
library(naturalsort)


####################################################################################
### 0. Load functions and colours ###
####################################################################################

source(paste0(func_dir, "6.fusion_call_heatmaps_functions.R"))

hm_cols <- c(
  pathology_detection = "#75EA3D",
  no_pathology_detection = "#58B9DB",
  detected = "#F4D30B",
  not_detected = "black",
  no_sample = "grey" )

annot_cols <- c(
  VAF = "#991425",
  sread = "#430F82" )

treatment_order <- c("naive1", "naive2", "naive3", "NACT1", "NACT2", "NACT3", 
  "resection", "relapse1", "relapse2")

####################################################################################
### 1. Load summary table and split into patient and cell line dfs ###
####################################################################################

summary_df <- read.table(
  file.path(VAF_dir, "tables/final_summary.tsv"),
  header = T,
  sep = "\t" )

# separate and order patients:
patient_df <- summary_df[
  summary_df$Patient_id != "A673" & summary_df$Patient_id != "HepG2", ]

patient_df$Treatment.dilution <- factor(patient_df$Treatment.dilution,
  levels = treatment_order )
patient_df <- patient_df[order(patient_df$Treatment.dilution),]

# separate and order cell line samples:
cl_df <- summary_df[
  summary_df$Patient_id == "A673" | summary_df$Patient_id == "HepG2", ]
cl_df <- cl_df[naturalorder(cl_df$Treatment.dilution, decreasing = TRUE),]
cl_df <- cl_df[order(cl_df$Patient_id),]

# multiply deletion VAFs by 100:
patient_df$Deletion_VAF <- as.numeric(patient_df$Deletion_VAF)*100
patient_df$Deletion_VAF[is.na(patient_df$Deletion_VAF)] <- "not_detected"
cl_df$Deletion_VAF <- as.numeric(cl_df$Deletion_VAF)*100
cl_df$Deletion_VAF[is.na(cl_df$Deletion_VAF)] <- "not_detected"


####################################################################################
### 2. Create longitudinal heatmaps ###
####################################################################################

patient_heatmaps <- mutation_heatmap(
  mutation_df = patient_df,
  type = "patient",
  hm_cols )

cl_heatmaps <- mutation_heatmap(
  mutation_df = cl_df[cl_df$Patient_id != "A673",],
  type = "dilution",
  hm_cols )

all_heatmaps <- c(patient_heatmaps, cl_heatmaps)
names(all_heatmaps) <- c(
  paste0("patient_heatmap_", names(all_heatmaps)[1:2]),
  paste0("cell_line_heatmap_", names(all_heatmaps)[3:4]) )

for (i in seq_along(all_heatmaps)) {
  
  # convert to grob:
  hm_grob <- grid.grabExpr(
    draw(all_heatmaps[[i]], gap = unit(6, "mm"))
  )
  dev.off()
  
  # write to png:
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 14.5,
      height = 9,
      unit = "in",
      res = 300
    )
    annot_y <-  0.45
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.798
    } else {
      annot_x <- 0.828
    }
    
  } else {
    
    png(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.png"),
      width = 20,
      height = 2,
      unit = "in",
      res = 300
    )
    annot_y <-  0.25
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.835
    } else {
      annot_x <- 0.845
    }
    
  }
  
    grid.newpage()
    
    pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                          just = c("left", "bottom")))
      grid.draw(hm_grob)
    popViewport()
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      
      pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
        just = c("left", "top")))
        grid.text(
          "x% = VAF", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "VAF"])
        )
      popViewport()
      
    } else {
      
      pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
        just = c("left", "top")))
        grid.text(
          "x = No. supporting reads", 
          gp=gpar(fontsize=10, fontface="bold", 
            col=annot_cols[names(annot_cols) == "sread"])
        )
      popViewport()
    }
  
  dev.off()
  
  # write to pdf:
  if (length(grep("patient", names(all_heatmaps)[i])) > 0) {
    
    pdf(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.pdf"),
      width = 14.5,
      height = 9,
    )
    annot_y <-  0.45
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.798
    } else {
      annot_x <- 0.828
    }
    
  } else {
    
    pdf(
      paste0(plot_dir, names(all_heatmaps)[i], "_annotated.pdf"),
      width = 26,
      height = 2,
    )
    annot_y <-  0.22
    
    if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
      annot_x <- 0.845
    } else {
      annot_x <- 0.855
    }
    
  }
  
  grid.newpage()
  
  pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
                        just = c("left", "bottom")))
  grid.draw(hm_grob)
  popViewport()
  
  if (length(grep("VAF", names(all_heatmaps)[i])) > 0) {
    
    pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
                          just = c("left", "top")))
    grid.text(
      "x% = VAF", 
      gp=gpar(fontsize=10, fontface="bold", 
              col=annot_cols[names(annot_cols) == "VAF"])
    )
    popViewport()
    
  } else {
    
    pushViewport(viewport(x = annot_x, y = annot_y, width = 0.2, height = 0.1, 
                          just = c("left", "top")))
    grid.text(
      "x = No. supporting reads", 
      gp=gpar(fontsize=10, fontface="bold", 
              col=annot_cols[names(annot_cols) == "sread"])
    )
    popViewport()
  }
  
  dev.off()
  
}


projectname <- "hepatoblastoma"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome"
project_dir <- file.path(home_dir, "projects", projectname)
func_dir <- file.path(project_dir, "scripts/functions")
result_dir <- file.path(project_dir, "results/")
ref_dir <- file.path(project_dir, "refs")
variant_dir <- file.path(result_dir, "smcounter2/tables")
in_path <- file.path(result_dir, "detection_and_VAF")

Robject_dir <- file.path(in_path, "Rdata")
table_dir <- file.path(in_path, "tables")
plot_dir <- file.path(in_path, "plots")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", plot_dir))

library(ggplot2)
library(cowplot)
library(ggrepel)

compare_VAF <- dget(file.path(func_dir, "compare_VAF.R"))


####################################################################################
### 1. Load VAFs ###
####################################################################################

summary_df <- read.table(
  file.path(variant_dir, "sample_summary_with_smcounter_SNV.tsv"),
  sep = "\t",
  header = TRUE )

VAFs <- lapply(summary_df$Library_id, function(x) {
  print(x)
  VAF <- as.data.frame(readRDS(file.path(in_path, x, "Rdata/VAF.rds")))
  if (!is.null(VAF)) {
    if (nrow(VAF) > 0) {
      VAF <- VAF[order(VAF$total_supp, decreasing = T),]
      VAF$Library_id <- rep(x, nrow(VAF))
      VAF$up_VAF <- round(as.numeric(VAF$up_VAF), 4)
      VAF$dn_VAF <- round(as.numeric(VAF$dn_VAF), 4)
      VAF$avg_VAF <- round(as.numeric(VAF$avg_VAF), 4)
    } else {
      VAF <- data.frame(
        seqnames = "not_detected",
        start = "not_detected",
        end = "not_detected",
        width = "not_detected",
        strand = "not_detected",
        bp_A_supp = "not_detected",
        bp_A_non_supp = "not_detected",
        bp_A_total = "not_detected",
        bp_B_supp = "not_detected",
        bp_B_non_supp = "not_detected",
        bp_B_total = "not_detected",
        total_supp = "not_detected",
        up_VAF = "not_detected",
        dn_VAF = "not_detected",
        avg_VAF = "not_detected",
        Library_id = x
      )
      rownames(VAF) <- x
    }
  }
  return(VAF)
})

VAF_df <- do.call("rbind", VAFs)

VAF_df$deletion <- paste0(VAF_df$seqnames, ":", VAF_df$start, "-", VAF_df$end)
VAF_df$deletion[grep("not_detected", VAF_df$deletion)] <- "not_detected"
VAF_df <- subset(VAF_df, select = c(Library_id, deletion, up_VAF, dn_VAF, avg_VAF, 
  bp_A_supp, bp_A_non_supp, bp_A_total, bp_B_supp, bp_B_non_supp, bp_B_total, total_supp ) )


####################################################################################
### 2. Add previous VAF estimates and write table ###
####################################################################################

# format and make values numeric:
all_VAFs <- merge(summary_df, VAF_df, by="Library_id")

all_VAFs <- subset(all_VAFs, select = c(
  Patient_id, Sample_id, Library_id, Treatment, Sanger_SNV, 
  smCounter2_SNV, ddPCR_SNV_VAF, GeneGlobe_SNV_VAF, smCounter2_VAF, 
  smCounter2_effect_size, Sanger_deletion, Andre_deletion, 
  Andre_deletion_confidence, deletion, 
  ddPCR_deletion_VAF, Andre_deletion_VAF, avg_VAF, 
  up_VAF, dn_VAF, bp_A_supp, 
  bp_A_non_supp, bp_A_total, bp_B_supp, 
  bp_B_non_supp, bp_B_total, total_supp ))

colnames(all_VAFs) <- c(
  "Patient_id", "Sample_id", "Library_id", "Treatment", "Sanger_SNV", 
  "smCounter2_SNV", "ddPCR_SNV_VAF", "GeneGlobe_SNV_VAF", "smCounter2_VAF", 
  "smCounter2_effect_size", "Sanger_deletion", "Andre_deletion", 
  "Andre_deletion_confidence", "Deletion",
  "ddPCR_deletion_VAF", "Andre_deletion_VAF", "Deletion_VAF", 
  "Upstream_deletion_VAF", "Downstream_deletion_VAF", "Upstream_supporting", 
  "Upstream_non_supporting", "Upstream_total", "Downstream_supporting", 
  "Downstream_non_supporting", "Downstream_total", "Total_supporting" )

# subset deletion columns only, order by patient and write:
deletion_VAFs <- subset(all_VAFs, select = -c(
  Sanger_SNV, smCounter2_SNV, ddPCR_SNV_VAF, GeneGlobe_SNV_VAF, smCounter2_VAF, 
  smCounter2_effect_size ))
deletion_VAFs$Library_id <- factor(
  deletion_VAFs$Library_id, levels = summary_df$Library_id )
deletion_VAFs <- deletion_VAFs[order(deletion_VAFs$Library_id),]

write.table(
  deletion_VAFs, 
  file.path(table_dir, "final_deletion_summary.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE )

# remove secondary deletion rows, order by patient and write:
all_VAFs <- all_VAFs[!duplicated(all_VAFs$Library_id),]
all_VAFs <- all_VAFs[match(summary_df$Library_id, all_VAFs$Library_id), ]
write.table(
  all_VAFs, 
  file.path(table_dir, "final_summary.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE )


####################################################################################
### 3. Plot VAF correlations ###
####################################################################################

# check distributions of VAFs:
plot(hist(all_VAFs$ddPCR[all_VAFs$ddPCR != 0]))
plot(hist(all_VAFs$Andre_deletion_VAF[all_VAFs$Andre_deletion_VAF != 0]))
plot(hist(all_VAFs$avg_VAF[all_VAFs$avg_VAF != 0]))
dev.off()

# create ddPCR vs Andre_deletion_VAF VAF plot:
VAF_df <- subset(all_VAFs, select = c(Sample_id, Treatment, ddPCR_deletion_VAF, Andre_deletion_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment <- gsub("dilution_.*$", "HepG2", VAF_df$treatment)
ddPCR_vs_Andre_avg <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "Andre_avg", lim = 0.4, cortype = "pearson" )

pdf(file.path(plot_dir, "ddPCR_vs_Andre_deletion_VAFs.pdf"),
    height = 3.5,
    width = 5)
  print(ddPCR_vs_Andre_avg)
dev.off()

png(file.path(plot_dir, "ddPCR_vs_Andre_deletion_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(ddPCR_vs_Andre_avg)
dev.off()

# create average VAF vs Andre's average VAF plot:
VAF_df <- subset(all_VAFs, select = c(Sample_id, Treatment, Andre_deletion_VAF, Deletion_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment <- gsub("dilution_.*$", "HepG2", VAF_df$treatment)
Andre_vs_new_avg <- compare_VAF(
  VAF_df, lab1 = "avg Andre", lab2 = "avg new", lim = 0.4, cortype = "pearson" )

png(file.path(plot_dir, "Andre_deletion_VAF_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(Andre_vs_new_avg)
dev.off()

pdf(file.path(plot_dir, "Andre_deletion_VAF_vs_new_VAFs.pdf"),
    height = 3.5,
    width = 5)
print(Andre_vs_new_avg)
dev.off()

# create ddPCR vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(Sample_id, Treatment, ddPCR_deletion_VAF, Deletion_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
VAF_df$treatment <- gsub("dilution_.*$", "HepG2", VAF_df$treatment)
ddPCR_vs_new <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "avg new", lim = 0.4, 
  cortype = "pearson" )

png(file.path(plot_dir, "ddPCR_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(ddPCR_vs_new)
dev.off()

pdf(file.path(plot_dir, "ddPCR_vs_new_VAFs.pdf"),
    height = 3.5,
    width = 5)
print(ddPCR_vs_new)
dev.off()


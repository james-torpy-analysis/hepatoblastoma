
#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/hepatoblastoma/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/5.vaf_scatterplots/")

detection_dir <- paste0(project_dir, "results/detection_heatmaps/")
Robject_dir <- paste0(detection_dir, "Rdata/")
plot_dir <- paste0(detection_dir, "plots/")

system(paste0("mkdir -p ", plot_dir))

library(ggplot2)
library(GenomicAlignments)
library(tibble)
library(cowplot)
library(ggrepel)
library(stats)

simplify_id <- dget(paste0(func_dir, "simplify_id.R"))


####################################################################################
### 0. Load mutation data ###
####################################################################################

# load SV VAFs:
SV_VAFs <- readRDS(paste0(Robject_dir, "all_VAFs.Rdata"))
SV_VAFs$VAF <- round(SV_VAFs$VAF, 1)
colnames(SV_VAFs) <- "SV_VAF"

# load and add ddPCR/Andre VAFs:
patient_meta <- read.table(
  paste0(ref_dir, "/hepatoblastoma_samples_by_patient.tsv"),
  header = T,
  fill = T
)
patient_meta$Sample <- simplify_id(patient_meta$Sample, divider = "-")
rownames(patient_meta) <- patient_meta$Sample
patient_VAFs <- merge(patient_meta, SV_VAFs, by=0, all = FALSE)


####################################################################################
### 1. Plot comparisons of VAFs ###
####################################################################################

comparison_df <- subset(
  patient_VAFs, select = c(Patient, Condition, ddPCR_VAF, Andre_VAF, SV_VAF))

# check distributions:
plot(hist(comparison_df$ddPCR_VAF[comparison_df$ddPCR_VAF != 0]))
plot(hist(comparison_df$Andre_VAF[comparison_df$Andre_VAF != 0]))
plot(hist(comparison_df$Andre_VAF[comparison_df$SV_VAF != 0]))

# loop through to plot each comparison:
l=1
for (i in 3:5) {
  for (j in 3:5) {
    
    if (i != j) {
      
      # # remove 0 values:
      # plot_df <- comparison_df[
      #   !is.na(comparison_df[,i]) & !is.na(comparison_df[,j]) & 
      #     comparison_df[,i] != 0 & comparison_df[,j] != 0, ]
      
      # name cols generically:
      plot_df <- comparison_df
      colnames(plot_df)[i] <- "var1"
      colnames(plot_df)[j] <- "var2"
      
      # # calculate correlation:
      # corr <- cor.test(
      #   x = plot_df[,i], 
      #   y = plot_df[,j],
      #   method = "kendall"
      # )
      # 
      # # Fit regression line:
      # reg <- lm(var1 ~ var2, data = plot_df)
      # coeff=coefficients(reg)
      
      p <- ggplot(plot_df, aes(x = var1, y = var2, color = Condition))
      p <- p + geom_point(size = 1.5)
      p <- p + theme_cowplot(12)
      p <- p + xlim(c(0, 50))
      p <- p + ylim(c(0, 50))
      p <- p + xlab(colnames(comparison_df)[i])
      p <- p + ylab(colnames(comparison_df)[j])
      p <- p + geom_text_repel(data=plot_df, aes(label=Patient))
      # p <- p + geom_abline(
      #   intercept = coeff[1], slope = coeff[2], color = "red")
      # p <- p + annotate(
      #   "text", x = 25, y = 45, 
      #   label = paste0(
      #     "R2=", round(corr$estimate, 2), ", p=", 
      #     formatC(corr$p.value, digits = 2, format = "e")),
      #   color='red', size = 4
      # )
      
      if (l == 1) {
        cor_plots <- list(p)
        names(cor_plots) <- paste0(
          colnames(comparison_df)[i], "_vs_", colnames(comparison_df)[j])
      } else {
        cor_plots[[l]] <- p
        names(cor_plots)[l] <- paste0(
          colnames(comparison_df)[i], "_vs_", colnames(comparison_df)[j])
      }
      
    } else {
      
      if (l == 1) {
        cor_plots <- list(NA)
        names(cor_plots)[l] <- paste0(
          colnames(comparison_df)[i], "_vs_", colnames(comparison_df)[j])
      } else {
        cor_plots[[l]] <- NA
        names(cor_plots)[l] <- paste0(
          colnames(comparison_df)[i], "_vs_", colnames(comparison_df)[j])
      }
      
    }
    
    l <<- l+1
    
  }
}

# plot:
png(paste0(plot_dir, "ddPCR_vs_Andre_VAF.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(cor_plots$ddPCR_VAF_vs_Andre_VAF)
dev.off()

png(paste0(plot_dir, "SV_vs_ddPCR_VAF.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(cor_plots$SV_VAF_vs_ddPCR_VAF)
dev.off()

png(paste0(plot_dir, "SV_vs_Andre_VAF.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(cor_plots$SV_VAF_vs_Andre_VAF)
dev.off()


####################################################################################
### 2. Plot SV vs GeneGlobe SNP VAFs ###
####################################################################################

# plot SV vs TP53 VAFs, leaving out samples with either = 0 or unknown:
SV_vs_SM_TP53 <- subset(
  patient_VAFs, select = c(SM_TP53_VAF, SV_VAF, Treatment, Patient)
)

SV_vs_SM_TP53 <- SV_vs_SM_TP53[
  SV_vs_SM_TP53$SM_TP53_VAF != 0 & 
    SV_vs_SM_TP53$SV_VAF != 0 &
    SV_vs_SM_TP53$SM_TP53_VAF != "unknown" & 
    SV_vs_SM_TP53$SV_VAF != "unknown",
]

# calculate correlation between SV and TP53 VAFs:
plot(hist(SV_vs_SM_TP53$SV_VAF))
plot(hist(SV_vs_SM_TP53$SM_TP53_VAF))
dev.off()

corr <- cor.test(
  x = SV_vs_SM_TP53$SM_TP53_VAF, 
  y = SV_vs_SM_TP53$SV_VAF,
  method = "pearson"
)

# Fit regression line
require(stats)
reg <- lm(SM_TP53_VAF ~ SV_VAF, data = SV_vs_SM_TP53)
coeff=coefficients(reg)

p <- ggplot(
  SV_vs_SM_TP53, 
  aes(x = SV_VAF, y = SM_TP53_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 TP53 VAF")
p <- p + xlab("SV VAF")
p <- p + geom_text_repel(data=SV_vs_SM_TP53, aes(label=Patient))
p <- p + annotate(
  "text", x = 80, y = 85, 
  label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 3)
  ), color='red', size = 4
)
p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")

png(paste0(plot_dir, "SV_vs_SM_TP53_VAFs.png"),
  height = 3.5,
  width = 5,
  res = 300,
  units = "in")
  print(p)
dev.off()

# plot SV vs STAG2 VAFs, leaving out samples with either = 0 or unknown:
SV_vs_SM_STAG2 <- subset(
  patient_VAFs, select = c(SM_STAG2_VAF, SV_VAF, Treatment, Patient)
)

SV_vs_SM_STAG2 <- SV_vs_SM_STAG2[
  SV_vs_SM_STAG2$SM_STAG2_VAF != 0 & 
    SV_vs_SM_STAG2$SV_VAF != 0 &
    SV_vs_SM_STAG2$SM_STAG2_VAF != "unknown" & 
    SV_vs_SM_STAG2$SV_VAF != "unknown",
]
# calculate correlation between SV and STAG2 VAFs:
plot(hist(SV_vs_SM_STAG2$SV_VAF))
plot(hist(SV_vs_SM_STAG2$SM_STAG2_VAF))
dev.off()

corr <- cor.test(
  x = SV_vs_SM_STAG2$SM_STAG2_VAF, 
  y = SV_vs_SM_STAG2$SV_VAF,
  method = "pearson"
)

# Fit regression line
require(stats)
reg <- lm(SM_STAG2_VAF ~ SV_VAF, data = SV_vs_SM_STAG2)
coeff=coefficients(reg)

p <- ggplot(
  SV_vs_SM_STAG2, 
  aes(x = SV_VAF, y = SM_STAG2_VAF, color = Treatment)
)
p <- p + geom_point(size = 3)
p <- p + theme_cowplot(12)
p <- p + ylim(c(0, 100))
p <- p + xlim(c(0, 100))
p <- p + ylab("smCounter2 STAG2 VAF")
p <- p + xlab("SV VAF")
p <- p + geom_text_repel(data=SV_vs_SM_STAG2, aes(label=Patient))
p <- p + annotate(
  "text", x = 80, y = 85, 
  label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 3)
  ), color='red', size = 4
)
p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")

png(paste0(plot_dir, "SV_vs_SM_STAG2_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(p)
dev.off()


####################################################################################
### 3. Plot SV vs smcounter2 SNP VAFs ###
####################################################################################



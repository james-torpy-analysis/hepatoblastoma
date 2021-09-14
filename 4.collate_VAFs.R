projectname <- "hepatoblastoma"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/torpor/clusterHome/"
project_dir <- paste0(home_dir, "projects/", projectname, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
result_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
in_path <- paste0(result_dir, "detection_and_VAF/")

Robject_dir <- paste0(in_path, "Rdata/")
table_dir <- paste0(in_path, "tables/")
plot_dir <- paste0(in_path, "plots/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", plot_dir))

library(ggplot2)
library(cowplot)
library(ggrepel)


####################################################################################
### 1. Load VAFs ###
####################################################################################

samplenames <- list.files(in_path, pattern = "324|409")
samplenames <- samplenames[!(samplenames %in% c(
  "324_001_DB674_TAAGGCGA-CTCTCTAT_L001", 
  "324_003_DB674_AGGCAGAA-CTCTCTAT_L001" ))]

VAFs <- lapply(samplenames, function(x) {
  VAF <- as.data.frame(readRDS(paste0(in_path, x, "/Rdata/selected_VAF.rds")))
  if (!is.null(VAF)) {
    if (nrow(VAF) > 0) {
      rownames(VAF) <- x
    } else {
      VAF <- NULL
    }
  }
  return(VAF)
})
VAFs <- VAFs[sapply(VAFs, function(x) !is.null(x))]
VAF_df <- do.call("rbind", VAFs)
VAF_df$id <- gsub("_D.*$", "", rownames(VAF_df))

VAF_df <- subset(VAF_df, select = c(id, up_VAF, dn_VAF, avg_VAF))


####################################################################################
### 2. Add previous VAF estimates ###
####################################################################################

prior_VAFs <- read.table(
  paste0(ref_dir, "previous_hepatoblastoma_VAFs.tsv"),
  sep = "\t",
  header = TRUE )

# format and make values numeric:
all_VAFs <- merge(VAF_df, prior_VAFs, by="id")
for (i in 2:8) {
  all_VAFs[,i] <- as.numeric(all_VAFs[,i])
}


# check distributions of VAFs:
plot(hist(all_VAFs$ddPCR[all_VAFs$ddPCR != 0]))
plot(hist(all_VAFs$up_Andre[all_VAFs$up_Andre != 0]))
plot(hist(all_VAFs$dn_Andre[all_VAFs$dn_Andre != 0]))
plot(hist(all_VAFs$avg_Andre[all_VAFs$avg_Andre != 0]))
plot(hist(all_VAFs$up_VAF[all_VAFs$up_VAF != 0]))
plot(hist(all_VAFs$dn_VAF[all_VAFs$dn_VAF != 0]))
plot(hist(all_VAFs$avg_VAF[all_VAFs$avg_VAF != 0]))
dev.off()


####################################################################################
### 3. Plot VAF correlations ###
####################################################################################

# function to create VAF correlation plots:
compare_VAF <- function(VAF_df, lab1, lab2, lim = 40, cortype = "kendall") {
  
  # remove samples without two values:
  VAF_df <- VAF_df[!is.na(VAF_df$VAF1) & !is.na(VAF_df$VAF1),]
  # remove samples above lim:
  VAF_df <- VAF_df[VAF_df$VAF1 <= lim & VAF_df$VAF2 <= lim,]
  
  # calculate correlation between ddPCR and Andre VAFs:
  corr <- cor.test(x = VAF_df$VAF1, y = VAF_df$VAF2, method = cortype )
  
  # Fit regression line
  require(stats)
  reg <- lm(VAF2 ~ VAF1, data = VAF_df)
  coeff=coefficients(reg)
  
  p <- ggplot(
    VAF_df, aes(x = VAF1, y = VAF2, color = treatment) )
  p <- p + geom_point()
  p <- p + theme_cowplot(12)
  p <- p + xlim(c(0, lim))
  p <- p + ylim(c(0, lim))
  p <- p + xlab(paste0(lab1, " VAF"))
  p <- p + ylab(paste0(lab2, " VAF"))
  p <- p + geom_text_repel(data=VAF_df, aes(label=id), size = 3)
  p <- p + annotate("text", x = lim-0.07, y = lim-0.2, label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", formatC(corr$p.value, format = "e", digits = 2) ), 
    color='red', size = 3.5 )
  p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")
  return(p)
}

# create ddPCR vs avg_Andre VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, ddPCR, avg_Andre))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
ddPCR_vs_Andre_avg <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "Andre_avg", lim = 0.4, cortype = "kendall" )

png(paste0(plot_dir, "ddPCR_vs_avg_Andre_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(ddPCR_vs_Andre_avg)
dev.off()

# create average VAF vs Andre's average VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, avg_Andre, avg_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
Andre_vs_new_avg <- compare_VAF(
  VAF_df, lab1 = "avg Andre", lab2 = "avg new", lim = 0.4, cortype = "kendall" )

png(paste0(plot_dir, "avg_Andre_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(Andre_vs_new_avg)
dev.off()

# extend to include outliers:
Andre_vs_new_avg_ext <- compare_VAF(
  VAF_df, lab1 = "avg Andre", lab2 = "avg new", lim = 1, cortype = "kendall" )

png(paste0(plot_dir, "avg_Andre_vs_new_VAFs_ext.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(Andre_vs_new_avg_ext)
dev.off()

# create up_Andre vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, up_Andre, up_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
up_Andre_vs_new <- compare_VAF(
  VAF_df, lab1 = "Andre upstream", lab2 = "new upstream", lim = 0.4, cortype = "kendall" )

png(paste0(plot_dir, "up_Andre_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(up_Andre_vs_new)
dev.off()

# create down_Andre vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, dn_Andre, dn_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
dn_Andre_vs_new <- compare_VAF(
  VAF_df, lab1 = "Andre downstream", lab2 = "new downstream", lim = 0.5, 
  cortype = "kendall" )

png(paste0(plot_dir, "dn_Andre_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(dn_Andre_vs_new)
dev.off()

# create ddPCR vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, ddPCR, avg_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
ddPCR_vs_new <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "avg new", lim = 0.5, 
  cortype = "kendall" )

png(paste0(plot_dir, "ddPCR_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(ddPCR_vs_new)
dev.off()

# extend to include outliers:
ddPCR_vs_new_ext <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "avg new", lim = 1, 
  cortype = "kendall" )

png(paste0(plot_dir, "ddPCR_vs_new_VAFs_ext.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
print(ddPCR_vs_new_ext)
dev.off()



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

VAFs <- lapply(samplenames, function(x) {
  print(x)
  VAF <- readRDS(paste0(in_path, x, "/Rdata/selected_VAF.rds"))
  if (!is.null(VAF)) {
    if (nrow(VAF) > 0) {
      colnames(VAF) <- gsub("upstream_", "", colnames(VAF))
      colnames(VAF) <- gsub("downstream_", "", colnames(VAF))
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


####################################################################################
### 2. Add previous VAF estimates ###
####################################################################################

prior_VAFs <- read.table(
  paste0(ref_dir, "previous_hepatoblastoma_VAFs.tsv"),
  sep = "\t",
  header = TRUE )
prior_VAFs$id <- gsub("-", "_", prior_VAFs$id)

# format and make values numeric:
all_VAFs <- merge(VAF_df, prior_VAFs, by="id")
all_VAFs$VAF <- all_VAFs$VAF*100
all_VAFs$ddPCR[all_VAFs$ddPCR == "not_identified"] <- 0
all_VAFs$ddPCR <- as.numeric(all_VAFs$ddPCR)
all_VAFs$Andre[all_VAFs$Andre == "not_identified"] <- 0
all_VAFs$Andre <- as.numeric(all_VAFs$Andre)
colnames(all_VAFs) <- gsub("VAF", "new_VAF", colnames(all_VAFs))

# check distributions of VAFs:
plot(hist(all_VAFs$ddPCR[all_VAFs$ddPCR != 0]))
plot(hist(all_VAFs$Andre[all_VAFs$Andre != 0]))
plot(hist(all_VAFs$new_VAF[all_VAFs$new_VAF != 0]))
dev.off()

# function to create VAF correlation plots:
compare_VAF <- function(VAF_df, lab1, lab2, cortype = "kendall" ) {
  # calculate correlation between ddPCR and Andre VAFs:
  corr <- cor.test(x = VAF_df$VAF1, y = VAF_df$VAF2, method = cortype )
  
  # Fit regression line
  require(stats)
  reg <- lm(VAF1 ~ VAF2, data = VAF_df)
  coeff=coefficients(reg)
  
  p <- ggplot(
    VAF_df, aes(x = VAF1, y = VAF2, color = treatment) )
  p <- p + geom_point()
  p <- p + theme_cowplot(12)
  p <- p + xlim(c(0, 50))
  p <- p + ylim(c(0, 50))
  p <- p + xlab(paste0(lab1, " VAF"))
  p <- p + ylab(paste0(lab2, " VAF"))
  p <- p + geom_text_repel(data=VAF_df, aes(label=id))
  p <- p + annotate("text", x = 30, y = 45, label = paste0(
    "R2=", round(corr$estimate, 2), ", p=", round(corr$p.value, 6) ), 
    color='red', size = 4 )
  p <- p + geom_abline(intercept = coeff[1], slope = coeff[2], color = "red")
  return(p)
}

# create ddPCR vs Andre VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, ddPCR, Andre))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
ddPCR_vs_Andre <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "Andre", cortype = "kendall" )

png(paste0(plot_dir, "ddPCR_vs_Andre_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(ddPCR_vs_Andre)
dev.off()

# create ddPCR vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, ddPCR, new_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
ddPCR_vs_new <- compare_VAF(
  VAF_df, lab1 = "ddPCR", lab2 = "new", cortype = "kendall" )

png(paste0(plot_dir, "ddPCR_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(ddPCR_vs_new)
dev.off()

# create Andre vs new VAF plot:
VAF_df <- subset(all_VAFs, select = c(id, treatment, Andre, new_VAF))
colnames(VAF_df) <- c("id", "treatment", "VAF1", "VAF2")
Andre_vs_new <- compare_VAF(
  VAF_df, lab1 = "Andre", lab2 = "new", cortype = "kendall" )

png(paste0(plot_dir, "Andre_vs_new_VAFs.png"),
    height = 3.5,
    width = 5,
    res = 300,
    units = "in")
  print(Andre_vs_new)
dev.off()

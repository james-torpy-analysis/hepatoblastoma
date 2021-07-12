
home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/ewing_ctDNA/")
result_dir <- paste0(project_dir, "results/")
in_path <- paste0(result_dir, "VAF_calculation/")
ref_dir <- paste0(project_dir, "refs/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(ggplot2)
library(cowplot)


####################################################################################
### 1. Load VAFs ###
####################################################################################

samplenames <- as.list(
  gsub(
    "/VAF.txt", "", 
    list.files(in_path, pattern = "VAF.txt", recursive = TRUE)
  )
)

# load metadata and split by sample:
dilution_df <- read.table(
  paste0(ref_dir, "ES_samples_by_dilution.tsv"),
  sep = "\t",
  header = T
)
split_dfs <- split(dilution_df, dilution_df$Original.sample)

sample_list <- list(
  sample1 = samplenames[unlist(samplenames) %in% split_dfs[[1]]$Sample],
  sample2 = samplenames[unlist(samplenames) %in% split_dfs[[2]]$Sample]
)

VAFs <- lapply(sample_list, function(x) {
  return(
    sapply(x, function(y) {
      return(
        read.table(
          paste0(in_path, y, "/VAF.txt"),
          sep = "\t",
          header = F,
          stringsAsFactors = F
        )
      )
    })
  )
})

for (i in 1:length(VAFs)) {
  
  df <- data.frame(
    id = unlist(sample_list[[i]]),
    VAF = unlist(VAFs[[i]])
  )
  m <- match(dilution_df$Sample, df$id)
  m <- m[!is.na(m)]
  df$dilution <- factor(
    dilution_df$Dilution[m],
    levels = sort(dilution_df$Dilution[m], decreasing = TRUE)
  )
  
  if (i==1) {
    plot_dfs <- list(df)
  } else {
    plot_dfs[[i]] = df
  }
  
}
names(plot_dfs) <- names(VAFs)


####################################################################################
### 2. Plot VAFs ###
####################################################################################

VAF_plots <- lapply(plot_dfs, function(x) {
  
  p <- ggplot(x, aes(x = dilution, y = VAF))
  p <- p + geom_bar(stat="identity")
  p <- p + scale_x_discrete()
  p <- p + theme_cowplot(12)
  p <- p + xlab("dilution (%)")
  p <- p + ylab("VAF (%)")
  p <- p + geom_text(
    aes(label=VAF), 
    vjust=-0.25,
    fontface = "bold"
  )
  return(p)
  
})

for (i in 1:length(VAF_plots)) {
  png(paste0(in_path, "serial_dilution_", names(VAF_plots)[i], "_VAFs.png"))
    print(VAF_plots[[i]])
  dev.off()
}



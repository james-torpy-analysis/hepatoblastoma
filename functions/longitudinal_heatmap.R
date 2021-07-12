longitudinal_heatmap <- function(
  fusion_df, 
  hm_title, 
  type,
  annotation = "false positives",
  hm_cols
) {
  
  library(tibble)
  library(ComplexHeatmap)
  library(naturalsort)
  
  if (type == "patient") {

    # split df for each patient:
    split_df <- split(fusion_df, fusion_df$Patient)

    condition_vec <- c(
      "tumour", "naive", "NACT1", "NACT2", 
      "resection", "ACT1", "targeted", "relapse", 
      "ACT2"
    )

  } else if (type == "dilution") {

    # split df for each cell line:
    split_df <- split(fusion_df, fusion_df$Original_sample)

    condition_vec <- c(
      "100", "50", "10", "4", 
      "3", "2", "1", "0.4", 
      "0.2", "0.1", "0.04"
    )

  }
  
  # for each patient, make complete df of all treatments vs detections:
  hm_split <- lapply(split_df, function(x) {

    if (type == "patient") {

      empty_df <- data.frame(
        Treatment = factor(condition_vec, levels = condition_vec)
      )
      
      temp_sub <- subset(
        x, 
        select = c(
          Treatment, Known_EWSR1_FLI1_fusion, Stringent_true_positives, 
          Less_stringent_true_positives, Stringent_false_positives,
          Less_stringent_false_positives, VAF, Supporting_read_pairs
        )
      )
      
      # preserve supporting read number column:
      temp_sub$Supporting_read_pair_nos <- temp_sub$Supporting_read_pairs
        
      # change non-NA/non-zero counts to 'yes' for true positives or supporting 
      # reads present:
      temp_sub$Stringent_true_positives[
        !is.na(temp_sub$Stringent_true_positives) & 
          as.numeric(temp_sub$Stringent_true_positives) > 0
      ] <- "stringent_detection"
      temp_sub$Less_stringent_true_positives[
        !is.na(temp_sub$Less_stringent_true_positives) & 
          as.numeric(temp_sub$Less_stringent_true_positives) > 0
      ] <- "less_stringent_detection"
      temp_sub$Supporting_read_pairs[
        !is.na(temp_sub$Supporting_read_pairs) & 
          temp_sub$Supporting_read_pairs > 1
      ] <- "supporting_reads"
      
      # change NA/non-zero counts to 'no' for true positives:
      temp_sub$Stringent_true_positives[
        is.na(temp_sub$Stringent_true_positives) | 
          temp_sub$Stringent_true_positives == 0
      ] <- "no_stringent_detection"
      temp_sub$Less_stringent_true_positives[
        is.na(temp_sub$Less_stringent_true_positives) | 
          temp_sub$Less_stringent_true_positives == 0
      ] <- "no_less_stringent_detection"
      temp_sub$Supporting_read_pairs[
        is.na(temp_sub$Supporting_read_pairs) | 
          temp_sub$Supporting_read_pairs <= 1
      ] <- "no_supporting_reads"
      
      # change NA counts to 0 for false positives:
      temp_sub$Stringent_false_positives[
        is.na(temp_sub$Stringent_false_positives)
      ] <- 0
      temp_sub$Less_stringent_false_positives[
        is.na(temp_sub$Less_stringent_false_positives)
      ] <- 0
      
      merged_df <- merge(
        empty_df,
        temp_sub,
        by = "Treatment",
        all = T
      )
      merged_df <- merged_df %>%
        column_to_rownames("Treatment")

    } else if (type == "dilution") {

      empty_df <- data.frame(
        Dilution = factor(condition_vec, levels = condition_vec)
      )
      
      temp_sub <- subset(
        x, 
        select = c(
          Dilution, Known_EWSR1_FLI1_fusion, Stringent_true_positives, 
          Less_stringent_true_positives, Stringent_false_positives,
          Less_stringent_false_positives, VAF, Supporting_read_pairs
        )
      )
      
      # preserve supporting read number column:
      temp_sub$Supporting_read_pair_nos <- temp_sub$Supporting_read_pairs
      
      # change non-NA/non-zero counts to 'yes' for true positives or supporting 
      # reads present:
      temp_sub$Stringent_true_positives[
        !is.na(temp_sub$Stringent_true_positives) & 
          as.numeric(temp_sub$Stringent_true_positives) > 0
      ] <- "stringent_detection"
      temp_sub$Less_stringent_true_positives[
        !is.na(temp_sub$Less_stringent_true_positives) & 
          as.numeric(temp_sub$Less_stringent_true_positives) > 0
      ] <- "less_stringent_detection"
      temp_sub$Supporting_read_pairs[
        !is.na(temp_sub$Supporting_read_pairs) & 
          temp_sub$Supporting_read_pairs > 1
      ] <- "supporting_reads"
      
      
      # change NA/non-zero counts to 'no' for true positives:
      temp_sub$Stringent_true_positives[
        is.na(temp_sub$Stringent_true_positives) | 
          temp_sub$Stringent_true_positives == 0
      ] <- "no_stringent_detection"
      temp_sub$Less_stringent_true_positives[
        is.na(temp_sub$Less_stringent_true_positives) | 
          temp_sub$Less_stringent_true_positives == 0
      ] <- "no_less_stringent_detection"
      temp_sub$Supporting_read_pairs[
        is.na(temp_sub$Supporting_read_pairs) | 
          temp_sub$Supporting_read_pairs <= 1
      ] <- "no_supporting_reads"
      
      # change NA counts to 0 for false positives:
      temp_sub$Stringent_false_positives[
        is.na(temp_sub$Stringent_false_positives)
      ] <- 0
      temp_sub$Less_stringent_false_positives[
        is.na(temp_sub$Less_stringent_false_positives)
      ] <- 0
      
      merged_df <- merge(
        empty_df,
        temp_sub,
        by = "Dilution",
        all = T
      )

      # add 'b' to any duplicated dilutions:
      temp_dilution <- as.character(merged_df$Dilution)
      temp_dilution[duplicated(temp_dilution)] <- paste0(
        temp_dilution[duplicated(temp_dilution)], "b"
      )
      merged_df$Dilution <- temp_dilution
      merged_df <- merged_df %>%
        column_to_rownames("Dilution")

    }

    # order merged_df:
    m <- match(condition_vec, rownames(merged_df))
    merged_df <- merged_df[m,]
    
    # record FISH results for first sample:
    FISH = merged_df$Known_EWSR1_FLI1_fusion[
      !is.na(merged_df$Known_EWSR1_FLI1_fusion)
    ][1]
    
    # make FISH entries distinct from others and merge with stringent calls:
    if (FISH == "yes") {
      FISH <- "FISH_detection"
    } else {
      FISH <- "no_FISH_detection"
    }

    detection_df <- subset(
      merged_df, 
      select = c(Stringent_true_positives, Less_stringent_true_positives, Supporting_read_pairs)
    )
    detection_df <- rbind(
      data.frame(
        row.names = "FISH", Stringent_true_positives = FISH, Less_stringent_true_positives = FISH,
        Supporting_read_pairs = FISH
      ),
      detection_df
    )
    
    # for those samples with no stringent calls, fetch non-stringent calls:
    temp_df <- detection_df[apply(detection_df, 1, function(x) any(!is.na(x))),]
    temp_df$Stringent_true_positives[
      temp_df$Stringent_true_positives == "no_stringent_detection"
    ] <- temp_df$Less_stringent_true_positives[
        temp_df$Stringent_true_positives == "no_stringent_detection"
      ]
    
    # for those samples with no less stringent calls, fetch supporting read pair no:
    temp_df$Stringent_true_positives[
      temp_df$Stringent_true_positives == "no_less_stringent_detection"
    ] <- temp_df$Supporting_read_pairs[
      temp_df$Stringent_true_positives == "no_less_stringent_detection"
    ]
    
    # merge with missing samples:
    detection_df <- rbind(
      temp_df,
      detection_df[apply(detection_df, 1, function(x) all(is.na(x))),]
    )
    
    # keep only detections column:
    detection_df <- subset(
      detection_df, select = "Stringent_true_positives"
    )
    colnames(detection_df) <- "Detections"
    
    # ensure correct row order:
    detection_df <- detection_df %>%
      rownames_to_column("type") %>%
      arrange(
        factor(
          type, levels = c(
            "FISH", "tumour", "naive", "NACT1", "NACT2", "resection", 
            "ACT1", "targeted", "relapse", "ACT2"
          )
        )
      ) %>%
      column_to_rownames("type")
    
    # create annotation df:
    if (annotation == "false positives") {
      
      # create a false positive df in parallel:
      annot_df <- subset(
        merged_df, 
        select = c(Stringent_false_positives, Less_stringent_false_positives)
      )
      temp_bind <- as.data.frame(
        t(data.frame(rep(0, ncol(annot_df))))
      )
      colnames(temp_bind) <- colnames(annot_df)
      rownames(temp_bind) <- "FISH"
      annot_df <- rbind(
        temp_bind,
        annot_df
      )
      
      # for samples with less stringent calls, update false positives column
      # with less stringent false positives:
      temp_df <- annot_df[apply(annot_df, 1, function(x) any(!is.na(x))),]
      detection_vec <- detection_df[
        rownames(detection_df) %in% rownames(temp_df),
      ]
      temp_df$Stringent_false_positives[
        detection_vec == "less_stringent_detection"
      ] <- temp_df$Less_stringent_false_positives[
        detection_vec == "less_stringent_detection"
      ]
      
      # merge with missing samples:
      annot_df <- rbind(
        temp_df,
        annot_df[apply(annot_df, 1, function(x) all(is.na(x))),]
      )
      
      # keep only detections column:
      annot_df <- subset(
        annot_df, select = "Stringent_false_positives"
      )
      colnames(annot_df) <- "False_positives"
      
      # ensure correct row order:
      annot_df <- annot_df %>%
        rownames_to_column("type") %>%
        arrange(
          factor(
            type, levels = c(
              "FISH", "tumour", "naive", "NACT1", "NACT2", "resection", 
              "ACT1", "targeted", "relapse", "ACT2"
            )
          )
        ) %>%
        column_to_rownames("type")
      
      return(
        list(
          detection_df = detection_df, 
          annot_df = annot_df
        )
      )
      
    } else {
      
      # create a VAF df:
      VAF_annot_df <- subset(
        merged_df, 
        select = VAF
      )
      temp_bind <- as.data.frame(
        t(data.frame(rep(0, ncol(VAF_annot_df))))
      )
      colnames(temp_bind) <- colnames(VAF_annot_df)
      rownames(temp_bind) <- "FISH"
      VAF_annot_df <- rbind(
        temp_bind,
        VAF_annot_df
      )
      
      # put into list:
      annot_dfs <- list(VAF = VAF_annot_df)
      
      # create sread_annot_df:
      sread_annot_df <- subset(
        merged_df, 
        select = Supporting_read_pair_nos
      )
      temp_bind <- as.data.frame(
        t(data.frame(rep(0, ncol(sread_annot_df))))
      )
      colnames(temp_bind) <- colnames(sread_annot_df)
      rownames(temp_bind) <- "FISH"
      sread_annot_df <- rbind(
        temp_bind,
        sread_annot_df
      )
      
      # put into list:
      annot_dfs$sread <- sread_annot_df
      
      # ensure correct row orders:
      annot_dfs <- lapply(annot_dfs, function(y) {
        y %>%
          rownames_to_column("type") %>%
          arrange(
            factor(
              type, levels = c(
                "FISH", "tumour", "naive", "NACT1", "NACT2", "resection", 
                "ACT1", "targeted", "relapse", "ACT2"
              )
            )
          ) %>%
          column_to_rownames("type")
        
        return(y)
        
      })
        
      return(
        list(
          detection_df = detection_df, 
          VAF_annot_df = annot_dfs$VAF,
          sread_annot_df = annot_dfs$sread
        )
      )

    }
      
  })
  
  if (annotation == "false positives") {
    
    both_split <- list(
      detection_df = lapply(hm_split, function(x) x$detection_df),
      annot_df = lapply(hm_split, function(x) x$annot_df)
    )
    
    # bind together:
    hm_dfs <- lapply(both_split, function(x) {
      
      hm_df <- as.data.frame(t(do.call("cbind", x)))
      
      if (type == "patient") {
        
        # format data to order by met and FISH status:
        hm_df$Patient <- names(hm_split)
        met_status <- data.frame(
          Patient = fusion_df$Patient,
          Site = fusion_df$Site
        )
        met_status <- met_status[!duplicated(met_status$Patient),]
        hm_df <- merge(hm_df, met_status, by = "Patient")
        
        # split by primary/met:
        temp_split <- split(hm_df, hm_df$Site)
        
        if (!all(is.na(hm_df$FISH))) {
          temp_split <- lapply(temp_split, function(y) {
            # change all NAs to "unknown":
            y[is.na(y)] <- "unknown"
            # make FISH column a factor with required order:
            y$FISH <- factor(y$FISH, levels = c("FISH_detection", "no_FISH_detection", "unknown"))
            # order each split by FISH column:
            y <- y[order(y$FISH),]
            return(y)
          })
        }
        names(temp_split) <- names(temp_split)
        
        # bind together in met status order:
        temp_split <- list(
          temp_split$primary,
          temp_split$met,
          temp_split$unknown
        )
        hm_df <- do.call("rbind", temp_split)
        # save local/met vector and move patient ids to rownames:
        met_order <- factor(hm_df$Site, levels = unique(hm_df$Site))
        hm_df <- subset(hm_df, select = -Site)
        rownames(hm_df) <- NULL
        hm_df <- hm_df %>%
          column_to_rownames("Patient")
        
        return(list(hm_df = hm_df, met_order = met_order))
        
      } else if (type == "dilution") {
        
        rownames(hm_df) <- paste0("Cell line sample ", names(hm_split))
        
        return(list(hm_df = hm_df))
        
      }
      
    })
    
    # adjust order of annot dataframe:
    hm_dfs$annot_df$hm_df <- hm_dfs$annot_df$hm_df[
      rownames(hm_dfs$detection_df$hm_df), colnames(hm_dfs$detection_df$hm_df)
    ]
    
    # change all NA or 0 values in annot_df to spaces:
    final_annot <- apply(hm_dfs$annot_df$hm_df, 2, function(x) {
      x[is.na(x)] <- " "
      x[x == 0] <- " "
      x[x == "unknown"] <- " "
      return(x)
    })
    
    # create treatment_split vector:
    treatment_split <- c("FISH", rep("NOT_FISH", ncol(hm_dfs$detection_df$hm_df) - 1))
    
    # create heatmaps:
    if (type == "patient") {
      
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, 
              gp = gpar(fontsize = 10, fontface = "bold", col = "#991425")
            )
          }
        )
      )
      
    } else if (type == "dilution") {
      
      return(
        Heatmap(
          as.matrix(hm_dfs$detection_df$hm_df), 
          name = "Fusion detections", 
          na_col = "grey",
          row_split = hm_dfs$detection_df$met_order,
          column_split = treatment_split,
          col = hm_cols,
          border = "black",
          rect_gp = gpar(col = "black", lwd = 1),
          column_title = hm_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
              final_annot[i, j], x, y, gp = gpar(
                fontsize = 10, fontface = "bold", col = "#991425"
              )
            )
          }
        )
      )
    }
    
  } else if (annotation == "VAF") {
    
    all_split <- list(
      detection_df = lapply(hm_split, function(x) x$detection_df),
      VAF_annot_df = lapply(hm_split, function(x) x$VAF_annot_df),
      sread_annot_df = lapply(hm_split, function(x) x$sread_annot_df)
    )
    
    # bind together:
    hm_dfs <- lapply(all_split, function(x) {
      
      hm_df <- as.data.frame(t(do.call("cbind", x)))
      
      if (type == "patient") {
        
        # format data to order by met and FISH status:
        hm_df$Patient <- names(hm_split)
        met_status <- data.frame(
          Patient = fusion_df$Patient,
          Site = fusion_df$Site
        )
        met_status <- met_status[!duplicated(met_status$Patient),]
        hm_df <- merge(hm_df, met_status, by = "Patient")
        
        # split by primary/met:
        temp_split <- split(hm_df, hm_df$Site)
        
        if (!all(is.na(hm_df$FISH))) {
          temp_split <- lapply(temp_split, function(y) {
            # change all NAs to "unknown":
            y[is.na(y)] <- "unknown"
            # make FISH column a factor with required order:
            y$FISH <- factor(y$FISH, levels = c("FISH_detection", "no_FISH_detection", "unknown"))
            # order each split by FISH column:
            y <- y[order(y$FISH),]
            return(y)
          })
        }
        names(temp_split) <- names(temp_split)
        
        # put in in met status order:
        temp_split <- list(
          temp_split$primary,
          temp_split$met,
          temp_split$unknown
        )
        hm_df <- do.call("rbind", temp_split)
        # save local/met vector and move patient ids to rownames:
        met_order <- factor(hm_df$Site, levels = unique(hm_df$Site))
        hm_df <- subset(hm_df, select = -Site)
        rownames(hm_df) <- NULL
        hm_df <- hm_df %>%
          column_to_rownames("Patient")
        
        return(list(hm_df = hm_df, met_order = met_order))
        
      } else if (type == "dilution") {
        
        rownames(hm_df) <- paste0("Cell line sample ", names(hm_split))
        
        return(list(hm_df = hm_df))
        
      }
      
    })
    
    # adjust order of annot dataframe:
    annot_dfs <- list(
      VAF = hm_dfs$VAF_annot_df$hm_df[
        rownames(hm_dfs$detection_df$hm_df), colnames(hm_dfs$detection_df$hm_df)
      ],
      sread = hm_dfs$sread_annot_df$hm_df[
        rownames(hm_dfs$detection_df$hm_df), colnames(hm_dfs$detection_df$hm_df)
      ]
    )
    
    # adjust order of annot dataframes and change all NA or 0 values in 
    # annot_df to spaces:
    final_annots <- lapply(annot_dfs, function(x) {
      x <- x[
        rownames(hm_dfs$detection_df$hm_df), 
        colnames(hm_dfs$detection_df$hm_df)
      ]
      apply(x, 2, function(y) {
        y[is.na(y)] <- " "
        y[y == "unknown"] <- " "
        return(y)
      })
    })
 
    # create treatment_split vector:
    treatment_split <- c("FISH", rep("NOT_FISH", ncol(hm_dfs$detection_df$hm_df) - 1))
    
    # create heatmaps:
    if (type == "patient") {
      
      return(
        list(
          VAF_annot = Heatmap(
            as.matrix(hm_dfs$detection_df$hm_df), 
            name = "Fusion detections", 
            row_split = hm_dfs$detection_df$met_order,
            column_split = treatment_split,
            col = hm_cols,
            border = "black",
            rect_gp = gpar(col = "black", lwd = 1),
            column_title = hm_title,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(
                final_annots$VAF[i, j], x, y, 
                gp = gpar(fontsize = 10, fontface = "bold", col = "#221699")
              )
            }
          ),
          sread_annot = Heatmap(
            as.matrix(hm_dfs$detection_df$hm_df), 
            name = "Fusion detections", 
            row_split = hm_dfs$detection_df$met_order,
            column_split = treatment_split,
            col = hm_cols,
            border = "black",
            rect_gp = gpar(col = "black", lwd = 1),
            column_title = hm_title,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(
                final_annots$sread[i, j], x, y, 
                gp = gpar(fontsize = 10, fontface = "bold", col = "#221699")
              )
            }
          )
        )
      )
      
    } else if (type == "dilution") {
      
      return(
        list(
          VAF_annot = Heatmap(
            as.matrix(hm_dfs$detection_df$hm_df), 
            name = "Fusion detections", 
            row_split = hm_dfs$detection_df$met_order,
            column_split = treatment_split,
            col = hm_cols,
            border = "black",
            rect_gp = gpar(col = "black", lwd = 1),
            column_title = hm_title,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(
                final_annots$VAF[i, j], x, y, 
                gp = gpar(fontsize = 10, fontface = "bold", col = "#221699")
              )
            }
          ),
          sread_annot = Heatmap(
            as.matrix(hm_dfs$detection_df$hm_df), 
            name = "Fusion detections", 
            row_split = hm_dfs$detection_df$met_order,
            column_split = treatment_split,
            col = hm_cols,
            border = "black",
            rect_gp = gpar(col = "black", lwd = 1),
            column_title = hm_title,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(
                final_annots$sread[i, j], x, y, 
                gp = gpar(fontsize = 10, fontface = "bold", col = "#115E0F")
              )
            }
          )
        )
      )
    }
    
  }
  
}

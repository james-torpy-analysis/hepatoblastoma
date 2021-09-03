longitudinal_patient_heatmap <- function(
  SV_df, 
  hm_title,
  annot = "false positives",
  condition_order,
  hm_cols,
  func_dir
) {
  
  library(tibble)
  library(ComplexHeatmap)
  library(naturalsort)
  
  annotate_detection <- dget(paste0(func_dir, "annotate_detection.R"))
  
  # split df for each patient:
  split_df <- split(SV_df, SV_df$Patient)

  # for each patient, make complete df of all treatments vs detections:
  hm_split <- lapply(split_df, function(x) {
    
    empty_df <- data.frame(
      Condition = factor(condition_order, levels = condition_order)
    )
    
    temp_sub <- subset(
      x, 
      select = c(
        Condition, Pathology_detection, Stringent_true_positives, 
        Less_stringent_true_positives, Stringent_false_positives,
        Less_stringent_false_positives, VAF
      )
    )

    # change non-NA/non-zero true positive counts to detection status:
    temp_sub$Stringent_true_positives <- annotate_detection(
      temp_sub$Stringent_true_positives, label = "stringent_detection"
    )
    temp_sub$Less_stringent_true_positives <- annotate_detection(
      temp_sub$Less_stringent_true_positives, label = "less_stringent_detection"
    )
      
    merged_df <- merge(
      empty_df,
      temp_sub,
      by = "Condition",
      all = T
    )
    merged_df <- merged_df %>%
      column_to_rownames("Condition")

    # order merged_df:
    m <- match(condition_order, rownames(merged_df))
    merged_df <- merged_df[m,]
    
    # record pathology results for first sample:
    pathology = merged_df$Pathology_detection[
      !is.na(merged_df$Pathology_detection)
    ][1]
    
    # make pathology entries distinct from others and merge with stringent calls:
    if (pathology == "yes") {
      pathology <- "pathology_detection"
    } else if (pathology == "no") {
      pathology <- "no_pathology_detection"
    } else {
      pathology <- "unknown"
    }

    detection_df <- subset(
      merged_df, 
      select = c(Stringent_true_positives, Less_stringent_true_positives)
    )
    detection_df <- rbind(
      data.frame(
        row.names = "pathology", Stringent_true_positives = pathology, 
        Less_stringent_true_positives = pathology
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
      arrange(factor(type, levels = c("pathology", condition_order))) %>%
      column_to_rownames("type")
    
    # create annotation df:
    if (annot == "false positives") {
      
      # create a false positive df in parallel:
      annot_df <- subset(
        merged_df, 
        select = c(Stringent_false_positives, Less_stringent_false_positives)
      )
      temp_bind <- as.data.frame(
        t(data.frame(rep(0, ncol(annot_df))))
      )
      colnames(temp_bind) <- colnames(annot_df)
      rownames(temp_bind) <- "pathology"
      annot_df <- rbind(temp_bind,annot_df)
      
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
        arrange(factor(type, levels = rownames(detection_df))) %>%
        column_to_rownames("type")
      
      return(
        list(
          detection_df = detection_df, 
          annot_df = annot_df
        )
      )
      
    } else {
      
      # create a VAF df:
      annot_df <- subset(
        merged_df, 
        select = VAF
      )
      temp_bind <- as.data.frame(
        t(data.frame(rep(0, ncol(annot_df))))
      )
      colnames(temp_bind) <- colnames(annot_df)
      rownames(temp_bind) <- "pathology"
      annot_df <- rbind(
        temp_bind,
        annot_df
      )
      
      # put into list:
      annot_dfs <- list(VAF = annot_df)
      
      # ensure correct row orders:
      annot_dfs <- lapply(annot_dfs, function(y) {
        y %>%
          rownames_to_column("type") %>%
          arrange(factor(type, levels = c("pathology", condition_order))) %>%
          column_to_rownames("type")
        
        return(y)
        
      })
        
      return(
        list(
          detection_df = detection_df, 
          annot_df = annot_dfs$VAF
        )
      )

    }
      
  })
  
  if (annot == "false positives") {
    
    both_split <- list(
      detection_df = lapply(hm_split, function(x) x$detection_df),
      annot_df = lapply(hm_split, function(x) x$annot_df)
    )
    
    # bind together:
    hm_dfs <- lapply(both_split, function(x) {
      
      hm_df <- as.data.frame(t(do.call("cbind", x)))
        
      # format data to order by met and pathology status:
      hm_df$Patient <- names(hm_split)
      met_status <- data.frame(
        Patient = SV_df$Patient,
        Site = SV_df$Site
      )
      met_status <- met_status[!duplicated(met_status$Patient),]
      hm_df <- merge(hm_df, met_status, by = "Patient")
      
      # split by primary/met:
      temp_split <- split(hm_df, hm_df$Site)
      
      if (!all(is.na(hm_df$pathology))) {
        temp_split <- lapply(temp_split, function(y) {
          # change all NAs to "unknown":
          y[is.na(y)] <- "unknown"
          # make pathology column a factor with required order:
          y$pathology <- factor(
            y$pathology, 
            levels = c("pathology_detection", "no_pathology_detection", 
            "unknown"))
          # order each split by pathology column:
          y <- y[order(y$pathology),]
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
    
    # change 'no_less_stringent_detection' to 'no_detections':
    hm_dfs$detection_df$hm_df <- apply(hm_dfs$detection_df$hm_df, 2, function(x) {
      x[x == "no_less_stringent_detection"] <- "no_detection"
      return(x)
    })
    
    # change df order to alpheabetical:
    hm_dfs$detection_df$hm_df <- hm_dfs$detection_df$hm_df[
      order(rownames(hm_dfs$detection_df$hm_df)),]
    final_annot <- final_annot[
      match(rownames(hm_dfs$detection_df$hm_df), rownames(final_annot)),]
    
    # create treatment_split vector:
    condition_split <- factor(
      c(
        "pathology", 
        rep("not_pathology", ncol(hm_dfs$detection_df$hm_df) - 1)
      ),
      levels = c("pathology", "not_pathology")
    )
    
    # create heatmap:
    return(
      Heatmap(
        as.matrix(hm_dfs$detection_df$hm_df), 
        name = "SV detections", 
        row_split = hm_dfs$detection_df$met_order,
        column_split = condition_split,
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
    
  } else if (annot == "VAF") {
    
    all_split <- list(
      detection_df = lapply(hm_split, function(x) x$detection_df),
      annot_df = lapply(hm_split, function(x) x$annot_df)
    )
    
    # bind together:
    hm_dfs <- lapply(all_split, function(x) {
      
      hm_df <- as.data.frame(t(do.call("cbind", x)))
      
      # format data to order by met and pathology status:
      hm_df$Patient <- names(hm_split)
      met_status <- data.frame(
        Patient = SV_df$Patient,
        Site = SV_df$Site
      )
      met_status <- met_status[!duplicated(met_status$Patient),]
      hm_df <- merge(hm_df, met_status, by = "Patient")
      
      # split by primary/met:
      temp_split <- split(hm_df, hm_df$Site)
      
      if (!all(is.na(hm_df$pathology))) {
        temp_split <- lapply(temp_split, function(y) {
          # change all NAs to "unknown":
          y[is.na(y)] <- "unknown"
          # make pathology column a factor with required order:
          y$pathology <- factor(y$pathology, levels = c("pathology_detection", "no_pathology_detection", "unknown"))
          # order each split by pathology column:
          y <- y[order(y$pathology),]
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
      
    })
    
    # change 'no_less_stringent_detection' to 'no_detections':
    hm_dfs$detection_df$hm_df <- apply(hm_dfs$detection_df$hm_df, 2, function(x) {
      x[x == "no_less_stringent_detection"] <- "no_detection"
      return(x)
    })
    
    # adjust order of annot dataframe:
    final_annot_df <- hm_dfs$annot_df$hm_df[
        rownames(hm_dfs$detection_df$hm_df), 
        colnames(hm_dfs$detection_df$hm_df)]
    
    # change all NA or 0 values in final_annot_df to spaces:
    final_annot_df <- apply(final_annot_df, 2, function(y) {
      y[is.na(y)] <- " "
      y[y == 0] <- " "
      y[y == "unknown"] <- " "
      y[y != " "] <- round(as.numeric(y[y != " "]), 1)
      y[y != " "] <- paste0(y[y != " "], "%")
      return(y)
    })

    # create condition_split vector:
    condition_split <- factor(
      c(
        "pathology", 
        rep("not_pathology", ncol(hm_dfs$detection_df$hm_df) - 1)
      ),
      levels = c("pathology", "not_pathology")
    )
    
    # create heatmaps:
    return(
      Heatmap(
        as.matrix(hm_dfs$detection_df$hm_df), 
        name = "SV detections", 
        row_split = hm_dfs$detection_df$met_order,
        column_split = condition_split,
        col = hm_cols,
        border = "black",
        rect_gp = gpar(col = "black", lwd = 1),
        column_title = hm_title,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(
            final_annot_df[i, j], x, y, 
            gp = gpar(fontsize = 9, fontface = "bold", col = "#221699")
          )
        }
      )
    )
      
  }
  
}

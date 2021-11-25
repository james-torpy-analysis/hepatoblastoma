mutation_heatmap <- function(
  mutation_df,
  type,
  hm_cols
) {
  
  library(naturalsort)
  library(reshape2)
  library(ComplexHeatmap)
  
  round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  
  # create VAF and detection dfs:
  VAF_df <- subset(mutation_df, select = c(Patient_id,  Treatment.dilution))
  VAF_df$mutation_VAF <- apply(mutation_df, 1, function(x) {
    x <- x[grepl("smCounter2|Deletion", names(x)) & grepl("VAF", names(x))]
    x <- x[grep("not_detected", x, invert = TRUE)]
    
    if (length(x) > 0) {
      # round VAFs:
      VAF <- round2(as.numeric(x), 1)
      names(VAF) <- names(x)
      for (i in seq_along(VAF)) {
        if (names(VAF)[i] == "smCounter2_VAF") {
          mname <- "SNV"
        } else {
          mname <- "deletion"
        }
        if (i==1) {
          annot <- paste0(mname, ": ", VAF[i], "%" )
        } else {
          annot <- paste0(annot, "\n", mname, ": ", VAF[i], "%" )
        }
      }
      return(annot)
    } else {
      return("")
    }
    
  })
  
  detect_df <- VAF_df
  
  # annotate detected/not_detected values:
  detect_df$mutation_VAF[detect_df$mutation_VAF == ""] <- "not_detected"
  detect_df$mutation_VAF[detect_df$mutation_VAF != "not_detected"] <- "detected"

  # cast to wide format:
  VAF_df <- dcast(
    VAF_df, Patient_id ~Treatment.dilution, value.var = "mutation_VAF" )
  detect_df <- dcast(
    detect_df, Patient_id ~Treatment.dilution, value.var = "mutation_VAF" )
  
  # make rownames patient ids and add pathology column:
  rownames(VAF_df) <- VAF_df$Patient_id
  VAF_df$pathology <- ""
  VAF_df <- VAF_df[
    ,colnames(VAF_df) %in% 
      c("pathology", unique(as.character(mutation_df$Treatment.dilution))) ]
  VAF_df <- VAF_df %>%
    select(pathology, everything())
  
  # merge pathology results with dfs:
  rownames(detect_df) <- detect_df$Patient_id
  detect_df <- subset(detect_df, select = -Patient_id)
  path_df <- subset(mutation_df, select = c(
    Patient_id, Sanger_point_mut, Deletion ))
  path_df <- path_df[!duplicated(path_df$Patient_id),]
  
  # if were are any mutations picked up, indicate "pathology_detection":
  path_df$pathology <- "no_pathology_detection"
  path_df$pathology[apply(
    subset(path_df, select = -c(Patient_id, pathology)), 1, function(x) {
      any(x != "not_detected" & x != "pending")
    })] <- "pathology_detection"
  path_df$pathology[apply(
    subset(path_df, select = -c(Patient_id, pathology)), 1, function(x) {
      all(is.na(x))
    })] <- "no_sample"
  
  # order as in detect_df:
  path_df <- path_df[match(rownames(detect_df), path_df$Patient_id),]
  rownames(path_df) <- path_df$Patient_id
  path_df <- subset(path_df, select = pathology)
  
  # bind to detect_df:
  detect_df <- cbind(path_df, detect_df)
  
  # order dfs by name and pathology status:
  if (type == "patient") {
    detect_df <- as.data.frame(
      detect_df[naturalorder(gsub("P|HB", "", rownames(detect_df))),] )
    detect_df$pathology <- factor(detect_df$pathology, 
      levels = c("pathology_detection", "no_pathology_detection", "no_sample") ) 
    detect_df <- detect_df[order(detect_df$pathology),]
  } else {
    detect_df <- detect_df[order(rownames(detect_df)),]
  }
  
  # define pathology(column) split vector:
  path_split <- factor(
    c("path", rep("non_path", length(unique(mutation_df$Treatment.dilution)))),
    levels = c("path", "non_path") )
  
  VAF_df <- VAF_df[match(rownames(detect_df), rownames(VAF_df)),]
  
  if (nrow(detect_df) > 1) {
    
    # make NAs 'no_sample':
    detect_df <- apply(detect_df, 2, function(x) {
      x[is.na(x)] <- "no_sample"
      return(x)
    })
    
    # remove underscores:
    detect_df <- gsub("_", " ", detect_df)
    names(hm_cols) <- gsub("_", " ", names(hm_cols))
    
    # make NAs blank:
    VAF_df <- apply(VAF_df, 2, function(x) {
      x[is.na(x)] <- ""
      return(x)
    })
    
  } else {
    
    # make NAs 'no_sample':
    detect_df[1,][is.na(detect_df[1,])] <- "no_sample"
    
    # remove underscores:
    detect_df[1,] <- gsub("_", " ", detect_df[1,])
    names(hm_cols) <- gsub("_", " ", names(hm_cols))
    
    # make NAs blank:
    VAF_df[1,][is.na(VAF_df[1,])] <- ""
    
  }
  
  if (type == "dilution") {
    # order by dilution level:
    detect_df <- detect_df[, 
      c("pathology", "100", "50", "10", "2", "1", "0.4", "0.2", "0.1", "0.04") ]
    VAF_df <- VAF_df[,
      c("pathology", "100", "50", "10", "2", "1", "0.4", "0.2", "0.1", "0.04") ]
    # add percentages to dilutions:
    colnames(detect_df)[colnames(detect_df) != "pathology"] <- paste0(
      colnames(detect_df)[colnames(detect_df) != "pathology"], "%" )
    colnames(VAF_df)[colnames(VAF_df) != "pathology"] <- paste0(
      colnames(VAF_df)[colnames(VAF_df) != "pathology"], "%" )
  }
  
  # create VAF heatmap:
  if (type == "patient") {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Mutation detections",
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient CTNNB1 mutation detections",
      column_title_gp = gpar(fontsize = 18, fontface = "bold")
      ,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          VAF_df[i, j], x, y, 
          gp = gpar(fontsize = 10.5, fontface = "bold", col = "#991425") )
      }, )
  } else {
    VAF_hm <- Heatmap(
      as.matrix(detect_df), 
      name = "Mutation detections",
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Cell line CTNNB1 mutation detections",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          VAF_df[i, j], x, y, 
          gp = gpar(fontsize = 14, fontface = "bold", col = "#991425") )
      }, )
  }
  
  # create supporting read annotation:
  supp_df <- subset(mutation_df, select = c(Patient_id, Treatment.dilution,
    Total_supporting ))
  
  # make NA or not detected values blank:
  supp_df$Total_supporting[
    is.na(supp_df$Total_supporting) | supp_df$Total_supporting == "not_detected"
  ] <- " "
  
  # cast to wide format:
  supp_df <- dcast(
    supp_df, Patient_id ~Treatment.dilution, value.var = "Total_supporting" )
  
  if (type == "dilution") {
    # add percentages to dilutions:
    colnames(supp_df)[!(colnames(supp_df) %in% c("Patient_id", "Site"))] <- paste0(
      colnames(supp_df)[!(colnames(supp_df) %in% c("Patient_id", "Site"))], "%" )
  }
  
  # make rownames patient ids and add pathology column:
  rownames(supp_df) <- supp_df$Patient_id
  supp_df$pathology <- " "
  supp_df <- subset(supp_df, select = -c(Patient_id))
  
  # remove underscores:
  rownames(supp_df) <- gsub("_", " ", rownames(supp_df))
  
  # order df:
  supp_df <- supp_df[match(rownames(detect_df), rownames(supp_df)), 
    match(colnames(detect_df), colnames(supp_df)) ]
  
  if (nrow(supp_df) > 1) {
    
    # make NAs blank:
    supp_df <- apply(supp_df, 2, function(x) {
      x[is.na(x)] <- " "
      return(x)
    })
    
  } else {
    
   supp_df[1,][is.na(supp_df[1,])] <- " "
    
  }
  
  # remake detection to include only fusions:
  supp_detect <- detect_df
  supp_detect[supp_df == " " & supp_detect == "detected"] <- "not detected"
  
  # create supporting reads heatmap:
  if (type == "patient") {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads",
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Patient CTNNB1 deletion-supporting reads",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 12, fontface = "bold", col = "#430F82") )
      },
    )
  } else {
    sread_hm <- Heatmap(
      as.matrix(supp_detect), 
      name = "Fusion-supporting reads", 
      column_split = path_split,
      col = hm_cols,
      border = "black",
      rect_gp = gpar(col = "black", lwd = 1),
      column_title = "Cell line CTNNB1 deletion-supporting reads",
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          supp_df[i, j], x, y, 
          gp = gpar(fontsize = 14, fontface = "bold", col = "#430F82") )
      },
    )
  }

  # return heatmaps:
  return(list(VAF=VAF_hm, sread=sread_hm))

}

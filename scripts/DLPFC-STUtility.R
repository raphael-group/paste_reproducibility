# setwd("path/to/paste_reproducibility/scripts/") change this to work on your machine
library("STutility")
library("dplyr")
library("glue")
library("imager") # needed for STutility...
library(data.table) # for fwrite
set.seed(123)

get_info_table <- function(pat) {
  fp <- glue("DLPFC/git-orig-files/{pat[1]}/")
  matnames <- c("_filtered_feature_bc_matrix.h5")
  s1 <- c(glue(fp, "{pat[2]}/{pat[2]}", matnames), 
          glue(fp, "{pat[2]}/tissue_positions_list.txt"),
          glue(fp, "{pat[2]}/tissue_hires_image.png"),
          glue(fp, "{pat[2]}/scalefactors_json.json"))
  s2 <- c(glue(fp, "{pat[3]}/{pat[3]}", matnames), 
          glue(fp, "{pat[3]}/tissue_positions_list.txt"),
          glue(fp, "{pat[3]}/tissue_hires_image.png"),
          glue(fp, "{pat[3]}/scalefactors_json.json"))
  s3 <- c(glue(fp, "{pat[4]}/{pat[4]}", matnames), 
          glue(fp, "{pat[4]}/tissue_positions_list.txt"),
          glue(fp, "{pat[4]}/tissue_hires_image.png"),
          glue(fp, "{pat[4]}/scalefactors_json.json"))
  s4 <- c(glue(fp, "{pat[5]}/{pat[5]}", matnames), 
          glue(fp, "{pat[5]}/tissue_positions_list.txt"),
          glue(fp, "{pat[5]}/tissue_hires_image.png"),
          glue(fp, "{pat[5]}/scalefactors_json.json"))
  df <- data.frame(t(cbind(s1,s2,s3,s4)))
  colnames(df) <- c("samples", "spotfiles", "imgs", "json")
  df 
}

align_se <- function(infotab, thres=TRUE) {
  se <- InputFromTable(infotable = infotab, 
                       minUMICountsPerGene = 100, 
                       minUMICountsPerSpot = 100,
                       platform =  "Visium")
  se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)
  se <- se %>% MaskImages() %>% AlignImages() # default thresholding is set to TRUE
  se
}

get_coords <- function(se){
  spots <- colnames(x = se)
  data <- FetchData(object = se, vars = c("MFGE8"), cells = spots, slot = "data") # arbitrarily pick "MFGE8"
  data[,  "sample"] <- se@tools[["Staffli"]][[spots, "sample", drop = TRUE]]
  
  # gsub here removes the last two chars of the barcodes; in seurat they append the slice number to the end as _1...
  data[, "barcode"] <- gsub('.{2}$', '', spots) 
  data[,"x"] <- se@tools[["Staffli"]]@meta.data[["original_x"]] # coords of spots translated to pixel space somehow
  data[,"y"] <- se@tools[["Staffli"]]@meta.data[["original_y"]]
  data[,"align_x"] <- se@tools[["Staffli"]]@meta.data[["warped_x"]] # coords of spots aligned by STUtil 
  data[,"align_y"] <- se@tools[["Staffli"]]@meta.data[["warped_y"]]
  data[,"orig_x"] <- se@tools[["Staffli"]]@meta.data[["x"]] # coords of spots originally before pixel scaling
  data[,"orig_y"] <- se@tools[["Staffli"]]@meta.data[["y"]]
  
  data <- subset(data, select = -c(MFGE8))
}

pat_1 <- c("pat-1", "151507", "151508", "151509", "151510")
pat_2 <- c("pat-2", "151669", "151670", "151671", "151672")
pat_3 <- c("pat-3", "151673", "151674", "151675", "151676")

info1 <- get_info_table(pat_1)
info2 <- get_info_table(pat_2)
info3 <- get_info_table(pat_3)

se1 <- align_se(info1)
se2 <- align_se(info2)
se3 <- align_se(info3)

data1 <- get_coords(se1)
data2 <- get_coords(se2)
data3 <- get_coords(se3)

fwrite(data1, "../data/DLPFC/saved_results/pat-1-all-slices-STutil-coords-default.csv")
fwrite(data2, "../data/DLPFC/saved_results/pat-2-all-slices-STutil-coords-default.csv")
fwrite(data3, "../data/DLPFC/saved_results/pat-3-all-slices-STutil-coords-default.csv")

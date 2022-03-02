# setwd("~/Desktop/Research/paste-extension/paste-experiments/DLPFC-Analysis/") # TODO: add correct path for this file
library("Seurat")
# library(SeuratData) 
# library(SeuratDisk)
library(Matrix)
library(data.table)
library(glue)
library(tidyverse)
set.seed(0)

# We found interoperating between Seurat and Scanpy objects to be unreliable,
# so as an intermediate we generated csv files that contain the original spot 
# by count matrices for all the DLPFC samples. 
# 
# The DLPFC data has three samples with four slices each:
# Sample I: Slice “151507”, Slice “151508", Slice “151509”, Slice “151510"
# Sample II: Slice “151669”, Slice “151670", Slice “151671”, Slice “151672" 
# Sample III: Slice “151673”, Slice “151674", Slice “151675”, Slice “151676"
# 
# In most of the analysis below, we consider a slice by the last three letters
# as it is shorter and still unique, i.e. slice 1 from sample 1 is 507 instead of 151507.

# list_of_slices(which_sample): Combine all slices from one sample into a list of seurat objects. 
# Input: Either "1", "2", or "3" to specify which sample. 
# Output: Returns an R list of all slices (each slice is a Seurat object) in sample specified by input, 
# i.e. l <- list_of_slices("1"); l[[1]] = seurat object containing slice 1 of sample 1 and length(l) = 4.
list_of_slices <- function(which_sample) {
  if (which_sample == "1") {
    slices <- c("507", "508", "509", "510")
  } else if (which_sample == "2") {
    slices <- c("669", "670", "671", "672")
  } else if ( which_sample == "3") {
    slices <- c("673", "674", "675", "676")
  } 
  base_fn <- "data"
  l <- vector("list", 4) # 4 slices for each sample
  i <- 1
  while (i < 5){
    s <- slices[i]
    X <- fread(glue("{base_fn}/X-{s}.csv"), header=FALSE, sep=",")
    snames <- fread(glue("{base_fn}/spotnames-{s}.csv"), header=FALSE)
    colnames(snames) <- "spot_ids"
    gnames <- fread(glue("{base_fn}/genenames-{s}.csv"), header=FALSE)
    colnames(gnames) <- "gene_ids"
    colnames(X) <- snames$spot_ids
    row.names(X) <- gnames$gene_ids
    sc <- CreateSeuratObject(counts = X, assay=glue("slice{i}"))#project=glue("slice-{i}"))
    l[[i]] <- sc
    i <- i + 1
  }
  l
}

process_list <- function(s_list) {
  s_list_processed <- lapply(X = s_list, FUN = function(x) {
    x <- NormalizeData(x) # log normalize by default
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  s_list_processed 
}

sample_1_slices <- c("507", "508", "509", "510")
sample_2_slices <- c("669", "670", "671", "672")
sample_3_slices <- c("673", "674", "675", "676")
snames <- c("slice-1", "slice-2", "slice-3", "slice-4")

sclist_1 <- process_list(list_of_slices("1"))
names(sclist_1) <- snames
sclist_2 <- process_list(list_of_slices("2"))
names(sclist_2) <- snames
sclist_3 <- process_list(list_of_slices("3"))
names(sclist_3) <- snames

# Store all slices of every sample in a single list called samples
samples <- vector("list", 3)
for (i in 1:3) {
  if (i == 1) {
    samples[[i]] <- sclist_1
  } else if (i == 2) {
    samples[[i]] <- sclist_2
  } else if (i == 3) {
    samples[[i]] <- sclist_3
  } 
}

# For each sample i, calculate the seurat anchors between all slices (our.anchors)
i <- 1
while (i < 4) {
  samp <- samples[[i]]
  seurat.anchors <- FindIntegrationAnchors(object.list = samp) # Compute Seurat anchors
  fwrite(x = seurat.anchors@anchors, file = glue("data/sample-{i}-anchors.csv"))
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors)
  sparse_mat <- seurat.integrated@assays[["integrated"]]@data
  for (j in 1:4) { # iterate over slices 
    M <- as.data.frame(as.matrix(sparse_mat[,seurat.integrated$orig.ident == snames[j]]))
    if (i == 1) {
      s <- sample_1_slices[j]
    } else if (i == 2) {
      s <- sample_2_slices[j]
    } else if (i == 3) {
      s <- sample_3_slices[j]
    } 
    fwrite(x = M, row.names = TRUE, file = glue("data/X-integrated-{s}.csv"))
  }
  i <- i + 1
}
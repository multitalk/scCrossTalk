#---------------------------------Functions----------------------------------
# read scRNA-seq data
read_matrix <- function(patient_id, sc_condition, filedir){
  get_cell_id <- function(x){
    x_char <- as.character(x)
    x_nchar <- nchar(x_char)
    x_zero <- 5-x_nchar
    x_zero_new <- NULL
    for (i in 1:length(x_zero)) {
      x_zero1 <- x_zero[i]
      if (x_zero1 == 0) {
        x_zero_new <- c(x_zero_new, "")
      }
      if (x_zero1 == 1) {
        x_zero_new <- c(x_zero_new, "0")
      }
      if (x_zero1 == 2) {
        x_zero_new <- c(x_zero_new, "00")
      }
      if (x_zero1 == 3) {
        x_zero_new <- c(x_zero_new, "000")
      }
      if (x_zero1 == 4) {
        x_zero_new <- c(x_zero_new, "0000")
      }
    }
    x_char <- paste0(x_zero_new, x_char)
    return(x_char)
  }
  library(Matrix)
  a1 <- Seurat::Read10X(filedir)
  barcode_meta <- data.frame(patient_id = patient_id, dir_name = filedir,
                             barcode_raw = colnames(a1),
                             barcode = paste0(patient_id,sc_condition),stringsAsFactors = F)
  cellid <- 1:ncol(a1)
  cellid <- get_cell_id(cellid)
  barcode_meta$barcode <- paste0(barcode_meta$barcode, cellid)
  colnames(a1) <- barcode_meta$barcode
  return(list(rawdata = a1, barcode_meta = barcode_meta))
}

#---------------------------------Pre-process--------------------------------
# [1] load data
# download SC-Raw.zip
setwd("/path/to/SC-Raw/")
sc_meta <- data.frame(patient_id = paste0("P0",1:7),
                      cp_dir = c(paste0("Con_",1:5), "NO", "NO"), 
                      pp_dir = paste0("IR_",1:7), stringsAsFactors = FALSE)
rawdata_list <- list()
sc_meta_barcode <- data.frame()
k <- 0
for (i in 1:nrow(sc_meta)) {
  patient_id <- sc_meta$patient_id[i]
  cp_dir <- sc_meta$cp_dir[i]
  if (cp_dir != "NO") {
    k <- k + 1
    sc_condition <- "CP"
    rawdata <- read_matrix(patient_id, sc_condition, cp_dir)
    rawdata_list[[k]] <- rawdata[[1]]
    names(rawdata_list)[k] <- paste0(patient_id, "CP")
    sc_meta_barcode <- rbind(sc_meta_barcode, rawdata[[2]])
  }
  pp_dir <- sc_meta$pp_dir[i]
  if (pp_dir != "NO") {
    k <- k + 1
    sc_condition <- "PP"
    rawdata <- read_matrix(patient_id, sc_condition, pp_dir)
    rawdata_list[[k]] <- rawdata[[1]]
    names(rawdata_list)[k] <- paste0(patient_id, "PP")
    sc_meta_barcode <- rbind(sc_meta_barcode, rawdata[[2]])
  }
}

# [2] quality control
library(Seurat)
rawdata_meta <- lapply(X = rawdata_list, FUN = function(x){
  x <- CreateSeuratObject(x)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- x@meta.data
  return(x)
})
d1 <- rawdata_meta[[1]]
for (i in 2:length(rawdata_meta)) {
  d1 <- rbind(d1, rawdata_meta[[i]])
}
all(rownames(d1) == sc_meta_barcode$barcode)

sc_meta_barcode <- cbind(sc_meta_barcode, d1)

# save rds data
saveRDS(rawdata_list,file = "rawdata_list.rds")
saveRDS(sc_meta,file = "sc_meta_human.rds")
saveRDS(sc_meta_barcode,file = "sc_meta_cell.rds")

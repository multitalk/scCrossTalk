#' @title Create scCrossTalk object
#'
#' @description create scCrossTalk object using single-cell transcriptomics data
#' @param sc_data A data.frame or matrix or dgCMatrix containing raw counts of single-cell RNA-seq data. see \code{\link{demo_sc_data}}
#' @param sc_celltype A character containing the cell type of the single-cell RNA-seq data
#' @param species A character meaning species of the single-cell transcriptomics data.\code{'Human'} or \code{'Mouse'}.
#' @param if_normalize Normalize sc_data with Seurat LogNormalize. Set it \code{FLASE} when sc_data has been normalized.
#' @return scCrossTalk object
#' @import Matrix methods
#' @importFrom utils installed.packages
#' @importFrom Seurat CreateSeuratObject NormalizeData
#' @export

create_scCrossTalk <- function(sc_data, sc_celltype, species, if_normalize = TRUE) {
    if (is(sc_data, "data.frame")) {
        sc_data <- .get_dgCMatrix(as.matrix(sc_data))
    }
    if (is(sc_data, "matrix")) {
        sc_data <- .get_dgCMatrix(as.matrix(sc_data))
    }
    if (!is(sc_data, "dgCMatrix")) {
        stop("sc_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (!is.character(sc_celltype)) {
        stop("sc_celltype is not a character!")
    }
    if (ncol(sc_data) != length(sc_celltype)) {
        stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse")) {
        stop("Please provide a correct species, i.e., 'Human' or 'Mouse'!")
    }
    sc_celltype_new <- .rename_chr(sc_celltype)
    warning_info <- .show_warning(sc_celltype, sc_celltype_new)
    if (!is.null(warning_info)) {
        warning(warning_info)
    }
    colnames(sc_data) <- .rename_chr(colnames(sc_data))
    sc_meta <- data.frame(cell = colnames(sc_data), celltype = sc_celltype_new, stringsAsFactors = FALSE)
    sc_data <- sc_data[rowSums(sc_data) > 0, ]
    # if_normalize
    if (if_normalize) {
        sc_data <- Seurat::CreateSeuratObject(sc_data)
        sc_data <- Seurat::NormalizeData(sc_data,verbose = FALSE)
        ver <- packageVersion("Seurat")
        ver <- substr(ver,1,1)
        if (ver >= 5) {
            genenames <- rownames(sc_data)
            cellnames <- colnames(sc_data)
            sc_data <- sc_data[["RNA"]]@layers$data
            rownames(sc_data) <- genenames
            colnames(sc_data) <- cellnames
        } else {
            sc_data <- sc_data[["RNA"]]@data
        }
    }
    # generate scCrossTalk object
    object <- new("scCrossTalk", data = list(data = sc_data), meta = sc_meta, species = species)
    return(object)
}

#' @title Find highly expressed ligand-receptor pairs
#'
#' @description Find highly expressed ligands and receptors between pairwise cell types using Z score
#' @param object scCrossTalk object after \code{\link{create_scCrossTalk}}
#' @param lrpairs A data.frame of the system data containing ligand-receptor pairs of \code{'Human'} and \code{'Mouse'}. see \code{\link{demo_lrpairs}}
#' @param min_cell_num Min cell number for each cell type. Default is \code{10}
#' @param cell_min_pct Include the ligand and receptor gene detected in at least this many cells in each cluster. Default is \code{0.1}
#' @param p_value Include the significantly highly expressed ligands and receptors with this cutoff of p value from Z score. Default is \code{0.05}
#' @return scCrossTalk object containing the enriched ligand-receptor pairs and cell-cell communications
#' @import Matrix progress
#' @importFrom crayon cyan
#' @importFrom stats pnorm sd
#' @export

find_lrpairs <- function(object, lrpairs, min_cell_num = 10, cell_min_pct = 0.1, p_value = 0.05) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    if (!all(c("ligand", "receptor", "species") %in% colnames(lrpairs))) {
        stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
    }
    sc_data <- object@data$data
    sc_meta <- object@meta
    species <- object@species
    # filter lrpairs
    lrpairs <- lrpairs[lrpairs$species == species,]
    if (nrow(lrpairs) == 0) {
        stop("No lrpairs matched for input species!")
    }
    lrpairs <- lrpairs[lrpairs$ligand %in% rownames(sc_data) & lrpairs$receptor %in% rownames(sc_data),]
    if (nrow(lrpairs) == 0) {
        stop("No lrpairs matched for input scRNA-seq data!")
    }
    object@lrdb <- lrpairs
    # select cell types with min_cell_num
    celltype_sender_meta <- as.data.frame(table(sc_meta$celltype), stringsAsFactors = FALSE)
    colnames(celltype_sender_meta) <- c("celltype_sender", "cell_number")
    rownames(celltype_sender_meta) <- celltype_sender_meta$celltype_sender
    celltype_sender_meta <- celltype_sender_meta[celltype_sender_meta$cell_number >= min_cell_num, ]
    celltype_sender_meta <- celltype_sender_meta[celltype_sender_meta$cell_number > 1, ]
    if (nrow(celltype_sender_meta) == 0) {
        stop("No cell types found with cells more than min_cell_num!")
    }
    # generate pair-wise cell types
    cellname <- celltype_sender_meta$celltype_sender
    celltype_pair <- NULL
    for (i in 1:length(cellname)) {
        d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), celltype_receiver = cellname, stringsAsFactors = FALSE)
        celltype_pair <- rbind(celltype_pair, d1)
    }
    cat(cyan(paste0("Finding highly expressed LR pairs"), "\n"))
    pb <- progress_bar$new(format = "[:bar] Finished::percent time::elapsedfull", total = nrow(celltype_pair), clear = FALSE, width = 60, complete = "+", incomplete = "-")
    # generating result file
    res_pairs <- data.frame()
    for (i in 1:nrow(celltype_pair)) {
        res_pairs_tmp <- lrpairs
        res_pairs_tmp$celltype_sender <- celltype_pair$celltype_sender[i]
        res_pairs_tmp$celltype_receiver <- celltype_pair$celltype_receiver[i]
        # celltype_sender
        celltype_sender <- celltype_pair$celltype_sender[i]
        cell_sender <- sc_meta[sc_meta$celltype == celltype_sender, ]
        sc_data_sender <- sc_data[,cell_sender$cell]
        ligand_info <- .get_zscore_sender(sc_data_sender)
        ligand_info <- ligand_info[ligand_info$ligand_pct > cell_min_pct & ligand_info$ligand_p_value < p_value, ]
        if (nrow(ligand_info) == 0) {
            next
        }
        res_pairs_tmp <- merge(ligand_info, res_pairs_tmp)
        # celltype_receiver
        celltype_receiver <- celltype_pair$celltype_receiver[i]
        cell_receiver <- sc_meta[sc_meta$celltype == celltype_receiver, ]
        sc_data_receiver <- sc_data[,cell_receiver$cell]
        receptor_info <- .get_zscore_receiver(sc_data_receiver)
        receptor_info <- receptor_info[receptor_info$receptor_pct > cell_min_pct & receptor_info$receptor_p_value < p_value, ]
        if (nrow(receptor_info) == 0) {
            next
        }
        res_pairs_tmp <- merge(receptor_info, res_pairs_tmp)
        if (nrow(res_pairs_tmp) > 0) {
            res_pairs <- rbind(res_pairs, res_pairs_tmp)
        }
        pb$tick()
    }
    if (nrow(res_pairs) > 0) {
        res_pairs <- res_pairs[,c(12,13,6,1,7:10,2:5)]
        res_pairs$score <- res_pairs$ligand_exp_avg * res_pairs$receptor_exp_avg
        object@cci <- res_pairs
    } else {
        warning("No LR pairs found!")
    }
    return(object)
}

#' @title Get top n LR-pairs
#'
#' @description Get top n LR-pairs for each cell pairs
#' @param object scCrossTalk object after \code{\link{create_scCrossTalk}}
#' @param top_n Number of top LR-pairs. Default is {10}
#' @return A data.frame containing the top n LR-pairs for each cell pairs
#' @import Matrix
#' @export

get_top_lrpairs <- function(object, top_n = 10){
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    cellpair <- unique(cci[,c("celltype_sender", "celltype_receiver")])
    cci_final <- data.frame()
    for (i in 1:nrow(cellpair)) {
        cci_tmp <- cci[cci$celltype_sender == cellpair$celltype_sender[i] & cci$celltype_receiver == cellpair$celltype_receiver[i], ]
        cci_tmp <- cci_tmp[order(-cci_tmp$score),]
        if (nrow(cci_tmp) > top_n) {
            cci_tmp <- cci_tmp[1:top_n, ]
        }
        cci_final <- rbind(cci_final, cci_tmp)
    }
    return(cci_final)
}
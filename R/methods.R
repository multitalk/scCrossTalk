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
    if (!is(object, "miRTalk")) {
        stop("Invalid class for object: must be 'miRTalk'!")
    }
    if (!all(c("ligand", "receptor", "species") %in% colnames(lrpairs))) {
        stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
    }
    # check miRNA
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
        res_pairs <- res_pairs[,c(12,13,6,1,7:10,2:4)]
        res_pairs$score <- res_pairs$ligand_exp_avg * res_pairs$receptor_exp_avg
        celltype_pair$LR_number <- 0
        celltype_pair$LR_score <- 0
        for (i in 1:nrow(celltype_pair)) {
            res_pairs1 <- res_pairs[res_pairs$celltype_sender == celltype_pair$celltype_sender[i] & res_pairs$celltype_receiver == celltype_pair$celltype_receiver[i],]
            if (nrow(celltype_pair) > 0) {
                celltype_pair$LR_number[i] <- nrow(res_pairs1)
                celltype_pair$LR_score[i] <- sum(res_pairs1$score)
            }
        }
        celltype_pair <- celltype_pair[celltype_pair$LR_number > 0,]
        object@result@cci <- celltype_pair
        object@result@lrpairs <- res_pairs
    } else {
        warning("No LR pairs found!")
    }
    return(object)
}

#' @title Get pathways
#'
#' @description Get pathways for target genes
#' @param target_genes Character of one or more target genes
#' @param pathways A data.frame of the system data containing gene-gene interactions and pathways from KEGG and Reactome for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_pathways}}
#' @param species A character meaning species of the target genes.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @return Pathways for one or more target genes
#' @import Matrix
#' @export

find_pathways <- function(target_genes, pathways, species) {
    if (!all(c("src", "dest", "pathway", "species") %in% colnames(pathways))) {
        stop("Please provide a correct pathways data.frame! See demo_pathways()!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    pathways <- pathways[pathways$species == species, ]
    pathways <- pathways[pathways$src %in% target_genes | pathways$dest %in% target_genes, ]
    return(pathways)
}

#' @title Get GO terms
#'
#' @description Get GO terms for target genes
#' @param target_genes Character of one or more target genes
#' @param gene2go A data.frame of the system data containing GO terms for \code{'Human'}, \code{'Mouse'} or \code{'Rat'}. see \code{\link{demo_gene2go}}
#' @param species A character meaning species of the target genes.\code{'Human'}, \code{'Mouse'} or \code{'Rat'}
#' @return GO terms for one or more target genes
#' @import Matrix
#' @export

find_gene2go <- function(target_genes, gene2go, species) {
    if (!all(c("symbol", "GO_term", "species") %in% colnames(gene2go))) {
        stop("Please provide a correct gene2go data.frame! See demo_gene2go()!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human', 'Mouse' or 'Rat'!")
    }
    gene2go <- gene2go[gene2go$species == species, ]
    gene2go <- gene2go[gene2go$symbol %in% target_genes, ]
    return(gene2go)
}

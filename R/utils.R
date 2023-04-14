.rename_chr <- function(x){
    x <- strsplit(x, split = " ")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    x <- strsplit(x_new, split = "-")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    return(x_new)
}

.show_warning <- function(celltype, celltype_new){
    sc_meta <- data.frame(celltype = celltype, celltype_new = celltype_new, stringsAsFactors = FALSE)
    sc_meta <- unique(sc_meta)
    sc_meta <- sc_meta[sc_meta$celltype != sc_meta$celltype_new, ]
    warning_info <- NULL
    if (nrow(sc_meta) > 0) {
        warning_info <- "celltype of "
        if (nrow(sc_meta) == 1) {
            warning_info <- paste0(warning_info, sc_meta$celltype[1], " has been replaced by ", sc_meta$celltype_new[1])
        } else {
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype[i], ", ")
                }
            }
            warning_info <- paste0(warning_info, " have been replaced by ")
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype_new[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype_new[i], ", ")
                }
            }
        }
    }
    return(warning_info)
}

.get_dgCMatrix <- function(sc_data){
    all_pac <- installed.packages()
    all_pac <- as.data.frame(all_pac, stringsAsFactors = FALSE)
    all_pac <- all_pac[all_pac$Package == "Matrix", ]
    mat_version <- all_pac$Version
    mat_version <- strsplit(mat_version, split = '-')
    mat_version <- mat_version[[1]][1]
    mat_version <- strsplit(mat_version, split = '\\.')
    mat_version <- mat_version[[1]][2]
    mat_version <- substr(mat_version,1,1)
    mat_version <- as.numeric(mat_version)
    if (mat_version >= 5) {
        sc_data <- as(object = sc_data, Class = "CsparseMatrix")
    } else {
        sc_data <- as(object = sc_data, Class = "dgCMatrix")
    }
}

.percent_cell <- function(x) {
    return(length(x[x > 0])/length(x))
}

.get_zscore_sender <- function(sc_data_sender){
    ligand_pct <- apply(sc_data_sender, 1, .percent_cell)
    ligand_mean <- apply(sc_data_sender, 1, mean)
    ligand_info <- data.frame(ligand = names(ligand_pct), ligand_pct = as.numeric(ligand_pct), ligand_exp_avg = as.numeric(ligand_mean), stringsAsFactors = FALSE)
    ligand_info$ligand_z_score <- as.numeric(scale(ligand_info$ligand_exp_avg))
    ligand_info$ligand_p_value <- pnorm(q = ligand_info$ligand_exp_avg, mean = mean(ligand_info$ligand_exp_avg), sd = sd(ligand_info$ligand_exp_avg), lower.tail = FALSE)
    return(ligand_info)
}

.get_zscore_receiver <- function(sc_data_receiver){
    receptor_pct <- apply(sc_data_receiver, 1, .percent_cell)
    receptor_mean <- apply(sc_data_receiver, 1, mean)
    receptor_info <- data.frame(receptor = names(receptor_pct), receptor_pct = as.numeric(receptor_pct), receptor_exp_avg = as.numeric(receptor_mean), stringsAsFactors = FALSE)
    receptor_info$receptor_z_score <- as.numeric(scale(receptor_info$receptor_exp_avg))
    receptor_info$receptor_p_value <- pnorm(q = receptor_info$receptor_exp_avg, mean = mean(receptor_info$receptor_exp_avg), sd = sd(receptor_info$receptor_exp_avg), lower.tail = FALSE)
    return(receptor_info)
}

.get_coord <- function(cci, cellpair, y_len, lrpair_name, x_len) {
    y <- 1:y_len
    x <- 1:x_len
    cellpair <- rep(cellpair, x_len)
    y <- rep(y, x_len)
    lrpair_name <- rep(lrpair_name, each = y_len)
    x <- rep(x, each = y_len)
    cci_lrpair_temp <- data.frame(cellpair = cellpair, lrpair = lrpair_name, x = x, y = y, stringsAsFactors = FALSE)
    cci_lrpair_temp$score <- 0
    for (i in 1:nrow(cci)) {
        cci_lrpair_temp[cci_lrpair_temp$cellpair == cci$cellpair[i] & cci_lrpair_temp$lrpair == cci$lrpair[i], ]$score <- cci$score[i]
    }
    return(cci_lrpair_temp)
}

#' @title Show scCrossTalk object
#'
#' @param object scCrossTalk object after \code{\link{create_scCrossTalk}}
#' @return scCrossTalk object
#' @import Matrix
#' @importFrom methods show
#'
#' @export

setMethod(
    f = 'show',
    signature = 'scCrossTalk',
    definition = function(object) {
        cat("An object of class scCrossTalk", "\n")
        cci <- object@cci
        cat(paste0(nrow(cci), " ligand-receptor interactions found!"), "\n")
        return(invisible(x = NULL))
    }
)

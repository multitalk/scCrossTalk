#' @title Pre-processing step: revising gene symbols
#'
#' @description Revise genes according to NCBI Gene symbols updated in June 19, 2022 for count matrix and user-custom lrpairs data.frame
#' @param data A data.frame or matrix or dgCMatrix
#' @param data_type A character to define the type of \code{data}, select \code{'count'} for the data matrix, \code{'lrpairs'} for the lrpairs data.frame
#' @param species Species of the data.\code{'Human'} or \code{'Mouse'}
#' @param geneinfo A data.frame of the system data containing gene symbols of \code{'Human'} or \code{'Mouse'} updated on June 19, 2022 for revising gene symbols
#' @return A new data.frame, matrix, or dgCMatrix.
#' @import Matrix
#' @importFrom crayon cyan
#' @export

rev_gene <- function(data = NULL, data_type = NULL, species = NULL, geneinfo = NULL) {
    if (is.null(data)) {
        stop("Please provide the data for revsing gene symbols!")
    }
    if (is.null(data_type) | !is.character(data_type)) {
        stop("Please provide a correct data_type, i.e., 'count' or 'lrpairs'!")
    }
    if (is.null(geneinfo)) {
        stop("Please provide geneinfo for revsing gene symbols, or use the system data like 'geneinfo = geneinfo'")
    }
    if (length(data_type) > 1 | !data_type %in% c("count", "lrpairs")) {
      stop("Please provide a correct data_type, i.e., 'count' or 'lrpairs'!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse", "Rat")) {
        stop("Please provide a correct species, i.e., 'Human' or 'Mouse'!")
    }
    # define the species
    if (species == "Human") {
        geneinfo <- geneinfo[geneinfo$species == "Human", ]
    }
    if (species == "Mouse") {
        geneinfo <- geneinfo[geneinfo$species == "Mouse", ]
    }
    # revise matrix
    if (data_type == "count") {
        if (is(data, "data.frame")) {
            data <- methods::as(as.matrix(data), "dgCMatrix")
        }
        if (is(data, "matrix")) {
            data <- methods::as(data, "dgCMatrix")
        }
        if (!is(data, "dgCMatrix")) {
            stop("st_data must be a data.frame or matrix or dgCMatrix!")
        }
        Sys.sleep(1)
        # revise gene symbols
        genename <- rownames(data)
        genename1 <- genename[genename %in% geneinfo$symbol]
        if (length(genename1) == 0) {
            stop("Please ensure the rownames of data are gene symbols! See demo_sc_data()!")
        }
        genename2 <- genename[!genename %in% geneinfo$symbol]
        if (length(genename2) > 0) {
            genename3 <- genename2[genename2 %in% geneinfo$synonyms]
            if (length(genename3) > 0) {
                genename4 <- rep("NA", length(genename3))
                for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                        genename4[i] <- d1
                    }
                }
                genename3 <- c(genename1, genename3)
                genename4 <- c(genename1, genename4)
                genedata <- data.frame(raw_name = genename3, new_name = genename4, stringsAsFactors = F)
                genedata <- genedata[!genedata$new_name == "NA", ]
                genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
                genedata1 <- genedata1[genedata1$Freq == 1, ]
                genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
                data <- data[genedata$raw_name, ]
                rownames(data) <- genedata$new_name
            }
        } else {
            data <- data[rownames(data) %in% geneinfo$symbol, ]
        }
    }
    # revise lrpairs
    if (data_type == "lrpairs") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'lrpairs'!")
        }
        cat(crayon::cyan("Revising gene symbols for lrpairs data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("ligand", "receptor", "species") %in% colnames(data))) {
            # ligand
            genename <- unique(data$ligand)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the ligand of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$ligand == genename3[i], ]$ligand <- d1
                        }
                    }
                }
            }
            data <- data[data$ligand %in% geneinfo$symbol, ]
            # receptor
            genename <- unique(data$receptor)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the ligand of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                    for (i in 1:length(genename3)) {
                        d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                        if (length(d1) == 1) {
                            data[data$receptor == genename3[i], ]$receptor <- d1
                        }
                    }
                }
            }
            data <- data[data$receptor %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
        }
    }
    return(data)
}

#' @title Demo data of sc_data
#'
#' @description Demo data of sc_data.
#' @details \code{sc_data} can be a data.frame, matrix, or dgCMatrix object, each column representing a cell, each row representing a gene.
#' @return A dgCMatrix object.
#' @import Matrix
#' @export
#' @examples sc_data_demo <- demo_sc_data()

demo_sc_data <- function() {
    cellname <- paste0("cell", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NATP")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = TRUE)
    sc_data <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(sc_data) <- genename
    colnames(sc_data) <- cellname
    sc_data <- as(object = sc_data, Class = "dgCMatrix")
    return(sc_data)
}

#' @title Demo data of geneinfo
#'
#' @description Demo data of geneinfo
#' @details \code{geneinfo} must be a \code{data.frame} object with three columns, namely \code{'symbol'}, \code{'synonyms'}, \code{'species'}.
#' @export
#' @examples geneinfo_demo <- demo_geneinfo()

demo_geneinfo <- function() {
    gene1 <- c("A1BG", "A1BG", "A2MP1", "Aco1")
    gene2 <- c("A1B", "ABG", "A2MP", "Aco")
    species <- c("Human", "Human", "Human", "Mouse")
    geneinfo_demo <- data.frame(symbol = gene1, synonyms = gene2, species = species, stringsAsFactors = FALSE)
    return(geneinfo_demo)
}

#' @title Demo data of lrpairs
#'
#' @description Demo data of lrpairs
#' @details \code{lrpairs} must be a \code{data.frame} object with three columns, namely \code{'ligand'}, \code{'receptor'}, \code{'species'}.
#' @return A data.frame.
#' @export
#' @examples lrpairs_demo <- demo_lrpairs()

demo_lrpairs <- function() {
  ligand <- c("CX3CL1", "TGFB1", "CCL2", "Sst")
  receptor <- c("CX3CR1", "TGFBR2", "CCR2", "Sstr1")
  species <- c("Human", "Human", "Human", "Mouse")
  lrpairs_demo <- data.frame(ligand = ligand, receptor = receptor, species = species,
                             stringsAsFactors = F)
  return(lrpairs_demo)
}

#' @title Demo data of pathways
#'
#' @description Demo data of pathways
#' @details \code{pathways} must be a \code{data.frame} object with four columns, namely \code{'src'}, \code{'dest'}, \code{'pathway'}, \code{'species'}
#' @export
#' @examples pathways_demo <- demo_pathways()

demo_pathways <- function() {
    src <- c("CDKN1A", "CDKN1A", "CDK2", "Akt1", "Tcirg1")
    dest <- c("CDK2", "CDK4", "TP53", "Atf2", "Ppa1")
    pathway <- c("p53 signaling pathway", "p53 signaling pathway", "p53 signaling pathway",
                 "PI3K-Akt signaling pathway", "Oxidative phosphorylation")
    species <- c("Human", "Human", "Human", "Mouse", "Rat")
    pathways_demo <- data.frame(src = src, dest = dest, pathway = pathway, species = species, stringsAsFactors = FALSE)
    return(pathways_demo)
}

#' @title Definition of 'scCrossTalk' class
#'
#' @description An S4 class containing the data, meta, and results of inferred cell-cell communications mediated by ligand-receptor interactions.
#' @slot data A list containing the data.
#' @slot meta A data.frame containing the meta data.
#' @slot species A character containing the species.
#' @slot lrdb A data.frame containing expressed ligand and receptor genes.
#' @slot cci A data.frame containing enriched ligand-receptor interactions.
#' @import methods
#' @name scCrossTalk
#' @rdname scCrossTalk
#' @aliases scCrossTalk-class
#' @exportClass scCrossTalk

setClass("scCrossTalk", representation(data = "list", meta = "data.frame", species = "character",
    lrdb = "data.frame", cci = "data.frame"), prototype(data = list(), meta = data.frame(), species = character(),
    lrdb = data.frame(), cci = data.frame()))
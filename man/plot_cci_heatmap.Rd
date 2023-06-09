% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_cci_heatmap}
\alias{plot_cci_heatmap}
\title{Heatmap plot of cell-cell communications}
\usage{
plot_cci_heatmap(
  object,
  celltype = NULL,
  ligand = NULL,
  receptor = NULL,
  show_type = "number",
  text_size = 10,
  viridis_option = "D",
  ...
)
}
\arguments{
\item{object}{scCrossTalk object after \code{\link{find_lrpairs}}}

\item{celltype}{which cell types to plot by order. Default is to plot all cell types}

\item{ligand}{which ligand to use. Default is to plot all inferred ligands}

\item{receptor}{which receptor to use. Default is to plot all inferred receptors}

\item{show_type}{which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}}

\item{text_size}{Size of text labels. Default is \code{10}}

\item{viridis_option}{option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".}

\item{...}{parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color}
}
\value{
Heatmap plot of cell-cell communications mediated by ligand-receptor interactions
}
\description{
Heatmap plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
}

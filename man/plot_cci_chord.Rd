% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_cci_chord}
\alias{plot_cci_chord}
\title{Chord plot of cell-cell communications}
\usage{
plot_cci_chord(
  object,
  celltype = NULL,
  celltype_color = NULL,
  ligand = NULL,
  receptor = NULL,
  edge_color = NULL,
  edge_type = "big.arrow",
  show_type = "number",
  if_show_autocrine = FALSE,
  text_size = 1.5,
  y_scale = 0.1,
  ...
)
}
\arguments{
\item{object}{scCrossTalk object after \code{\link{find_lrpairs}}}

\item{celltype}{which cell types to plot by order. Default is to plot all cell types}

\item{celltype_color}{Colors for the cell types, whose length must be equal to \code{celltype}}

\item{ligand}{which ligand to use. Default is to plot all inferred ligands}

\item{receptor}{which receptor to use. Default is to plot all inferred receptors}

\item{edge_color}{Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}}

\item{edge_type}{Types for the edges from the sender cell type. Default is \code{"big.arrow"}. \code{"ellipse"} for ellipse, "triangle" for triangle, "curved" for curved. Details see \code{\link[circlize]{chordDiagram}}}

\item{show_type}{which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}}

\item{if_show_autocrine}{Whether to show autocrine. Default is \code{FALSE}}

\item{text_size}{Size of text labels. Default is \code{1.5}}

\item{y_scale}{y_scale to adjust the text. Default is \code{0.1}}

\item{...}{parameters pass to \code{\link[circlize]{chordDiagram}}, e.g., link.arr.width, link.arr.length, link.arr.col}
}
\value{
Chord plot of cell-cell communications mediated by ligand-receptor interactions
}
\description{
Chord plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
}

#' @title Chord plot of cell-cell communications
#'
#' @description Chord plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot by order. Default is to plot all cell types
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}
#' @param edge_type Types for the edges from the sender cell type. Default is \code{"big.arrow"}. \code{"ellipse"} for ellipse, "triangle" for triangle, "curved" for curved. Details see \code{\link[circlize]{chordDiagram}}
#' @param show_type which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param text_size Size of text labels. Default is \code{1.5}
#' @param y_scale y_scale to adjust the text. Default is \code{0.1}
#' @param ... parameters pass to \code{\link[circlize]{chordDiagram}}, e.g., link.arr.width, link.arr.length, link.arr.col
#' @import ggplot2 circlize Matrix
#' @importFrom scales hue_pal
#' @importFrom graphics text
#' @return Chord plot of cell-cell communications mediated by ligand-receptor interactions
#' @export

plot_cci_chord <- function(object, celltype = NULL, celltype_color = NULL, ligand = NULL, receptor = NULL, edge_color = NULL, edge_type = "big.arrow", show_type = "number", if_show_autocrine = FALSE, text_size = 1.5, y_scale = 0.1, ...) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype!")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
    } else {
        if (length(edge_color) != length(celltype)) {
            stop("The length of edge_color must be equal to celltype!")
        }
        link_color <- edge_color
        names(link_color) <- celltype
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "ligand", "receptor")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    chordDiagram(x = cci_pair, grid.col = clu_col[unique(cci_pair$from)], col = link_color[cci_pair$from], preAllocateTracks = 1, transparency = 0.25, directional = 1,
        direction.type = c("arrows", "diffHeight"), diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + y_scale, sector.name, facing = "inside", niceFacing = FALSE, adj = c(0.5, 0), cex = text_size)
    }, bg.border = NA)
    graphics::text(0, 1, paste0("Show_type: ", show_type_new), cex = text_size)
}

#' @title Circle plot of cell-cell communications
#'
#' @description Circle plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}
#' @param edge_type Types for the edges. \code{"fan"} by default, \code{"link"}, \code{"hive"}
#' @param show_type which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param edge_alpha Transparency of edge. Default is \code{0.5}
#' @param node_size Size of node. Default is \code{10}
#' @param text_size Size of text. Default is \code{5}
#' @import ggplot2 Matrix ggraph
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @return Circle plot of cell-cell communications mediated by ligand-receptor interactions
#' @export

plot_cci_circle <- function(object, celltype = NULL, ligand = NULL, receptor = NULL, celltype_color = NULL, edge_color = NULL, edge_type = "fan",
    show_type = "number", if_show_autocrine = FALSE, edge_alpha = 0.5, node_size = 10, text_size = 5) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
    } else {
        if (length(edge_color) != length(celltype)) {
            stop("The length of edge_color must be equal to celltype!")
        }
        link_color <- edge_color
        names(link_color) <- celltype
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "ligand", "receptor")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    cci_pair$sender <- as.character(cci_pair$from)
    celltype_node <- data.frame(name = celltype, celltype = celltype, id = 1:length(celltype), stringsAsFactors = FALSE)
    # angle
    angle <- 360 * (celltype_node$id - 0.5)/nrow(celltype_node)
    celltype_node$angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
    celltype_node$hjust <- ifelse(angle > 180, 1, 0)
    celltype_ordered <- celltype[order(celltype)]
    clu_col <- clu_col[celltype_ordered]
    celltype_ordered <- celltype_ordered[celltype_ordered %in% cci_pair$sender]
    link_color <- link_color[celltype_ordered]
    cci_pair_new <- data.frame()
    for (i in 1:length(celltype_ordered)) {
        cci_pair1 <- cci_pair[cci_pair$sender == celltype_ordered[i],]
        cci_pair_new <- rbind(cci_pair_new, cci_pair1)
    }
    mygraph <- graph_from_data_frame(cci_pair_new, vertices = celltype_node, directed = FALSE)
    p <- ggraph(mygraph, layout = "linear", circular = TRUE)
    if (edge_type == "fan") {
        p <- p + geom_edge_fan(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "link") {
        p <- p + geom_edge_link(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (edge_type == "hive") {
        p <- p + geom_edge_hive(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    if (if_show_autocrine) {
        p <- p + geom_edge_loop(aes(edge_colour = sender, edge_width = value), edge_alpha = edge_alpha)
    }
    p + geom_node_point(aes(color = celltype), size = node_size) + scale_color_manual(values = clu_col) +
        geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust), size = text_size) +
        scale_edge_color_manual(values = link_color) + expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
        coord_fixed() + theme_minimal() + theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"))
}

#' @title Sankey plot of cell-cell communications
#'
#' @description Sankey plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param celltype_color Colors for the cell types, whose length must be equal to \code{celltype}
#' @param edge_color Colors for the edges from the sender cell type, whose length must be equal to \code{celltype}, Or use \code{"NO"} to cancel it
#' @param show_type which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param edge_alpha Transparency of edge. Default is \code{0.5}
#' @param node_size Size of node. Default is \code{40}
#' @param text_size Size of text. Default is \code{15}
#' @param node_pad Size of node padding. Numeric essentially influences the width height. Default is \code{20}
#' @param ... parameters pass to \code{\link[networkD3]{sankeyNetwork}}
#' @import networkD3 Matrix
#' @importFrom scales hue_pal
#' @return Sankey plot of cell-cell communications mediated by ligand-receptor interactions
#' @export

plot_cci_sankey <- function(object, celltype = NULL, ligand = NULL, receptor = NULL, celltype_color = NULL, edge_color = NULL, show_type = "number", if_show_autocrine = FALSE,
    edge_alpha = 0.5, node_size = 40, text_size = 15, node_pad = 20, ...) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
    } else {
        if (length(celltype_color) != length(celltype)) {
            stop("The length of celltype_color must be equal to celltype")
        }
        clu_col <- celltype_color
    }
    names(clu_col) <- celltype
    if (is.null(edge_color[1])) {
        link_color <- clu_col
        names(link_color) <- celltype
        edge_color_new <- TRUE
    } else {
        if (edge_color[1] == "NO") {
            edge_color_new <- FALSE
            link_color <- clu_col
            names(link_color) <- celltype
        } else {
            if (length(edge_color) != length(celltype)) {
                stop("The length of edge_color must be equal to celltype!")
            }
            edge_color_new <- TRUE
            link_color <- edge_color
            names(link_color) <- celltype
        }
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "ligand", "receptor")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("from", "to", "value")
    cci_pair$signal <- cci_pair$from
    cci_pair$fromID <- 0
    cci_pair$toID <- 0
    celltype_sender <- data.frame(name = unique(cci_pair$from), stringsAsFactors = FALSE)
    celltype_sender$id <- (1:nrow(celltype_sender)) - 1
    celltype_receiver <- data.frame(name = unique(cci_pair$to), stringsAsFactors = FALSE)
    celltype_receiver$id <- (1:nrow(celltype_receiver)) - 1 + nrow(celltype_sender)
    celltype_node <- rbind(celltype_sender, celltype_receiver)
    # sender
    for (i in 1:nrow(celltype_sender)) {
        cci_pair[cci_pair$from == celltype_sender$name[i], ]$fromID <- celltype_sender$id[i]
    }
    # receiver
    for (i in 1:nrow(celltype_receiver)) {
        cci_pair[cci_pair$to == celltype_receiver$name[i], ]$toID <- celltype_receiver$id[i]
    }
    # link color
    celltype_sender <- celltype[celltype %in% cci_pair$from]
    cci_pair$link_type <- paste0("type", 1:nrow(cci_pair))
    cci_pair$link_color <- link_color[cci_pair$from]
    # node color
    celltype_node$node_color <- clu_col[celltype_node$name]
    ColourScal <- 'd3.scaleOrdinal() .domain(['
    cci_pair_link <- unique(cci_pair[,c("link_type","link_color")])
    cci_pair_link_type <- cci_pair_link$link_type
    cci_pair_link_color <- cci_pair_link$link_color
    celltype_node_new <- unique(celltype_node[,c("name","node_color")])
    celltype_node_type <- celltype_node_new$name
    celltype_node_color <- celltype_node_new$node_color
    cci_pair_link_type <- c(cci_pair_link_type, celltype_node_type)
    for (i in 1:length(cci_pair_link_type)) {
        if (i != length(cci_pair_link_type)) {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_type[i], '",')
        } else {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_type[i], '"])')
        }
    }
    ColourScal <- paste0(ColourScal, " .range([")
    cci_pair_link_color <- c(cci_pair_link_color, celltype_node_color)
    for (i in 1:length(cci_pair_link_color)) {
        if (i != length(cci_pair_link_color)) {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_color[i], '",')
        } else {
            ColourScal <- paste0(ColourScal, '"',cci_pair_link_color[i], '"])')
        }
    }
    if (edge_color_new) {
        sankeyNetwork(Links = cci_pair, Nodes = celltype_node, Source = "fromID", Target = "toID", Value = "value", NodeID = "name",colourScale = ColourScal,
            nodeWidth = node_size, fontSize = text_size, nodePadding = node_pad, fontFamily = "Arial", sinksRight = FALSE, LinkGroup = "link_type", NodeGroup = "name", ...)
    } else {
        sankeyNetwork(Links = cci_pair, Nodes = celltype_node, Source = "fromID", Target = "toID", Value = "value", NodeID = "name",colourScale = ColourScal,
            nodeWidth = node_size, fontSize = text_size, nodePadding = node_pad, fontFamily = "Arial", sinksRight = FALSE, NodeGroup = "name", ...)
    }
}

#' @title Heatmap plot of cell-cell communications
#'
#' @description Heatmap plot of cell-cell communications from senders to receivers with the sum of inferred number or score of ligand-receptor interactions
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot by order. Default is to plot all cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param show_type which type to show, \code{"number"} and \code{"score"} for sum of inferred LR number and score, respectively. Default is \code{"number"}
#' @param text_size Size of text labels. Default is \code{10}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color
#' @import heatmaply viridis reshape2 Matrix
#' @return Heatmap plot of cell-cell communications mediated by ligand-receptor interactions
#' @export

plot_cci_heatmap <- function(object, celltype = NULL, ligand = NULL, receptor = NULL, show_type = "number", text_size = 10, viridis_option = "D", ...) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_raw <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (is.null(celltype[1])) {
        celltype <- celltype_raw
    } else {
        if (!all(celltype %in% celltype_raw)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype & cci$celltype_receiver %in% celltype, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these cell types!")
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    cci_pair <- unique(cci[, c("celltype_sender", "celltype_receiver")])
    if (show_type == "number") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "ligand", "receptor")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- nrow(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ])
        }
        show_type_new <- "number"
    }
    if (show_type == "score") {
        cci <- unique(cci[, c("celltype_sender", "celltype_receiver", "score")])
        cci_pair$value <- 0
        for (i in 1:nrow(cci_pair)) {
            cci_pair$value[i] <- sum(cci[cci$celltype_sender == cci_pair$celltype_sender[i] & cci$celltype_receiver == cci_pair$celltype_receiver[i], ]$score)
        }
        show_type_new <- "score"
    }
    colnames(cci_pair) <- c("sender", "receiver", "value")
    cci_pair <- reshape2::dcast(data = cci_pair, formula = sender ~ receiver, value.var = "value", fill = 0)
    rownames(cci_pair) <- cci_pair$sender
    cci_pair <- cci_pair[, -1]
    heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = viridis_option)
    heatmaply::heatmaply(x = as.matrix(cci_pair), colors = heat_col, dendrogram = "none", xlab = "Receiver", ylab = "Sender", main = paste0("Show type: ", show_type_new), margins = c(60,100,40,20),
        grid_color = "white", grid_width = 0.00001, branches_lwd = 0.1, fontsize_row = text_size, fontsize_col = text_size, labCol = colnames(cci_pair), labRow = rownames(cci_pair),
        heatmap_layers = theme(axis.line=element_blank()),label_names = c("from","to",show_type_new), ...)
}

#' @title Heatmap plot of all lrpairs among cell pairs
#'
#' @description heatmap plot of inferred lrpairs. Rows for cell pairs, and columns for LR pairs by default.
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param text_size Size of text labels. Default is \code{10}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color
#' @import heatmaply viridis Matrix
#' @importFrom reshape2 dcast
#' @return Heatmap plot of all lrpairs among cell pairs
#' @export

plot_cci_lrpairs_heatmap <- function(object, celltype = NULL, ligand = NULL, receptor = NULL, if_show_autocrine = FALSE, if_horizontal = TRUE, text_size = 10, viridis_option = "D", ...) {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_sender <- unique(cci$celltype_sender)
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (!is.null(celltype[1])) {
        if (all(celltype %in% celltype_sender)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype, ]
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci$cellpair <- paste0(cci$celltype_sender, " | ", cci$celltype_receiver)
    cci$lrpair <- paste0(cci$ligand, " | ", cci$receptor)
    cci <- unique(cci[, c("cellpair", "lrpair", "score")])
    if (length(unique(cci$cellpair)) < 2 | length(unique(cci$lrpair)) < 2) {
        stop("Limited cellpair and lrpair to plot!")
    }
    cci_pair <- reshape2::dcast(data = cci, formula = cellpair ~ lrpair, value.var = "score", fill = 0)
    rownames(cci_pair) <- cci_pair$cellpair
    cci_pair <- cci_pair[, -1]
    if (!if_horizontal) {
        cci <- as.data.frame(t(cci))
    }
    heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = viridis_option)
    if (if_horizontal) {
        heatmaply::heatmaply(x = as.matrix(cci_pair), colors = heat_col, dendrogram = "none", xlab = "Ligands | Receptors", ylab = "Senders | Receivers", main = "score",
            fontsize_row = text_size, fontsize_col = text_size, labCol = colnames(cci_pair), labRow = rownames(cci_pair),
            heatmap_layers = theme(axis.line=element_blank()),label_names = c("cellpair","lrpair","score"), ...)
    } else {
        heatmaply::heatmaply(x = as.matrix(cci_pair), colors = heat_col, dendrogram = "none", xlab = "Senders | Receivers", ylab = "Ligands | Receptors", main = "score",
            fontsize_row = text_size, fontsize_col = text_size, labCol = colnames(cci_pair), labRow = rownames(cci_pair),
            heatmap_layers = theme(axis.line=element_blank()),label_names = c("cellpair","lrpair","score"), ...)
    }
}

#' @title Bubble plot of inferred ligand-receptor pairs
#'
#' @description Bubble plot of inferred ligand-receptor pairs from senders top receivers. Rows for cell pairs, and columns for LR pairs by default.
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype which cell types to plot. Default is to plot all cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param if_show_autocrine Whether to show autocrine. Default is \code{FALSE}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param size Size of the bubble. Default is \code{2}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @import ggplot2 viridis Matrix
#' @return Bubble plot of inferred ligand-receptor pairs
#' @export

plot_cci_lrpairs_bubble <- function(object, celltype = NULL, ligand = NULL, receptor = NULL, if_show_autocrine = FALSE, if_horizontal = TRUE, size = 2, viridis_option = "D") {
    # check input
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    celltype_sender <- unique(cci$celltype_sender)
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (!is.null(celltype[1])) {
        if (all(celltype %in% celltype_sender)) {
            stop("Please input the right celltype name!")
        }
        cci <- cci[cci$celltype_sender %in% celltype, ]
    }
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    if (!if_show_autocrine) {
        cci <- cci[cci$celltype_sender != cci$celltype_receiver, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these paracrine!")
    }
    cci$cellpair <- paste0(cci$celltype_sender, " | ", cci$celltype_receiver)
    cci$lrpair <- paste0(cci$ligand, " | ", cci$receptor)
    cci <- unique(cci[, c("cellpair", "lrpair", "score")])
    if (length(unique(cci$cellpair)) < 2 | length(unique(cci$lrpair)) < 2) {
        stop("Limited cellpair and lrpair to plot!")
    }
    cellpair <- unique(cci$cellpair)
    y_len <- length(cellpair)
    lrpair_name <- unique(cci$lrpair)
    x_len <- length(lrpair_name)
    cci <- .get_coord(cci, cellpair, y_len, lrpair_name, x_len)
    if (if_horizontal) {
        p <- ggplot2::ggplot(data = cci) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = score), size = size) + 
            viridis::scale_colour_viridis(option = viridis_option) + ggplot2::scale_y_continuous(name = "Senders | Receivers",
            breaks = 1:y_len, labels = cellpair, limits = c(1, y_len)) + ggplot2::scale_x_continuous(name = "Ligands | Receptors",
            breaks = 1:x_len, labels = lrpair_name, limits = c(1, x_len)) + ggplot2::labs(color = "score") + 
            ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else {
        p <- ggplot2::ggplot(data = cci) + ggplot2::geom_point(ggplot2::aes(x = y, y = x, color = score), size = size) + 
            viridis::scale_colour_viridis(option = viridis_option) + ggplot2::scale_x_continuous(name = "Senders | Receivers",
            breaks = 1:y_len, labels = cellpair, limits = c(1, y_len)) + ggplot2::scale_y_continuous(name = "Ligands | Receptors",
            breaks = 1:x_len, labels = lrpair_name, limits = c(1, x_len)) + ggplot2::labs(color = "score") + 
            ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
    return(p)
}

#' @title Chord plot of ligand-receptor pairs
#'
#' @description Chord plot of ligand-receptor pairs from senders to receivers.
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype_sender Name of celltype_sender. One or more cell types
#' @param celltype_receiver Name of celltype_receiver. One or more cell types
#' @param celltype_color Colors for the celltype_sender nodes and celltype_receiver nodes, or use \code{"NO"} to make it simple
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param edge_color Colors for the edges from the sender cell type
#' @param edge_type Types for the edges from the sender cell type. Default is \code{"circle"}. \code{"big.arrow"} for big arrow, "triangle" for triangle, "ellipse" for ellipse, "curved" for curved. Details see \code{\link[circlize]{chordDiagram}}
#' @param text_size Size of text labels. Default is \code{0.5}
#' @param y_scale y_scale to adjust the text. Default is \code{1}
#' @param ... parameters pass to \code{\link[circlize]{chordDiagram}}, e.g., link.arr.width, link.arr.length, link.arr.col
#' @import ggplot2 Matrix circlize
#' @importFrom scales hue_pal
#' @importFrom graphics text
#' @return Chord plot of ligand-receptor pairs from senders to receivers
#' @export

plot_lrpairs_chord <- function(object, celltype_sender, celltype_receiver, celltype_color = NULL, ligand = NULL, receptor = NULL, edge_color = NULL, edge_type = "circle", text_size = 0.5, y_scale = 1, ...) {
    # check object
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    celltype <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    # check celltype_sender
    if (!all(celltype_sender %in% cci$celltype_sender)) {
        stop("Please provide the correct name of celltype_sender!")
    }
    # check celltype_receiver
    if (!all(celltype_receiver %in% cci$celltype_receiver)) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    cci <- cci[cci$celltype_sender %in% celltype_sender & cci$celltype_receiver %in% celltype_receiver, ]
    if (nrow(cci) == 0) {
        stop("No lrpairs found from senders to receivers!")
    }
    # check lrpairs
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    celltype_sender_new <- unique(cci$celltype_sender)
    celltype_receiver_new <- unique(cci$celltype_receiver)
    if (length(celltype_sender) != length(celltype_sender_new)) {
        celltype_sender_no <- celltype_sender[!celltype_sender %in% celltype_sender_new]
        if (length(celltype_sender_no) == 1) {
            warning(paste0("celltype_sender of ", celltype_sender_no, " has been removed!"))
        }
        else {
            celltype_sender_no1 <- celltype_sender_no[1]
            for (i in 2:length(celltype_sender_no)) {
                celltype_sender_no1 <- paste0(celltype_sender_no1, ", ", celltype_sender_no[i])
            }
            warning(paste0("celltype_sender of ", celltype_sender_no1, " have been removed!"))
        }
        celltype_sender <- celltype_sender_new
    }
    if (length(celltype_receiver) != length(celltype_receiver_new)) {
        celltype_receiver_no <- celltype_receiver[!celltype_receiver %in% celltype_receiver_new]
        if (length(celltype_sender_no) == 1) {
            warning(paste0("celltype_receiver of ", celltype_receiver_no, " has been removed!"))
        }
        else {
            celltype_receiver_no1 <- celltype_receiver_no[1]
            for (i in 2:length(celltype_receiver_no)) {
                celltype_receiver_no1 <- paste0(celltype_receiver_no1, ", ", celltype_receiver_no[i])
            }
            warning(paste0("celltype_receiver of ", celltype_receiver_no1, " have been removed!"))
        }
        celltype_receiver <- celltype_receiver_new
    }
    # color
    if (is.null(celltype_color[1])) {
        clu_col <- scales::hue_pal()(length(celltype))
        names(clu_col) <- celltype
        clu_col <- clu_col[unique(c(celltype_sender, celltype_receiver))]
        celltype_color_new <- FALSE
    } else {
        if (celltype_color[1] == "NO") {
            clu_col <- scales::hue_pal()(length(celltype))
            names(clu_col) <- celltype
            clu_col <- clu_col[unique(c(celltype_sender, celltype_receiver))]
            celltype_color_new <- TRUE
        } else {
            if (length(celltype_color) != length(unique(c(celltype_sender, celltype_receiver)))) {
                stop("The length of celltype_color must be equal to sum of celltype_sender and celltype_receiver!")
            }
            celltype_color_new <- TRUE
            clu_col <- celltype_color
            names(clu_col) <- unique(c(celltype_sender, celltype_receiver))
        }
    }
    if (is.null(edge_color[1])) {
        link_color <- clu_col
        link_color <- link_color[celltype_sender]
    } else {
        if (length(edge_color) != length(celltype_sender)) {
            stop("The length of edge_color must be equal to celltype_sender!")
        }
        link_color <- edge_color
        names(link_color) <- celltype_sender
    }
    cci$ligand <- paste0(cci$ligand, " (", cci$celltype_sender, ")")
    cci$receptor <- paste0(cci$receptor, " (", cci$celltype_receiver, ")")
    celltype <- c(cci$celltype_sender, cci$celltype_receiver)
    grid_col <- clu_col[celltype]
    names(grid_col) <- c(cci$ligand, cci$receptor)
    cci <- unique(cci[ ,c("ligand", "receptor", "score")])
    if (celltype_color_new) {
        chordDiagram(x = cci, grid.col = grid_col, col = link_color, preAllocateTracks = 1, transparency = 0.25, directional = 1, direction.type = c("arrows", "diffHeight"),
            diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    } else {
        chordDiagram(x = cci, grid.col = grid_col, col = link_color, preAllocateTracks = 1, transparency = 0.25, directional = 1, direction.type = c("arrows", "diffHeight"),
            diffHeight = -0.04, annotationTrack = "grid", link.arr.type = edge_type, ...)
    }
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1]+y_scale, sector.name, facing = "clockwise", niceFacing = TRUE, cex = text_size)
    }, bg.border = NA)
}

#' @title Heatmap plot of ligand-receptor pairs
#'
#' @description Heatmap plot of ligand-receptor pairs from senders to receivers.
#' @param object scCrossTalk object after \code{\link{find_lrpairs}}
#' @param celltype_sender Name of celltype_sender. One or more cell types
#' @param celltype_receiver Name of celltype_receiver. One or more cell types
#' @param ligand which ligand to use. Default is to plot all inferred ligands
#' @param receptor which receptor to use. Default is to plot all inferred receptors
#' @param text_size Size of text labels. Default is \code{3}
#' @param if_horizontal Whether to plot with the horizontal direction. Default is \code{TRUE}
#' @param viridis_option option in \code{\link[viridis]{scale_color_viridis}}, can be "A", "B", "C", "D", "E", "F", "G", "H". Default is "D".
#' @param ... parameters pass to \code{\link[heatmaply]{heatmaply}}, e.g., grid_color
#' @import heatmaply Matrix
#' @importFrom reshape2 dcast
#' @return Heatmap plot of ligand-receptor pairs
#' @export

plot_lrpairs_heatmap <- function(object, celltype_sender, celltype_receiver, ligand = NULL, receptor = NULL, text_size = 5, if_horizontal = TRUE, viridis_option = "D", ...) {
    # check object
    if (!is(object, "scCrossTalk")) {
        stop("Invalid class for object: must be 'scCrossTalk'!")
    }
    cci <- object@cci
    celltype <- unique(c(cci$celltype_sender, cci$celltype_receiver))
    if (nrow(cci) == 0) {
        stop("No cci found in object!")
    }
    # check celltype_sender
    if (!celltype_sender %in% cci$celltype_sender) {
        stop("Please provide the correct name of celltype_sender!")
    }
    # check celltype_receiver
    if (!celltype_receiver %in% cci$celltype_receiver) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    cci <- cci[cci$celltype_sender == celltype_sender & cci$celltype_receiver == celltype_receiver, ]
    if (nrow(cci) == 0) {
        stop("No lrpairs found from senders to receivers!")
    }
    # check lrpairs
    ligand_name <- unique(cci$ligand)
    receptor_name <- unique(cci$receptor)
    if (!is.null(ligand[1])) {
        if (!all(ligand %in% ligand_name)) {
            stop("Please input the right ligand name!")
        }
        cci <- cci[cci$ligand %in% ligand, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these ligands!")
    }
    if (!is.null(receptor[1])) {
        if (!all(receptor %in% receptor_name)) {
            stop("Please input the right receptor name!")
        }
        cci <- cci[cci$receptor %in% receptor, ]
    }
    if (nrow(cci) == 0) {
        stop("No cci found for these receptors!")
    }
    cci <- unique(cci[ ,c("ligand", "receptor", "score")])
    cci <- reshape2::dcast(data = cci, formula = ligand ~ receptor, value.var = "score", fill = NA)
    rownames(cci) <- cci$ligand
    cci <- cci[, -1]
    if (!if_horizontal) {
        cci <- as.data.frame(t(cci))
    }
    if (if_horizontal) {
        heatmaply::heatmaply(x = as.matrix(cci), dendrogram = "none", xlab = paste0(celltype_receiver, ": receptor"), ylab = paste0(celltype_sender, ": ligand"), main = "LR score",
            grid_color = "black",fontsize_row = text_size, na.value = "white", fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci),
            heatmap_layers = theme(axis.line=element_blank()), label_names = c("ligand","receptor","score"), ...)
    } else{
        heatmaply::heatmaply(x = as.matrix(cci), dendrogram = "none", xlab = paste0(celltype_sender, ": ligand"), ylab = paste0(celltype_receiver, ": receptor"), main = "LR score",
            grid_color = "black",fontsize_row = text_size, na.value = "white", fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci),
            heatmap_layers = theme(axis.line=element_blank()), label_names = c("receptor","ligand","score"), ...)
      
    }
}

library(scCrossTalk)

# load example data
load(paste0(system.file(package = 'scCrossTalk'), "/extdata/example.rda"))

sc_data <- rev_gene(data = sc_data,data_type = "count",species = "Human",geneinfo = geneinfo)

obj <- create_scCrossTalk(sc_data = sc_data,sc_celltype = sc_celltype,species = "Human",if_normalize = FALSE)
obj <- find_lrpairs(object = obj, lrpairs = lrpairs, min_cell_num = 10, cell_min_pct = 0.1, p_value = 0.05)

plot_cci_chord(obj)
plot_cci_circle(obj)
plot_cci_heatmap(obj)
plot_cci_sankey(obj)
plot_cci_lrpairs_bubble(obj)
plot_cci_lrpairs_heatmap(obj)
plot_lrpairs_chord(obj,celltype_sender = "alpha_cell",celltype_receiver = "beta_cell",link.arr.width = 0.5, link.arr.length = 0.2)
plot_lrpairs_heatmap(obj,celltype_sender = "alpha_cell",celltype_receiver = "beta_cell",if_horizontal = F)

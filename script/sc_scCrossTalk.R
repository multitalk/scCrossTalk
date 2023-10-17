library(scCrossTalk)
obj_sc_nead <- create_scCrossTalk(sc_data = obj_data_nead, sc_celltype = obj_meta_nead$subtype, species = "Human", if_normalize = F)
obj_sc_nead <- find_lrpairs(obj_sc_nead, lrpairs = lrpairs)
obj_sc_ead <- create_scCrossTalk(sc_data = obj_data_ead, sc_celltype = obj_meta_ead$subtype, species = "Human", if_normalize = F)
obj_sc_ead <- find_lrpairs(obj_sc_ead, lrpairs = lrpairs)

# [a]
cci_nead <- obj_sc_nead@cci
cci_nead$cellpair <- paste0(cci_nead$celltype_sender,cci_nead$celltype_receiver)
cci_nead$lrpair <- paste0(cci_nead$ligand,cci_nead$receptor)
max(table(cci_nead$cellpair)) # 54
min(table(cci_nead$cellpair)) # 4

cci_ead <- obj_sc_ead@cci
cci_ead$cellpair <- paste0(cci_ead$celltype_sender,cci_ead$celltype_receiver)
cci_ead$lrpair <- paste0(cci_ead$ligand,cci_ead$receptor)
max(table(cci_ead$cellpair)) # 78
min(table(cci_ead$cellpair)) # 6

celltypecolor <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
p <- plot_cci_circle(object = obj_sc_ead,celltype_color = celltypecolor)+
  scale_edge_width_continuous(range = c(0,4))
p
ggsave(filename = "D:/circle_ead.pdf",plot = p,device = "pdf",width = 8,height = 8)


p <- plot_cci_chord(object = obj_sc_ead,celltype_color = celltypecolor)
p
ggsave(filename = "D:/chord_ead.pdf",plot = p,device = "pdf",width = 8,height = 8)


p <- plot_cci_circle(object = obj_sc_nead,celltype_color = celltypecolor)+
  scale_edge_width_continuous(range = c(0,4*(54/78)))
p
ggsave(filename = "D:/circle_nead.pdf",plot = p,device = "pdf",width = 8,height = 8)

p <- plot_cci_chord(object = obj_sc_nead,celltype_color = celltypecolor)
p
ggsave(filename = "D:/chord_nead.pdf",plot = p,device = "pdf",width = 8,height = 8)


# [b] compare

d1 <- as.data.frame(table(cci_nead$cellpair),stringsAsFactors = F)
d2 <- as.data.frame(table(cci_ead$cellpair),stringsAsFactors = F)
all(d1$Var1 == d2$Var1)
d3 <- cbind(d1,d2)
write.csv(x = d3,file = "D:/num_lr_all.csv")

# [d] score

d1 <- cbind(cci_pair_nead,cci_pair)
write.csv(x = d1,file = "D:/score_lr_all.csv")

# [f] overlap
d1 <- obj_sc_nead@cci
#d1 <- d1[d1$celltype_sender == "Hep",]
d1 <- d1[d1$celltype_sender == "Endo",]
d1 <- d1[d1$celltype_sender != d1$celltype_receiver,]
d1$celltype_sender <- paste0(d1$celltype_sender,"1")
d1$celltype_receiver <- paste0(d1$celltype_receiver,"1")
d1$cellpair <- paste0(d1$celltype_sender,d1$celltype_receiver)
d1$lrpair <- paste0(d1$ligand,d1$receptor)

d2 <- obj_sc_ead@cci
#d2 <- d2[d2$celltype_sender == "Hep",]
d2 <- d2[d2$celltype_sender == "Endo",]
d2 <- d2[d2$celltype_sender != d2$celltype_receiver,]
d2$celltype_sender <- paste0(d2$celltype_sender,"2")
d2$celltype_receiver <- paste0(d2$celltype_receiver,"2")
d2$cellpair <- paste0(d2$celltype_sender,d2$celltype_receiver)
d2$lrpair <- paste0(d2$ligand,d2$receptor)

res_plot <- as.data.frame(matrix(0, nrow = 8, ncol = 8))
colnames(res_plot) <- c("Endo1Mait1","Endo1NK_GZMB_GZMK1","Endo1S100A12+1","Endo1Hep1","Endo2Mait2","Endo2NK_GZMB_GZMK2","Endo2S100A12+2","Endo2Hep2")
rownames(res_plot) <- c("Endo1Mait1","Endo1NK_GZMB_GZMK1","Endo1S100A12+1","Endo1Hep1","Endo2Mait2","Endo2NK_GZMB_GZMK2","Endo2S100A12+2","Endo2Hep2")
d3 <- rbind(d1,d2)
for (i in 1:nrow(res_plot)) {
  a1_row <- rownames(res_plot)[i]
  a1_row <- d3[d3$cellpair == a1_row,]
  for (j in 1:ncol(res_plot)) {
    a1_col <- colnames(res_plot)[j]
    a1_col <- d3[d3$cellpair == a1_col,]
    res_plot[i,j] <- length(intersect(a1_row$lrpair, a1_col$lrpair))
  }
}
p <- pheatmap::pheatmap(as.matrix(res_plot),cluster_cols = F,cluster_rows = F,
                        color = colorRampPalette(c("white", "#ff7f0e"))(50),display_numbers = T,border_color = "black")
ggsave(filename = "D:/heat_overlap_endo.pdf",plot = p,device = "pdf",width = 8,height = 8)

# [f]
d1 <- obj_sc_nead@cci
d1$cellpair <- paste0(d1$celltype_sender,d1$celltype_receiver)
d1$lrpair <- paste0(d1$ligand,d1$receptor)
d1 <- d1[d1$celltype_receiver == "Hep",]
d1 <- d1[d1$celltype_sender != d1$celltype_receiver,]
table(d1$cellpair)
# 4,4,5,13

d1 <- obj_sc_ead@cci
d1$cellpair <- paste0(d1$celltype_sender,d1$celltype_receiver)
d1$lrpair <- paste0(d1$ligand,d1$receptor)
d1 <- d1[d1$celltype_receiver == "Hep",]
d1 <- d1[d1$celltype_sender != d1$celltype_receiver,]
table(d1$cellpair)
# 8,6,7,19

# [h] sender:hep and ec
obj_filter_nead <- obj_sc_nead
obj_filter_nead@cci <- get_top_lrpairs(object = obj_filter_nead, top_n = 10)
d1 <- obj_filter_nead@cci
d1 <- d1[d1$celltype_sender %in% c("Hep", "Endo"),]
d1 <- d1[d1$celltype_sender != d1$celltype_receiver,]
cci <- d1
cci$cellpair <- paste0(cci$celltype_sender, " | ", cci$celltype_receiver)
cci$lrpair <- paste0(cci$ligand, " | ", cci$receptor)
cci <- unique(cci[, c("cellpair", "lrpair", "score")])
if (length(unique(cci$cellpair)) < 2 | length(unique(cci$lrpair)) < 2) {
  stop("Limited cellpair and lrpair to plot!")
}
cci1 <- cci
# ead
obj_filter_ead <- obj_sc_ead
obj_filter_ead@cci <- get_top_lrpairs(object = obj_filter_ead, top_n = 10)
d1 <- obj_filter_ead@cci
d1 <- d1[d1$celltype_sender %in% c("Hep", "Endo"),]
d1 <- d1[d1$celltype_sender != d1$celltype_receiver,]
d1$celltype_sender <- paste0(d1$celltype_sender, "1")
d1$celltype_receiver <- paste0(d1$celltype_receiver, "1")
cci <- d1
cci$cellpair <- paste0(cci$celltype_sender, " | ", cci$celltype_receiver)
cci$lrpair <- paste0(cci$ligand, " | ", cci$receptor)
cci <- unique(cci[, c("cellpair", "lrpair", "score")])
if (length(unique(cci$cellpair)) < 2 | length(unique(cci$lrpair)) < 2) {
  stop("Limited cellpair and lrpair to plot!")
}

cci <- rbind(cci1,cci)
cellpair <- unique(cci$cellpair)
y_len <- length(cellpair)
lrpair_name <- unique(cci$lrpair)
x_len <- length(lrpair_name)
cci <- .get_coord(cci, cellpair, y_len, lrpair_name, x_len)

p <- ggplot2::ggplot(data = cci) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = score, size = score)) + 
  viridis::scale_colour_viridis(option = viridis_option) + ggplot2::scale_y_continuous(name = "Senders | Receivers",
  breaks = 1:y_len, labels = cellpair, limits = c(1, y_len)) + ggplot2::scale_x_continuous(name = "Ligands | Receptors",
  breaks = 1:x_len, labels = lrpair_name, limits = c(1, x_len)) + ggplot2::labs(color = "score") + 
  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggsave(filename = "D:/heat_all.pdf",plot = p,device = "pdf",width = 14,height = 4)

# [i] SAA1 in Hep, FPR1 in S100A12
d1 <- obj_data_nead["SAA1",rownames(obj_meta_nead[obj_meta_nead$subtype == "Hep",])]
d2 <- obj_data_ead["SAA1",rownames(obj_meta_ead[obj_meta_ead$subtype == "Hep",])]
d1 <- data.frame(con = "NO", value = as.numeric(d1),stringsAsFactors = F)
d2 <- data.frame(con = "YES", value = as.numeric(d2),stringsAsFactors = F)
t.test(d1$value, d2$value,alternative = "less")
# p-value < 2.2e-16
d3 <- rbind(d1,d2)
write.csv(d3, file = "D:/SAA1_hep.csv")


d1 <- obj_data_nead["FPR1",rownames(obj_meta_nead[obj_meta_nead$subtype == "S100A12+",])]
d2 <- obj_data_ead["FPR1",rownames(obj_meta_ead[obj_meta_ead$subtype == "S100A12+",])]
d1 <- data.frame(con = "NO", value = as.numeric(d1),stringsAsFactors = F)
d2 <- data.frame(con = "YES", value = as.numeric(d2),stringsAsFactors = F)
t.test(d1$value, d2$value,alternative = "less")
# p-value = 5.483e-15
d3 <- rbind(d1,d2)
write.csv(d3, file = "D:/FPR1_neu.csv")


d1 <- obj_data_nead["ITGB2",rownames(obj_meta_nead[obj_meta_nead$subtype == "S100A12+",])]
d2 <- obj_data_ead["ITGB2",rownames(obj_meta_ead[obj_meta_ead$subtype == "S100A12+",])]
d1 <- data.frame(con = "NO", value = as.numeric(d1),stringsAsFactors = F)
d2 <- data.frame(con = "YES", value = as.numeric(d2),stringsAsFactors = F)
t.test(d1$value, d2$value,alternative = "less")
# p-value = 2.62e-09
d1 <- obj_data_nead["ITGB2",rownames(obj_meta_nead[obj_meta_nead$subtype == "NK_GZMB_GZMK",])]
d2 <- obj_data_ead["ITGB2",rownames(obj_meta_ead[obj_meta_ead$subtype == "NK_GZMB_GZMK",])]
d1 <- data.frame(con = "NO", value = as.numeric(d1),stringsAsFactors = F)
d2 <- data.frame(con = "YES", value = as.numeric(d2),stringsAsFactors = F)
t.test(d1$value, d2$value,alternative = "less")
# p-value = 2.749e-09

# [j] SAA1 in Hep, FPR1 in S100A12
d1 <- obj_data_nead["SAA1",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["SAA1",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/SAA1_mean.csv")

# FPR1
d1 <- obj_data_nead["FPR1",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["FPR1",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/FPR1_mean.csv")

# ITGB2
d1 <- obj_data_nead["ITGB2",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["ITGB2",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/ITGB2_mean.csv")

# FGA
d1 <- obj_data_nead["FGA",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["FGA",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/FGA_mean.csv")

# FGB
d1 <- obj_data_nead["FGB",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["FGB",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/FGB_mean.csv")

# FGG
d1 <- obj_data_nead["FGG",]
names(d1) <- obj_meta_nead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_1 <- d2

d1 <- obj_data_ead["FGG",]
names(d1) <- obj_meta_ead$subtype
d2 <- NULL
d2 <- c(d2, mean(d1[names(d1) == "Endo"]))
d2 <- c(d2, mean(d1[names(d1) == "Hep"]))
d2 <- c(d2, mean(d1[names(d1) == "Mait"]))
d2 <- c(d2, mean(d1[names(d1) == "NK_GZMB_GZMK"]))
d2 <- c(d2, mean(d1[names(d1) == "S100A12+"]))
exp_2 <- d2

exp12 <- data.frame(celltype = c("Endo","Hep","Mait","NK_GZMB_GZMK","S100A12+"),nead = exp_1, ead = exp_2,stringsAsFactors = F)
write.csv(exp12, file = "D:/FGG_mean.csv")


# [k] DEGs for hep
obj_markers_hep$cluster <- as.character(obj_markers_hep$cluster)
d1 <- obj_markers_hep[obj_markers_hep$cluster == "YES",]
d1 <- d1[order(-d1$avg_log2FC),]
deg_up <- d1$gene[1:100]
d1 <- obj_markers_hep[obj_markers_hep$cluster == "NO",]
d1 <- d1[order(-d1$avg_log2FC),]
deg_down <- d1$gene[1:100]
obj_hep <- ScaleData(obj_hep)
p <- DoHeatmap(obj_hep,features = c(c(deg_up, deg_down)))
p
ggsave(filename = "D:/heat_degs.pdf",plot = p,device = "pdf",width = 8,height = 8)

degs <- data.frame(nead = deg_down, ead = deg_up,stringsAsFactors = F)
write.csv(degs,file = "D:/degs.csv")

# GSEA analysis # C2, C5, H
library(msigdbr)
gmt <- msigdbr(species = "Homo sapiens", category = "H")
gmt <- unique(gmt[,c("gs_name","gene_symbol")])
colnames(gmt) <- c("term","gene")
gmt$term <- tolower(gmt$term)
sc_data <- obj_hep[["RNA"]]@data
colnames(sc_data) <- obj_hep$ead
geneList <- apply(sc_data,1,function(x) {
  x1 <- mean(x[names(x) == "NO"])
  x2 <- mean(x[names(x) == "YES"])
  return(x2-x1)
})
geneList <- geneList[order(-geneList)]

res<- clusterProfiler::GSEA(geneList = geneList,TERM2GENE = gmt)
res1 <- clusterProfiler::gseKEGG(geneList = geneList,organism = "hsa")
enrichplot::gseaplot2(res,1,color="red",pvalue_table = T)
sc_data <- obj_hep[["RNA"]]@counts
sc_data <- log2(sc_data+1)
sc_data <- as.matrix(sc_data)
colnames(sc_data) <- obj_hep$ead
geneList <- apply(sc_data,1,function(x) {
  x1 <- mean(x[names(x) == "NO"])
  x2 <- mean(x[names(x) == "YES"])
  return(x2-x1)
})
geneList <- geneList[order(-geneList)]

p <- enrichplot::gseaplot2(res,1,color="red",pvalue_table = T)
ggsave(filename = "D:/gsea1.png",plot = p,device = "png",width = 12,height = 8,dpi = 300)

p <- enrichplot::gseaplot2(res,3,color="red",pvalue_table = T)
ggsave(filename = "D:/gsea2.png",plot = p,device = "png",width = 12,height = 8,dpi = 300)

# [m]
celltypecolor <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "S100A12+", celltype_receiver = "Mait",
                   celltype_color = c("#9467bd","#2ca02c"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "S100A12+", celltype_receiver = "Mait",
                   celltype_color = c("#9467bd","#2ca02c"),edge_type = "big.arrow")

plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "S100A12+", celltype_receiver = "NK_GZMB_GZMK",
                   celltype_color = c("#9467bd","#d62728"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "S100A12+", celltype_receiver = "NK_GZMB_GZMK",
                   celltype_color = c("#9467bd","#d62728"),edge_type = "big.arrow")


plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "Mait", celltype_receiver = "NK_GZMB_GZMK",
                   celltype_color = c("#2ca02c","#d62728"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "Mait", celltype_receiver = "NK_GZMB_GZMK",
                   celltype_color = c("#2ca02c","#d62728"),edge_type = "big.arrow")

plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "Mait", celltype_receiver = "S100A12+",
                   celltype_color = c("#2ca02c","#9467bd"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "Mait", celltype_receiver = "S100A12+",
                   celltype_color = c("#2ca02c","#9467bd"),edge_type = "big.arrow")


plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "NK_GZMB_GZMK", celltype_receiver = "Mait",
                   celltype_color = c("#d62728","#2ca02c"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "NK_GZMB_GZMK", celltype_receiver = "Mait",
                   celltype_color = c("#d62728","#2ca02c"),edge_type = "big.arrow")

plot_lrpairs_chord(object = obj_sc_nead,celltype_sender = "NK_GZMB_GZMK", celltype_receiver = "S100A12+",
                   celltype_color = c("#d62728","#9467bd"),edge_type = "big.arrow")
plot_lrpairs_chord(object = obj_sc_ead,celltype_sender = "NK_GZMB_GZMK", celltype_receiver = "S100A12+",
                   celltype_color = c("#d62728","#9467bd"),edge_type = "big.arrow")


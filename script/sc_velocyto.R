library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(SeuratDisk)

revise_gene <- function(obj_loom, geneinfo){
  rawdata <- obj_loom[["spliced"]]@counts
  genename<- rownames(rawdata)
  genename1<- genename[genename %in% geneinfo$symbol]
  genename2<- genename[!genename %in% geneinfo$symbol]
  genename3<- genename2[genename2 %in% geneinfo$synonyms]
  genename4<- rep('NA',length(genename3))
  for (i in 1:length(genename3)) {
    d1<- geneinfo[geneinfo$synonyms == genename3[i],]$symbol
    if(length(d1) == 1){
      genename4[i]<- d1
    }
  }
  genename3<- c(genename1,genename3)
  genename4<- c(genename1,genename4)
  genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
  genedata<- genedata[!genedata$new_name == 'NA',]
  genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
  genedata1<- genedata1[genedata1$Freq == 1,]
  genedata<- genedata[genedata$new_name %in% genedata1$Var1,]
  rawdata<- rawdata[genedata$raw_name,]
  rownames(rawdata)<- genedata$new_name
  obj_loom[["spliced"]]@counts <- rawdata
  obj_loom[["spliced"]]@data <- rawdata
  # unspliced
  rawdata <- obj_loom[["unspliced"]]@counts
  genename<- rownames(rawdata)
  genename1<- genename[genename %in% geneinfo$symbol]
  genename2<- genename[!genename %in% geneinfo$symbol]
  genename3<- genename2[genename2 %in% geneinfo$synonyms]
  genename4<- rep('NA',length(genename3))
  for (i in 1:length(genename3)) {
    d1<- geneinfo[geneinfo$synonyms == genename3[i],]$symbol
    if(length(d1) == 1){
      genename4[i]<- d1
    }
  }
  genename3<- c(genename1,genename3)
  genename4<- c(genename1,genename4)
  genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
  genedata<- genedata[!genedata$new_name == 'NA',]
  genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
  genedata1<- genedata1[genedata1$Freq == 1,]
  genedata<- genedata[genedata$new_name %in% genedata1$Var1,]
  rawdata<- rawdata[genedata$raw_name,]
  rownames(rawdata)<- genedata$new_name
  obj_loom[["unspliced"]]@counts <- rawdata
  obj_loom[["unspliced"]]@data <- rawdata
  # ambiguous
  rawdata <- obj_loom[["ambiguous"]]@counts
  genename<- rownames(rawdata)
  genename1<- genename[genename %in% geneinfo$symbol]
  genename2<- genename[!genename %in% geneinfo$symbol]
  genename3<- genename2[genename2 %in% geneinfo$synonyms]
  genename4<- rep('NA',length(genename3))
  for (i in 1:length(genename3)) {
    d1<- geneinfo[geneinfo$synonyms == genename3[i],]$symbol
    if(length(d1) == 1){
      genename4[i]<- d1
    }
  }
  genename3<- c(genename1,genename3)
  genename4<- c(genename1,genename4)
  genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
  genedata<- genedata[!genedata$new_name == 'NA',]
  genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
  genedata1<- genedata1[genedata1$Freq == 1,]
  genedata<- genedata[genedata$new_name %in% genedata1$Var1,]
  rawdata<- rawdata[genedata$raw_name,]
  rownames(rawdata)<- genedata$new_name
  obj_loom[["ambiguous"]]@counts <- rawdata
  obj_loom[["ambiguous"]]@data <- rawdata
  return(obj_loom)
}
integrated_data <- function(obj_loom){
  obj_meta <- obj_loom@meta.data
  obj_loom$sample <- substr(rownames(obj_meta),1,5)
  sc_list <- SplitObject(object = obj_loom,split.by = "sample")
  sc_list <- lapply(X = sc_list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE, assay = "spliced")
    x <- FindVariableFeatures(x, verbose = FALSE, assay = "spliced")
    return(x)
  })
  features <- SelectIntegrationFeatures(object.list = sc_list,verbose = F)
  rawdata_seurat <- FindIntegrationAnchors(object.list = sc_list, anchor.features = features, verbose = F)
  gc()
  rawdata_seurat <- IntegrateData(anchorset = rawdata_seurat,verbose = F)
  gc()
  DefaultAssay(rawdata_seurat) <- "integrated"
  rawdata_seurat <- ScaleData(rawdata_seurat, verbose = FALSE)
  rawdata_seurat <- RunPCA(rawdata_seurat, verbose = FALSE)
  rawdata_seurat <- RunUMAP(rawdata_seurat, reduction = "pca", dims = 1:30,verbose = FALSE)
  return(rawdata_seurat)
}

# read tnk data
# obj_tnk <- subset(obj_seurat, subset = celltype == "T_NK")
obj_seurat <- readRDS("/path/to/obj_tnk.rds")
obj_tnk_umap <- obj_tnk@reductions[["umap"]]@cell.embeddings
obj_tnk_pca <- obj_tnk@reductions[["pca"]]@cell.embeddings
gene_info <- readRDS("/path/to/data/geneinfo.rds")
gene_info <- gene_info[gene_info$species == "Human",]
obj_tnk_loom <- readRDS("obj_tnk_loom.rds")
obj_tnk_loom <- revise_gene(obj_tnk_loom, geneinfo)
obj_tnk_loom[["RNA"]] <- obj_tnk_loom[["spliced"]]
obj_tnk_loom <- integrated_data(obj_tnk_loom)
obj_tnk_loom_umap <- obj_tnk_loom@reductions[["umap"]]@cell.embeddings
obj_tnk_loom_pca <- obj_tnk_loom@reductions[["pca"]]@cell.embeddings
obj_tnk_umap <- obj_tnk_umap[rownames(obj_tnk_loom_umap),]
obj_tnk_pca <- obj_tnk_umap[rownames(obj_tnk_loom_pca),]
obj_tnk_loom@reductions[["umap"]]@cell.embeddings <- obj_tnk_umap
obj_tnk_loom@reductions[["pca"]]@cell.embeddings <- obj_tnk_pca
Idents(obj_tnk_loom) <- obj_tnk_loom$subtype
DefaultAssay(obj_tnk_loom) <- "RNA"
SaveH5Seurat(obj_tnk_loom, filename = "velocyto/obj_tnk_loom.h5Seurat")
Convert("velocyto/obj_tnk_loom.h5Seurat", dest = "h5ad")

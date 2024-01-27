# scCrossTalk
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen?logo=github)](https://github.com/multitalk/scCrossTalk/actions) [![CellTalkDB v1.0](https://img.shields.io/badge/CellTalkDB-v1.0-yellow)](http://tcm.zju.edu.cn/celltalkdb/)

### A cell-cell communication inference method for single-cell transcriptomic data

<img src='https://github.com/multitalk/scCrossTalk/blob/main/img/github.png'>

Cell-cell communications in multi-cellular organisms generally involve secreted ligand-receptor (LR) interactions, which is vital for various biological phenomena. Recent advancements in single-cell RNA sequencing (scRNA-seq) have effectively resolved cellular phenotypic heterogeneity and the cell-type composition of complex tissues, facilitating the systematic investigation of [cell-cell communications at single-cell resolution](https://pubmed.ncbi.nlm.nih.gov/32435978/). Here, we introduce scCrossTalk, a  cell-cell communication inference approach for single-cell transcriptomic data based on [CellTalkDB](https://pubmed.ncbi.nlm.nih.gov/33147626/) by enriching the highly expressed ligand-receptor pairs with the Z-test statistical method. scCrossTalk is an effective method that can help scientists analyze and visualize cell-cell communications for single-cell transcriptomic data.

# Install

```
# install dependent packages `devtools` and install
> install.packages(pkgs = 'devtools')
> devtools::install_github('multitalk/scCrossTalk')

# or download the repository as ZIP
> devtools::install_local("/path/to/scCrossTalk-main.zip")
```

# Usage
scCrossTalk method consists of two parts, wherein the first is to enrich the LR pairs mediating cell-cell communications and the second is to visualize the cell-cell communications and the underlying LR interactions. Classification and description of scCrossTalk functions are shown in the __[document](https://github.com/multitalk/scCrossTalk/blob/main/vignettes/scCrossTalk.pdf)__ and __[tutorial](https://raw.githack.com/multitalk/scCrossTalk/main/vignettes/tutorial.html)__

- ### Enrich the LR pairs mediating cell-cell communications
```
# sc_data: A matrix containing counts of scRNA-seq data
# sc_celltype: A character containing the cell types for scRNA-seq data

> obj <- create_scCrossTalk(sc_data, sc_celltype, species)
> 
> obj
An object of class scCrossTalk
0 ligand-receptor interactions found!

# object: scCrossTalk object containg scRNA-seq data
# lrpairs: A data.frame of the system data containing ligand-receptor pairs

> obj <- find_lrpairs(object, lrpairs = lrpairs)
Finding highly expressed LR pairs
[++++++++++++++++++++++++++++++] Finished:100% time:00:00:55
>
> obj
An object of class scCrossTalk 
679 ligand-receptor interactions found!
```

- ### Visualize the cell-cell communications and the underlying LR interactions
```
# [a] cell-cell communication analysis

# object: scCrossTalk object after performing find_lrpairs()

> plot_cci_chord(object = obj)
> plot_cci_circle(object = obj)
> plot_cci_heatmap(object = obj)
> plot_cci_sankey(object = obj)

# [b] ligand-receptor interactions analysis

# object: scCrossTalk object after performing find_lrpairs()
# celltype_sender: sender cell type
# celltype_receiver: receiver cell type

> plot_cci_lrpairs_bubble(object = obj)
> plot_cci_lrpairs_heatmap(object = obj)
> plot_lrpairs_chord(object = obj, celltype_sender, celltype_receiver)
> plot_lrpairs_heatmap(object = obj, celltype_sender, celltype_receiver)
```

# About
scCrossTalk was designed in our manuscript to dissect the different mechanism underlying cell-cell communications within the transplanted livers of EAD and non-EAD patients. 

<img src='https://github.com/multitalk/scCrossTalk/blob/main/img/ccc.png'>

__(a, c)__ Cell-cell communications mediated by LR interactions in non-EAD and EAD patients, where the circle and heat plots displayed the number and score of LR pairs, respectively. __(b, d)__ Comparison of the number and score of LR pairs between non-EAD and EAD patients with one-sided Welch¡¯s test. __(e)__ Cell-cell communications represented by the number of LR pairs with ligands sent from hepatocytes and received by endothelial cells, MAIT cells, GZMB+ GZMK+ NK cells, and S100A12+ neutrophils. (f) Shared LR pairs underlying cell-cell communications among hepatocytes, endothelial cells, MAIT cells, GZMB+ GZMK+ NK cells, and S100A12+ neutrophils in non-EAD and EAD patients. __(g)__ Cell-cell communications represented by the number of LR pairs with ligands sent from endothelial cells, MAIT cells, GZMB+ GZMK+ NK cells, and S100A12+ neutrophils and received by hepatocytes. __(h)__ Top 10 LR pairs mediating cell-cell communications from hepatocytes and endothelial cells to endothelial cells, hepatocytes, MAIT cells, GZMB+ GZMK+ NK cells, and S100A12+ neutrophils in non-EAD and EAD patients. __(i)__ Comparison of the expression of SAA1 in hepatocytes and FPR1 in S100A12+ neutrophils between non-EAD and EAD patients. __(j)__ Expression of SAA1 and FPR1 among endothelial cells, hepatocytes, MAIT cells, GZMB+ GZMK+ NK cells, and S100A12+ neutrophils in non-EAD and EAD patients. __(k)__ Top 100 DEGs of hepatocytes and the enriched pathways and biological processes between non-EAD and EAD patients. __(l)__ Enriched hallmarks of TNF¦Á signaling via NFKB and bile acid metabolism by comparing EAD patients to non-EAD patients with GSEA.

# Cite
Please cite us as: X. Shao, Z. Wang, K. Wang, X. Lu, P. Zhang, R. Guo, J. Liao, P. Yang, X. Xu, X. Fan,
A Single-Cell Landscape of Human Liver Transplantation Reveals a Pathogenic Immune Niche Associated with
Early Allograft Dysfunction, Engineering (2024), doi: https://doi.org/10.1016/j.eng.2023.12.004

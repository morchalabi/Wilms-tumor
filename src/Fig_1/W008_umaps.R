# This script plots Fig. 2d and f (UMAP of clusters and histology of P10). P10 refers to patient PN14 in the Supplementary Table 1.

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Reading in snRNA-seq data of P10 (PN14) ####

load(file = '../data/untreated/W008/data_QC_COMPRT.RData')
DefaultAssay(s_obj) = 'RNA'

# Annotating cancer clusters ####

s_obj$histology = factor('healthy', levels = c('cancer','healthy','mixed'))
s_obj$histology[!s_obj$compartment %in% c(5,6)] = as.factor('cancer')
s_obj$histology[s_obj$compartment %in% 4] = as.factor('mixed')

# Plotting UMAPs ####

p_ = DimPlot(s_obj, group.by = c('compartment','histology'), label = T, pt.size = .5, label.size = 4) & NoLegend()
ggsave(plot = p_, filename = 'W008_umaps.pdf', device = 'pdf', width = 14, height = 7)

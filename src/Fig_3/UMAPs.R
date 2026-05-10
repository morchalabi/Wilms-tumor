# This script plots Fig. 3a (UMAP projection of cancer compartment).

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)

# Reading in cancer compartment ####

load(file = '../data/integration/healthy/immune/lymphoid/integrated_compartment_lymphoid_processed.RData'); lymph_ = s_objs
load(file = '../data/integration/healthy/immune/myeloid/integrated_compartment_myeloid_processed.RData'); mye_ = s_objs
load(file = '../data/integration/healthy/vasculature/integrated_compartment_vasculature_processed.RData'); vasc_ = s_objs

# Plotting UMAP ####

DimPlot(lymph_, group.by = 'cell_type', label = T, label.size = 5, pt.size = 1.5, repel = T) & NoLegend()
ggsave(filename = 'lymphpid.pdf', device = 'pdf', width = 7, height = 7)


DimPlot(mye_, group.by = 'cell_type', label = T, label.size = 5, pt.size = 1.5, repel = T) & NoLegend()
ggsave(filename = 'myeloid.pdf', device = 'pdf', width = 7, height = 7)


DimPlot(vasc_, group.by = 'cell_type', label = T, label.size = 5, pt.size = 1.5, repel = T) & NoLegend()
ggsave(filename = 'vasculature.pdf', device = 'pdf', width = 7, height = 7)


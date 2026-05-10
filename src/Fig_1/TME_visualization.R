# This script plots UMAP presentation of TME for several options

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratWrappers)
library(grid)
library(gridExtra)
library(gridBase)
library(viridis)


load('../data/integration/healthy_cancer/TME.RData')

pdf(file = 'TME.pdf', width = 8, height = 8)

  # patients
  p_ = DimPlot(tme_, group.by = 'orig.ident', raster = F, label = F, shuffle = T, repel = T)+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
        labs(title = NULL, color = 'Patients', x = 'UMAP', y = 'UMAP')
  plot(p_)

  # treatment condition
  p_ = DimPlot(tme_, group.by = 'condition', raster = F, label = F, repel = T, pt.size = .7, shuffle = T, cols = c('red','blue'))+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Treatment', x = 'UMAP', y = 'UMAP')
  plot(p_)

  # compartment
  p_ = DimPlot(tme_, group.by = 'compartment_name', raster = F, label = F, repel = T, pt.size = .7)+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Compartment', x = 'UMAP', y = 'UMAP')
  plot(p_)

  # subcompartment
  p_ = DimPlot(tme_, group.by = 'subcompartment', raster = F, label = F, repel = F, pt.size = .6)+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Subcompartment', x = 'UMAP', y = 'UMAP')
  plot(p_)

  # cell type
  p_ = DimPlot(tme_, group.by = 'cell_type_2', raster = F, label = T, repel = T, pt.size = .7)+
  theme(legend.position = 'none',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 7, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = 'Cell type', x = 'UMAP', y = 'UMAP')
  plot(p_)

  # cell cycle phase
  p_ = DimPlot(tme_, group.by = 'Phase', raster = F, label = F, repel = T)+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Cell cycle', x = 'UMAP', y = 'UMAP')
  plot(p_)

graphics.off()

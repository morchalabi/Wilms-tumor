# This scripts plot cancer compartment UMAP

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

load('../data/integration/cancer/integrated_compartment_cancer_processed.RData')

DimPlot(s_objs, group.by = 'cell_type', pt.size = .4, repel = T, shuffle = T, label = T, raster = F, label.size = 4,
        cols = c("blastema 1" = 'red4',
                 "S blastema 1" = '#460080',
                 "G2M blastema 1"= '#C55300',
                 
                 "blastema 2" = 'red3',
                 "G2M blastema 2" = '#E65300',
                 
                 "blastema 3" = 'red',
                 
                 "smooth myocyte-like precursors" = 'lightgreen',
                 "G2M smooth myocyte-like precursors" = 'green2',
                 
                 "nascent smooth myocyte-like" = 'green3',
                 
                 "smooth myocyte-like" = 'green4',
                 "S smooth myocyte-like" = '#004680',
                 "G2M smooth myocyte-like" = '#809800',
                 
                 "PAX3+ myogenic precursors" = '#464600',
                 "myoblast-like" = '#E68000',
                 "striated myocyte-like" = '#DFB25F',
                 
                 "tubules" = 'blue',
                 "S tubules" = 'lightblue3',
                 "G2M tubules" = 'lightblue4',
                 
                 "podocyte-like" = 'skyblue'))+
theme(legend.position = 'none',
      axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())+
labs(title = 'Cancer')+
coord_equal()
ggsave(filename = 'cancer.pdf', device = 'pdf', width = 7.5, height = 7.5)

# This script generates a ternary diagram (x,y,z in 2D) of cancer cells where x, y and z are fetal signature scores for blastema, epithelium
# and stroma, and superimposes density contour on them. In other words, each cell is assigned three coordinates
# (x = blastema, y = epithelium, z = stroma). These coordinates dictate their positions.
# The script uses ggtern package that uses a 2D kernel density estimation using kde2d and display the results with contours.

library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(ggplot2)
library(ggtern)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(parallel)
set.seed(42)

# loading cancer data ####

load('../data/integration/cancer/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'

# Subsetting ####
mitotic_cells = grepl(x = s_objs$cell_type_2, pattern = '(mitotic)', perl = T)
s_objs = subset(s_objs, cells = colnames(s_objs)[!mitotic_cells])

e_ = colnames(s_objs)[s_objs$subcompartment %in% 'epithelium']      # smallest cluster's size
b_ = sample(colnames(s_objs)[s_objs$subcompartment %in% 'blastema'], size = length(e_))
s_ = sample(colnames(s_objs)[s_objs$subcompartment %in% 'stroma'], size = length(e_))
s_objs = subset(s_objs, cells = c(e_,b_,s_)); gc()
s_objs$cell_type_2 = droplevels(s_objs$cell_type_2)
s_objs$cell_type = droplevels(s_objs$cell_type)

DimPlot(s_objs, group.by = c('subcompartment','cell_type_2'), label = T) & NoLegend()
#####################################################################################

# Loading fetal kidney signatures ####

markers_1 = read.delim(file = '../out/benchmark/fetal/kidneycellatlas(PMID-31604275)/celltypes_markers.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
markers_1 = markers_1[trimws(markers_1$gene) %in% rownames(s_objs),]

markers_2 = read.delim(file = '../out/benchmark/fetal/descartes(PMID-33184181)/celltypes_markers.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
markers_2 = markers_2[trimws(markers_2$gene) %in% rownames(s_objs),]

markers_ = rbind(markers_1, markers_2)

# blastema:
blastema = c('Cap mesenchyme',
             'Distal renal vesicle','Proximal renal vesicle',
             'Proximal S shaped body','Medial S shaped body','Distal S shaped body')
blast_markers = unique(unlist(lapply(X = unique(blastema), FUN = function(c_){ trimws(markers_$gene[markers_$seurat_clusters %in% c_ & 0 < markers_$avg_log2FC])[1:20] })))

# epithelium (tubules and podocyte):
epithelium = c('Proximal tubule','Loop of Henle',
               'Podocyte')
epi_markers = unique(unlist(lapply(X = unique(epithelium), FUN = function(c_){ trimws(markers_$gene[markers_$seurat_clusters %in% c_ & 0 < markers_$avg_log2FC])[1:20] })))

# stroma (fibroblast, myofibroblast and mesangium):
stroma = c('Kidney-Mesangial cells','Mesangial',
           'Kidney-Stromal cells', 'Stroma progenitor',
           'Fibroblast','Myofibroblast',
           'Smooth muscle')
stroma_markers = unique(unlist(lapply(X = unique(stroma), FUN = function(c_){ trimws(markers_$gene[markers_$seurat_clusters %in% c_ & 0 < markers_$avg_log2FC])[1:20] })))

# Signature scores ####

signatures_ = list(blastema = blast_markers, epithelium = epi_markers, stroma = stroma_markers)     # all signatures

# average expression as signature score
s_objs$blastema = colMeans(s_objs[['RNA']]@data[signatures_$blastema,])
s_objs$epithelium = colMeans(s_objs[['RNA']]@data[signatures_$epithelium,])
s_objs$stroma = colMeans(s_objs[['RNA']]@data[signatures_$stroma,])

# Plotting ####

# scaling
vals_ = s_objs@meta.data[,names(signatures_)]
colnames(vals_) = c("tumor blastema","tumor epithelium","tumor stroma")
vals_ = (vals_ - min(vals_))/ diff(range(vals_))
vals_ = vals_[rowSums(vals_) > 0,]
vals_ = sqrt(vals_)
vals_ = cbind(s_objs@meta.data[rownames(vals_),!colnames(s_objs@meta.data) %in% names(signatures_)], vals_)

# plotting

# for treatment condition

pdf(file = 'cancer_density_plasticity.pdf', width = 14.3, height = 12)
breaks_ = seq(0,1,length.out = 5)
p_ = ggtern(data = vals_, mapping = aes(x = `tumor blastema`, y = `tumor epithelium`, z = `tumor stroma`))+
     labs(title = "Tumor blastemal plasticity", subtitle = 'Fetal signature expressoion', x = 'fetal blastema', z = 'fetal stroma', y = 'fetal epithelium')+
     geom_point(size = .2, alpha = .2)+
     stat_density_tern(geom = 'polygon', aes(fill  = after_stat(level)), alpha = .5,
                       n = 400,        # smoothing
                       h = .5,         # smoothing
                       bins = 15,      # number of contours
                       color = 'grey50', linewidth = .1)
  
lvl_rng = range(ggplot_build(p_)$data[[2]]$level)
p_ = p_ + scale_fill_gradient(low = "black", high = "red", name = "Density level", breaks = seq(lvl_rng[1],lvl_rng[2], length.out = 3),      # breaks is based on density levels computed by stat_density_tern. To see levels: p_ = ggplot_build(p_); range(p_$data[[2]]$level)
          labels = c("low", "mid", "high"))+
          scale_R_continuous(breaks = breaks_, labels = breaks_)+
          scale_L_continuous(breaks = breaks_, labels = breaks_)+
          scale_T_continuous(breaks = breaks_, labels = breaks_)+
          guides(fill = guide_colorbar(order = 1), alpha = guide_none()) +
          theme_bvbw(base_size = 15, base_family = 'Helvetica')+
          theme_ticklength(major = unit(5,'mm'))+
          theme(plot.title = element_text(hjust = .5, face = 'bold', family = 'Helvetica'), plot.subtitle = element_text(hjust = .5, face = 'bold', family = 'Helvetica'),
                axis.title = element_text(size = 15, face = 'bold', family = 'Helvetica'), axis.text  = element_text(size = 15, face = 'bold', family = 'Helvetica'),
                legend.text = element_text(size = 15, face = 'bold', family = 'Helvetica'), legend.title = element_text(size = 15, face = 'bold', family = 'Helvetica'),
                # Axis lines
                tern.axis.line.L = element_line(color = "red"),
                tern.axis.line.T = element_line(color = "blue"),
                tern.axis.line.R = element_line(color = "green3"),
                # Axis tick text
                tern.axis.text.L = element_text(color = "red"),
                tern.axis.text.T = element_text(color = "blue"),
                tern.axis.text.R = element_text(color = "green3"),
                # Axis titles
                tern.axis.title.L = element_text(color = "red"),
                tern.axis.title.T = element_text(color = "blue"),
                tern.axis.title.R = element_text(color = "green3"),
                # Axis arrows
                tern.axis.arrow.L = element_line(color = "red"),
                tern.axis.arrow.T = element_line(color = "blue"),
                tern.axis.arrow.R = element_line(color = "green3"),
                # Axis arrow text
                tern.axis.arrow.text.L = element_text(color = "red"),
                tern.axis.arrow.text.T = element_text(color = "blue"),
                tern.axis.arrow.text.R = element_text(color = "green3"))

plot(p_)
graphics.off()

# for tumor subcompartments

pdf(file = 'cancer_density_plasticity_subcompartments.pdf', width = 14.3, height = 12)

vals_$subcompartment = factor(x = vals_$subcompartment, levels = c('blastema','stroma','epithelium'))
p_ = ggtern(data = vals_, mapping = aes(x = `tumor blastema`, y = `tumor epithelium`, z = `tumor stroma`, color = subcompartment))+
            labs(title = "Tumor blastemal plasticity", subtitle = 'Fetal signature expressoion', x = 'fetal blastema', z = 'fetal stroma', y = 'fetal epithelium', color = 'Tumor')+
            geom_point(size = 1, alpha = 1)+
            guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
            theme_bvbw(base_size = 15, base_family = 'Helvetica')+
            theme_ticklength(major = unit(5,'mm'))+
            theme(plot.title = element_text(hjust = .5, face = 'bold', family = 'Helvetica'), plot.subtitle = element_text(hjust = .5, face = 'bold', family = 'Helvetica'),
                  axis.title = element_text(size = 15, face = 'bold', family = 'Helvetica'), axis.text  = element_text(size = 15, face = 'bold', family = 'Helvetica'),
                  legend.text = element_text(size = 15, face = 'bold', family = 'Helvetica'), legend.title = element_text(size = 15, face = 'bold', family = 'Helvetica'),
                  legend.key = element_blank())
plot(p_)
graphics.off()


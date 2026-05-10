# This script was updated March 30, 2025

# This scripts takes infercnv.references.txt and infercnv.observations.txt from inferCNV output and generates CNV heatmap per patient
#
# in inferCNV expression values are cut in 15 equal intervals with those in the 8th interval as neutral/wildtype CNV.
# The centroid in the 8th interval is 1, meaning all gene expressions within (1 - epi, 1 + epi) are neutral. Those expressions
# from 1st through 7th intervals are shown by shades of blue representing deletion. Those expression values from 9t through
# 15th intervals are amplification shown by shades of red

# Att. In case of using integrated data, an extra id was added to cells during integration by Seurat to avoid duplicates.
# Thus, inferCNV's cells have that id, perform this:
# s_obj = RenameCells(object = s_obj, new.names = annot$cells[annot$orig.ident %in% smpl_])

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(gridExtra)
library(pheatmap)
library(grid)
library(viridis)
set.seed(42)


# integrated cancer data
load('../data/integration/cancer/integrated_compartment_cancer_processed.RData')

# ####

all_genes = NULL
obs_mat = annot_ = list()
for(smpl_ in unique(s_objs$orig.ident))
{
  message(smpl_)
  
  # reading current sample's cells
  cancer_cells = colnames(s_objs)[s_objs$orig.ident %in% smpl_]     # cancer cells of current sample
  cond_ = s_objs$condition[cancer_cells[1]]                         # treatment condition
  cell_types = s_objs$cell_type[cancer_cells]                       # cancer cell subtypes
  names(cell_types) = cancer_cells = sub(pattern = '_[0-9]+', replacement = '', x = cancer_cells)                                           # as inferCNV was run per patient, cell barcodes do not have _[1-9]+ ids added by Seurat integration
  
  # observation mat
  obs_mat_tmp = t(as.matrix(read.delim(file = paste0('../out/',cond_,'/',smpl_,'/inferCNV/inferCNV_TCells/infercnv.observations.txt'),      # all non-reference cells, including nonmalignant cells
                                       header = T, row.names = 1, sep = ' ', quote = '"', as.is = T, check.names = F)))
  
  # not all cancer cells were analyzed by inferCNV?
  cancer_cells = cancer_cells[cancer_cells %in% rownames(obs_mat_tmp)]                                                                      # which cancer cells are present at observation matrix?
  obs_mat_tmp = obs_mat_tmp[cancer_cells,]                                                                                                  # retaining cancer cells in observation mat
  rownames(obs_mat_tmp) = paste0(smpl_,'_',rownames(obs_mat_tmp))                                                                           # adding sample name to cell ids to make them unique across cohort, required for heatmap
  
  # updating all cells and all genes variable
  all_genes = if(is.null(all_genes)){ colnames(obs_mat_tmp) }else{ intersect(all_genes, colnames(obs_mat_tmp))}                             # updating list of genes shared by all samples
  annot_[[smpl_]] = data.frame(cells = rownames(obs_mat_tmp),                                                              # annotation table required for pheatmap
                               sample = smpl_,
                               codition = cond_,
                               cell_type = cell_types[cancer_cells],
                               row.names = NULL)
  obs_mat[[smpl_]] = obs_mat_tmp
}
rm(s_objs); gc()
annot_ = do.call(annot_, what = rbind)
rownames(annot_) = annot_$cells

# Ordering the genes used by inferCNV in the cohort based on their genomic locus ####

genes_chr = read.delim(file = '../data/datasets/gene_ordering_file.txt', header = F, sep = '\t', quote = "", as.is = T, check.names = F, col.names = c('gene','chr','str','end'))     # gene-chromosome file; MUST BE THE ONE USED BY inferCNV
rownames(genes_chr) = genes_chr$gene
all_genes = genes_chr$gene[genes_chr$gene %in% all_genes]
names(all_genes) = genes_chr[all_genes,"chr"]

# forming universal obs mat ####

mat_ = matrix(data = 1, nrow = length(annot_$cells), ncol = length(all_genes), dimnames = list(annot_$cells, all_genes))
for(m_ in obs_mat)
{
  m_ = m_[,all_genes]
  mat_[ rownames(m_), all_genes] = m_
}
rm(obs_mat); gc()

# plotting ####

# preparing annotation table for pheatmap()


# grouping cancer cell subtypes into subcompartments
cell_types = c("blastema 1","S blastema 1","G2M blastema 1","blastema 2","G2M blastema 2","blastema 3",
               "smooth myocyte-like precursors","G2M smooth myocyte-like precursors","nascent smooth myocyte-like","smooth myocyte-like","S smooth myocyte-like","G2M smooth myocyte-like",
               "PAX3+ myogenic precursors","myoblast-like","striated myocyte-like",
               "tubules","S tubules","G2M tubules","podocyte-like")

annot_$codition = factor(annot_$codition, levels = c('untreated','neoadjuvant'))
annot_$compartment[annot_$cell_type %in% c("blastema 1","S blastema 1","G2M blastema 1","blastema 2","G2M blastema 2","blastema 3")] = 'blastema'
annot_$compartment[annot_$cell_type %in% c("smooth myocyte-like precursors","G2M smooth myocyte-like precursors","nascent smooth myocyte-like","smooth myocyte-like","S smooth myocyte-like","G2M smooth myocyte-like",
                                           "PAX3+ myogenic precursors","myoblast-like","striated myocyte-like")] = 'rhabdoid'
annot_$compartment[annot_$cell_type %in% c("tubules","S tubules","G2M tubules")] = 'tubules'
annot_$compartment[annot_$cell_type %in% c("podocyte-like")] = 'podocyte-like'
annot_$compartment = factor(x = annot_$compartment, levels = c('blastema','rhabdoid','tubules','podocyte-like'))

# setting up color annotation
comp_cols = c('red','green','orange','purple')
names(comp_cols) = levels(annot_$compartment)

cond_cols = c('blue','red')
names(cond_cols) = levels(annot_$codition)

annoCol = list(compartment = comp_cols,
               codition = cond_cols)

# ordering cancer cell types for every sample in the merged meta data table (annot_)

annot_tmp = list()
for(s_ in unique(annot_$sample))
{
  tmp_ = annot_[annot_$sample %in% s_,]
  tmp_ = tmp_[order(tmp_$compartment, tmp_$cell_type),]
  annot_tmp[[s_]] = tmp_
}
annot_ = do.call(annot_tmp, what = rbind)
rownames(annot_) = annot_$cells
mat_ = mat_[annot_$cells,]

# de-noising?

# cnv_lvls = seq(min(mat_), max(mat_), length.out = 16)
# mat_[ mat_ < cnv_lvls[6]] = min(mat_)
# mat_[ cnv_lvls[11] < mat_] = max(mat_)
# mat_[ cnv_lvls[6] <= mat_ & mat_ <= cnv_lvls[11] ] = 1

cnv_lvls = seq(0, 2, length.out = 16)
neutral_intvl = cnv_lvls[8:9]
mat_[ cnv_lvls[11] < mat_] = max(mat_)
mat_[ mat_ < cnv_lvls[6]] = min(mat_)
mat_[ neutral_intvl[1] <= mat_ & mat_ <= neutral_intvl[2] ] = 1

# breaks

# den_ = density(mat_)
# neutral_peak_ = den_$x[which.max(den_$y)]                                             # neutral CNV value
# 
# breaks_ = seq(min(mat_), max(mat_), length.out = 16)                                  # inferCNV has 15 intervals
# lower_ = max(which(breaks_ < neutral_peak_))                                          # neutral interval is here
# dels_ = colorRampPalette(colors =  c('blue4','white'))(8)[max(7-(lower_-2),1):7]      # based on which interval in neutral, some of the 8 blue shades are selected
# neutral_ = colorRampPalette(colors =  c('white'))(1)
# amps_ = colorRampPalette(colors =  c('white','red4'))(8)[2:(1+min(7,15-lower_))]      # based on which interval in neutral, some of the 8 red shades are selected
# cols_ = c(dels_, neutral_, amps_)

# Or

# neutral_peak_ = 1
# breaks_ = cnv_lvls
# lower_ = 8
# dels_ = colorRampPalette(colors =  c('blue4','white'))(8)[max(7-(lower_-2),1):7]      # based on which interval in neutral, some of the 7 blue shades are selected
# neutral_ = colorRampPalette(colors =  c('white'))(1)
# amps_ = colorRampPalette(colors =  c('white','red4'))(8)[2:(1+min(7,15-lower_))]      # based on which interval in neutral, some of the 8 red shades are selected
# cols_ = c(dels_, neutral_, amps_)
  
# row labels and gaps

feature_ = 'sample'
row_gaps = sapply(X = unique(annot_[,feature_]), FUN = function(s_){ max(which(annot_[,feature_] %in% s_))})                   # end position of each sample on the heatmap
row_labs = rep('',nrow(mat_))
row_labs[row_gaps] = as.character(unique(annot_[,feature_]))

# col labels and gaps

col_gaps = sapply(X = unique(names(all_genes)), FUN = function(chr_){ max(which(names(all_genes) %in% chr_)) })                  # gaps: end position of chromosomes on the heatmap
col_labs = rep('', ncol(mat_))
chr_index = round(sapply(X = unique(names(all_genes)), FUN = function(chr_){ mean(which(names(all_genes) %in% chr_)) }))      # where to put chromosome labels on x-axis
col_labs[chr_index] = names(chr_index)                                                                                         # column labels

# heatmap

p_ = pheatmap(mat = mat_,
              annotation_row = annot_[,c("codition","compartment")], annotation_colors = annoCol,
              gaps_col = col_gaps, labels_col = col_labs, labels_row = row_labs, gaps_row = row_gaps,
              legend_labels = c('loss','neutral','gain'), legend_breaks = c(min(mat_),1,max(mat_)),
              angle_col = 90,
              color = colorRampPalette(colors =  c('blue4','white','red4'))(15),
              # color = cols_,
              border_color = NA,
              cluster_rows = F, cluster_cols = F, treeheight_row = 0, treeheight_col = 0, clustering_method = 'complete',
              legend = T,
              show_rownames = T, show_colnames = T, silent = T)
p_$gtable$grobs[[2]]$gp = gpar(col = "white", fontsize = 10)      # col labels: p_$gtable
p_$gtable$grobs[[3]]$gp = gpar(col = "white", fontsize = 10)      # row lables
p_$gtable$grobs[[5]]$gp = gpar(col = 'white', fontsize = 10)      # row_annotation_names
p_$gtable$grobs[[6]]$gp = gpar(col = 'white', fontsize = 10)      # annotation_legend
p_$gtable$grobs[[7]]$gp = gpar(col = "white", fontsize = 10)      # legend labels

png(file = 'cnv_heatmap.png', width = 14, height = 13, units = 'in', res = 400, bg = 'black')
  grid::grid.newpage()
  grid::grid.draw(p_$gtable)
graphics.off()


# This scripts takes infercnv.references.txt and infercnv.observations.txt from inferCNV output and generates CNV heatmap per cluster and CNV UMAPs per chromosome
#
# in inferCNV expression values are cut in 15 equal intervals with those in the 8th interval as neutral/wildtype CNV.
# The centroid in the 8th interval is 1, meaning all gene expressions within (1 - epi, 1 + epi) are neutral. Those expressions
# from 1st through 7th intervals are shown by shades of blue representing deletion. Those expression values from 9t through
# 15th intervals are amplification shown by shades of red

# Att. In case of using integrated data, an extra id was added to cells during integration by Seurat to avoid duplicates.
# Thus, inferCNV's cells have that id, perform this:
# s_obj = RenameCells(object = s_obj, new.names = annot$cells[annot$orig.ident %in% smpl_])

library(Seurat)
library(infercnv)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(viridis)
set.seed(41)

in_dir = '../data/untreated/'
out_dir = '../out/untreated/'

# Reading in data ####

# W008 sample data
load(file = '../data/untreated/W008/data_QC_COMPRT.RData')      # original sample

# gene ordering file
genes_chr = read.delim(file = '../data/datasets/gene_ordering_file.txt',
                       header = F, sep = '\t', quote = "", as.is = T, check.names = F, col.names = c('gene','chr','str','end'))     # gene-chromosome file; MUST BE THE ONE USED BY inferCNV
rownames(genes_chr) = genes_chr$gene

# reference matrix
ref_mat = t(as.matrix(read.delim(file = '../out/untreated/W008/inferCNV/inferCNV_TCells/infercnv.references.txt',
                                 header = T, row.names = 1, sep = ' ', quote = '"', as.is = T, check.names = F)))                   # reference matrix normalized and returned by inferCNV; has the same cells for all samples but different genes

# extracting T reference cells related to this sample (if any)
ref_cells = which(grepl(x = rownames(ref_mat), pattern = paste0('^(W008)_')))                                                       # reference cells in this sample
if(length(ref_cells) != 0)
{
  ref_mat = ref_mat[ref_cells,]
  rownames(ref_mat) = gsub(pattern = 'W[0-9]+_', replacement = '', x = rownames(ref_mat))                                           # removing WXXXX_ from the beginning of reference cell names
}else
{
  ref_mat = NULL
}

# observation matrix
obs_mat = t(as.matrix(read.delim(file = '../out/untreated/W008/inferCNV/inferCNV_TCells/infercnv.observations.txt',
                                 header = T, row.names = 1, sep = ' ', quote = '"', as.is = T, check.names = F)))                   # all non-reference cells normalized and returned by inferCNV including immune cells

# Forming cnv matrix (observation + reference) ####

cnv_mat = do.call(list(ref_mat, obs_mat), what = rbind); rm(ref_mat, obs_mat); gc()
s_obj = subset(s_obj, cells = colnames(s_obj)[colnames(s_obj) %in% row.names(cnv_mat)])                                             # low-quality cells were not passed on to inferCNV
cnv_mat = cnv_mat[colnames(s_obj),]                                                                                                 # reordering cells in cnv_mat to match original count matrix
genes_ = colnames(cnv_mat)                                                                                                          # all genes in this samples
chrs_ = genes_chr[genes_,"chr"]                                                                                                     # chromosome of each gene
names(genes_) = chrs_

# Plotting CNV heatmap ####

clusters_ = s_obj$compartment                                                           # clusters/compartments to show
clusters_ = droplevels(clusters_)
clusters_= sort(clusters_)                                                              # sorting clusters in ascending order

mat_ = cnv_mat[names(clusters_),]                                                       # rearranging CNV matrix by cluster order
breaks_ = seq(min(mat_), max(mat_), length.out = 16)                                    # 15 intervals for heatmap; in inferCNV there are 15 intervals/color shades.
                                                                                        # In pheatmap break points (ticks) should be one element larger than color vector

gaps_cols = sapply(X = unique(chrs_), FUN = function(c_){ max(which(chrs_ == c_))})     # vertical gaps in heatmap, one for each chromosome
plot_list = list()                                                                      # list to store plots
for(clst_ in levels(clusters_))
{
  # main heatmap
  message(clst_)
  tmp_ = mat_[clusters_ %in% clst_,]
  tmp_ = tmp_[sample(x = 1:nrow(tmp_), size = min(1000, nrow(tmp_))),]      # subsampling to shorten processing time!
  p_ = pheatmap(mat = tmp_,
                color = colorRampPalette(colors =  c('blue4','white','red4'))(15),
                border_color = NA,
                gaps_col = gaps_cols,
                cluster_cols = F, cluster_rows = T, clustering_method = 'ward.D', treeheight_row = 0,
                show_rownames = F, show_colnames = F,
                cellheight = 32/nrow(tmp_),
                breaks = breaks_,
                legend = F,
                silent = T)
  plot_list[[clst_]] = p_[[4]]

  # right heatmap
  labs_row = rep('',nrow(tmp_))
  labs_row[round(length(labs_row)/2)] = clst_
  p_ = pheatmap(mat = matrix(0, ncol = 1, nrow = nrow(tmp_)),
                color = 'white',                                            # color of rugs in front of cluster labels
                border_color = NA,
                cluster_rows = F, cluster_cols = F,
                show_rownames = T,labels_row = labs_row,
                show_colnames = F,
                fontsize_row = 10,
                cellheight = 32/nrow(tmp_), cellwidth = .75,
                breaks = 0,
                legend = F,
                silent = T)
  p_[[4]]$grobs[[2]]$gp$col = 'white'                                       # color of cluster labels
  plot_list[[length(plot_list)+1]] = p_[[4]]
}

# adding chromosome names to the bottom of heatmap
labels_col = chrs_
labels_col[-gaps_cols] = ''                                                             # column names should contain only one chromosome name for each column
p_ = pheatmap(mat = matrix(0, nrow = 1, ncol = ncol(mat_)),
              color = 'white',
              border_color = NA,
              gaps_col = gaps_cols,
              cluster_rows = F, cluster_cols = F,
              show_rownames = F, show_colnames = T, labels_col = labels_col,
              cellheight = 0.75,
              angle_col = 90, fontsize_col = 9.5,
              breaks = 0,
              legend = F,
              silent = T)
p_[[4]]$grobs[[2]]$gp$col = 'white'
plot_list[[length(plot_list)+1]] = p_[[4]]

# plotting
h_ = .6*length(levels(clusters_))
jpeg(file = 'inferCNV_W008.jpg', width = 3.3*h_, height = h_, units = 'in', res = 400, bg = 'black')
  grid.arrange(grobs= plot_list, nrow = length(levels(clusters_))+1, ncol = 2, as.table = T, widths = c(1, 0.01))
graphics.off()


# Calculating average CNV per chromosome ####

vals_ = seq(0,2, length.out = 16)
neutral_range = vals_[8:9]                                                  # the 8th interval in inferCNV is considered wildtype CNV
cnv_mat[ neutral_range[1] <= cnv_mat & cnv_mat <= neutral_range[2]] = 1     # all values falling in the 8th interval are set to 1 (the centroid)
cnv_mat = abs(cnv_mat-1)                                                    # to make sure neutral values are assigned black color; amplification or deletion contribute synergistically

s_obj$avg_cnv = 0                                                           # average CNV across all chromosomes
i_ = 0                                                                      # chromosome counter
for(chr_ in unique(chrs_))
{
  s_obj[[paste0('avg_cnv_',chr_)]] = rowMeans(cnv_mat[, names(genes_) %in% chr_])
  s_obj$avg_cnv = s_obj$avg_cnv + s_obj[[paste0('avg_cnv_',chr_)]]
  i_ = i_+1
}
s_obj$avg_cnv = (s_obj$avg_cnv/i_)

# Plotting average CNV for each cell per chromosome ####

pdf(file = 'inferCNV_avg_W008.pdf', width = 40, height = 40)

## average CNV per chromosome ####

p_ = list()     # plot list
for(chr_ in unique(chrs_))
{
  col_ = paste0('avg_cnv_',chr_)      # column of meta.data containing average cnv for chr
  breaks_ = c(min(s_obj[[col_]]), mean(range(s_obj[[col_]])), max(s_obj[[col_]]))
  p_[[chr_]] =  FeaturePlot(s_obj, features = col_, pt.size = 1, order = T, coord.fixed = F)+
                theme(legend.spacing.y = unit(0.5, 'cm'), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), text = element_text(size = 22, face = 'bold'), plot.title = element_text(size = 25, face = 'bold'))+
                labs(title = chr_, color = 'CNV level')+
                scale_color_viridis(option = 'A', breaks = breaks_, labels = round(breaks_,2))
}

## average CNV across aberrant chromosomes ####

s_obj$avg_cnv_comb = rowMeans(s_obj@meta.data[,c('avg_cnv_chr1','avg_cnv_chr6','avg_cnv_chr7','avg_cnv_chr8','avg_cnv_chr12','avg_cnv_chr13')])
cutoff_ = quantile(x = s_obj$avg_cnv_comb, probs = c(.2,.99))
breaks_ = c(cutoff_, mean(range(cutoff_)))
p_[[length(p_)+1]] = FeaturePlot(s_obj, features = 'avg_cnv_comb', pt.size = 1, order = T, coord.fixed = F,
                                 min.cutoff = cutoff_[1], max.cutoff = cutoff_[2])+
                     theme(plot.title = element_text(size = 25, face = 'bold'),
                           axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                           text = element_text(size = 22, face = 'bold'))+
                     labs(title = 'Aberrant chromosomes', color = 'CNV level')+
                     scale_color_viridis(option = 'A', breaks = breaks_, labels = round(breaks_,2))

do.call(what = grid.arrange, args = c(grobs = p_, nrow = 5))

graphics.off()


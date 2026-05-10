# This script plots a cleaned Jaccard similarity matrix of clustered factors where each cluster is called a program.

# Input: It needs Jaccard_similarity_matrix_programs.tsv and manually annotated programs_annotation_processed.tsv from NMF_metaprogram_heatmap.R.
# Output: It outputs a cleaned heatmap (NMF_programs_cleaned.pdf) of clustered factors where each cluster is called a program.
#         The file programs_annotation.tsv MUST BE ANNOTATED (now named programs_annotation_processed.tsv) USIGING THE OUTPUT
#         HEATMAP (programs.pdf) AND NMF_GSEA_1e-3.tsv.

# NMF-specific terminology:

# Feature: genes
# Sample: cells (in single cell data)
# Factor: basis component (similar to principle component) or a metagene
# Rank: max number of factors (similar to rank in prcomp() for PCA)
# V: data matrix: genes × cells (unlike prcomp that input is cells x genes)
# W: basis matrix (similar to rotation/loadings matrix in PCA): genes × factors
# H: Coefficient Matrix (similar to x matrix in PCA but always non-negative): factors × samples (cells). Unlike PCA, each value is a metagene expression not a coordinate!

library(pheatmap)
library(gridExtra)
library(viridis)
set.seed(42)

# Reading in data ####

# reading in raw Jaccard similarity matrix
mat_ = as.matrix(read.delim(file = 'Jaccard_similarity_matrix_programs.tsv', header = T, sep = '\t', as.is = T, check.names = F))

# reading in modified annotation file
# Att.: FIRST ANNOTATE open programs_annotation.tsv using NMF_programs.pdf and NMF_GSEA_1e-3.tsv files

annot_ = read.delim(file = 'programs_annotation_processed.tsv', header = T, sep = '\t', as.is = T, check.names = T)     # manually annotated programs_annotation.tsv

# programs identified in programs_annotation_processed.tsv for ordering levels
progs_ = c('blastema 1','S blastema 1','G2M blastema 1',
           'blastema 2',
           'PAX3+ myogenic precursors','striated myocyte',
           'smooth-myocyte precursors','nascent smooth myocyte','smooth myocyte',
           'podocyte',
           'tubular','connecting tubule')

# reformatting annotation file for ggpolt2
ls_ = list()
for(p_ in progs_)
{
  ls_[[p_]] = annot_[annot_$program %in% p_,]
  ls_[[p_]] = ls_[[p_]][order(ls_[[p_]]$condition),]
}
annot_ = do.call(ls_, what = rbind)
rownames(annot_) = annot_$factor_sample
annot_$program = factor(x = annot_$program, levels = progs_)

# plotting final Jaccard heatmap

pdf(file = 'NMF_programs_cleaned.pdf', width = 10, height = 10)

cols_1 = sample(viridis(n = length(progs_), option = 'turbo'))                    # color for each program
names(cols_1) = levels(annot_$program)

cols_2 = sample(viridis(n = length(unique(annot_$sample)), option = 'turbo'))     # color for each patient sample
names(cols_2) = unique(annot_$sample)

cols_3 = c(untreated = 'blue', neoadjuvant = 'red')                               # color for treatment condition
annot_cols = list(program = cols_1, sample = cols_2, condition = cols_3)

m_ = mat_[rownames(annot_), rownames(annot_)]
m_[m_ <= 0.05] = 0
p_ = pheatmap(mat = m_, annotation_row = annot_[,-1, drop = F], annotation_col = annot_[,-1, drop = F],
              cluster_cols = F, cluster_rows = F, clustering_method = 'ward.D2', treeheight_row = 0, treeheight_col = 0,
              border_color = NA, scale = 'none', color = rev(viridis(n = 100, option = 'A')),
              show_rownames = F, show_colnames = F, fontsize = 2,
              annotation_colors = annot_cols,
              legend_breaks = c(0,.5,1),
              silent = F)
p_
graphics.off()

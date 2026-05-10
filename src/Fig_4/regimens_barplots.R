library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)
library(dplyr)
library(viridis)

# Reading in data ####

load('../data/integration/healthy_cancer/TME.RData')

tme_$regimen = NA
tme_$regimen[tme_$orig.ident %in% c('W010','W009','W003','W441')] = 'EE-4A'
tme_$regimen[tme_$orig.ident %in% c('W184','W180','W370')] = 'DD-4A'
tme_$regimen[tme_$orig.ident %in% c('W002')] = 'I'
tme_$regimen[tme_$orig.ident %in% c('W050','W944')] = 'V+D'
tme_$regimen[tme_$orig.ident %in% c('W006','W001')] = 'Unk'

pdf(file = 'regimens_TME.pdf', width = 7.5, height = 7.5)

#_______________________________________________________________
# TME ####

cells_ = colnames(tme_)[tme_$condition %in% 'neoadjuvant']
dt_tme = tme_@meta.data[cells_,]
DimPlot(tme_, cells = cells_, group.by = 'subcompartment', label = T, repel = T, raster = F)

## plotting ####

# formatting data.frame for ggplot2
dt_ = list()
for(p_ in unique(dt_tme$orig.ident))
{
  if(p_ %in% c('W006','W001')){ next() }
  tmp_ = dt_tme[dt_tme$orig.ident %in% p_,]
  portions_ = round(table(tmp_$subcompartment)/nrow(tmp_)*100,2)
  dt_[[p_]] = data.frame(Patient = p_, portions_[c('blastema','lymphoid','myeloid','vasculature')], Regimen = tmp_$regimen[1])
}
dt_tme = do.call(dt_, what = rbind)
colnames(dt_tme)[2] = 'Component'
dt_tme = dt_tme[!is.na(dt_tme$Component),]     # not all compartments are present in all patients
dt_tme$Patient = factor(x = dt_tme$Patient, level = unique(dt_tme$Patient[order(dt_tme$Regimen)]))
dt_tme$Component = factor(x = dt_tme$Component, levels = c('vasculature','lymphoid','myeloid', 'blastema'))

# plotting
ggplot(dt_tme, aes(x = Patient, y = Freq, fill = Component))+
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(x = "Patient ID", y = "Percentage (%)", fill = "Compartment")+
geom_bar(stat = "identity")+
geom_text_repel(aes(label = sprintf("%.2f%%", Freq)),
                position = position_stack(vjust = .5),
                size = 2.5,
                color = "yellow", bg.colour = 'black', fontface = 'bold',
                min.segment.length = 0,
                box.padding = 0.1,
                point.padding = 0.1,
                segment.color = "transpblackarent", seed = 42)+
scale_fill_manual(values = c("lymphoid" = "green3",
                             "vasculature" = "blue",
                             "myeloid" = "green4",
                             "blastema" = 'red'))
#_______________________________________________________________
# Cancer compartment ####

cells_ = colnames(tme_)[tme_$compartment_name %in% 'cancer' & tme_$condition %in% 'neoadjuvant']
dt_cancer = tme_@meta.data[cells_,]
DimPlot(tme_, cells = cells_,group.by = 'subcompartment', label = T, repel = T, raster = F)

## plotting ####

# formatting data.frame for ggplot2
dt_ = list()
for(p_ in unique(dt_cancer$orig.ident))
{
  if(p_ %in% c('W006','W001')){ next() }
  tmp_ = dt_cancer[dt_cancer$orig.ident %in% p_,]
  portions_ = round(table(tmp_$subcompartment)/nrow(tmp_)*100,2)
  dt_[[p_]] = data.frame(Patient = p_, portions_, Regimen = tmp_$regimen[1])
}
dt_cancer = do.call(dt_, what = rbind)
colnames(dt_cancer)[2] = 'Component'
dt_cancer$Patient = factor(x = dt_cancer$Patient, level = unique(dt_cancer$Patient[order(dt_cancer$Regimen)]))
dt_cancer$Component = factor(x = dt_cancer$Component, levels = c('epithelium','stroma','blastema'))

# plotting
ggplot(dt_cancer, aes(x = Patient, y = Freq, fill = Component))+
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(x = "Patient ID", y = "Percentage (%)", fill = "Tumor Component")+
geom_bar(stat = "identity")+
geom_text_repel(aes(label = sprintf("%.2f%%", Freq)),
                position = position_stack(vjust = .5),
                size = 2.5,
                color = "yellow", bg.colour = 'black', fontface = 'bold',
                min.segment.length = 0,
                box.padding = 0.1,
                point.padding = 0.1,
                segment.color = "transpblackarent", seed = 42)+
scale_fill_manual(values = c("blastema" = "red",
                             "epithelium" = "blue",
                             "stroma" = "green4"))

#_______________________________________________________________
# Myeloid compartment ####

cells_ = colnames(tme_)[tme_$subcompartment %in% 'myeloid' & tme_$condition %in% 'neoadjuvant']
dt_myeloid = tme_@meta.data[cells_,]
DimPlot(tme_, cells = cells_,group.by = 'cell_type', label = T, repel = T, raster = F)

## plotting ####

# formatting data.frame for ggplot2
dt_ = list()
for(p_ in unique(dt_myeloid$orig.ident))
{
  if(p_ %in% c('W006','W001')){ next() }
  tmp_ = dt_myeloid[dt_myeloid$orig.ident %in% p_,]
  portions_ = round(table(tmp_$cell_type)/nrow(tmp_)*100,2)
  dt_[[p_]] = data.frame(Patient = p_, portions_, Regimen = tmp_$regimen[1])
}
dt_myeloid = do.call(dt_, what = rbind)
colnames(dt_myeloid)[2] = 'Cell_Type'
dt_myeloid$Patient = factor(x = dt_myeloid$Patient, level = unique(dt_myeloid$Patient[order(dt_myeloid$Regimen)]))
dt_myeloid$Cell_Type = factor(x = dt_myeloid$Cell_Type, levels = c("tDC", "cDC2", "pDC", "mast", "cDC1", "monocyte", "M2d-like", "M2a-like", "M2c-like"))

# plotting
ggplot(dt_myeloid, aes(x = Patient, y = Freq, fill = Cell_Type))+
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(x = "Patient ID", y = "Percentage (%)", fill = "Cell_Type")+
geom_bar(stat = "identity")+
geom_text_repel(aes(label = sprintf("%.2f%%", Freq)),
                position = position_stack(vjust = .5),
                size = 2.5,
                color = "yellow", bg.colour = 'black', fontface = 'bold',
                min.segment.length = 0,
                box.padding = 0.1,
                point.padding = 0.1,
                segment.color = "transpblackarent", seed = 42)+
# scale_fill_manual(values = )
scale_fill_viridis(option = 'H', discrete = T)

#_______________________________________________________________
# Lymphoid compartment ####

cells_ = colnames(tme_)[tme_$subcompartment %in% 'lymphoid' & tme_$condition %in% 'neoadjuvant']
dt_lymphoid = tme_@meta.data[cells_,]
DimPlot(tme_, cells = cells_,group.by = 'cell_type', label = T, repel = T, raster = F)

## plotting ####

# formatting data.frame for ggplot2
dt_ = list()
for(p_ in unique(dt_lymphoid$orig.ident))
{
  if(p_ %in% c('W006','W001')){ next() }
  tmp_ = dt_lymphoid[dt_lymphoid$orig.ident %in% p_,]
  portions_ = round(table(tmp_$cell_type)/nrow(tmp_)*100,2)
  dt_[[p_]] = data.frame(Patient = p_, portions_, Regimen = tmp_$regimen[1])
}
dt_lymphoid = do.call(dt_, what = rbind)
colnames(dt_lymphoid)[2] = 'Cell_Type'
dt_lymphoid$Patient = factor(x = dt_lymphoid$Patient, level = unique(dt_lymphoid$Patient[order(dt_lymphoid$Regimen)]))
dt_lymphoid$Cell_Type = factor(x = dt_lymphoid$Cell_Type, levels = c("plasma", "CD8 Tem", "Treg", "B", "CD8+ Tcm", "NKT", "CD4+ Tcm", "NK", "CD8+ Teff", "CD8+ Tex", "CD4+ Tex"))

# plotting
ggplot(dt_lymphoid, aes(x = Patient, y = Freq, fill = Cell_Type))+
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(x = "Patient ID", y = "Percentage (%)", fill = "Cell_Type")+
geom_bar(stat = "identity")+
geom_text_repel(aes(label = sprintf("%.2f%%", Freq)),
                position = position_stack(vjust = .5),
                size = 2.5,
                color = "yellow", bg.colour = 'black', fontface = 'bold',
                min.segment.length = 0,
                box.padding = 0.1,
                point.padding = 0.1,
                segment.color = "transpblackarent", seed = 42)+
# scale_fill_manual(values = )
scale_fill_viridis(option = 'H', discrete = T)

graphics.off()

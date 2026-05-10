# This script computes FGA using the outpout of ichorCNA applied to lpWGS.
# Each layer is called track (gender, age, etc.) and each section (a patient) of a track is called sector

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(SeuratWrappers)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(gridBase)
library(viridis)
library(tidyverse)


# Reading in data ####

# cohort data

cohort_ = read.delim(file = '../doc/cohort_study/cohort.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
rownames(cohort_) = cohort_$patient

# FGA and TMB data

fg_tmb = read.delim(file = '../out/FGA_TMB.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
rownames(fg_tmb) = fg_tmb$Patient
for(p_ in rownames(cohort_))
{
  fga_tmb = fg_tmb[fg_tmb$Patient %in% p_, c("FGA","TMB")]
  if(nrow(fga_tmb) == 1){ cohort_[p_, c('FGA', 'TMB')] = fga_tmb }else{ cohort_[p_, c('FGA', 'TMB')] = c(NA,NA) }
}
cohort_$neoadjuvant[cohort_$treated %in% 'untreated'] = 'unk'

# tumor cellularity data

load('../data/integration/healthy_cancer/TME.RData')
tme_ = tme_@meta.data; gc()
tme_ = cbind(blastema   = round(table(tme_$orig.ident[tme_$subcompartment %in% 'blastema'])/table(tme_$orig.ident)*100, 2),
             epithelium = round(table(tme_$orig.ident[tme_$subcompartment %in% 'epithelium'])/table(tme_$orig.ident)*100, 2),
             stroma     = round(table(tme_$orig.ident[tme_$subcompartment %in% 'stroma'])/table(tme_$orig.ident)*100, 2))
tme_ = as.data.frame(tme_[rownames(cohort_),])
tme_$patient = rownames(tme_)

# converting mutant data to long format using pivot_longer()

tme_ = as.data.frame(pivot_longer(tme_,
                                  cols = c(blastema,epithelium,stroma),     # columns to be reshaped
                                  names_to = "subcompartment",              # name of the new column for variable names
                                  values_to = "frequency"))                 # name of the new column for values
tme_$subcompartment = factor(x = tme_$subcompartment, levels = c('blastema','stroma','epithelium'))

# Plotting circlize ####

graphics.off()

pdf(file = 'cohort_summary.pdf', width = 15.71, height = 14.81)

# initialize circos plot

par(mar = c(1, 1, 1, 1))

circos.par(start.degree = 90,                 # Start at 45 degrees
           gap.degree = c(rep(1, 23),90),     # Add a 90-degree gap for the unused part
           cell.padding = c(0, 0, 0, 0),      # No padding
           track.height = 0.1,              # default sector height
           track.margin = c(0.0, 0.008))

# create the circos plot
circos.initialize(factors = cohort_$patient, xlim = c(0, 3))

#########################################################################################


# genomic events _____________________________________________________

# 1a: tumor type content
col_fun = list(blastema    = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T),    max(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T)),    c("pink", "red3")),
               stroma      = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T),      max(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T)),      c("lightgreen", "green4")),
               epithelium  = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "epithelium"] , na.rm = T), max(tme_$frequency[tme_$subcompartment == "epithelium"] , na.rm = T)), c("lightblue","blue3")))

circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, max(tme_$frequency)*1.2),
                        track.height = 0.15,
                        bg.border = 'black',
                        bg.col = NA,
                        panel.fun = function(x, y)
                        {
                          sector_index = get.cell.meta.data("sector.index")
                          
                          # filter data for the current sector (patient)
                          patient_data = tme_[tme_$patient == sector_index, ]
                          
                          # define bar positions and colors
                          bar_positions = c(.5, 1.5, 2.5)  # Positions for 3 bars
                          
                          for (i_ in 1:length(bar_positions))
                          {
                            subcompartment = levels(patient_data$subcompartment)[i_]
                            feature_data = patient_data[patient_data$subcompartment == subcompartment, ]
                            
                            # get the color based on the frequency using the gradient
                            bar_color = col_fun[[subcompartment]](feature_data$frequency)
                            
                            # draw bars for the current sector
                            circos.barplot(value = feature_data$frequency,
                                           bar_width = .85,
                                           pos = bar_positions[i_],
                                           col = bar_color,
                                           border = NA)
                            
                            # grid lines
                            # breaks = seq(0, CELL_META$ylim[2], length.out = 6)
                            # for(b in breaks){ circos.lines(x = CELL_META$cell.xlim, y = rep(b, 2), lty = 3, col = "black") }
                            breaks_ = round(seq(min(tme_$frequency, na.rm = T), max(tme_$frequency, na.rm = T), length.out = 3),2)
                            for(b_ in breaks_)
                            {
                              # grid lines
                              circos.lines(x = CELL_META$cell.xlim, y = rep(b_, 2), lty = 3, col = "black")
                              
                              # adding y-axis ticks
                              if(CELL_META$sector.numeric.index == 1){ circos.yaxis(side = 'left', at = b_, labels = T) }                                                     # first sector
                              if(1 < CELL_META$sector.numeric.index & CELL_META$sector.numeric.index < nrow(cohort_)){ circos.yaxis(side = 'left', at = b_, labels = F) }     # all other sectors
                              if(CELL_META$sector.numeric.index == nrow(cohort_)){ circos.yaxis(side = 'right', at = b_, labels = F) }                                        # last sector
                            }
                          }
                        })

# 1b: FGA
col_fun = colorRamp2(breaks = c(min(cohort_$FGA, na.rm = T), max(cohort_$FGA, na.rm = T)), colors = c("pink", "purple4"))
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, max(cohort_$FGA, na.rm = T)*1.2),
                        bg.border = "black",  # Box around the track
                        bg.col = NA,
                        panel.fun = function(x, y)
                        {
                          # drawing bars for the current sector
                          val_ = cohort_$FGA[CELL_META$sector.numeric.index]
                          circos.barplot(value = val_,
                                         bar_width = .85,
                                         pos = 1.5,
                                         col = col_fun(val_),
                                         border = NA)

                          # adding ticks and grid lines
                          breaks_ = round(seq(min(cohort_$FGA, na.rm = T), max(cohort_$FGA, na.rm = T), length.out = 3),2)
                          for(b_ in breaks_)
                          {
                            # grid lines
                            circos.lines(x = CELL_META$cell.xlim, y = rep(b_, 2), lty = 3, col = "black")
                            
                            # adding y-axis ticks
                            if(CELL_META$sector.numeric.index == 1){ circos.yaxis(side = 'left', at = b_, labels = T) }                                                     # first sector
                            if(1 < CELL_META$sector.numeric.index & CELL_META$sector.numeric.index < nrow(cohort_)){ circos.yaxis(side = 'left', at = b_, labels = F) }     # all other sectors
                            if(CELL_META$sector.numeric.index == nrow(cohort_)){ circos.yaxis(side = 'right', at = b_, labels = F) }                                        # last sector
                          }
                        })

# 1c: TMB
col_fun = colorRamp2(breaks = c(min(cohort_$TMB, na.rm = T), max(cohort_$TMB, na.rm = T)), colors = c("khaki", "brown4"))
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, max(cohort_$TMB, na.rm = T)*1.2),
                        bg.border = "black",  # Box around the track
                        bg.col = NA,
                        panel.fun = function(x, y)
                        {
                          # drawing bars for the current sector
                          val_ = cohort_$TMB[CELL_META$sector.numeric.index]
                          circos.barplot(value = val_,
                                         bar_width = .85,
                                         pos = 1.5,
                                         col = col_fun(val_),
                                         border = NA)
                          
                          # adding ticks and grid lines
                          breaks_ = round(seq(min(cohort_$TMB, na.rm = T), max(cohort_$TMB, na.rm = T), length.out = 3),2)
                          for(b_ in breaks_)
                          {
                            # grid lines
                            circos.lines(x = CELL_META$cell.xlim, y = rep(b_, 2), lty = 3, col = "black")
                            
                            # adding y-axis ticks
                            if(CELL_META$sector.numeric.index == 1){ circos.yaxis(side = 'left', at = b_, labels = T) }                                                     # first sector
                            if(1 < CELL_META$sector.numeric.index & CELL_META$sector.numeric.index < nrow(cohort_)){ circos.yaxis(side = 'left', at = b_, labels = F) }     # all other sectors
                            if(CELL_META$sector.numeric.index == nrow(cohort_)){ circos.yaxis(side = 'right', at = b_, labels = F) }                                        # last sector
                          }
                        })

# demographics ___________________________________________________________

# 2: race
cols_ = c('Black,Hispanic/Latino' = 'yellow4', "White" = "pink", "Black" = "violet", "Asian" = "purple2", 'Hispanic/Latino' = 'purple4', 'unk' = NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.027,
                        bg.col = cols_[cohort_$race],
                        bg.border = NA,  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA, #cohort_$race[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 3: sex
cols_ = c("female" = "lightgreen", "male" = "darkgreen", 'unk' = NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.027,
                        bg.col = cols_[cohort_$sex],
                        bg.border = NA,  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,#cohort_$sex[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 4: age
col_fun = colorRamp2(breaks = c(min(cohort_$age, na.rm = T), max(cohort_$age, na.rm = T)), colors = c("yellow","darkorange4"))
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = range(cohort_$age),
                        track.height = 0.027,
                        bg.col = col_fun(cohort_$age),
                        bg.border = NA,  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = mean(cohort_$age),
                                      labels = NA,#cohort_$age[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })


# 5: patient
cols_ = c("treated" = "red", "untreated" = "blue")
fontcols_ = rep(c('black','white'), each = 12)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        bg.col = cols_[cohort_$treatment],
                        bg.border = NA,  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = cohort_$patient[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = T,
                                      adj = c(0.5, 0.5),
                                      cex = 1,
                                      col = fontcols_[CELL_META$sector.numeric.index],
                                      font = 2)
                        })


# pathology info ________________________________________________________________________________________

# 6: survival
cols_ = c('dead' = 'black', 'alive' = 'grey75', 'unk' =  NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$survival],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 7: cancer stage (I-V)
col_fun = colorRamp2(breaks = c(min(cohort_$stage, na.rm = T), max(cohort_$stage, na.rm = T)), colors = c("pink","purple"))
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = range(cohort_$stage),
                        track.height = 0.026,
                        bg.col = col_fun(cohort_$stage),
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = mean(cohort_$stage),
                                      labels = NA,#cohort_$stage[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 8: recurrence
cols_ = c( 'yes' = 'black', 'no' = 'yellow')
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0,1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$recurrence],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 9: tumor laterality
cols_ = c('left' = '#D2B48C', 'right' = 'skyblue', 'bilateral' = '#4A2F26', 'unk' =  NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$laterality],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 10: tumor size (<=10, >=10)
cols_ = c( '<=10' = 'orange', '>=10' = 'black')
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0,1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$size],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 11: tumor weight
cohort_$weight = as.numeric(cohort_$weight)      # as.numeric() converts 'unk' to NA
col_fun = colorRamp2(breaks = c(min(cohort_$weight, na.rm = T), max(cohort_$weight, na.rm = T)), colors = c("lightblue","darkblue"))
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = range(cohort_$weight, na.rm = T),
                        track.height = 0.026,
                        bg.col = col_fun(as.numeric(cohort_$weight)),
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = mean(cohort_$weight, na.rm = T),
                                      labels = NA,#cohort_$weight[CELL_META$sector.numeric.index],
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "red",
                                      font = 2)
                        })

# 12: tumor histology (favorable, anaplastic)
cols_ = c( 'favorable' = 'green', 'anaplastic' = 'black')
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0,1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$histology],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 13: tumor focality (unifocal, multifocal)
cols_ = c( 'unifocal' = '#602CA1', 'multifocal' = '#BF5842')
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0,1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$focality],
                        bg.border = 'white',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 14: neoadjuvant
cols_ = c('EE-4A' = 'orange', 'DD-4A' = 'red', 'Vincristine+Doxorubicin' = '#AA9922','I' = 'grey35', 'unk' =  NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$neoadjuvant],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# 15: response to neoadjuvant
cols_ = c('decreased' = 'green2', 'increased' = 'red', 'unchanged' = 'grey50', 'unk' =  NA)
circos.trackPlotRegion( sectors = cohort_$patient,
                        ylim = c(0, 1),
                        track.height = 0.026,
                        bg.col = cols_[cohort_$response],
                        bg.border = 'black',  # Box around the track
                        panel.fun = function(x, y)
                        {
                          circos.text(x = CELL_META$xcenter, y = 0.5,
                                      labels = NA,
                                      facing = "clockwise",
                                      niceFacing = TRUE,
                                      adj = c(0.5, 0),
                                      cex = 0.6,
                                      col = "black",
                                      font = 2)
                        })

# Legends ###################################################################################################

# genomic events ______________________________________________________________

# 1a: tumor type cellularity
col_fun = list(blastema    = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T),    max(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T)),    c("pink", "red3")),
               stroma      = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T),      max(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T)),      c("lightgreen", "green4")),
               epithelium  = colorRamp2( c(min(tme_$frequency[tme_$subcompartment == "epithelium"] , na.rm = T), max(tme_$frequency[tme_$subcompartment == "epithelium"] , na.rm = T)), c("lightblue","blue3")))

blastema_ = Legend(at = round(c(min(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T), sum(range(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T))/2, max(tme_$frequency[tme_$subcompartment == "blastema"], na.rm = T)),2),
                              col_fun = col_fun$blastema,
                              title_position = "topleft",
                              title = "Blastema content (%)",
                              direction = 'horizontal')
stroma_ = Legend(at = round(c(min(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T), sum(range(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T))/2, max(tme_$frequency[tme_$subcompartment == "stroma"], na.rm = T)),2),
                            col_fun = col_fun$stroma,
                            title_position = "topleft",
                            title = "Stroma content (%)",
                            direction = 'horizontal')
epi_    = Legend(at = round(c(min(tme_$frequency[tme_$subcompartment == "epithelium"], na.rm = T), sum(range(tme_$frequency[tme_$subcompartment == "epithelium"], na.rm = T))/2, max(tme_$frequency[tme_$subcompartment == "epithelium"], na.rm = T)),2),
                            col_fun = col_fun$epithelium,
                            title_position = "topleft",
                            title = "Epithelium content (%)",
                            direction = 'horizontal')

# 1b: FGA

col_fun = colorRamp2(breaks = c(min(cohort_$FGA, na.rm = T), max(cohort_$FGA, na.rm = T)), colors = c("pink", "purple4"))
fga_ = Legend(at = round(c(min(cohort_$FGA, na.rm = T), sum(range(cohort_$FGA, na.rm = T))/2, max(cohort_$FGA, na.rm = T)),2),
              col_fun = col_fun,
              title_position = "topleft",
              title = "FGA (%)",
              direction = 'horizontal')

# 1c: TMB
col_fun = colorRamp2(breaks = c(min(cohort_$TMB, na.rm = T), max(cohort_$TMB, na.rm = T)), colors = c("khaki", "brown4"))
TMB_ = Legend(at = round(c(min(cohort_$TMB, na.rm = T), sum(range(cohort_$TMB, na.rm = T))/2, max(cohort_$TMB, na.rm = T)),2),
              col_fun = col_fun,
              title_position = "topleft",
              title = "TMB (Mut/Mb)",
              direction = 'horizontal')


# demographics ______________________________________________________________

# 2: race
race_ = Legend(at = c("Black,Hispanic/Latino", "White", "Black", "Asian", "Hispanic/Latino"), type = "points",
               legend_gp = gpar(col = c("yellow4","pink", "violet", "purple2", "purple4")),
               pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Race", title_position = "topleft",
               background = NULL,
               border = NA,
               ncol = 1)

# 3: sex
gender_ = Legend(at = c("Female", "Male"),
                 legend_gp = gpar(col = c('lightgreen','darkgreen')),
                 type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
                 title = "Sex", title_position = "topleft",
                 background = NULL,
                 border = NA,
                 nrow = 1)

# 4: age
col_fun = colorRamp2(breaks = c(min(cohort_$age), max(cohort_$age, na.rm = T)), colors = c("yellow","darkorange4"))
age_ = Legend(at = round(c(min(cohort_$age), sum(range(cohort_$age))/2, max(cohort_$age))),
              col_fun = col_fun,
              title_position = "topleft",
              title = "Age (months)",
              border = NA,
              direction = 'horizontal')

# 5: patient
patient_ = Legend(at = c("Treated", "Untreated"),
                  legend_gp = gpar(col = c('red','blue')),
                  type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
                  title = "Patient", title_position = "topleft",
                  background = NULL,
                  border = 'black',
                  nrow = 1)

# pathology ______________________________________________________________

# 6: survival
surv_ = Legend(at = c( 'dead','alive'),
               legend_gp = gpar(col = c('black', 'grey75')),
               type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Survival", title_position = "topleft",
               background = NULL,
               nrow = 1)

# 7: cancer stage
col_fun = colorRamp2(breaks = c(min(cohort_$stage, na.rm = T), max(cohort_$stage, na.rm = T)), colors = c("pink","purple"))
stage_ = Legend(at = unique(sort(cohort_$stage)),
                col_fun = col_fun,
                title_position = "topleft",
                title = "Cancer stage",
                direction = 'horizontal')

# 8: recurrence
rec_ = Legend(at = c( 'Yes','No'),
               legend_gp = gpar(col = c('black','yellow')),
               type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Recurrence", title_position = "topleft",
               background = 'black',
               nrow = 1)

# 9: tumor laterality
late_ = Legend(at = c( 'left','right','bilateral'),
               legend_gp = gpar(col = c('#D2B48C', 'skyblue','#4A2F26')),
               type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Laterality", title_position = "topleft",
               background = NULL,
               nrow = 1)

# 10: tumor size
size_ = Legend(at = c( '<=10','>=10'),
               legend_gp = gpar(col = c('orange', 'black')),
               type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Tumor size (cm)", title_position = "topleft",
               background = NULL,
               nrow = 1)

# 11: tumor weight
col_fun = colorRamp2(breaks = c(min(cohort_$weight, na.rm = T), max(cohort_$weight, na.rm = T)), colors = c("lightblue","darkblue"))
wgt_ = Legend(at = c(min(cohort_$weight, na.rm = T), sum(range(cohort_$weight, na.rm = T))/2, max(cohort_$weight, na.rm = T)),
              col_fun = col_fun,
              title_position = "topleft",
              title = "Tumor weight (g)",
              direction = 'horizontal')

# 12: histology
hist_ = Legend(at = c( 'Favorable','Anaplastic'),
               legend_gp = gpar(col = c('green','black')),
               type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
               title = "Histology", title_position = "topleft",
               background = NULL,
               nrow = 1)

# 13: tumor focality
fcl_ = Legend(at = c( 'Unifocal','Multifocal'),
              legend_gp = gpar(col = c('#602CA1','#BF5842')),
              type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
              title = "Tumor focality", title_position = "topleft",
              background = 'black',
              nrow = 1)

# 14: neoadjuvant
neoadj_ = Legend(at = c('EE-4A', 'DD-4A', 'Vincristine+Doxorubicin','I'),
                 legend_gp = gpar(col = c('orange','red','#AA9922','grey35')),
                 type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
                 title = "Neoadjuvant", title_position = "topleft",
                 background = NULL,
                 border = 'black',
                 ncol = 1)

# 15: response to neoadjuvant
rep_ = Legend(at = c('decreased', 'increased', 'unchanged'),
              legend_gp = gpar(col = c('green2','red','grey50')),
              type = "points", pch = 15, size = unit(x = 5, units = 'mm'),
              title = "Response", title_position = "topleft",
              background = NULL,
              border = 'black',
              ncol = 1)

# drawing legends _____________________________________________________________________

# 1st column of legends
lgd_list = packLegend(direction = 'vertical',
                      surv_,
                      stage_,
                      rec_,
                      late_,
                      size_,
                      wgt_,
                      hist_,
                      fcl_,
                      neoadj_,
                      rep_)
draw(lgd_list, x = unit(30, "mm"), y = unit(360, "mm"), just = c('top'))

# 2nd column of legends
lgd_list = packLegend(direction = 'vertical',
                      race_,
                      gender_,
                      age_,
                      patient_)
draw(lgd_list, x = unit(85, "mm"), y = unit(360, "mm"), just = c('top'))

# 3rd column of legends
lgd_list = packLegend(direction = 'vertical',
                      blastema_, stroma_, epi_,
                      fga_,
                      TMB_)
draw(lgd_list, x = unit(135, "mm"), y = unit(360, "mm"), just = c('top'))

# clear the circos plot
circos.clear()

graphics.off()



rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

# Load Libraries                    ####
library(openxlsx)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

## Supplementary Figure 4 ##

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4])

names <- c("Human1",
           "Human1-O1", 
           "Human1-D1", 
           "Human1-T1")

gMCSs.list <- list()
for (model in names){
  gMCSs.list[[model]] <- as.data.frame(read.xlsx('../calculated_gMCSs_of_each_model/gMCSs_list.xlsx', sheet = model))
  gMCSs.list[[model]] <- unique(apply(gMCSs.list[[model]],1,function(y){paste(sort(y[y!=""]),collapse = "--")}))
}

nameVector <- unlist(mapply(function(x,y){ rep(y, length(x)) }, gMCSs.list, names(gMCSs.list)))
resultDF <- cbind.data.frame(unlist(gMCSs.list), nameVector)

colnames(resultDF) <- c("ENSEMBL", "layer")

data <- resultDF %>% 
  pivot_wider(id_cols = ENSEMBL,
              names_from = layer, 
              values_from = layer, 
              values_fn = list(x = length), 
              values_fill = list(x = 0))

layer <- c("Human1", "Human1-O1",
           "Human1-D1", 
           "Human1-T1")

data <- t(data)
colnames(data) <- data[1,]

data <- data[-1,]
data[!is.na(data)] <- TRUE
data[is.na(data)] <- FALSE

data <- as.data.frame(t(data))

plot_final_gmcs <- ComplexUpset::upset(
  data,
  rev(layer), sort_sets = F,  stripes='white', name='Model',
  base_annotations=list(
    'gMCSs Intersection'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
      
    ) + scale_fill_manual(values=c('bars_color'='black'), guide = 'none')
  ),matrix=(
    intersection_matrix(geom=geom_point(shape='circle filled', size=3))
    + scale_color_manual(
      values=c('Human1'=colors_final[1], 'Human1-O1'= colors_final[2], 
               'Human1-D1'=colors_final[3], 
               'Human1-T1'=colors_final[4]),
      guide='none'
    )
  ),queries=list(
    upset_query(set='Human1', fill=colors_final[1]),
    upset_query(set='Human1-O1', fill= colors_final[2]),
    upset_query(set='Human1-D1', fill=colors_final[3]),
    upset_query(set='Human1-T1', fill=colors_final[4])
    
  ),
  width_ratio=0.2, wrap=TRUE,  
  set_sizes=FALSE,
  themes=upset_modify_themes(
    list(
      'gMCSs Intersection'=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")))) 
)

plot_final_gmcs
ggsave(paste0("./plots/SuppFig4.pdf"),
       plot = plot_final_gmcs, device = "pdf", width = 7, height = 5, units = "in", dpi = 300, bg = "white")

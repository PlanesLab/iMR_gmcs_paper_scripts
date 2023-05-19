rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

## Supp Figure 11 ##

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3")

# database <- 'DepMap'

length <- 5
ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults_adaptation/DepMap_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
AllThresholdAnalysis <- do.call(rbind,ResultList)

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")

colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[3],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[3],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[3],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2])

colnames(AllThresholdAnalysis)[1]<- "Cell line"
colnames(AllThresholdAnalysis)[20]<- "Model"

Threshold <- 5
ThresholdVector <- 5

Results_t <- AllThresholdAnalysis %>% filter(thMethod == paste0("gMCS_T", Threshold))

Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                      id.vars = c("Model", "Cell line"),
                                      measure.vars = c("True Positives", "False Positives",
                                                       "False Negatives", "True Negatives",
                                                       "Accuracy", "Sensitivity", "Specificity",
                                                       "Positive Predicted Value", "Matthew's Cor. Coef."))

TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) + theme(axis.text.x = element_text(angle = 45, hjust=1))+ #scale_fill_grey()+
  theme_bw() +  
  scale_fill_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs")

FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) +  #scale_fill_grey()+
  theme_bw() +  
  scale_fill_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) + #scale_fill_grey()+
   theme_bw() +  
  scale_fill_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
  ggtitle("Positive Predictive Value (PPV)")

plot_final_final <- ggarrange(plotlist = list(TP, FP, PPV), ncol = 1, nrow = 3, common.legend = T, legend = "bottom")

plot_final_final

ggsave(paste0("./plots/SuppFig11.pdf"),
       plot = plot_final_final, device = "pdf", width = 8, height = 12, units = "in", dpi = 300, bg = "white")


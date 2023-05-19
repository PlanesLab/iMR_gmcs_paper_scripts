rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

## Supp Figure 5 ##
length <- 5
database <- "Depmap"

names <- c("Human1",
           "Human1-O1", 
           "Human1-D1", 
           "Human1-T1")

ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults/", database, "_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
AllThresholdAnalysis <- do.call(rbind,ResultList)

colnames(AllThresholdAnalysis)[1]<- "Cell line"
colnames(AllThresholdAnalysis)[20]<- "Model"

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")


colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4]
)
Threshold <- 5
Results_t <- AllThresholdAnalysis %>% filter(thMethod == paste0("gMCS_T", Threshold))

# All Statistics #
Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                      id.vars = c("Model", "Cell line"),
                                      measure.vars = c("True Positives", "False Positives",
                                                       "False Negatives", "True Negatives",
                                                       "Accuracy", "Sensitivity", "Specificity",
                                                       "Positive Predicted Value", "Matthew's Cor. Coef."))

TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) +  #scale_fill_grey()+
  theme_bw() +  
  scale_fill_manual(values = colors_final) + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs")

FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) +  #scale_fill_grey()+
  theme_bw() +  
  scale_fill_manual(values = colors_final) + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
  geom_boxplot(aes(x = Model, y = value, fill = Model)) +  #scale_fill_grey()+
  theme_bw() +  
  scale_fill_manual(values = colors_final) + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
  ggtitle("Positive Predictive Value (PPV)")

plot_final_final <- ggarrange(ggarrange(TP, FP, ncol = 2, legend = "none"), ggarrange(NULL, PPV, NULL, widths = c(0.33, 0.66, 0.33),ncol = 3, common.legend = T, legend = "bottom"), nrow = 2)

ggsave(paste0("./plots/SuppFig5.pdf"),
       plot = plot_final_final, device = "pdf", width = 8, height = 8, units = "in", dpi = 300, bg = "white")

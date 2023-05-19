rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Supp Figure 10 ##
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3")

database_array <- c("Hart", "DepMap")
length <- 5

plot_final_final <- list()

for (database in database_array){
  ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults_adaptation/", database, "_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
  AllThresholdAnalysis <- do.call(rbind,ResultList)
  
  AllThresholdAnalysis$adaptation <- "yes"
  
  ResultList_no <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults/", database, "_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
  AllThresholdAnalysis_no <- do.call(rbind,ResultList_no)
  AllThresholdAnalysis_no$adaptation <- "no"
  
  AllThresholdAnalysis <- rbind(AllThresholdAnalysis, AllThresholdAnalysis_no)
  AllThresholdAnalysis$adaptation <- factor(AllThresholdAnalysis$adaptation)
  
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
  
  # All Statistics #
  Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                        id.vars = c("Model", "Cell line", "adaptation"),
                                        measure.vars = c("True Positives", "False Positives",
                                                         
                                                         "Positive Predicted Value"))
  
  TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
    geom_boxplot(aes(x = Model, y = value, fill = Model, alpha = adaptation), position = position_dodge(0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))+ scale_alpha_manual(values = c(1, 0.5), guide = "none")+
    scale_fill_manual(values = colors_final) +
    theme(plot.title = element_text(size = 14, 
                                    hjust = 0.5),
          strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs") #+

  FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
    geom_boxplot(aes(x = Model, y = value, fill = Model, alpha = adaptation), position = position_dodge(0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))+ scale_alpha_manual(values = c(1, 0.5), guide = "none")+
    scale_fill_manual(values = colors_final) + 
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")
  
  PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
    theme_bw() + 
    geom_boxplot(aes(x = Model, y = value, fill = Model, alpha = adaptation), position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+ scale_alpha_manual(values = c(1, 0.5), guide = "none")+
    scale_fill_manual(values = colors_final) +
    theme(plot.title = element_text(size = 14,  hjust = 0.5),
          strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
    ggtitle("Positive Predictive Value (PPV)")
  plot_final_final[[database]] <- annotate_figure(ggarrange(plotlist = list(TP, FP, PPV), ncol = 1, nrow = 3, common.legend = T, legend = F), top = text_grob(database, face = "bold", size = 16))

  
}
ggarrange(plotlist = plot_final_final, ncol = 2, common.legend = T, labels = c("A", "B"))
ggsave(paste0("./plots/SuppFig10.pdf"),
              device = "pdf", width = 12, height = 14, units = "in", dpi = 300, bg = "white")

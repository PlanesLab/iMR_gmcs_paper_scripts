rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Figure 5 ##
library(dplyr)
library(reshape2)
library(ggpubr)

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3")

# database <- 'DepMap'

length <- 5
ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults_adaptation/Hart_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
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

# All Statistics #
Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                      id.vars = c("Model", "Cell line"),
                                      measure.vars = c("True Positives", "False Positives",
                                                       "False Negatives", "True Negatives",
                                                       "Accuracy", "Sensitivity", "Specificity",
                                                       "Positive Predicted Value", "Matthew's Cor. Coef."))

TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +   scale_shape_manual(values=c(16, 17, 15, 11, 7)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs")

FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +   scale_shape_manual(values=c(16, 17, 15, 11, 7))  + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +   scale_shape_manual(values=c(16, 17, 15, 11, 7))  + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
  ggtitle("Positive Predictive Value (PPV)")

plot_final_final <- ggarrange(plotlist = list(TP, FP, PPV), ncol = 1, nrow = 3, common.legend = T, legend = "right")
plot_final_final
ggsave(paste0("./plots/Figure5.pdf"),
       plot = plot_final_final, device = "pdf", width = 8, height = 14, units = "in", dpi = 300, bg = "white")

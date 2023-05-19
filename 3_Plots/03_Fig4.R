## Figure 4 ##

rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

names <- c("Human1",
           "Human1-O1", 
           "Human1-D1", 
           "Human1-T1")

length <- 5
ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults/Hart_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
AllThresholdAnalysis <- do.call(rbind,ResultList)

colnames(AllThresholdAnalysis)[1]<- "Cell line"
colnames(AllThresholdAnalysis)[20]<- "Model"

Threshold <- 5
ThresholdVector <- 5

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4])


Results_t <- AllThresholdAnalysis %>% filter(thMethod == paste0("gMCS_T", Threshold))

# All Statistics #
Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                      id.vars = c("Model", "Cell line"),
                                      measure.vars = c("True Positives", "False Positives",
                                                       "False Negatives", "True Negatives",
                                                       "Accuracy", "Sensitivity", "Specificity",
                                                       "Positive Predicted Value", "Matthew's Cor. Coef."))

# mean <- Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",]%>% group_by(Model)%>%summarise(mean_val=mean(value))

TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +  #
  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs")

FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
  ggtitle("Positive Predictive Value (PPV)")

TP_PPV <- ggplot(Results_t, aes(`True Positives`, `Positive Predicted Value`, color = Model)) +
  geom_point(aes(color=Model, shape = `Cell line`), size = 5) + #scale_color_grey(start = 0, end = 0.85,) +
  theme_bw()  +  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 18, face = "bold")) + ggtitle("TPs vs PPV") + ylab("PPV") + xlab("TPs")

plot_final_final <- ggarrange(plotlist = list(TP, FP, PPV, TP_PPV), ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
plot_final_final


ggsave(paste0("./plots/Figure4.pdf"),
       plot = plot_final_final, device = "pdf", width = 9.5, height = 9, units = "in", dpi = 300, bg = "white")



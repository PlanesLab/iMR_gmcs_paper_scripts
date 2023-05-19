rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Supp Figure 7 ##

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

names <- c("Human1",
           "Human1-O1", 
           "Human1-D1", 
           "Human1-T1", 
           "Human1-D1∩O1", 
           "Human1-D1∩T1", 
           "Human1-O1∩T1", 
           "Human1-D1∩O1∩T1",
           "Human1-D1UO1",
           "Human1-D1UT1",
           "Human1-O1UT1",
           "Human1-D1UO1UT1")

# database <- 'DepMap'

length <- 5
ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("../2_Essentiality_prediction/RDSResults_adaptation/Hart_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
AllThresholdAnalysis <- do.call(rbind,ResultList)

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")

colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2],
                  colorRampPalette(c(colors[7], colors[8], "black"))(6)[2],
                  "magenta4",
                  "grey",
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[7], colors[8], "black"))(6)[4],
                  "maroon4",
                  "black"
)

colnames(AllThresholdAnalysis)[1]<- "Cell line"
colnames(AllThresholdAnalysis)[20]<- "Model"

Threshold <- 5
ThresholdVector <- 5

levels(AllThresholdAnalysis$Model) <- c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", 
                                        "Human1-O∩D1", "Human1-D∩T1",
                                        "Human1-O∩T1","Human1-O∩D∩T1",
                                        
                                        "Human1-OUD1",
                                        "Human1-DUT1", 
                                        
                                        "Human1-OUT1",
                                        "Human1-OUDUT1", "NA"
)

Results_t <- AllThresholdAnalysis %>% filter(thMethod == paste0("gMCS_T", Threshold))

Melted_gMCS_All_Statistics_th <- melt(data = Results_t,
                                      id.vars = c("Model", "Cell line"),
                                      measure.vars = c("True Positives", "False Positives",
                                                       "False Negatives", "True Negatives",
                                                       "Accuracy", "Sensitivity", "Specificity",
                                                       "Positive Predicted Value", "Matthew's Cor. Coef."))
mean <- Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",]%>% group_by(Model)%>%summarise(mean_val=mean(value))

TP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "True Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  scale_color_manual(values = colors_final) +
  theme_bw() + 
  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + 
  scale_x_discrete(labels=c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA"='',
                            "Human1-O∩D1", "Human1-D∩T1",
                            "Human1-O∩T1","Human1-O∩D∩T1",
                            "NA"= '',
                            "Human1-OUD1",
                            "Human1-DUT1", 
                            
                            "Human1-OUT1",
                            "Human1-OUDUT1"),
                   limits =c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA",
                             "Human1-O∩D1", "Human1-D∩T1",
                             "Human1-O∩T1","Human1-O∩D∩T1",
                             "NA",
                             "Human1-OUD1",
                             "Human1-DUT1", 
                             
                             "Human1-OUT1",
                             "Human1-OUDUT1"
                   )) +
  geom_vline(xintercept = 5, size = 1, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 10, size = 1, color = "black", linetype = "dashed") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("True Positives (TPs)") + ylab("TPs")

FP <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "False Positives",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  scale_color_manual(values = colors_final) +
  theme_bw() +  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + 
  scale_x_discrete(labels=c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA"='',
                            "Human1-O∩D1", "Human1-D∩T1",
                            "Human1-O∩T1","Human1-O∩D∩T1",
                            "NA"= '',
                            "Human1-OUD1",
                            "Human1-DUT1", 
                            
                            "Human1-OUT1",
                            "Human1-OUDUT1"),
                   limits =c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA",
                             "Human1-O∩D1", "Human1-D∩T1",
                             "Human1-O∩T1","Human1-O∩D∩T1",
                             "NA",
                             "Human1-OUD1",
                             "Human1-DUT1", 
                             
                             "Human1-OUT1",
                             "Human1-OUDUT1"
                   )) +
  geom_vline(xintercept = 5, size = 1, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 10, size = 1, color = "black", linetype = "dashed") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank()) + ggtitle("False Positives (FPs)") + ylab("FPs")

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "Positive Predicted Value",])+
  geom_point(aes(x = Model, y = value, shape = `Cell line`, color = Model), size = 5) + 
  stat_summary(aes(x = Model, y = value), geom="point",fun="mean",shape="*",color="black",size=10) +
  theme_bw() +  scale_shape_manual(values=c(16, 17, 15, 11, 7)) + 
  scale_x_discrete(labels=c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA"='',
                            "Human1-O∩D1", "Human1-D∩T1",
                            "Human1-O∩T1","Human1-O∩D∩T1",
                            "NA"= '',
                            "Human1-OUD1",
                            "Human1-DUT1", 
                            
                            "Human1-OUT1",
                            "Human1-OUDUT1"),
                   limits =c("Human1", "Human1-O1", "Human1-D1", "Human1-T1", "NA",
                             "Human1-O∩D1", "Human1-D∩T1",
                             "Human1-O∩T1","Human1-O∩D∩T1",
                             "NA",
                             "Human1-OUD1",
                             "Human1-DUT1", 
                             
                             "Human1-OUT1",
                             "Human1-OUDUT1"
                   )) +
  geom_vline(xintercept = 5, size = 1, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 10, size = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = colors_final) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x =element_blank(), axis.title.x = element_blank())+ ylab("PPV")+
  ggtitle("Positive Predictive Value (PPV)")

plot_final_final <- ggarrange(plotlist = list(TP, FP, PPV), ncol = 1, nrow = 3, common.legend = T, legend = "right")

plot_final_final
# ggsave(paste0("plots/SuppFig7.png"),
#        plot = plot_final_final, device = "png", width =9, height = 14, units = "in", dpi = 300, bg = "white")

ggsave(paste0("plots/SuppFig7.pdf"),
       plot = plot_final_final, device = "pdf", width =9, height = 14, units = "in", dpi = 300, bg = "white")

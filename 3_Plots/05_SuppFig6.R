rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Supp Figure 6 ##
library(openxlsx)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(ComplexUpset)

source('../2_Essentiality_prediction/ComputeBinaryMatrix.R')

## Supp Figure 3##  
names <- c("Human1",
           "Human1-O1", 
           "Human1-D1", 
           "Human1-T1")

max_length <- 5
database <- 'Hart'

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")


colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4]
)

Hart_Effect <- read.xlsx("../DATA/Hart2015_BayesFactors.xlsx")[,-c(1,8,9,10)]
Hart_Exp    <- read.table("../DATA/Hart2015_RNAseq.txt", header = TRUE)

rownames(Hart_Exp) <- Hart_Exp[,1]
Hart_Exp <- Hart_Exp[,-1]
toENSEMBL <- clusterProfiler::bitr(Hart_Effect[,1], "SYMBOL", "ENSEMBL", org.Hs.eg.db)
colnames(Hart_Effect)[1] <- "SYMBOL"

Hart_Effect_ENSEMBL <- merge(toENSEMBL, Hart_Effect, by = "SYMBOL")

Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$SYMBOL),]

Effect_Tresholds <-  as.data.frame(rbind(c('DLD1', 3.57), c('GBM', 3.20), c('HCT116', 1.57), c('HELA', 15.47), c('RPE1', 6.84)))
Effect_Tresholds$V2 <- as.numeric(Effect_Tresholds$V2)

Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$ENSEMBL),]
Hart_Effect_SYMBOL  <- Hart_Effect_ENSEMBL[,1]
rownames(Hart_Effect_ENSEMBL) <- Hart_Effect_ENSEMBL[,2]
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[,-c(1,2)]




Threshold <- 5
ThresholdVector <- 5
# pos <- 1
Binary_Matrix_Length_List <- list()
gMCS_order_one <- gMCS_order_two <- list()
for (model in names){
  gMCS.info <- new.env()
  
  load(paste0('../2_Essentiality_prediction/Data_essentiality/gMCS_', model, '_', database, '_length', max_length, '.RData'), envir = gMCS.info)
  gMCS.info <- as.list(gMCS.info)

  HumanGEM_Genes <- as.data.frame(fread(paste0("../2_Essentiality_prediction/Gene_lists/Genes_", model, ".txt")))
  table.genes.HumanGEM <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  Table_HumanGEM_Genes <- table.genes.HumanGEM
  
  
  gMCS.info$table.genes.HumanGEM = table.genes.HumanGEM
  
  RoLResults <- readRDS(paste0("../2_Essentiality_prediction/RDSResults/", database, "_RoLAllThresholds", "_", model, "_length_", max_length, ".RDS"))
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  
  
  Binary_Matrix_Length_List[[model]] <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                        gMCSs.ENSEMBL.list, unlist(HumanGEM_Genes))
  
  
}

sdf <- reshape2::melt(Binary_Matrix_Length_List)
colnames(sdf) <- c("ENSEMBL", "CellLine", "Prediction", "layer")
sdf$Prediction <- factor(as.character(sdf$Prediction), levels = c("1", "0"), labels = c("Essential", "Not essential"))


colnames(Hart_Effect_ENSEMBL) <- c("HCT116", "HELA", "GBM", "RPE1", "DLD1")
Hart_Effect_ENSEMBL2 <- Hart_Effect_ENSEMBL
Hart_Effect_ENSEMBL2 <- Hart_Effect_ENSEMBL2 %>% as.matrix() %>% reshape2::melt() 
colnames(Hart_Effect_ENSEMBL2) <- c("ENSEMBL", "CellLine", "Value")

Effect_Tresholds <-  as.data.frame(rbind(c('DLD1', 3.57), c('GBM', 3.20), c('HCT116', 1.57), c('HELA', 15.47), c('RPE1', 6.84)))
colnames(Effect_Tresholds) <- c("CellLine", "Threshold")

final <- merge(sdf, Hart_Effect_ENSEMBL2, by = c("ENSEMBL", "CellLine"))
final2 <- merge(final, Effect_Tresholds, by = "CellLine")

final2 <- na.omit(final2)
final2$Essential <- as.factor(ifelse(final2$Value >= as.numeric(final2$Threshold), "Essential", "Not essential"))
HumanGEM_Genes <- as.data.frame(fread(paste0("../2_Essentiality_prediction/Gene_lists/Genes_Human1.txt")))
final2$metabolic <- as.factor(ifelse(final2$ENSEMBL %in% HumanGEM_Genes$Var1, "Metabolic", "Not metabolic"))

final_final <- final2[final2$Prediction == "Essential",]

essential <- unique(final_final[, c("Essential", "ENSEMBL", "CellLine")])

BM <- bitr(essential$ENSEMBL, "ENSEMBL", "SYMBOL", org.Hs.eg.db)

essential2 <- merge(essential, BM)

X<-split(final_final, final_final$CellLine)

layer <- c("Human1", "Human1-O1", 
           "Human1-D1",
           "Human1-T1")

plot_list <- list()

for (i in 1:length(X)){
  a <- X[[i]] %>%
    pivot_wider(id_cols = ENSEMBL,
                names_from = layer,
                values_from = layer,
                values_fn = list(x = length),
                values_fill = list(x = 0))
  
  a <- t(a)
  colnames(a) <- a[1,]
  
  a <- a[-1,]
  a[!is.na(a)] <- TRUE
  a[is.na(a)] <- FALSE
  
  b <- as.data.frame(t(a))
  
  b$`Gene type` <- as.factor(ifelse(rownames(b) %in% HumanGEM_Genes$Var1, "Metabolic", "Regulatory"))
  b$ENSEMBL <- rownames(b)
  essential_GBM <- essential[essential$CellLine %in% names(X)[i], c("ENSEMBL", "Essential")]
  b <- as.data.frame(merge(b, essential_GBM, by = "ENSEMBL"))
  rownames(b) <- b$ENSEMBL
  b <- b[,-1]
  
  ####################################
  plot_list[[i]] <- ComplexUpset::upset(
    b,
    rev(layer), sort_sets = F,  stripes='white', name='Model',
    base_annotations=list(
      'Positives Intersection'=intersection_size(
        counts=TRUE, text_colors=c(
          on_background='black', on_bar='black'
        ),
        mapping=aes(fill=`Gene type`, color = "black")
      ) + scale_fill_manual(values=c(
        'Metabolic'= 'white', 'Regulatory'='black'),  guide = 'none'
      )),   matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=3))
        + scale_color_manual(
          values=c('Human1'=colors_final[1], 'Human1 + Omnipath'= colors_final[2], 'Human1 + Dorothea'=colors_final[3],  'Human1 + TRRUST'=colors_final[4]),
          guide='none'
        )
      ),queries=list(
        upset_query(set='Human1', fill=colors_final[1]),
        upset_query(set='Human1-O1', fill= colors_final[2]),
        upset_query(set='Human1-D1', fill=colors_final[3]),
        upset_query(set='Human1-T1', fill=colors_final[4])
      ),
    width_ratio=0.2, height_ratio = 0.25,
    wrap=TRUE, set_sizes = FALSE,
    themes=upset_modify_themes(
      list(
        'Positives Intersection'=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
                                       axis.title.x = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))))
  ) + ggtitle(paste0(names(X)[i], " cell line")) +  theme(plot.title = element_text(size = 16, face = "bold"))
  
  plot_list[[i]]
  ggsave(paste0("./plots/SuppFig6/", names(X)[i], ".pdf"),
         plot = plot_list[[i]], device = "pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "white")
  
}

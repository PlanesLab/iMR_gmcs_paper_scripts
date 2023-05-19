## ################################ ##
##                                  ##
##        Analysis of gMCSs         ##
##                                  ##
## ################################ ##


rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/2_Essentiality_prediction/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

# Load Libraries                    ####
library(data.table)
library(doSNOW)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(openxlsx)
library(org.Hs.eg.db)
library(pheatmap)
library(pROC)
library(reshape2)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)

# Prepare Data                      ####
# Prepare Dummy Sample Class Data #
Sample_Class_Dummy <- data.frame("Sample_Class" = "Cells")
levels(Sample_Class_Dummy)  <- "Cells"

# Get HumanGEM information #

# Store function to compute Essential Genes #
source('fun_CalculateHartEssentialGenes_gmcsTH.R')
source('ComputeBinaryMatrix.R')
source('fun_TweakedHartgMCS_Results.R')
source('fun_TweakedDepMapgMCS_Results.R')


database <- "Hart"
# database <- "DepMap"

length <- 3
# length <- 5

if (database == "Hart"){
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
  
} else if (database == "DepMap"){
  
  load("../DATA/CCLE_DepMap_Data.RData")
  rownames(CCLE_Exp) <- CCLE_Exp[,1]
  CCLE_Exp <- CCLE_Exp[,-1]
  
  rownames(Achilles) <- Achilles[,1]
  Achilles <- Achilles[,-1]
  
  suffix <- sub(".*\\..", "", colnames(CCLE_Exp))
  suffix <- sub("\\.", "", suffix)
  colnames(CCLE_Exp) <- paste0("E", suffix)
  
  suffix <- sub("\\..*.", "", colnames(Achilles))
  colnames(Achilles) <- suffix
  
  Achilles <- as.data.frame(t(Achilles))
  Achilles$ALIAS <- rownames(Achilles)
  
  toENSEMBL <- clusterProfiler::bitr(Achilles[,"ALIAS"], "ALIAS", "ENSEMBL", org.Hs.eg.db)
  
  Hart_Effect_ENSEMBL <- merge(toENSEMBL, Achilles, by = "ALIAS")
  
  Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$ALIAS),]
 
  Hart_Exp <- t(CCLE_Exp)
}


Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$ENSEMBL),]
Hart_Effect_SYMBOL  <- Hart_Effect_ENSEMBL[,1]
rownames(Hart_Effect_ENSEMBL) <- Hart_Effect_ENSEMBL[,2]
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[,-c(1,2)]

if (database == "Hart"){
  colnames(Hart_Effect_ENSEMBL) <- c("HCT116", "HELA", "GBM", "RPE1", "DLD1")
  Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[, c("DLD1", "GBM", "HCT116", "HELA", "RPE1")]
}



names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3",
           "Human1-D1∩O1", "Human1-D1UO1",
           "Human1-D1∩T1", "Human1-D1UT1",
           "Human1-O1∩T1", "Human1-O1UT1",
           "Human1-D1∩O1∩T1", "Human1-D1UO1UT1")

names <- c(#"Human1",
           # "Human1-O1", "Human1-O2", "Human1-O3",
           # "Human1-D1", "Human1-D2", "Human1-D3",
           # "Human1-T1", "Human1-T2", "Human1-T3",
           # "Human1-D1∩O1", "Human1-D1∩T1",
           # "Human1-O1∩T1","Human1-D1∩O1∩T1",
           # "Human1-D1UO1",
           # "Human1-D1UT1", 
           # "Human1-O1UT1",
           # "Human1-D1UO1UT1",
           
           "Human1-D2∩O2", "Human1-D2∩T2",
           "Human1-O2∩T2","Human1-D2∩O2∩T2",
           "Human1-D2UO2",
           "Human1-D2UT2", 
           "Human1-O2UT2",
           "Human1-D2UO2UT2",
           
           "Human1-D3∩O3", "Human1-D3∩T3",
           "Human1-O3∩T3","Human1-D3∩O3∩T3",
           "Human1-D3UO3",
           "Human1-D3UT3", 
           "Human1-O3UT3",
           "Human1-D3UO3UT3")


pos <- 1
for (model in names){
  
  print(paste0(model))
  
  load(paste0('./Data_essentiality/gMCS_', model, '_', database, '_length', length, '.RData'))

  HumanGEM_Genes <- as.data.frame(fread(paste0("Gene_lists/Genes_", model, ".txt")))
  table.genes.HumanGEM <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  Table_HumanGEM_Genes <- table.genes.HumanGEM
  
  # })
  
  gMCS.info <- list(gMCSs.ENSEMBL = gMCSs.ENSEMBL, gMCSs.ENSEMBL.list = gMCSs.ENSEMBL.list, 
                    gMCSs.ENSEMBL.mat = gMCSs.ENSEMBL.mat, table.gMCSs = table.gMCSs, genes.gMCSs.ENSEMBL = genes.gMCSs.ENSEMBL,
                    gMCSs.ENSEMBL.txt = gMCSs.ENSEMBL.txt, gMCSs.ENSEMBL.length = gMCSs.ENSEMBL.length, table.genes.HumanGEM = table.genes.HumanGEM)
  

  Threshold <- 5
  ThresholdVector <- 5
  RoLResults <- vector("list", length(ThresholdVector))
  names(RoLResults) <- paste0("Th",ThresholdVector)
  

    # for(Threshold in ThresholdVector){
      if (model == "Human1"){
        ratio_threshold <- NA
      } else{
        ratio_threshold <- readRDS(paste0("./RDSResults/", database, "_RoLAllThresholds", "_Human1_length_", length, ".RDS"))
        ratio_threshold <- ratio_threshold[[paste0("Th",Threshold)]]$ratio_threshold
      }
      
      RoLResults[[paste0("Th",Threshold)]] <- CalculateEssentialGenes_gmcsTH(gene.exp = Hart_Exp,
                                                                             gMCS_Analysis = gMCS.info,
                                                                             sample.class = Sample_Class_Dummy,
                                                                             gmcsTH_perc = Threshold/100,
                                                                             nWorkers = 8,
                                                                             gMCS_order = length,
                                                                             order1 = TRUE,
                                                                             ratio_threshold = ratio_threshold)
    # }
    saveRDS(RoLResults, paste0("./RDSResults/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
    ResultList <- vector("list", length(ThresholdVector))
    names(ResultList) <- paste0("Th",ThresholdVector)
    
    # for(Threshold in ThresholdVector){
      Binary_Matrix_Length_List <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                       gMCSs.ENSEMBL.list, unlist(HumanGEM_Genes))
      if (database == "Hart"){
        ResultsInList <- TweakedHartgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List,
                                                 Effect = Hart_Effect_ENSEMBL,
                                                 Hart_Thresholds = Effect_Tresholds,
                                                 Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                 gMCS_order =  length,
                                                 Threshold_Info =  paste0("gMCS_T",Threshold))
      } else if(database == "DepMap"){
        ResultsInList <- TweakedDepMapgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List,
                                                   Effect = Hart_Effect_ENSEMBL,
                                                   Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                   gMCS_order = length,
                                                   Threshold_Info =  paste0("gMCS_T",Threshold))
      }
      
      
      # ResultsInList <- do.call(rbind, ResultsInList)
      ResultsInList$Order <- factor(model,
                                    levels=c(model))
      ResultsInList <-  ResultsInList %>% dplyr::rename("True Positives" = TP, "False Positives" = FP,
                                                 "True Negatives" = TN, "False Negatives" = FN,
                                                 "Positive Predicted Value" =  PPV,
                                                 "Sensitivity" = sensitivity, "Specificity" = specificity,
                                                 "Accuracy" = accuracy , "Matthew's Cor. Coef." = MCC)
      
      ResultList[[paste0("Th",Threshold)]] <- ResultsInList
    # }
    saveRDS(ResultList, file = paste0("./RDSResults/", database, "_AllThresholdsResults_", model, "_length_", length, ".RDS"))
}





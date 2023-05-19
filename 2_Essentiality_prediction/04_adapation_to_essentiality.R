## ################################ ##
##                                  ##
##  Analysis of adaptation in gMCSs         ##
##                                  ##
## ################################ ##

rm (list = ls())
gc()
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


source('ComputeBinaryMatrix.R')
source('fun_TweakedHartgMCS_Results.R')
source('fun_TweakedDepMapgMCS_Results.R')

database <- "Hart"
# database <- "DepMap"

# length <- 3
length <- 5

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
  "Human1-D1∩O1", "Human1-D1∩T1",
  "Human1-O1∩T1","Human1-D1∩O1∩T1",
  "Human1-D1UO1",
  "Human1-D1UT1",
  "Human1-O1UT1",
  "Human1-D1UO1UT1",
  
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

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3",
           "Human1-D1∩O1", "Human1-D1∩T1",
           "Human1-O1∩T1","Human1-D1∩O1∩T1",
           "Human1-D1UO1",
           "Human1-D1UT1",
           "Human1-O1UT1",
           "Human1-D1UO1UT1")

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3")

names <- c("Human1-D1∩O1", "Human1-D1∩T1",
           "Human1-O1∩T1","Human1-D1∩O1∩T1",
           "Human1-D1UO1",
           "Human1-D1UT1",
           "Human1-O1UT1",
           "Human1-D1UO1UT1")

names <- c("Human1-T2")

################################### ####
#           All Analysis            ####
################################### ####
# Compute Many Thresholds           ####
# ThresholdVector <- c(1,2,2.5,3.5,5,10,20)
ThresholdVector <- c(5)
# number <- 3

Binary_Matrix_Length_List <- list()
gMCS_order_one <- gMCS_order_two <- list()
for (model in names){
  gMCS.info <- new.env()
  
  load(paste0('./Data_essentiality/gMCS_', model, '_', database, '_length', length, '.RData'),  envir = gMCS.info)
  gMCS.info <- as.list(gMCS.info)
  
  HumanGEM_Genes <- as.data.frame(fread(paste0("Gene_lists/Genes_", model, ".txt")))
  table.genes.HumanGEM <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
  Table_HumanGEM_Genes <- table.genes.HumanGEM
  
  
  gMCS.info$table.genes.HumanGEM = table.genes.HumanGEM
  
  RoLResults <- readRDS(paste0("./RDSResults/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
  ResultList <- vector("list", length(ThresholdVector))
  names(ResultList) <- paste0("Th",ThresholdVector)
  Threshold <- 5
  
  Binary_Matrix_Length_List[[model]] <- ComputeBinaryMatrix(RoLResults[[paste0("Th",Threshold)]],
                                                            gMCSs.ENSEMBL.list, unlist(HumanGEM_Genes))
}


sdf <- reshape2::melt(Binary_Matrix_Length_List)
colnames(sdf) <- c("ENSEMBL", "CellLine", "Prediction", "layer")
sdf$Prediction <- factor(as.character(sdf$Prediction), levels = c("1", "0"), labels = c("Essential", "Not essential"))

Hart_Effect_ENSEMBL2 <- Hart_Effect_ENSEMBL
Hart_Effect_ENSEMBL2 <- Hart_Effect_ENSEMBL2 %>% as.matrix() %>% reshape2::melt() 
colnames(Hart_Effect_ENSEMBL2) <- c("ENSEMBL", "CellLine", "Value")


Effect_Tresholds <-  as.data.frame(rbind(c('DLD1', 3.57), c('GBM', 3.20), c('HCT116', 1.57), c('HELA', 15.47), c('RPE1', 6.84)))
Effect_Tresholds$V2 <- as.numeric(Effect_Tresholds$V2)
colnames(Effect_Tresholds) <- c("CellLine", "Threshold")

final <- merge(sdf, Hart_Effect_ENSEMBL2, by = c("ENSEMBL", "CellLine"))

if (database == "DepMap"){
  final2 <- na.omit(final)
  final2$Essential <- as.factor(ifelse(final2$Value <= as.numeric(-0.6), "Essential", "Not essential"))
  
}else if (database == "Hart"){
  final2 <- merge(final, Effect_Tresholds, by = "CellLine")
  colnames(Effect_Tresholds) <- c("V1", "V2")
  final2 <- na.omit(final2)
  final2$Essential <- as.factor(ifelse(final2$Value >= as.numeric(final2$Threshold), "Essential", "Not essential"))
  
}


HumanGEM_Genes <- as.data.frame(fread(paste0("Gene_lists/Genes_Human1.txt")))

final2$metabolic <- as.factor(ifelse(final2$ENSEMBL %in% HumanGEM_Genes$total_genes, "Metabolic", "Not metabolic"))

BM <- bitr(final2$ENSEMBL, "ENSEMBL", "SYMBOL", org.Hs.eg.db)

final3 <- merge(final2, BM)

ResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("./RDSResults/", database, "_AllThresholdsResults_", model, "_length_", length, ".RDS")))})
AllThresholdAnalysis <- do.call(rbind,ResultList)

AllThresholdAnalysis_2 <- AllThresholdAnalysis

final_list <- list()

for (model in names){
    RoLResults <- readRDS(paste0("./RDSResults/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
    ResultList <- vector("list", length(ThresholdVector))
    names(ResultList) <- paste0("Th",ThresholdVector)
    RoLResults_new <- RoLResults
    
    ResultList <- readRDS(paste0("./RDSResults/",  database, "_AllThresholdsResults_", model, "_length_", length, ".RDS"))
    
    load(paste0('./Data_essentiality/gMCS_', model, '_', database, '_length', length, '.RData'))
    
    if (file.exists(paste0("Results_txt/Results_KO_", model, "_3.txt"))){
      KO_genes <- as.data.frame(fread(paste0("adaptation_txt/Results_KO_", model, "_", length, ".txt"), header = F))
      gMCSs_problem_genes <- fread(paste0("adaptation_txt/Results_gMCSs_", model, "_", length, ".txt"), header = F)
      final_final <- final3[final3$layer %in% model,]
      
      final_final <- final_final[final_final$ENSEMBL %in% KO_genes$V1,]
      
      library(dplyr)
      
      gMCSs.ENSEMBL.4.aux <- apply(gMCSs_problem_genes,1,function(x){paste0(x[x!=""],collapse = "--")})
      pos_pos <- which(gMCSs.ENSEMBL.txt %in% gMCSs.ENSEMBL.4.aux)
      
      KO_genes_final <- KO_genes[which(gMCSs.ENSEMBL.4.aux %in% gMCSs.ENSEMBL.txt),]
      
      final_final <- final_final[final_final$ENSEMBL %in% KO_genes_final,]
      final_final_essential <- final_final[final_final$Prediction == "Essential",]
      
      essential_genes <- unique(final_final[final_final$Prediction == "Essential","ENSEMBL"])
      if (!isEmpty(essential_genes)){
        for (i in 1:length(essential_genes)){
          pos <- which(names(RoLResults$Th5$list.gMCS.essential.mat) %in% essential_genes[i])
          list_pos <- apply(RoLResults$Th5$list.gMCS.essential.mat[[pos]],2, which) 
          
          for (j in 1:length(list_pos)){
            if (!identical(list_pos[[j]], integer(0))){
              pos_pos_final <- list_pos[[j]] %in% pos_pos
              if (sum(list_pos[[j]] %in% pos_pos)==length(list_pos[[j]])){
                
                
                final_final_essential$Prediction[final_final_essential$CellLine %in% names(list_pos)[j] & final_final_essential$ENSEMBL %in% essential_genes[i]] <- "Not essential - NB"
                RoLResults_new$Th5$list.gMCS.essential.mat[[pos]][list_pos[[j]], j] <- FALSE
                if (essential_genes[i] %in% rownames(RoLResults_new$Th5$mat.essential.gene)){
                  pos_new <- which( rownames(RoLResults_new$Th5$mat.essential.gene) %in% essential_genes[i])
                  RoLResults_new$Th5$mat.essential.gene[pos, j] <- 0
                }
              }
            }
          }
          if (length(which(RoLResults_new$Th5$list.gMCS.essential.mat[[pos]]))==0){
            RoLResults_new <- RoLResults_new[-pos]
          }
          
        }
        saveRDS(object = RoLResults_new, paste0("./RDSResults_adaptation/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
        HumanGEM_Genes <- as.data.frame(fread(paste0("Gene_lists/Genes_", model, ".txt")))
        table.genes.HumanGEM <- clusterProfiler::bitr(unlist(HumanGEM_Genes), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
        Table_HumanGEM_Genes <- table.genes.HumanGEM
        
        Binary_Matrix_Length_List_2 <- ComputeBinaryMatrix(RoLResults_new[[paste0("Th",Threshold)]],
                                                           gMCSs.ENSEMBL.list, unlist(HumanGEM_Genes))
        if (database == "Hart"){
          ResultsInList <- TweakedHartgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List_2,
                                                   Effect = Hart_Effect_ENSEMBL,
                                                   Hart_Thresholds = Effect_Tresholds,
                                                   Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                   gMCS_order =  length,
                                                   Threshold_Info =  paste0("gMCS_T",Threshold))
        } else if(database == "DepMap"){
          ResultsInList <- TweakedDepMapgMCS_Results(Binary_Matrix = Binary_Matrix_Length_List_2,
                                                     Effect = Hart_Effect_ENSEMBL,
                                                     Table_HumanGEM_Genes = Table_HumanGEM_Genes,
                                                     gMCS_order = length,
                                                     Threshold_Info =  paste0("gMCS_T",Threshold))
        }
        ResultsInList$Order <- factor(model,
                                      levels=c(model))
        ResultsInList <-  ResultsInList %>% dplyr::rename("True Positives" = TP, "False Positives" = FP,
                                                          "True Negatives" = TN, "False Negatives" = FN,
                                                          "Positive Predicted Value" =  PPV,
                                                          "Sensitivity" = sensitivity, "Specificity" = specificity,
                                                          "Accuracy" = accuracy , "Matthew's Cor. Coef." = MCC)
        
        ResultList[[paste0("Th",Threshold)]] <- ResultsInList
        saveRDS(ResultList, file =paste0("./RDSResults_adaptation/", database, "_AllThresholdsResults_", model, "_length_", length, ".RDS"))
      }
      else{
        saveRDS(object = RoLResults, paste0("./RDSResults_adaptation/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
        saveRDS(ResultList, file =paste0("./RDSResults_adaptation/", database, "_AllThresholdsResults_", model, "_length_", length, ".RDS"))
        
      }
      
      
    }else{
      saveRDS(object = RoLResults, paste0("./RDSResults_adaptation/", database, "_RoLAllThresholds", "_", model, "_length_", length, ".RDS"))
      saveRDS(ResultList, file =paste0("./RDSResults_adaptation/", database, "_AllThresholdsResults_", model, "_length_", length, ".RDS"))
      
    }
}

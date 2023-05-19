#############################################
######                                 ######
######        Analysis of gMCSs        ######
######        and gene expression      ######
######                                 ######
#############################################

rm (list = ls())

library(foreach)
library(doParallel)
library(parallel)
library(data.table)
library(biomaRt)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/2_Essentiality_prediction/")

database <- "Hart"
# database <- "DepMap"

if (database == "Hart"){
  Hart_Effect <- read.xlsx("../DATA/Hart2015_BayesFactors.xlsx")[,-c(1,8,9,10)]
  genes <- bitr(Hart_Effect$Gene, "SYMBOL", "ENSEMBL", OrgDb = org.Hs.eg.db)
  genes_ALIAS <- bitr(Hart_Effect$Gene, "ALIAS", "ENSEMBL", OrgDb = org.Hs.eg.db)
  CCLE_genes <- unique(genes$ENSEMBL, genes_ALIAS$ENSEMBL)
  
} else if (database == "DepMap"){
  load("../DATA/CCLE_DepMap_Data.RData")
  rownames(CCLE_Exp) <- CCLE_Exp[,1]
  CCLE_Exp <- CCLE_Exp[,-1]
  
  suffix <- sub(".*\\..", "", colnames(CCLE_Exp))
  suffix <- sub("\\.", "", suffix)
  colnames(CCLE_Exp) <- paste0("E", suffix)
  
  CCLE_genes <- colnames(CCLE_Exp)
}


length <- 3

names <- c("Human1",
           "Human1-O1", "Human1-O2", "Human1-O3",
           "Human1-D1", "Human1-D2", "Human1-D3",
           "Human1-T1", "Human1-T2", "Human1-T3",
           "Human1-D1∩O1", "Human1-D1UO1",
           "Human1-D1∩T1", "Human1-D1UT1",
           "Human1-O1∩T1", "Human1-O1UT1",
           "Human1-D1∩O1∩T1", "Human1-D1UO1UT1")

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


## Load the gene information
for (model in names){
  gMCSs.ENSEMBL <- as.matrix(read.xlsx('../calculated_gMCSs_of_each_model/gMCSs_list.xlsx', sheet = model))
  gMCSs.ENSEMBL <- unique(gMCSs.ENSEMBL)
  gMCSs.ENSEMBL[is.na(gMCSs.ENSEMBL)] <- ""
  
  gMCSs.ENSEMBL <- list(gMCSs.ENSEMBL)
  names(gMCSs.ENSEMBL) <- 57
  
  lapply(gMCSs.ENSEMBL, dim)
  
  # generate info for Cell Lines with culture media ####
  lengths <- array()
  lengths[length+1] <- nrow(gMCSs.ENSEMBL[[1]])

  gMCSs.ENSEMBL_new <- gMCSs.ENSEMBL
  
  lapply(gMCSs.ENSEMBL, dim)
  
  for (i in length:1){
    lengths[i] <- length(which(gMCSs.ENSEMBL_new[[1]][,i] == ""))
    gMCSs.ENSEMBL_new[[1]] <- gMCSs.ENSEMBL_new[[1]][gMCSs.ENSEMBL_new[[1]][,i] == "",1:i]
  }
  
  a <- apply(gMCSs.ENSEMBL[[1]], 1, function(x){
    x %in% CCLE_genes
  })
  
  a <- colSums(a) 
  pos <- a>0
  
  pos <- array()
  
  for (i in 1:(length)){
    pos[(lengths[i]+1):(lengths[i+1])] <- a[(lengths[i]+1):(lengths[i+1])]>i-1
  }
  
  gMCSs.ENSEMBL[[1]] <- gMCSs.ENSEMBL[[1]][pos,]
  max_ncol <- max(unlist(lapply(gMCSs.ENSEMBL, ncol)))
  
  gMCSs.ENSEMBL.2 <- gMCSs.ENSEMBL
  for (t in names(gMCSs.ENSEMBL)){
    try({
      if (dim(gMCSs.ENSEMBL.2[[t]])[2] < max_ncol){
        # t = 1
        gMCSs.ENSEMBL.2[[t]] <- as.data.frame(gMCSs.ENSEMBL.2[[t]])
        aux <- dim(gMCSs.ENSEMBL.2[[t]])[2]+1
        aux <- aux:max_ncol
        gMCSs.ENSEMBL.2[[t]][,aux] <- ""
      }
    })
  }
  
  gMCSs.ENSEMBL.2 <- lapply(gMCSs.ENSEMBL.2, as.matrix)
  for (t in names(gMCSs.ENSEMBL.2)){
    if (nrow(gMCSs.ENSEMBL.2[[t]])>0){
      rownames(gMCSs.ENSEMBL.2[[t]]) <- paste0(t, '--', 1:dim(gMCSs.ENSEMBL.2[[t]])[1])
    }
  }
  
  gMCSs.ENSEMBL.3 <- do.call(rbind, gMCSs.ENSEMBL.2[unlist(lapply(gMCSs.ENSEMBL.2, function(x){prod(dim(x))}))>0])

  table.gMCSs <- data.frame(task = unlist(lapply(as.list(rownames(gMCSs.ENSEMBL.3)), function(x){strsplit(x, '--')[[1]][1]})), 
                            gMCS = unlist(lapply(as.list(rownames(gMCSs.ENSEMBL.3)), function(x){strsplit(x, '--')[[1]][2]})),
                            idx = NA)
  
  table.gMCSs$task <- as.character(table.gMCSs$task)
  table.gMCSs$gMCS <- as.character(table.gMCSs$gMCS)
  
  gMCSs.ENSEMBL.4 <- unique(gMCSs.ENSEMBL.3)
  idx <- order(apply(gMCSs.ENSEMBL.4=="",1,sum), decreasing = T)
  gMCSs.ENSEMBL.4 <- gMCSs.ENSEMBL.4[idx,]
  head(gMCSs.ENSEMBL.4)
  
  gMCSs.ENSEMBL.3.aux <- apply(gMCSs.ENSEMBL.3,1,function(x){paste0(x[x!=""],collapse = "--")})
  gMCSs.ENSEMBL.4.aux <- apply(gMCSs.ENSEMBL.4,1,function(x){paste0(x[x!=""],collapse = "--")})
  
  table.gMCSs$idx <- match(gMCSs.ENSEMBL.3.aux, gMCSs.ENSEMBL.4.aux)
  
  genes.ENSEMBL.list <- list()
  for( n in names(gMCSs.ENSEMBL)){
    genes.ENSEMBL.list[[n]] <- unique(as.character(as.matrix(gMCSs.ENSEMBL[[n]])))
    genes.ENSEMBL.list[[n]] <- genes.ENSEMBL.list[[n]][genes.ENSEMBL.list[[n]]!=""]
    genes.ENSEMBL.list[[n]] <- genes.ENSEMBL.list[[n]][!is.na(genes.ENSEMBL.list[[n]])]
    genes.ENSEMBL.list[[n]]
  }
  
  genes.ENSEMBL.essential.list <- list()
  for( n in names(gMCSs.ENSEMBL)){
    genes.ENSEMBL.essential.list[[n]] <- gMCSs.ENSEMBL[[n]][apply(gMCSs.ENSEMBL[[n]],1,function(x){sum(x!="")})==1,]
    genes.ENSEMBL.essential.list[[n]] <- unique(as.character(as.matrix(genes.ENSEMBL.essential.list[[n]])))
    genes.ENSEMBL.essential.list[[n]] <- genes.ENSEMBL.essential.list[[n]][genes.ENSEMBL.essential.list[[n]]!=""]
    genes.ENSEMBL.essential.list[[n]] <- genes.ENSEMBL.essential.list[[n]][!is.na(genes.ENSEMBL.essential.list[[n]])]
    genes.ENSEMBL.essential.list[[n]]
  }
  
  genes.ENSEMBL <- unique(unlist(genes.ENSEMBL.list))
  
  gMCSs.ENSEMBL.mat <- matrix(0,nrow = nrow(gMCSs.ENSEMBL.4), ncol = length(genes.ENSEMBL))
  rownames(gMCSs.ENSEMBL.mat) <- 1:nrow(gMCSs.ENSEMBL.4)
  colnames(gMCSs.ENSEMBL.mat) <- genes.ENSEMBL
  
  for (i in 1:nrow(gMCSs.ENSEMBL.4)){
    if (i %in% seq(0,1e6,1e3)) {print(paste(i,nrow(gMCSs.ENSEMBL.4)))}
    gmcs <- gMCSs.ENSEMBL.4[i,]
    gmcs <- gmcs[gmcs!=""]
    gMCSs.ENSEMBL.mat[i,gmcs] <- 1
  }
  
  gMCSs.ENSEMBL.mat <- Matrix(gMCSs.ENSEMBL.mat, sparse = T)
  gMCSs.ENSEMBL.txt <- gMCSs.ENSEMBL.4.aux
  genes.gMCSs.ENSEMBL <- genes.ENSEMBL

  gMCSs.ENSEMBL.length <- apply(gMCSs.ENSEMBL.mat,1,sum)
  
  gMCSs.ENSEMBL.list <- strsplit(gMCSs.ENSEMBL.txt, '--')

  save(gMCSs.ENSEMBL, 
       gMCSs.ENSEMBL.mat, table.gMCSs,
       genes.gMCSs.ENSEMBL, gMCSs.ENSEMBL.txt, 
       gMCSs.ENSEMBL.length, gMCSs.ENSEMBL.list,
       file = paste0('Data_essentiality/gMCS_', model, '_', database, '_length', length, '.rdata'))
  
}







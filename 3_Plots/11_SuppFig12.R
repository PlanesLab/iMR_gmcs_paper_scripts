rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Supp Figure 12 ##
library(pheatmap)

#Heatmap 
RoLResults <- readRDS(paste0("../2_Essentiality_prediction/RDSResults/Hart_RoLAllThresholds", "_Human1_length_5.RDS"))
gene.ratio <- RoLResults$Th5$gene.ratio
gene_exp <- log2(as.matrix(gene.ratio[c("ENSG00000118260", "ENSG00000076555"),]))
rownames(gene_exp) <- c("CREB1", "ACACB" )
ha <- data.frame(class = factor(rep("Hart.responder", 5) , levels = c("Hart.non.responder", "Hart.responder")), row.names = colnames(gene_exp))
colors_annotation <- list(scales::brewer_pal(palette = "Paired"))[[1]](2)
names(colors_annotation) <- levels(ha$class)
pheatmap_plot <- pheatmap(gene_exp,
                           cluster_rows = F, cluster_cols = FALSE, 
                           show_rownames = TRUE, show_colnames = T,  
                           color = colorRampPalette(unlist(strsplit(c("royalblue1_white_firebrick1", "#902EBB_white_#F4831B")[1],'_')))(100), 
                           breaks = seq(-1,+1, length.out = 100),
                           annotation_col = ha,
                           legend = T,
                           annotation_legend = T,
                           annotation_colors = list(class = colors_annotation),
                           main = paste0("Hart"), 
                           scale = "none",
                           gp=grid::gpar(fontfamily = "sans-serif"),
                           silent = T)[[4]]
#Correlation plot
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
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[!duplicated(Hart_Effect_ENSEMBL$ENSEMBL),]
Hart_Effect_SYMBOL  <- Hart_Effect_ENSEMBL[,1]
rownames(Hart_Effect_ENSEMBL) <- Hart_Effect_ENSEMBL[,2]
Hart_Effect_ENSEMBL <- Hart_Effect_ENSEMBL[,-c(1,2)]

cells <- intersect(rownames(CCLE_Exp), colnames(Hart_Effect_ENSEMBL))

#Expression of ACACB (ENSG00000076555) gene and essentiality of CREB1 (ENSG00000118260):
data <- data.frame(exp = CCLE_Exp[cells,"ENSG00000076555"], ess = t(Hart_Effect_ENSEMBL["ENSG00000118260",cells]))
correlation_plot <- ggscatter(data, x = "exp", y = 'ENSG00000118260', 
                              add = "reg.line", conf.int = TRUE, 
                              cor.coef = TRUE, cor.method = "pearson",
                              xlab = "ACACB expression [ log2(TPM+1) ]", ylab = "Essentiality of CREB1 [ Achilles score]")


ggarrange(pheatmap_plot, correlation_plot, labels = c("A", "B"))

ggsave(paste0("./plots/SuppFig12.pdf"),
       device = "pdf", width = 14, height = 6, units = "in", dpi = 300, bg = "white")

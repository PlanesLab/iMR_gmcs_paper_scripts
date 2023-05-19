rm (list = ls())
setwd("C:/Users/nbarrenac/OneDrive - Tecnun/2022_Barrena_signaling/CODE/3_Plots/")
Day <-  format(Sys.Date(), format="%Y-%m-%d")

## Supp Figure 2 ##
library(ggvenn)
library(data.table)

TRRUST <- fread("../DATA/TRRUST.txt")
Omnipath <- fread("../DATA/Omnipath.txt")
Dorothea <- fread("../DATA/Dorothea.txt")

genes_T <- unique(c(TRRUST$source_ENSEMBL, TRRUST$target_ENSEMBL))
genes_D <- unique(c(Dorothea$source_ENSEMBL, Dorothea$target_ENSEMBL))
genes_O <- unique(c(Omnipath$source_ENSEMBL, Omnipath$target_ENSEMBL))

TRRUST_all <- list(gsub(" ", "", paste(TRRUST$source_ENSEMBL, '--', TRRUST$target_ENSEMBL, '--', TRRUST$interaction), fixed = TRUE))
names(TRRUST_all) <- "TRRUST"
Omnipath_all <- list(gsub(" ", "", paste(Omnipath$source_ENSEMBL, '--', TRRUST$target_ENSEMBL, '--', TRRUST$interaction), fixed = TRUE))
names(Omnipath_all) <- "Omnipath"
Dorothea_all <- list(gsub(" ", "", paste(Dorothea$source_ENSEMBL, '--', Dorothea$target_ENSEMBL, '--', Dorothea$interaction), fixed = TRUE))
names(Dorothea_all) <- "Dorothea"

gMCSs.list.2 <- lapply(c(TRRUST_all, Omnipath_all, Dorothea_all), function(x) split(x, lengths(strsplit(x, '--'))))

gMCSs.list.2 <-reshape2::melt(gMCSs.list.2)
gMCSs.list.2 <- split(gMCSs.list.2, gMCSs.list.2$L2)
gMCSs.list.2 <- lapply(gMCSs.list.2, function(x) lapply(split(x, x$L1), function(y) y$value))

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")


colors_final <- c(colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4]
)

# colors <- scales::hue_pal()(5)

b <- ggarrange(plotlist = lapply(as.character(sort(as.numeric(names(gMCSs.list.2)))), function(NAME){
  ggvenn(gMCSs.list.2[[NAME]],c("Omnipath", "Dorothea", "TRRUST"), fill_color = colors_final) + ggtitle(paste0("All interactions"))
}))
b

all_genes <- list(Omnipath = genes_O, Dorothea = genes_D, TRRUST = genes_T)

a <- ggvenn(all_genes,names(all_genes), fill_color = colors_final) + ggtitle(paste0("Genes"))


ggarrange(b,a)

ggsave(filename = "plots/SuppFig2.pdf", width = 10, height = 5, bg = "white", dpi = 300, device = "pdf")
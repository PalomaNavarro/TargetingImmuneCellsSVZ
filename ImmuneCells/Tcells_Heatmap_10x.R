## Make heatmap of T cells from young and old SVZ of combined object

library(Seurat)
library(dplyr)
library(ggthemes)
library(ggplot2)

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code")
svz <- readRDS("data/seurat/HARMONY_dulken_buckley_consolidated_April24.rds")


tc <- subset (svz, Celltype=="T_cell")
table(tc$AGE)
rm(svz)

tc$AGE <- factor(tc$AGE, levels = c("YOUNG", "OLD"))

#Genes for heatmap
genelist <-c("Gzmk","Cxcr6","Pdcd1","Tigit","Lag3","Ctla4","Cd160","Nr4a2","Tox","Tcf7","Eomes","Ifng","Tnf", "Il2", "Gzmb","Cd69","Xcl1","Ccr7","S1pr1", "Lef1","Il7r","Sell")

DefaultAssay(tc) <- "RNA"
data <- as.matrix(GetAssayData(object = tc))
meta <- tc@meta.data

input <- data[rownames(data) %in% genelist,] # SUBSET
input <- input[rowSums(input)>1, ] # 13 10755
input <- data[genelist,]

orig <- c("#188ad0","#C70039")
color <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = -1, option = "magma")

#ordered heatmap
pdf(file = paste0("plots/Figure1/Order_SVZ2_Tcells", Sys.Date(), "_.pdf"), width = 6, height = 7)
heat <- heatmap.2(x = input,
                  Rowv = NULL,
                  Colv = NULL,
                  key = F,
                  dendrogram = "none",
                  labCol = F,
                  ColSideColors = orig[factor(meta$AGE)],
                  trace = "none",
                  col=color
)

dev.off()

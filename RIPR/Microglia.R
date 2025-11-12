
library(Seurat)
library(MAST)

setwd("~/Dropbox/Code/")
svz <- readRDS("data/seurat/RIPR_svz_tcr_seurat.rds")


Microglia <- subset(svz, Celltype == "Microglia")
rm(svz)

Microglia <- FindNeighbors(Microglia, dims = 1:20)
Microglia <- FindClusters(Microglia, resolution = 0.25)
Microglia <- RunUMAP(Microglia, dims = 1:20,  min.dist = 0.7, spread=2)
DimPlot(Microglia,  pt.size =1)
ggsave("plots/Microglia_tsne.pdf", width=7, height=6)

cluster_id <- 6
clusterMarkers6 <- FindMarkers(object = Microglia, ident.1 = cluster_id)

## Remove small clusters of doublets
Microglia <- subset(Microglia, seurat_clusters!=6)
Microglia <- subset(Microglia, seurat_clusters!=5)
DimPlot(Microglia,  pt.size =1)


## Calculate DEGs

    DefaultAssay(Microglia) <- "RNA"
    Idents(Microglia) <- "orig.ident"
    
    Microglia_de <-     FindMarkers(object = Microglia,
                             ident.1 = "OldRIPR",
                             ident.2 = "OldControl",
                             test.use = "MAST",
                             max.cells.per.ident = 1000,
                             min.pct = 0.0,
                             logfc.threshold = log(0),
                             random.seed = 3)
    Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
    
    write.table(Microglia_de, file = "data/DE/Microglia_Old_RIPR_control_DE_April24.txt", sep = "\t", quote = F)

## Top 100 up or downregulated genes from this list were imputed into enrichR


## For plotting percentage of PD1 positive cells

    table(Microglia$orig.ident)
    #OldControl      OldRIPR YoungControl    YoungRIPR 
    #1875         1842         1175         1087 

    pd1 <- subset(Microglia, Pdcd1>0)
    table(pd1$orig.ident)
    #OldControl      OldRIPR YoungControl    YoungRIPR 
    #320          284           41           31 

    pd1 <- subset(svz, Pdcd1>0)
    pd1_cell <- table(pd1$Celltype)
    #Microglia Oligodendrocytes      Neuroblasts      Endothelial Astrocytes_qNSCs       aNSCs_NPCs      Mural_cells      Macrophages          T_cells        Ependymal 
    #680               11                1                7                3                2                0                4              163                0 
    #OPC          Neurons 
    #0                0 
    
    table(svz$Celltype)
    
    #Microglia Oligodendrocytes      Neuroblasts      Endothelial Astrocytes_qNSCs       aNSCs_NPCs      Mural_cells      Macrophages          T_cells        Ependymal 
    #6087             6176             4451             3229             2176             2113             1069              559              417              127 
    #OPC          Neurons 
    #124               62 
    
    

# RIPR SVZ Seurat Processing
# Paloma Navarro


library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(tidyr)
library(ggplot2)
sessionInfo()


# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code/")

t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)
#======================================================================================
# Create Seurat object and basic processing
#======================================================================================

    old.control <- Read10X("data/RIPR/G1/filtered_feature_bc_matrix") # 31053  6625
    old.ripr <- Read10X("data/RIPR/G2/filtered_feature_bc_matrix") # 31053  7656
    young.control <- Read10X("data/RIPR/G3/filtered_feature_bc_matrix") # 31053  7276
    young.ripr <- Read10X("data/RIPR/G4/filtered_feature_bc_matrix") # 31053  6758

    cosvz <- CreateSeuratObject(counts = old.control, project = "OldControl")
    rosvz <- CreateSeuratObject(counts = old.ripr, project = "OldRIPR")
    cysvz <- CreateSeuratObject(counts = young.control, project = "YoungControl")
    rysvz <- CreateSeuratObject(counts = young.ripr, project = "YoungRIPR")

    osvz <- merge(cosvz, y = rosvz, add.cell.ids = c("OC", "OR")) # 14281 cells
    ysvz <- merge(cysvz, y = rysvz, add.cell.ids = c("YC", "YR")) # 14034 cells

    svz <- merge(cosvz, y = c(rosvz,cysvz,rysvz), add.cell.ids = c("OC", "OR","YC", "YR")) # 28315 cells
#======================================================================================
# Filter QC and normalize
#======================================================================================
    svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")
    VlnPlot(svz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size= .025)
    
    svz <- subset(svz, subset = nFeature_RNA > 500 & percent.mt < 10) # 26590 cells left
    svz <- NormalizeData(svz)
    svz <- FindVariableFeatures(svz)
    svz <- ScaleData(svz)

#======================================================================================
# Scale Normalize based on Negative Binomial Model and add cell cycle states
#======================================================================================
    
    svz <- SCTransform(svz)
    
    svz <- RunPCA(svz, verbose = FALSE)
    
    svz <- FindNeighbors(svz, dims = 1:20)
    svz <- FindClusters(svz, resolution = 0.15)
    svz <- RunUMAP(svz, dims = 1:20,  min.dist = 0.4, spread=1, seed=42)
    DimPlot(svz, reduction = "umap", pt.size = .9, cols = t20)
    
    #adds state of cell cycle to metadata
    s.genes <- stringr::str_to_title(cc.genes$s.genes)
    g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)
    svz <- CellCycleScoring(svz, s.features = s.genes, g2m.features = g2m.genes)

#======================================================================================
# Add TCR Data
#======================================================================================
    
    tcr_folder <- "data/RIPR/T1/outs/"
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
    tcr <- tcr[!duplicated(tcr$barcode), ]
    tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","j_gene", "cdr3")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "frequency")])
    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3:8)]
    rownames(tcr) <- paste0("OC-", tcr[,1])
    tcr[,1] <- NULL
    con_old_tcr <- tcr
    
    #do the same for other samples samples
    tcr_folder <- "data/RIPR/T2/outs/"
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
    tcr <- tcr[!duplicated(tcr$barcode), ]
    tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","j_gene", "cdr3")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "frequency")])
    tcr <- tcr[, c(2,1,3:8)]
    rownames(tcr) <- paste0("OR-", tcr[,1])
    tcr[,1] <- NULL
    ripr_old_tcr <- tcr
    
    tcr_folder <- "data/RIPR/T3/outs/"
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
    tcr <- tcr[!duplicated(tcr$barcode), ]
    tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","j_gene", "cdr3")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "frequency")])
    tcr <- tcr[, c(2,1,3:8)]
    rownames(tcr) <- paste0("YC-", tcr[,1])
    tcr[,1] <- NULL
    con_young_tcr <- tcr
    
    tcr_folder <- "data/RIPR/T4/outs/"
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
    tcr <- tcr[!duplicated(tcr$barcode), ]
    tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","j_gene", "cdr3")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "frequency")])
    tcr <- tcr[, c(2,1,3:8)]
    rownames(tcr) <- paste0("YR-", tcr[,1])
    tcr[,1] <- NULL
    ripr_young_tcr <- tcr

    
    # Add to the Seurat object
    all_tcrs <- rbind(con_old_tcr, ripr_old_tcr,con_young_tcr,ripr_young_tcr)
    rownames(all_tcrs) <- gsub("-", "_", rownames(all_tcrs))
    rownames(all_tcrs) <- gsub("_1", "-1", rownames(all_tcrs))
    
    svz <- AddMetaData(object=svz, metadata=all_tcrs)
    
    # ISSUE: Resolve clonotype names
    cdr3s <- svz@meta.data$cdr3s_aa
    cdr3s_set <- unique(cdr3s)
    cdr3s_df <- as.data.frame(cdr3s_set)
    cdr3s_df <- tibble::rownames_to_column(cdr3s_df, "ID")
    cdr3s <- as.character(cdr3s)
    x <- plyr::mapvalues(cdr3s, cdr3s_df$cdr3s_set,cdr3s_df$ID)
    x[x==1] <- NA
    c <- paste0("clone_", x)
    c[c=="clone_NA"] <- "NA"
    svz$CloneID <- c
    

#======================================================================================
    # Assign celltypes
#======================================================================================  
    
    svz <- FindNeighbors(svz, dims = 1:20)
    svz <- FindClusters(svz, resolution = 0.1)
    svz <- RunUMAP(svz, dims = 1:20,  min.dist = 0.4, spread=1, seed=42)
    DimPlot(svz, reduction = "umap", pt.size = .9, cols = t20)
    

    
    #Find markers
    clusterMarkers <- FindAllMarkers(object=svz)
    saveRDS(clusterMarkers, file = paste0("data/RIPR_svz_tcr_markers_", Sys.Date(), ".rds"))
    
    #Heatmap
    top10 <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    DoHeatmap(object = svz, features = top10$gene) + NoLegend()
    ggsave("plots/RIPR_markers.pdf", width=20, height=20)
    
    cluster_id <- 13
    clusterMarkers13 <- FindMarkers(object = svz, ident.1 = cluster_id)
    
    #assing celltypes
    new.cluster.ids <- c("Microglia",
                         "Oligodendrocytes",
                         "Neuroblasts",
                         "Endothelial",
                         "Astrocytes_qNSCs",
                         "aNSCs_NPCs",
                         "Mural_cells",
                         "Oligodendrocytes",
                         "Macrophages",
                         "T_cells",
                         "Ependymal",
                         "OPC",
                         "Neurons", 
                         "Microglia"
                         )
    
    names(new.cluster.ids) <- levels(svz)
    svz <- RenameIdents(svz, new.cluster.ids)
    svz[["Celltype"]] <- Idents(svz)
    unique(svz@meta.data$Celltype)
    
    DimPlot(svz, group.by = "Celltype",label = TRUE, label.size = 4,pt.size = 0.5 ) + scale_color_tableau(palette = "Tableau 20")
    
#====================================================================
#Downsample to lowest number of cells per treatment condition
#====================================================================
    
    table(svz$orig.ident)
    #OldControl      OldRIPR YoungControl    YoungRIPR 
   # 6236         7181         6738         6435 
    
    Idents(svz) <- 'orig.ident'
    downsampled <- subset(svz, downsample = 6236 )
    
    old <- subset(downsampled, orig.ident== "OldControl" |orig.ident== "OldRIPR")

    orig <- c('#343a40','#9d4edd')
    
    DimPlot(old,  pt.size =1.5, group.by = "orig.ident", cols=alpha(orig, 1)) +
        theme(plot.title = element_blank(), legend.position = "none")
    ggsave("plots/umap_downsample_condition_old_Fig2d.pdf", width=10, height=10)

    
    CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
               "OPC", "Oligodendrocytes", "Ependymal","Endothelial", "Mural_cells",
               "Microglia", "Macrophages", "T_cells")
    
    old$Celltype <- factor(old$Celltype,  levels=CELLS, ordered=T)
    
    # Adjust colors
    new_colors <- c("#499894","#59A14F","#8CD17D", "#004643", 
                             "#4E79A7" ,"#A0CBE8", "#ffe45e","#FFBE7D","#B6992D",
                             "#b56576","#FF9D9A","#E15759" )
    
    DimPlot(old,  pt.size =1.5, group.by = "Celltype", cols=new_colors) +
        theme(plot.title = element_blank())
    
    ggsave("plots/umap_downsample_celltype_Fig2c.pdf", width=12, height=10)
    
    
    ### get markers for cell types
    unique(svz$Celltype)
    svz[["Celltype"]] <- Idents(svz)
    clusterMarkers <- FindAllMarkers(object=svz)
    saveRDS(clusterMarkers, file = paste0("data/svz_RIPR_markers_NEW_", Sys.Date(), ".rds"))
    
    
    top5 <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    plot <- DoHeatmap(object = old, features = top5$gene, group.by= "Celltype", angle = 45,size = 3, group.colors = new_colors)
    ggplot2::ggsave(filename = "plots/heatmap_old_RIPR_celltypes_5genes_legend.pdf", plot = plot, width = 24, height = 12)
    
    
    

#====================================================================    
### Save final Labeled seurat object and plotting data
#====================================================================
    
    saveRDS(svz, paste0("data/RIPR_svz_tcr_seurat.rds"))
    

    meta <- svz[[]]
    pca <- Embeddings(object = svz, reduction = "pca")[, 1:20]
    umap <- Embeddings(object = svz, reduction = "umap")
    rna <- t(as.matrix(GetAssayData(object = svz[["RNA"]], slot = "data")))
    d <- as.data.frame(cbind(meta, umap, pca, rna))
    saveRDS(d, file = paste0("data/RIPR_svzPlotDf.rda"))

#====================================================================
## Cell proportions
#====================================================================

library(tidyverse)
library(Seurat)
library(ggthemes)
library(ggpubr)

svz <- readRDS("data/RIPR_svz_tcr_seurat.rds")

DefaultAssay(svz) <- "SCT"
meta <- svz[[]]

# Forget obj to free memory
rm(svz)

# Convert to tidy format
meta <- meta %>% rownames_to_column("Cell") %>% tibble::as_tibble()


# Modify the group_by arguments to treatment, celltype columns
Celltype.Freq <- meta %>% 
    group_by(orig.ident, Celltype) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    mutate(percent = paste0(round(100 * n/sum(n), 2), "%"))


# calculate change in proportion as log2 RIPR/control in frequency for each celltype
Celltype.Freq <- Celltype.Freq %>%
    filter(orig.ident %in% c("OldControl", "OldRIPR")) %>%
    group_by(Celltype) %>%
    mutate(proportion = log2(freq[orig.ident == "OldRIPR"] / freq[orig.ident == "OldControl"])) %>%
    ungroup()


CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts","Neurons",
            "Oligodendrocytes", "OPC","Ependymal","Endothelial", "Mural_cells",
           "Microglia", "Macrophages","T_cells")


Celltype.Freq$Celltype <- factor(Celltype.Freq$Celltype,  levels=CELLS, ordered=T)

Celltype.Freq <- subset(Celltype.Freq, orig.ident=="OldControl")

#colors to match tsne
new_colors <- c("#499894","#59A14F","#8CD17D", "#86BCB6", 
                         "#4E79A7" ,"#A0CBE8", "#ffe45e","#FFBE7D","#B6992D",
                         "#b56576","#FF9D9A","#E15759")

p3 <- ggplot(Celltype.Freq, aes(x = Celltype, y = proportion)) +
    geom_col(aes(x = Celltype, y = proportion ,fill = Celltype), show.legend = FALSE )+
    theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
    ggtitle("") +
    theme(axis.title.y = element_text(size = 14, face = "plain")) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(axis.title.y = element_text(size = 14, face = "plain", margin = margin(r = 20))) +
    theme(axis.title.x = element_blank()) +
    theme(plot.title = element_text(size=16, face = "plain")) +
    labs(y = "Log2 Enrichment (RIPR/Control)") +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values = new_colors) +
    #scale_color_manual(values = new_colors) +
    ggtitle("Changes in cell number with RIPR in old mice") +
    geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")+
    theme(axis.line = element_line(colour = "black"))

p3

ggsave("plots/Proportions_RIPR_control.pdf", p3, width = 8, height = 4)


    

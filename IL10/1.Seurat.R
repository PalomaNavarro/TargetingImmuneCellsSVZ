# GEX IL-10 concentration SVZ Seurat Processing
# Paloma Navarro

rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
library(tidyr)
sessionInfo()

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================

    Control <- Read10X("data/IL10/C1/filtered_feature_bc_matrix_C1") # 31053 
    WT1 <- Read10X("data/IL10/C2/filtered_feature_bc_matrix_C2") # 31053  
    Mut1 <- Read10X("data/IL10/C3/filtered_feature_bc_matrix_C3") # 31053  
    WT5 <- Read10X("data/IL10/C4/filtered_feature_bc_matrix_C4") # 31053  
    Mut5 <- Read10X("data/IL10/C5/filtered_feature_bc_matrix_C5") # 31053  
    WT50 <- Read10X("data/IL10/C6/filtered_feature_bc_matrix_C6") # 31053  
    Mut50 <- Read10X("data/IL10/C7/filtered_feature_bc_matrix_C7") # 31053  
    Young <- Read10X("data/IL10/C8/filtered_feature_bc_matrix_C8") # 31053  
    
    Control <- CreateSeuratObject(counts = Control, project = "Control")
    WT1 <- CreateSeuratObject(counts = WT1, project = "WT1")
    Mut1 <- CreateSeuratObject(counts = Mut1, project = "Mut1")
    WT5 <- CreateSeuratObject(counts = WT5, project = "WT5")
    Mut5 <- CreateSeuratObject(counts = Mut5, project = "Mut5")
    WT50 <- CreateSeuratObject(counts = WT50, project = "WT50")
    Mut50 <- CreateSeuratObject(counts = Mut50, project = "Mut50")
    Young <- CreateSeuratObject(counts = Young, project = "Young")

    svz <- merge(Control, y = c(WT1,Mut1,WT5,Mut5,WT50,Mut50,Young), add.cell.ids = c("CO", "WT1", "MU1", "WT5", "MU5", "WT50", "MU50", "Y")) 
    
    rm(Control, WT1,Mut1,WT5,Mut5,WT50,Mut50,Young)
    
#======================================================================================
# QC Metrics, normalize and SCT transform
#======================================================================================
    svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")
    svz <- subset(svz, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 
    
    svz <- NormalizeData(svz)
    svz <- FindVariableFeatures(svz)
    svz <- ScaleData(svz)
    
    svz <- SCTransform(svz)
    svz <- RunPCA(svz, verbose = FALSE)
    ElbowPlot(svz, ndims = 50)
 
    svz <- FindNeighbors(svz, dims = 1:20)
    svz <- FindClusters(svz, resolution = 0.15)
    svz <- RunUMAP(svz, dims = 1:20,  min.dist = 0.8)
    DimPlot( svz, reduction = "umap")

#=======================================================================================
# Run Seurat's Cell Cycle Analysis
#=======================================================================================
    
    #adds state of cell cycle to metadata
    s.genes <- stringr::str_to_title(cc.genes$s.genes)
    g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)
    svz <- CellCycleScoring(svz, s.features = s.genes, g2m.features = g2m.genes)
    
#======================================================================================
# Set parameters for UMAP and Visualize Seurat Style
#======================================================================================
    t20 <- tableau_color_pal(palette = "Tableau 20")(20)
    
    svz <- FindNeighbors(svz, dims = 1:20)
    svz <- FindClusters(svz, resolution = 0.1)
    svz <- RunUMAP(svz, dims = 1:20,  min.dist = 0.3, spread = 1.5, seed.use = 42)
    
    DimPlot(svz,reduction = "umap", cols = t20, pt.size = .3)
    ggsave("plots/UMAP.pdf", width = 6, height = 5)
    
  
    ## Find markers for each cluster
    table(svz[[c("orig.ident", "seurat_clusters")]])
    
    #Find markers
    clusterMarkers <- FindAllMarkers(object=svz)
    saveRDS(clusterMarkers, file = paste0("data/svz_markers_", Sys.Date(), ".rds"))
    
    #Heatmap
    top10 <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    png("plots/marker_heatmap.png", width = 18, height = 20, units = "in", res = 300)
    DoHeatmap(object = svz, features = top10$gene) + NoLegend()
    dev.off()
    
    # Assign celltypes

    #label clusters with first set (11 clusters)  
    new.cluster.ids <- c("Microglia",
                         "Endothelial",
                         "Oligodendrocytes",
                         "Neuroblasts",
                         "Astrocytes_qNSCs",
                         "Mural_cells",
                         "aNSCs_NPCs",
                         "Oligodendrocytes",
                         "Ependymal/Mix",
                         "Macrophages",
                         "T_cells",
                         "Neurons",
                         "Remove",
                         "OPC"
    )
    
    names(new.cluster.ids) <- levels(svz)
    svz <- RenameIdents(svz, new.cluster.ids)
    svz[["Celltype.low"]] <- Idents(svz)
    unique(svz@meta.data$Celltype.low)
    
    ## Remove doublet cluster and save
    svz2<-subset(svz,Celltype.low!="Remove" )
    
    
   
#======================================================================================
## simplify clusters
#=====================================================================================   
    svz <- FindNeighbors(svz, dims = 1:20)
    svz <- FindClusters(svz, resolution = 0.2)
    svz <- RunUMAP(svz, dims = 1:20,  min.dist = 0.3, spread = 1.5, seed.use = 42)
    
    DimPlot(svz,reduction = "umap", cols = t20, pt.size = .3)
    
    #rename new clusters
    
    new.cluster.ids <- c("Microglia",
                         "Endothelial",
                         "Oligodendrocytes",
                         "Neuroblasts",
                         "Astrocytes_qNSCs",
                         "Mural_cells",
                         "aNSCs_NPCs",
                         "Oligodendrocytes",
                         "Remove",
                         "Macrophages",
                         "T_cells",
                         "Oligodendrocytes",
                         "Ependymal",
                         "Neurons",
                         "Remove",
                         "OPC",
                         "Microglia",
                         "Remove",
                         "Remove"
    )
    
    names(new.cluster.ids) <- levels(svz)
    svz <- RenameIdents(svz, new.cluster.ids)
    svz[["Celltype"]] <- Idents(svz)
    unique(svz@meta.data$Celltype)
    
    svz <- subset (svz, Celltype != "Remove")
    
    tplot= DimPlot(svz,  pt.size =.3, group.by = "Celltype",label = TRUE, label.size = 4, cols = t20, combine = TRUE) & NoLegend()
    tplot[[1]]$layers[[1]]$aes_params$alpha = .3
    tplot
    
    ### get markers for cell types
    
    unique(svz$Celltype)
    svz[["Celltype"]] <- Idents(svz)
    
    
    ## Get top genes for each cluster
    clusterMarkers <- FindAllMarkers(object=svz)
    saveRDS(clusterMarkers, file = paste0("data/svz_markers_NEW", Sys.Date(), ".rds"))
    
    svz2 <- subset(svz, orig.ident=="Control"|orig.ident=="WT5"|orig.ident=="Mut5")
    
    top5 <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    plot <- DoHeatmap(object = svz2, features = top5$gene, angle = 45,size = 3) + NoLegend()
    ggplot2::ggsave(filename = "plots/heatmap_control_Wt_Mut_celltypes_5genes.pdf", plot = plot, width = 20, height = 12)
    
    
    #Downsample to lowest number of cells per treatment condition
    
    il10 <- subset(svz,orig.ident=="Control"| orig.ident=="WT5"|orig.ident=="Mut5")
    table(il10$orig.ident)
    #Control    Mut5     WT5 
    #3365    6162    7760
    
    Idents(il10) <- 'orig.ident'
    downsampled <- subset(il10, downsample = 3365)
    
    table(downsampled$orig.ident)
    
    orig <- c('#adb5bd','#00afb9','#023047')
    tplot= DimPlot(downsampled,  pt.size =0.3, group.by = "orig.ident", cols= orig, order=T)
    tplot[[1]]$layers[[1]]$aes_params$alpha = .7
    tplot
    ggsave("plots/umap_downsample_WT5_Mut5_2.pdf",width=5.8, height=5)


    CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts", "Neurons",
               "OPC", "Oligodendrocytes", "Ependymal","Endothelial", "Mural_cells",
               "Microglia", "Macrophages", "T_cells")
    
    downsampled$Celltype <- factor(downsampled$Celltype,  levels=CELLS, ordered=T)
    
    # Adjust colors
    new_colors <- c("#499894","#59A14F","#8CD17D", "#86BCB6", 
                             "#4E79A7" ,"#A0CBE8", "#ffe45e","#FFBE7D","#B6992D",
                             "#b56576","#FF9D9A","#E15759" )
                             
    
    DimPlot(downsampled,  pt.size =0.7, group.by = "Celltype", cols=new_colors) +
        theme(plot.title = element_blank())
    
    tplot= DimPlot(downsampled,  pt.size =1, group.by = "Celltype", cols = new_colors, combine = TRUE)
    tplot[[1]]$layers[[1]]$aes_params$alpha = .5
    tplot
    
    ggsave("plots/umap_downsample_celltype_paper_il10.pdf", width=7, height=5)
    
    
#======================================================================================
   
     # Save Labelled Seurat Object 
    saveRDS(svz, file = paste0("data/IL10_svz_seurat.rds"))  
    
    # Save plotting data
    meta <- svz[[]]
    pca <- Embeddings(object = svz, reduction = "pca")[, 1:20]
    umap <- Embeddings(object = svz, reduction = "umap")
    rna <- t(as.matrix(GetAssayData(object = svz[["RNA"]], slot = "data")))
    d <- as.data.frame(cbind(meta, umap, pca, rna))
    
    saveRDS(d, file = paste0("data/IL10_svzPlotDf.rda"))
    
    
    
#======================================================================================
    # Gene Expression Vertical Violin Plots by cell type > Figure S6E
    
    setwd("~/Dropbox/Code")
    d <- readRDS("data/IL10_svzPlotDf.rda")
    
    d <- subset(d, orig.ident != "WT50")
    d <- subset(d, orig.ident != "Mut50")
    
    # Randomize rows to mix up plotting order
    d <- d[sample(nrow(d)),]
    
    # Adjust factors
    CELLS <- c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts","Neurons",
               "Oligodendrocytes", "OPC","Ependymal","Endothelial", "Mural_cells",
               "Microglia","Macrophages","T_cells")
    d$celltype_factor <- factor(d$Celltype,  levels=CELLS, ordered=T)
    
    d$orig_factor <- factor(d$orig.ident,  levels=c("Young","Control","WT1","WT5","Mut1", "Mut5"), ordered=T)
    orig<- c('#4cc9f0','#f84b35','#144552','#065a60','#00afb9', '#90e0ef')
    
    #######################
    # Individual PLOTTING #
    #######################
    
    gene <- "Ifnb1"
    p_ifna <- ggplot(d, aes(x=celltype_factor, y= Ifnb1))
    p_ifna <- p_ifna + geom_violin(trim=T,scale="width",alpha=0.8,width=0.5)
    p_ifna <- p_ifna + geom_jitter(aes(color=orig_factor), size = 0.7, position = position_jitter(width=.3), alpha=0.8)
    p_ifna <- p_ifna + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
    p_ifna <- p_ifna + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    p_ifna <- p_ifna + theme(axis.title.y=element_blank())
    p_ifna <- p_ifna + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = 5)))
    p_ifna <- p_ifna + labs(title=gene)
    p_ifna <- p_ifna + theme(legend.position="none")
    p_ifna <- p_ifna + scale_color_manual(values = orig)
    p_ifna <- p_ifna + theme(plot.margin=unit(c(0.5,0.5,0.5,2),"cm"))
    p_ifna
    
    gene <- "Ifng"
    p_ifng <- ggplot(d, aes(x=`celltype_factor`, y= Ifng))
    p_ifng <- p_ifng + geom_violin(trim=T,scale="width",alpha=0.8,width=0.5)
    p_ifng <- p_ifng + geom_jitter(aes(color=orig_factor), size = 0.7, position = position_jitter(width=.3), alpha=0.8)
    p_ifng <- p_ifng + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
    p_ifng <- p_ifng + theme(axis.text.x = element_text(angle=45, hjust=1, size = 12), axis.title.x = element_blank())
    p_ifng <- p_ifng + theme(axis.title.y=element_blank())
    p_ifng <- p_ifng + theme(plot.title = element_text(hjust = 0.5,size=15, margin = margin(b = 5)))
    p_ifng <- p_ifng + labs(title=gene)
    p_ifng <- p_ifng + theme(legend.position="none")
    p_ifng <- p_ifng + scale_color_manual(values = orig)
    p_ifng <- p_ifng + theme(plot.margin=unit(c(0.5,0.5,0.5,2),"cm"))
    p_ifng
    
    pq <- plot_grid(p_ifna, p_ifng, ncol = 1, rel_heights = c(3,5), align = "v")
    pq
    
    ggsave("plots/IL10_conc/IFN_violins.pdf", height = 5, width = 6)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        

    
    
    
    
    

library(tidyverse)
library(reshape2)
library(viridis)
library(Seurat)
library(gplots)
library(corrplot)
library(gtools)
library(patchwork)
library(dplyr)
library(cowplot)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(ggplot2)

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code/")
data <- readRDS("data/Brain_lung_Tcells_seurat.rds")

t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)
orig <- c('#702670','#0e9594')

    #Initial clustering
    DimPlot(data, reduction = "umap", group.by = "orig.ident", cols= orig, pt.size = 1)
    
    #Find markers for initial clusters
    clusterMarkersAll <- FindAllMarkers(object=data)
    top10 <- clusterMarkersAll %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    png("plots/marker_heatmap.png", width = 18, height = 20, units = "in", res = 300)
    DoHeatmap(object = data, features = top10$gene) + NoLegend()
    dev.off()
    
    #Label clusters
    new.cluster.ids <- c("Central_memory_Tcells",
                         "Exhausted_Tcells",
                         "B_cells",
                         "Effector_Tcells",
                         "Alveolar_Macrophages",
                         "Tcells",
                         "Microglia",
                         "Monocytes",
                         "Neutrophils",
                         "B_cells2",
                         "Endothelial")
    
    names(new.cluster.ids) <- levels(data)
    data <- RenameIdents(data, new.cluster.ids)
    data[["Celltype"]] <- Idents(data)
    unique(data@meta.data$Celltype)
    
    DimPlot(data, reduction = "umap", label = TRUE, label.size = 5, pt.size = .7, cols = t20)
    ggsave("plots/all_cells/umap.celltype.pdf", width=13, height=10)  
    
    
    ### Subset T cells only 
    tc <- subset(data, Celltype == "Tcells"|Celltype == "Central_memory_Tcells" |Celltype == "Exhausted_Tcells"| Celltype == "Effector_Tcells")
    tc <- FindNeighbors(tc, dims = 1:20)
    tc <- FindClusters(tc, resolution = 0.4)
    tc <- RunUMAP(tc, dims = 1:20,  min.dist = 0.3, spread = 2, seed=40)
    DimPlot(tc,cols= t20,  pt.size =1)
    
    
    clusterMarkersAll <- FindAllMarkers(object=tc)
    top10 <- clusterMarkersAll %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    png("plots/marker_heatmap_8clusters.png", width = 18, height = 20, units = "in", res = 300)
    DoHeatmap(object = tc, features = top10$gene) + NoLegend()
    dev.off()
    write.csv(clusterMarkersAll, "clusterMarkers_Tcells.csv", row.names = FALSE)
    
    
    ## Exhaustion markers
    col.red <- CustomPalette(low = "#dcdcdd", high = "#720026", mid ="#ce4251", k = 50)
    ex <- c("Pdcd1","Tigit", "Lag3","Havcr2","Tox", "Nr4a2","Entpd1",  "Gzmb", "Eomes","Tcf7", "Slamf6", "Il7r")
    
    FeaturePlot(tc, features= ex, order=T,cols=alpha(col.red,1),pt.size=0.3)& NoLegend() & NoAxes()
    ggsave("plots/T_cell_exhaustion_markers_all_June24.pdf", width=8, height=8)
    
    
    new.cluster.ids <- c("Exhausted-like_Tcells",
                         "Central_Memory_Tcells",
                         "Central_Memory_Tcells",
                         "Central_Memory_Tcells",
                         "Early_effector_TR_Tcells",
                         "Cytotoxic_Effector_Tcells",
                         "Central_Memory_Tcells",
                         "Exhausted_like_IFN_responsive",
                         "Tissue_resident_effector_memory"
                         )

    
    
    names(new.cluster.ids) <- levels(tc)
    tc <- RenameIdents(tc, new.cluster.ids)
    tc[["Celltype.Final"]] <- Idents(tc)
    unique(tc@meta.data$Celltype.Final)
    
    
    DimPlot(tc,  pt.size =0.3, group.by = "Celltype.Final", cols= t20, order=F)
    
    

##====================================================================================================================
    
    ## Subset effector and exhausted t cells
    
    tc <- subset(tc, Celltype.Final!="Central_Memory_Tcells")
    tc <-subset(tc, Celltype.Final!="Cytotoxic_Effector_Tcells")
    
    tc <- FindNeighbors(tc, dims = 1:20)
    tc <- FindClusters(tc, resolution = 0.3)
    tc <- RunUMAP(tc, dims = 1:20,  min.dist = 0.1, spread = 2, seed=40)
    DimPlot(tc,cols= t20,  pt.size =1)
    
    col.red <- CustomPalette(low = "#dcdcdd", high = "#720026", mid ="#ce4251", k = 50)
    ex <- c("Gzmk", "Pdcd1","Tigit", "Tox","Havcr2","Entpd1", "Lag3",  "Gzmb", "Eomes","Tcf7", "Slamf6", "Il7r")
    
    FeaturePlot(tc, features= ex, order=T,cols=alpha(col.red,1),pt.size=0.5)& NoLegend() & NoAxes()
    ggsave("plots/T_cell_exhaustion_markers.pdf", width=7, height=8)
    
 ##====================================================================================================================
    
    
    tc <- readRDS("data/Brain_lung_Tcells_seurat.rds")
   
     ## Get top genes for each cluster and make heatmaps (Fig1f)
     clusterMarkers <- FindAllMarkers(object=tc)
     saveRDS(clusterMarkers, "data/Brain_lung_Tcells_markers.rds")
     write.csv(clusterMarkers, "data/Brain_lung_Tcells_markers.csv")
     
     
     top10 <- clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
     plot <- DoHeatmap(object = tc, features = top10$gene, angle = 0,size = 3)
     ggplot2::ggsave(filename = "All_T_cells_heatmap_legend.pdf", plot = plot, width = 10, height = 10)


     brain <- subset (tc, orig.ident == "BrainTcell")
     clusterbrain <- FindAllMarkers(object=brain)
     tcells <- c("Exhausted_Tcells","Central_Memory_Tcells","Resident_Memory_Tcells","Effector_Tcells")
     brain@active.ident <- factor(brain@active.ident, levels = tcells)
     
     top10 <- clusterbrain %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
     top10 <- top10 %>% arrange(factor(cluster, levels = tcells), desc(avg_log2FC))
     
     plot <- DoHeatmap(object = brain, features = top10$gene, angle = 45,size = 3) +  scale_fill_viridis() + NoLegend()
     plot
     ggplot2::ggsave(filename = "Plots/Brain_Tcells_heatmap_legend_10_2.pdf", plot = plot, width = 5, height = 6)
     

     
##====================================================================================================================     
     
    ##DOWNSAMPLE to same number of T cells per condition and plot UMAPS (Fig1 b,c)
     
    Idents(tc) <- 'orig.ident'
    downsampled <- subset(tc, downsample = 3724)
    
    
    tplot= DimPlot(downsampled,  pt.size =0.3, group.by = "orig.ident", cols= orig, order=F)
    tplot
    ggsave("plots/UMAP_downsampled_brain_lung.pdf",width=5.5, height=5)
    
      
    tplot= DimPlot(downsampled,  pt.size =0.3, group.by = "Celltype.Final", cols= c('#0e9594',"#491E54",'#7A428F',"#59A14F"), order=F)
    tplot
    ggsave("plots/UMAP_brain_lung_types_downsampled.pdf",width=6.5, height=5)

    

##====================================================================================================================     
    #Proportions (Fig1d)
    
      DefaultAssay(tc) <- "SCT"
      meta <- tc[[]]
    
      # Convert to tidy format
      meta <- meta %>% rownames_to_column("Cell") %>% tibble::as_tibble()
      
      Celltype.Freq <- meta %>% 
          group_by(orig.ident, Celltype.Final) %>%
          summarise (n = n()) %>%
          mutate(freq = n / sum(n)) %>%
          mutate(percent = paste0(round(100 * n/sum(n), 2), "%"))
      
      write.table(Celltype.Freq, file = "Freq_brain_lung_final.txt", sep = "\t", quote = F)
      
      
      # reorganize table in excell to calculate enrichment
      enrich <- read.csv("enrich_final.csv")
      
      CELLS <- c("Exhausted_Tcells",
                 "Resident_Memory_Tcells",
                 "Central_Memory_Tcells",
                 "Effector_Tcells")
      enrich$Celltype.Final <- factor( enrich$Celltype.Final,  levels=CELLS, ordered=T)
      
      
      p3 <- ggplot(enrich, aes(x = Celltype.Final, y = freq)) +
          geom_col(aes(x = Celltype.Final, y = freq ,fill = Celltype.Final), show.legend = FALSE )+
          theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
          ggtitle("") +
          theme(axis.text.y = element_text(size = 10)) +
          theme(axis.text.x = element_text(size = 10)) +
          theme(axis.text.x=element_text(colour="black"))+
          theme(axis.title.y = element_text(size = 10, face = "plain", margin = margin(r = 10))) +
          theme(axis.title.x = element_blank()) +
          labs(y = "Percent of brain T cells") +
          theme(panel.background = element_blank()) +
          scale_fill_manual(values = c('#0e9594',"#59A14F","#491E54",'#7A428F')) +
          geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")+
          theme(axis.line = element_line(colour = "black"))
      p3
      
      ggsave("plots/Proportions_Tcells_brain_lung_enrich_new.pdf", p3, width = 3.5, height = 5)
      
      
    
    

      

####
##Harmony Integration of Dulken et al., Liu et al. SVZ aging datasets.
##DR 2024

#===================================================================================================================================
# 1. Merge seurat objects
#===================================================================================================================================
        library(Seurat)
        
        setwd("~/Dropbox/Code_checking/")
        svz <- readRDS("data/DBN_seuratObject_UpdatedV3.rds")
        
        library(sctransform)
        library(glmGamPoi)
        
        svz <- SCTransform(svz, vars.to.regress = "percent.mito", verbose = FALSE, vst.flavor="v2")
        
        svz_meta <- svz@meta.data
        svz_meta$AGE <- "OLD"
        svz_meta$AGE[grepl("y", svz_meta$orig.ident)] <- "YOUNG"
        
        svz@meta.data <- svz_meta
        
        ##load liu et al dataset
        buck_exercise <- readRDS("data/ex_seurat.SVZ.annotated.2020-04-27.rds")
        
        ex_meta <- buck_exercise@meta.data
        ex_meta$AGE <- "OLD"
        ex_meta$AGE[grepl("Y", ex_meta$AgeCond)] <- "YOUNG"
        ##remove exercise groups
        to_remove_ex <- grepl("Exercise", ex_meta$AgeCond)
        ex_meta_filt <- ex_meta[!to_remove_ex,]
        buck_exercise_filt <- buck_exercise[ , colnames(buck_exercise) %in% rownames(ex_meta_filt)]
        buck_exercise_filt@meta.data <- ex_meta_filt
        rm(buck_exercise)
        buck_exercise_filt <- SCTransform(buck_exercise_filt, vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor="v2")
        
        svz@meta.data$TYPE <- "Dulken"
        buck_exercise_filt@meta.data$TYPE <- "buck_exercise"
        
        ##it's unlikely, but in the event that there are overlapping UTI index names for integrated cell sets.
        
        rownames(svz@meta.data) <- paste0("DU_", rownames(svz@meta.data))
        svz <- RenameCells(svz, paste0("DU_", colnames(svz)))
        
        rownames(buck_exercise_filt@meta.data) <- paste0("EX_", rownames(buck_exercise_filt@meta.data))
        buck_exercise_filt <- RenameCells(buck_exercise_filt, paste0("EX_", colnames(buck_exercise_filt)))
        
        seurat_set_UPD <- list(svz, buck_exercise_filt)
        
        ##let's homogenize cell type names
        
        homo_names <- function(harmonized_seurat_filt) {
            harmonized_seurat_filt@meta.data$Celltype <- as.character(harmonized_seurat_filt@meta.data$Celltype)
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSCs_NPCs"] <- "aNSCs_NPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSCs_NPC_1"] <- "aNSCs_NPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSCs_NPC_2"] <- "aNSCs_NPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSC_NPC_1"] <- "aNSCs_NPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSC_NPC_2"] <- "aNSCs_NPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Astrocytes_qNSCs"] <- "Astrocyte_qNSC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Macrophages"] <- "Macrophage"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Mural_cells"] <- "Mural"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblasts"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Oligodendrocytes"] <- "Oligodendro"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Oligodendro_1"] <- "Oligodendro"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Oligodendro_2"] <- "Oligodendro"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Oligodendro_3"] <- "Oligodendro"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Oligodendro_4"] <- "Oligodendro"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neurons"] <- "Neuron"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuronal"] <- "Neuron"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblasts"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblast_1"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblast_2"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblasts_1"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Neuroblasts_2"] <- "Neuroblast"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "OPC_1"] <- "OPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "OPC_2"] <- "OPC"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "T_cells"] <- "T_cell"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "T_Cell"] <- "T_cell"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Microglia_1"] <- "Microglia"
            harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "Microglia_2"] <- "Microglia"
            return(harmonized_seurat_filt)
        }
        
        ##not run
        
        ##let's also output some files for safekeeping.
        
        seurat_set_output <- list(svz,  buck_exercise_filt)
        names(seurat_set_output) <- c("DULKEN", "EXERCISE")
        seurat_set_output <- lapply(seurat_set_output, homo_names)
        lapply(seurat_set_output, function(x) table(x@meta.data$Celltype))
        lapply(seurat_set_output, function(x) length(which(is.na(x@meta.data$Celltype))))
        
        seurat_set_output <- seurat_set_output[names(seurat_set_output) %in% c("DULKEN", "EXERCISE")]
        
        ##Merging objects with Seurat prior to Harmony integration:
        #https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
        big_object <- merge(seurat_set_output[["DULKEN"]], y = seurat_set_output[["EXERCISE"]], add.cell.ids =c("DULKEN", "BUCKLEY"), project = "AGE")
        
        library(dplyr)
        ##Integrate datasets.
        ##first, find variable features using SCT-normalized seurat data
        variable_features <- SelectIntegrationFeatures(seurat_set_output)
        #https://github.com/satijalab/seurat/issues/5761
        VariableFeatures(big_object) <- variable_features
        big_object <- big_object %>% Seurat::NormalizeData(verbose = FALSE)  %>% ScaleData(verbose = FALSE) %>% 
            RunPCA(pc.genes = big_object@var.genes, npcs = 20, verbose = FALSE)
        
        library(harmony)
        
        ##now run harmony with the experiment as a control variable.
        pdf("harmony_outputs.pdf", width = 12, height = 12)
        big_harmony <- big_object %>% RunHarmony("TYPE", plot_convergence = T)
        dev.off()
        
        harmony_embeddings <- Embeddings(big_harmony, 'harmony')
        library(ggplot2)
        library(cowplot)
        pdf("harmony_outputs2.pdf", width = 12, height = 12)
        options(repr.plot.height = 5, repr.plot.width = 12)
        p1 <- DimPlot(object = big_harmony, reduction = "harmony", pt.size = .1, group.by = "TYPE")
        p2 <- VlnPlot(object = big_harmony, features = "harmony_1", group.by = "TYPE", pt.size = .1)
        plot_grid(p1,p2)
        dev.off()
        
        ##visually this looks good.
        
        big_harmony_proc <- big_harmony %>% RunUMAP(reduction = "harmony", dims = 1:20) %>%
            FindNeighbors(reduction = "harmony", dims = 1:20, graph.name = "test") %>%
            FindClusters(resolution = 0.5, graph.name = "test") %>% identity()
        
        pdf("harmony_outputs3.pdf", width = 12, height = 12)
        options(repr.plot.height = 4, repr.plot.width = 10)
        DimPlot(big_harmony_proc, reduction = "umap", group.by = "TYPE", pt.size = .1, split.by = 'TYPE')
        DimPlot(big_harmony_proc, reduction = "umap", group.by = "AGE", pt.size = .1, split.by = 'AGE')
        options(repr.plot.height = 4, repr.plot.width = 6)
        DimPlot(big_harmony_proc, reduction = "umap", label = TRUE, pt.size = .1)
        DimPlot(big_harmony_proc, reduction = "umap", label = T, group.by = "Celltype")
        dev.off()
        
        saveRDS(big_harmony_proc, "data/HARMONY_dulken_buckley_consolidated_April24.rds")


#===================================================================================================================================
## 2. Run sub-clustering to define Immune cell subsets
#===================================================================================================================================
           

        svz <- readRDS("data/seurat/HARMONY_dulken_buckley_consolidated_April24.rds")
        
        
        #remove doublets
        Idents(svz) <- svz$Celltype
        cells_to_remove <- grep("Doublet_", Idents(svz))
        svz <- subset(svz, cells = -cells_to_remove)
        
        table(svz$AGE)
        #OLD YOUNG 
        #22721 29775
        
        
        table(svz$orig.ident)
        
        
        
        DefaultAssay(svz) <- "RNA"
        
        # Tcells
        #subset for cd3 expression
        tc<-subset(svz, Cd3e>0 |Cd3d>0|Cd3g>0)
        table(tc@meta.data$AGE)
        #OLD YOUNG 
        #368    48
        # Get number by mouse
        table(tc@meta.data$orig.ident)
        
        
        
        ##CD4+ T cells
        cd4 <- subset(tc, Cd4 >0)
        table(cd4@meta.data$AGE)
        #OLD YOUNG 
        # 7     4 
        table(cd4@meta.data$orig.ident)
       
        ##CD8+ T cells
        cd8<-subset(tc, Cd8a >0|Cd8b1>0)
        table(cd8@meta.data$AGE)
        # OLD YOUNG 
        # 261    10 
        table(cd8@meta.data$orig.ident)
        
        ##gamma delta T cells
        tcr.genes <- grep(pattern = "^Tcr", x = rownames(x = svz$RNA@data), value = TRUE)
        Tcrg <-subset(tc, `Tcrg-C2` >0|`Tcrg-C4` >0)
        table(Tcrg@meta.data$AGE)
        #OLD 
        #43
        table(Tcrg@meta.data$orig.ident)
        
        nkt <- subset(tc, Zbtb16>0)
        table(nkt@meta.data$AGE)
        #OLD 
        #10
        table(nkt@meta.data$orig.ident)
        
        
        
        ## NK cells
        Klrb1c<-subset(svz, Klrb1c >0&Ncr1>0)
        table(Klrb1c@meta.data$AGE)
        #OLD YOUNG 
        #1     2 
        table(Klrb1c@meta.data$orig.ident)
        
        ## Foxp3+ cells
        Foxp3 <- subset(svz, Foxp3>0)
        table(Foxp3@meta.data$AGE)
        #OLD YOUNG 
        #1     2 
        
        ##B cells
        bcell <-subset(svz, Cd79a>0)
        table(bcell@meta.data$AGE)
        #OLD YOUNG 
        #13  13
        
        #Microglia
        Microglia <- subset(svz, Celltype=="Microglia")
        table(Microglia$AGE)
        #O_Control Y_Control 
        #6287  5763 
        
        #Macrophages
        Macrophages <- subset(svz, Celltype=="Macrophage")
        #remove potential b cells
        Macrophages <- subset(Macrophages, Cd79a==0)
        Macrophages <- subset(Macrophages, Cd19==0)
        table(Macrophages$AGE)
        #O_Control Y_Control 
        #226  264
        
        
    
        
        ### These numbers were used in excel to calculate the percentage of each cell type and change in percentage with age
        
#=================================================================================================================================
#DE MAST and manhattan plot
#=================================================================================================================================
         library(dplyr)
         library(MAST)
         library(tidyverse)
         library(ggthemes)
         library(ggrepel)
         library(ggpubr)

        
         ## Microglia
         DefaultAssay(Microglia) <- "RNA"
         Idents(Microglia) <- "AGE"
         
         Microglia_de <-     FindMarkers(object = Microglia,
                                         ident.1 = "OLD",
                                         ident.2 = "YOUNG",
                                         test.use = "MAST",
                                         max.cells.per.ident = 1000,
                                         min.pct = 0.0,
                                         logfc.threshold = log(0),
                                         random.seed = 3)
         
         Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
         write.table(Microglia_de, file = "data/Microglia_Old_vs_young_combined.txt", sep = "\t", quote = F)
         
         
         ## Macrophages
         
         DefaultAssay(Macrophages) <- "RNA"
         Idents(Macrophages) <- "AGE"
         
         Macrophages_de <-     FindMarkers(object = Macrophages,
                                           ident.1 = "OLD",
                                           ident.2 = "YOUNG",
                                           test.use = "MAST",
                                           max.cells.per.ident = 1000,
                                           min.pct = 0.0,
                                           logfc.threshold = log(0),
                                           random.seed = 3)
         
         Macrophages_de$fdr <- p.adjust(Macrophages_de$p_val, method = "fdr", n = 31053)
         write.table(Macrophages_de, file = "data/Macrophages_Old_vs_young_combined.txt", sep = "\t", quote = F)
         
         
         ## CD8
         DefaultAssay(cd8) <- "RNA"
         Idents(cd8) <- "AGE"
         
         cd8_de <-     FindMarkers(object = cd8,
                                   ident.1 = "OLD",
                                   ident.2 = "YOUNG",
                                   test.use = "MAST",
                                   max.cells.per.ident = 1000,
                                   min.pct = 0.0,
                                   logfc.threshold = log(0),
                                   random.seed = 3)
         
         cd8_de$fdr <- p.adjust(cd8_de$p_val, method = "fdr", n = 31053)
         write.table(cd8_de, file = "data/cd8_Old_vs_young_combined.txt", sep = "\t", quote = F)
         
         ## CD4
         DefaultAssay(cd4) <- "RNA"
         Idents(cd4) <- "AGE"
         
         cd4_de <-     FindMarkers(object = cd4,
                                   ident.1 = "OLD",
                                   ident.2 = "YOUNG",
                                   test.use = "MAST",
                                   max.cells.per.ident = 1000,
                                   min.pct = 0.0,
                                   logfc.threshold = log(0),
                                   random.seed = 3)
         
         cd4_de$fdr <- p.adjust(cd4_de$p_val, method = "fdr", n = 31053)
         view(cd4_de)
         write.table(cd4_de, file = "data/cd4_Old_vs_young_combined.txt", sep = "\t", quote = F)
         
         ##B cell
         DefaultAssay(bcell) <- "RNA"
         Idents(bcell) <- "AGE"
         
         bc_de <-     FindMarkers(object = bcell,
                                  ident.1 = "OLD",
                                  ident.2 = "YOUNG",
                                  test.use = "MAST",
                                  max.cells.per.ident = 1000,
                                  min.pct = 0.0,
                                  logfc.threshold = log(0),
                                  random.seed = 3)
         
         bc_de$fdr <- p.adjust(bc_de$p_val, method = "fdr", n = 31053)
         view(bc_de)
         write.table(bc_de, file = "data/bcell_Old_vs_young_combined.txt", sep = "\t", quote = F)
         
        
         
         ## Manhattan plot
         
         library(NCmisc)
         
         Microglia_de <-read.table("data/DE/Fig1/Microglia_Old_vs_young_combined.txt")
         Microglia_de <- tibble::rownames_to_column(Microglia_de, "gene")
         Microglia_de$Celltype <- "Microglia"
         
         Macrophages_de <- read.table("data/DE/Fig1/Macrophages_Old_vs_young_combined.txt")
         Macrophages_de <- tibble::rownames_to_column(Macrophages_de, "gene")
         Macrophages_de$Celltype <- "Macrophage"
         
         cd8_de <- read.table("data/DE/Fig1/cd8_Old_vs_young_combined.txt")
         cd8_de <- tibble::rownames_to_column(cd8_de, "gene")
         cd8_de$Celltype <- "CD8_Tcell"
         
         cd4_de <- read.table("data/DE/Fig1/cd4_Old_vs_young_combined.txt")
         cd4_de <- tibble::rownames_to_column(cd4_de, "gene")
         cd4_de$Celltype <- "CD4_Tcell"
         
         bc_de <- read.table("data/DE/Fig1/bcell_Old_vs_young_combined.txt")
         bc_de <- tibble::rownames_to_column(bc_de, "gene")
         bc_de$Celltype <- "B_cell"
         
         df <- rbind(Microglia_de,Macrophages_de,cd8_de,cd4_de,bc_de)
         
         
         # Add z-score based on two sided null hypothesis.
         df$z <- p.to.Z(df$p_val) * sign(df$avg_log2FC)
         df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_log2FC)
         
         # Randomize rows to reduce overplotting issues.
         df <- df[sample(nrow(df)), ]
         
         # Adjust factors
         CELLS <- c("Microglia","Macrophage","CD8_Tcell","CD4_Tcell", "B_cell")
         df$Celltype <- factor(df$Celltype,  levels=CELLS, ordered=T)
        
         new_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
         names(new_colors) <- levels(df$Celltype)
         
         # Make custom color column to facilitate grey coloring by threshold.
         col <- new_colors[df$Celltype]
         col[df$p_val_adj > 0.05] <- "#edede9" # grey
             df$col <- as.factor(col)
             
             q <- ggplot(df, aes(x = Celltype, y = z, color = col)) +
                 geom_jitter(shape = 16, width = 0.40, alpha = .9, size = 2) +
                 theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 14), axis.title.x = element_blank()) +
                 #ggtitle("Changes in gene expression with age (old/young)") +
                 theme(axis.title.y = element_text(size = 20, face = "plain", color= "black")) +
                 theme(axis.text.y = element_text(size = 20, color= "black")) +
                 theme(axis.text.x = element_text(size = 18,color= "black")) +
                 theme(plot.title = element_text(size=20, face = "plain")) +
                 labs(y = "Z-score") +
                 theme(legend.position="none") +
                 scale_color_manual(values = levels(df$col)) +
                 geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")+
                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
                 theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
             q
             
             ggsave("Manhattan_immune_cells_Jun25.pdf", q, width = 10, height = 5)
             
             
             
             
    
    ## Select only genes that are significantly changed
             
      df <- subset(df, df$p_val_adj < 0.05)      
      
      microglia <- subset(df, Celltype== "Microglia") 
      nrow(microglia)
      # 3447
      
      macrophage <- subset(df, Celltype== "Macrophage") 
      nrow(macrophage)
      #162
      
      
#=================================================================================================================================
      #Make heatmap of T cells from young and old SVZ of combined object > Figure 1G
#=================================================================================================================================      
      
      library(Seurat)
      library(dplyr)
      library(ggthemes)
      library(ggplot2)
      
      
      svz <- readRDS("data/seurat/HARMONY_dulken_buckley_consolidated_April24.rds")
      
      
      tc <- subset (svz, Celltype=="T_cell")
      table(tc$AGE)
      #OLD YOUNG 
      #323    53 
      rm(svz)
      
      tc$AGE <- factor(tc$AGE, levels = c("YOUNG", "OLD"))
     
      
      #Gene list
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
      pdf(file = paste0("plots/Figure1/Heatmap_Tcells", Sys.Date(), "_.pdf"), width = 6, height = 7)
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
      
      

      
      
    
             
             


library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(scales)
library(ggthemes)
library(ggpubr)

    t20 <- tableau_color_pal(palette = "Tableau 20")(20)
    t10 <- tableau_color_pal(palette = "Tableau 10")(10)
    
# Modify path below. All subsequent paths are relative.
    setwd("~/Dropbox/Code/")
    svz <- readRDS("data/IL10_svz_seurat_IL10.rds")
    
    
    Microglia <- subset(svz, Celltype== "Microglia")
    rm(svz)
    
    Microglia <- FindNeighbors(Microglia, dims = 1:20)
    Microglia <- FindClusters(Microglia, resolution = 0.2)
    Microglia <- RunUMAP(Microglia, dims = 1:20,  min.dist = 0.4, spread=1,seed=40)
    DimPlot(Microglia,cols= t20,  pt.size =1)
    
    
    clusterMarkerstc <- FindAllMarkers(object=Microglia)
    saveRDS(clusterMarkerstc, file = paste0("Microglia.markers_", Sys.Date(), ".rds"))
    #subset > remove cluster 5, cluster 4 is dividing microglia, cluster 3 inflammed, cluster 2 old, cluster 1, cluster 0 young
    
    new.cluster.ids <- c("Young_Microglia","Inter_Microglia", "Old_Microglia", "Inflammed_microglia", "Dividing_Microglia", "Remove")
    
    names(new.cluster.ids) <- levels(Microglia)
    Microglia <- RenameIdents(Microglia, new.cluster.ids)
    Microglia[["Micro_types"]] <- Idents(Microglia)
    unique(Microglia@meta.data$Micro_types)
    Microglia <- subset (Microglia, Micro_types != "Remove" )

#========================================================================================
# Data frame and violins
#=========================================================================================

    meta <- Microglia[[]]
    pca <- Embeddings(object = Microglia, reduction = "pca")[, 1:20]
    umap <- Embeddings(object = Microglia, reduction = "umap")
    rna <- t(as.matrix(GetAssayData(object = Microglia[["RNA"]], slot = "data")))
    d <- as.data.frame(cbind(meta, umap, pca, rna))
    
    d1<- subset(d,orig.ident=="Young"| orig.ident=="Control"|orig.ident=="Mut5"|orig.ident=="WT5")
    d1$orig_factor <- factor(d1$orig.ident,  levels=c("Young","Control","WT5", "Mut5"), ordered=T)
    
    orig <- c('#4cc9f0','#f84b35','#023047','#00afb9')
    
    a1 = 0.5
    a2 = 0.7
    s = 0.5
    my_comparisons <- list( c("Mut5", "Control"), c("WT5", "Control"), c("Young", "Control"), c("WT5", "Mut5") )
    
    
    ##Figure 3E
    
    gene <- "Socs3"
    p <- ggplot(d1, aes(x=orig_factor, y= Socs3)) +
        geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
        geom_violin(aes(fill=orig_factor), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
        scale_fill_manual(values=orig) +
        scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
        theme_classic()+
        theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
        theme(axis.text.y = element_text(size = 8,color="black")) +
        theme(strip.text.x = element_text(size = 10)) +
        labs(title=gene)+
        ylab(gene)+
        xlab(NULL) +
        theme(legend.position="none")+
        stat_compare_means(method = "wilcox", comparisons = my_comparisons, size=3, tip.length = 0.01, step.increase = 0.05)
    p
    
    ggsave("plots/IL10_conc/Microglia/Socs3_Microglia_WT5_Mut5_young.pdf", height=3.5, width=4.5)
    
    ##Figure 3F

    gene <- "Tnf"
    p <- ggplot(d1, aes(x=orig_factor, y= Tnf)) +
        geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
        geom_violin(aes(fill=orig_factor), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
        scale_fill_manual(values=orig) +
        scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
        theme_classic()+
        theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
        theme(axis.text.y = element_text(size = 8,color="black")) +
        theme(strip.text.x = element_text(size = 10)) +
        labs(title=gene)+
        ylab(gene)+
        xlab(NULL) +
        theme(legend.position="none")+
        stat_compare_means(comparisons = my_comparisons, size=3, tip.length = 0.01, step.increase = 0.06)+
        stat_compare_means(label.y = 4.2,size=3)
    p
    
    ggsave("plots/IL10_conc/Microglia/Tnf_stats_Microglia_WT5_Mut5.pdf", height=3.5, width=4.5)





## Plots with other concentrations

        d3<- subset(d,orig.ident=="Young"| orig.ident=="Control"|orig.ident=="Mut1"|orig.ident=="WT1")
        d3$orig_factor <- factor(d3$orig.ident,  levels=c("Young","Control","WT1","Mut1"), ordered=T)
        
        orig <- c('#4cc9f0','#f84b35','#065a60', '#90e0ef')
        
        a1 = 0.5
        a2 = 0.7
        s = 0.5
        
        my_comparisons <- list( c("Mut1", "WT1"), c("Young", "Control"),c("WT1", "Control"),c("Mut1", "Control"))
        
        
        ## Figure S7A
        genes <- c("Socs3","Il6","Tnf")
        
        
        for (gene in genes){
            print(
                ggplot(d3, aes(x=orig_factor, y= d3[,gene])) +
                    geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
                    geom_violin(aes(fill=orig_factor), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
                    scale_fill_manual(values=orig) +
                    scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
                    theme_classic()+
                    theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
                    theme(axis.text.y = element_text(size = 8,color="black")) +
                    theme(strip.text.x = element_text(size = 10)) +
                    labs(title=gene)+
                    ylab(gene)+
                    xlab(NULL) +
                    theme(legend.position="none")+
                    stat_compare_means(comparisons = my_comparisons, size=3, tip.length = 0.01, step.increase = 0.06)+
                    stat_compare_means(label.y = 4.2,size=3)
            )
            
            ggsave(paste0("plots/IL10_conc/Microglia/",gene,"_Microglia_ALL.pdf"), height=3.5, width=4.5)
        }



#====================================================================================================================================
#### Interferon gamma sum pathway violins 
#====================================================================================================================================

load(file = "mouse.gene.hallmark.kegg.reactome.rda")
gly <- hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE

##Figure 4C

    ifn_data <- d1[, colnames(d1) %in% gly] # dim: 565 308
    ifn_data$ifng_response <- rowSums(ifn_data)
    meta <- d1[, colnames(d1) %in% c("orig.ident")]
    d3 <- cbind(meta, ifn_data)
    
    d3$meta <- factor(d3$meta,  levels=c("Young","Control","WT5","Mut5"), ordered=T)
    
    my_comparisons <- list( c("Mut5", "Control"), c("WT5", "Control"), c("Young", "Control") )
    
    gene <- "Interferon gamma Response"
    p <- ggplot(d3, aes(x=meta, y= ifng_response)) +
        geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
        geom_violin(aes(fill=meta), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
        scale_fill_manual(values=orig) +
        scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
        theme_classic()+
        theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
        theme(axis.text.y = element_text(size = 8,color="black"), axis.title=element_text(size=10) ) +
        labs(title=gene)+
        ylab("Sum of IFN gamma response genes")+
        xlab(NULL) +
        theme(legend.position="none")+
        stat_compare_means(method = "wilcox", comparisons = my_comparisons, size=3, tip.length = 0.01, step.increase = 0.05)
    p
    
    ggsave("plots/IFNg_Response_microglia_KEGG_fig4C.pdf", p, height=3.5, width=4.5)



##Figure S7C
    ifn_data <- d3[, colnames(d3) %in% gly] # dim: 565 308
    ifn_data$ifng_response <- rowSums(ifn_data)
    meta <- d3[, colnames(d3) %in% c("orig.ident")]
    d4 <- cbind(meta, ifn_data)
    
    d4$meta <- factor(d4$meta,  levels=c("Young","Control","WT1","Mut1"), ordered=T)
    
    my_comparisons <- list( c("Mut1", "Control"), c("WT1", "Control"), c("Young", "Control"),c("WT1", "Mut1") )
    
    gene <- "Interferon gamma Response"
    p <- ggplot(d4, aes(x=meta, y= ifng_response)) +
        geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
        geom_violin(aes(fill=meta), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
        scale_fill_manual(values=orig) +
        scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
        theme_classic()+
        theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
        theme(axis.text.y = element_text(size = 8,color="black"), axis.title=element_text(size=10) ) +
        labs(title=gene)+
        ylab("Sum of IFN gamma response genes")+
        xlab(NULL) +
        theme(legend.position="none")+
        stat_compare_means(method = "wilcox", comparisons = my_comparisons, size=3, tip.length = 0.01, step.increase = 0.05)
    p
    
    ggsave("plots/IFNg_Response_microglia_KEGG_figS7C.pdf", p, height=3.5, width=4.5)



## more conditions


#=================================================================================================================================
#differential expression with MAST
#=================================================================================================================================

    DefaultAssay(Microglia) <- "RNA"
    Idents(Microglia) <- "orig.ident"
    
    Microglia_de <-     FindMarkers(object = Microglia,
                                    ident.1 = "WT5",
                                    ident.2 = "Control",
                                    test.use = "MAST",
                                    max.cells.per.ident = 1000,
                                    min.pct = 0.0,
                                    logfc.threshold = log(0),
                                    random.seed = 3)
    
    Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
    view(Microglia_de)
    write.table(Microglia_de, file = "data/Microglia_WT5_vs_control_DE.txt", sep = "\t", quote = F)
    
    
    Microglia_de <-     FindMarkers(object = Microglia,
                                    ident.1 = "Mut5",
                                    ident.2 = "Control",
                                    test.use = "MAST",
                                    max.cells.per.ident = 1000,
                                    min.pct = 0.0,
                                    logfc.threshold = log(0),
                                    random.seed = 3)
    
    Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
    view(Microglia_de)
    write.table(Microglia_de, file = "data/Microglia_Mut5_vs_control_DE.txt", sep = "\t", quote = F)

    
    DefaultAssay(Microglia) <- "RNA"
    Idents(Microglia) <- "orig.ident"
    
    Microglia_de <-     FindMarkers(object = Microglia,
                                    ident.1 = "WT1",
                                    ident.2 = "Control",
                                    test.use = "MAST",
                                    max.cells.per.ident = 1000,
                                    min.pct = 0.0,
                                    logfc.threshold = log(0),
                                    random.seed = 3)
    
    Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
    view(Microglia_de)
    write.table(Microglia_de, file = "data/Microglia_WT1_vs_control_DE.txt", sep = "\t", quote = F)
    
    
    Microglia_de <-     FindMarkers(object = Microglia,
                                    ident.1 = "Mut1",
                                    ident.2 = "Control",
                                    test.use = "MAST",
                                    max.cells.per.ident = 1000,
                                    min.pct = 0.0,
                                    logfc.threshold = log(0),
                                    random.seed = 3)
    
    Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
    view(Microglia_de)
    write.table(Microglia_de, file = "data/Microglia_Mut1_vs_control_DE.txt", sep = "\t", quote = F)
    



#=================================================================================================================
## DF dot plot
#=================================================================================================================


## FC plot with IL-10 pathway genes > Figure 3G
    
        il10 <- c("Il10ra","Il10rb","Jak1", "Tyk2","Stat3","Socs3",
                 "Hmox1", "Blvra","Blvrb", "Bcl3","Etv3","Nfil3","Sbno2","Cdkn2d",
                 "Il1a","Il1b","Il6","Tnf","Tgfb1","Nfkbia","Ccl3", "Ccl4",
                 "H2-Aa","H2-Eb1","H2-DMa","H2-Oa", "H2-Ob","Cd86")
        
        
        df4 <- read.table("data/DE/IL10/Microglia_Mut5_vs_control_DE.txt")
        df4 <- tibble::rownames_to_column(df4, "gene")
        df4 <- df4%>% subset(gene%in%il10)
        df4 <- df4 %>% mutate(celltype = "Microglia")
        df4 <- df4 %>% mutate(rep = "Mut")
        
        
        df3 <- read.table("data/DE/IL10/Microglia_WT5_vs_control_DE.txt")
        df3 <- tibble::rownames_to_column(df3, "gene")
        df3 <- df3 %>% subset(gene%in%il10)
        df3 <- df3 %>% mutate(celltype = "Microglia")
        df3 <- df3 %>% mutate(rep = "WT")
        
        df <- rbind (df3,df4)
        df$gene_factor <- factor(df$gene,  levels=il10, ordered=F)
        
        
        t <- ggplot(df, aes(x = fct_rev(rep), y = gene_factor, fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_distiller(palette="RdBu", limits = c(-0.5,0.5), oob=squish) +
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 10, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.key.size = unit(0.4, 'cm'),legend.title = element_text(size=8))
        
        t
        
        ggsave("plots/FC_microglia_Mut_WT_IL10_Fig3G.pdf", width = 3, height =5.5)



## FC plot with IFNg pathway genes > Figure 4D

        ifng <- c("Stat1","Bst2","Isg15", "Ifit1","Ifit2","Ifit3","Ifitm3","Ifitm1","Ifi27l2a","Ifi27", "Irf2", "Irf5", "Irf7","Irf9", "B2m", "H2-K1", "H2-D1","Ccl5", "Gbp1", "Adar")
        
        df4 <- read.table("data/DE/IL10/Microglia_Mut5_vs_control_DE.txt")
        df4 <- tibble::rownames_to_column(df4, "gene")
        df4 <- df4%>% subset(gene%in%ifng)
        df4 <- df4 %>% mutate(celltype = "Microglia")
        df4 <- df4 %>% mutate(rep = "Mut")
        
        df3 <- read.table("data/DE/IL10/Microglia_WT5_vs_control_DE.txt")
        df3 <- tibble::rownames_to_column(df3, "gene")
        df3 <- df3 %>% subset(gene%in%ifng)
        df3 <- df3 %>% mutate(celltype = "Microglia")
        df3 <- df3 %>% mutate(rep = "WT")
        
        
        df <- rbind (df3,df4)
        df$gene_factor <- factor(df$gene,  levels=ifng, ordered=F)
        
        
        t <- ggplot(df, aes(x = fct_rev(rep), y = gene_factor, fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_distiller(palette="RdBu", limits = c(-1.5,1.5), oob=squish) +
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 10, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.key.size = unit(0.4, 'cm'),legend.title = element_text(size=8))
        
        t
        ggsave("plots/FC_microglia_Mut_WT_Ifng_Fig4D.pdf", width = 3, height =3)
        

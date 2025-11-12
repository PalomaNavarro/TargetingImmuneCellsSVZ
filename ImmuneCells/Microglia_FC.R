
# Make heatmap of Log2FC changes in microglia

rm(list = ls())
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)
library(ggthemes)

setwd("/Dropbox/Code_checking")
        
 ## Load MAST results from microglia only

        #Liu et al 2023 dataset
        ex <- readRDS("data/Microglia_Old_vs_young_mast_EXERCISE.rds")
        ex$cond <- "Liu"
        ex <- tibble::rownames_to_column(ex, "gene")
        
        #Dulken et al 2019 dataset
        dulken <- readRDS("data/Celltype_DE/Microglia_Old_vs_young_mast_Dulken.rds")
        dulken$cond <-"Dulken"
        dulken <- tibble::rownames_to_column(dulken, "gene")
        
        #Buckley et al 2023 dataset
        bu <- readRDS("data/Microglia_Old_vs_young_mast_Buckley.rds")
        bu$cond <-"Buckley"
        bu <- tibble::rownames_to_column(bu, "gene")
        
        #combine all 3 datasets together
        data <- rbind(ex, dulken,bu)
        

        ##Activated/homeostatic microglia genes
        genes <- c("Apoe","Axl","Cst7","Ctsb","Itgax","Il1b", "Lgals3","Csf1","Ccr5","Cx3cr1","Glul","Gpr34","Serinc3", "Siglech","Spi1", "Tmem119", "P2ry12")
    
        input <- data[data$gene %in% genes,] # SUBSET GENES
        input$gene <- factor(input$gene,  levels= genes, ordered=T)
        input$cond <- factor(input$cond,  levels= c("Liu","Dulken","Buckley"), ordered=T)# order 
        
        t <- ggplot(data = input, mapping = aes(x = gene, y = cond, fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_gradient2( low = "#48bfe3",
                                  mid = "white",
                                  high = "#e5383b")+
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 8, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.position = "right") +theme(legend.key.size = unit(0.15, "in"))
        t
        ggsave("plots/microglia_liu_dulken_buckley_activated_homeostatic.pdf", width = 10, height =1.5)
        
        
        ##Pro and anti-inflammatory genes
        genes <- c("Stat1","Bst2","Ifi27l2a", "Ifitm3","H2-K1", "H2-D1","Cd74","Tnf","Il1b","Il10ra","Il10rb","Il10","Socs3", "Stat3","Tgfb1","Cd86")
        
        input <- data[data$gene %in% genes,] # SUBSET GENES
        input$gene <- factor(input$gene,  levels= genes, ordered=T)
        input$cond <- factor(input$cond,  levels= c("Liu","Dulken","Buckley"), ordered=T)# order 
        
        t <- ggplot(data = input, mapping = aes(x = gene, y = cond, fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_gradient2( low = "#48bfe3",
                                  mid = "white",
                                  high = "#e5383b")+
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 8, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.position = "right") +theme(legend.key.size = unit(0.15, "in"))
        t
        
        
        ggsave("plots/microglia_liu_dulken_buckley_il10_ifn.pdf", width = 10, height =1.5)
        

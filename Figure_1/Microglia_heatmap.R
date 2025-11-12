library(Seurat)
library(dplyr)
library(ggthemes)
library(ggplot2)

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code")
svz <- readRDS("data/DBN_seuratObject_UpdatedV3.rds")

micro <- subset (svz, Celltype=="Microglia")
rm(svz)

micro$Age <- factor(micro$Age, levels = c("y", "o"))
table(micro$Age)


orig <- c( "skyblue","darkred")
color <- magma(n = 20, direction = -1)

genelist <-c("Lgals3","Cst7","Id2","Ccl4","Stat1","Bst2","Ifitm3","Ifi27l2a","Il6","Il1b","Il1a",
             "Il10rb","Il10ra","Socs3", "Stat3", "Cd86")


mapal <- c("#fefae0","#FCFDBFFF", "#FB8861FF","#E85362FF","#D6456CFF", "#B63679FF" ,"#802582FF")
DoHeatmap(object = micro, features = genelist, group.by = "Age",group.colors = orig) + scale_fill_gradientn(colours = mapal) + 
    theme(axis.text.y = element_text(color = "black", size = 10))

ggsave("plots/heatmap2.pdf", width = 5, height = 5)

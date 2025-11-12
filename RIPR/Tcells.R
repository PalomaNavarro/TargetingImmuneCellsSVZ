## Analyze T cells only from RIPR repeat 1 experiment

library(Seurat)
library(dplyr)
library(tidyverse)
library(scales)
library(ggthemes)
library(cowplot)
library(ggpubr)
library(rstatix)
library(ggplot2)
sessionInfo()

        t20 <- tableau_color_pal(palette = "Tableau 20")(20)
        t10 <- tableau_color_pal(palette = "Tableau 10")(10)
        col.red <- CustomPalette(low = "#dcdcdd", high = "#720026", mid ="#ce4251", k = 50)
        orig <- c('#343a40','#7209b7')

        
setwd("~/Dropbox/Code/")
svz <- readRDS("data/seurat/RIPR_svz_tcr_seurat.rds")
tc <- subset(svz, Celltype =="T_cells")

# Total number of cells
table(svz$orig.ident)
#OldControl      OldRIPR 
# 6236         7181
rm(svz)

#Total number of T cells
table(tc$orig.ident)
#OldControl      OldRIPR
#   97          243     

#Total number of T cells that express Ki67
ki67 <- subset(tc, Mki67>0)
table(ki67$orig.ident)
#OldControl      OldRIPR
# 5           33   

cd8 <- subset(tc, Cd8a>0|Cd8b1>0)
table(cd8$orig.ident)
#OldControl      OldRIPR 
#69          180           

cd8ki67 <- subset(cd8, Mki67>0)
table(cd8ki67$orig.ident)
#OldControl      OldRIPR 
#3           25           



tc <- FindNeighbors(tc, dims = 1:20)
tc <- FindClusters(tc, resolution = 0.7)


tc <- RunUMAP(tc, dims = 1:20,  min.dist = 0.5, spread = 1.4, seed=1)
DimPlot(tc,cols= t20,  pt.size =2)

#Conditions UMAP
tc <- subset(tc, orig.ident=="OldControl"|orig.ident=="OldRIPR")
DimPlot(tc, group.by = "orig.ident", cols = orig, pt.size =2, order=F)+ NoLegend() +  labs(title = NULL)
ggsave("plots/Tcells_UMAP_condition_RIPR.pdf", width=4, height=4)


#Mki67 UMAP
FeaturePlot(tc, features = c("Mki67"),pt.size =2,cols = col.red, order=T)+  labs(title = NULL)+ NoLegend()
ggsave("plots/Tcells_UMAP_Mki67.pdf", width=4, height=4)


###Calculating T-cell state UCell scores
##DR - Feb 2024

t_ref_dataset <- tc

ripr <- subset(t_ref_dataset, subset = orig.ident == "OldRIPR")
control <- subset(t_ref_dataset, subset = orig.ident == "OldControl")

library(SignatuR)
programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)
ripr$Celltype <- "T_cell"
control$Celltype <- "T_cell"

ripr.score <- AddModuleScore_UCell(ripr, features = programs, assay = "RNA", name = NULL)
control.score <- AddModuleScore_UCell(control, features = programs, assay = "RNA", name = NULL)

ripr_frame <- do.call("rbind", lapply(unique(names(programs)), 
                                      function(x) data.frame(score = x, VALUE = as.numeric(ripr.score[[x]][,1]), TYPE = "RIPR")))
control_frame <- do.call("rbind", lapply(unique(names(programs)), 
                                         function(x) data.frame(score = x, VALUE = as.numeric(control.score[[x]][,1]), TYPE = "Control")))

combined <- do.call("rbind", list(ripr_frame, control_frame))


##Clean program names
combined <- combined[combined$score %in% c( "Tcell.stemness", "Tcell.cytotoxicity", "Tcell.exhaustion"),]
combined$score[combined$score == "Tcell.stemness"] <- "Stemness"
combined$score[combined$score == "Tcell.cytotoxicity"] <- "Cytotoxicity"
combined$score[combined$score == "Tcell.exhaustion"] <- "Exhaustion"
combined$TYPE <- factor(combined$TYPE, levels = c("Control", "RIPR"))


bxp <- ggboxplot(combined, x = "score", y = "VALUE", 
                 color = "TYPE", palette = c('#adb5bd','#7209b7'))
bxp

stat.test <- combined %>%
    group_by(score) %>%
    wilcox_test(VALUE ~ TYPE)

stat.test
#score        .y.   group1 group2    n1    n2 statistic         p
#1 Cytotoxicity VALUE Control RIPR      98   244    12496. 0.501     
#2 Exhaustion   VALUE Control RIPR      98   244     8234  0.00000527
#3 Stemness     VALUE Control RIPR      98   244    13324  0.0696    



##annoying alphabetic-requirements.
stat.test <- stat.test[unlist(lapply(unique(combined$score), function(x) which(stat.test$score == x))),]

stat.test <- stat.test %>%
    add_xy_position(x = "score", dodge = 0.8)

stata <- stat.test[, 1:which(grepl("x", colnames(stat.test)))[1] - 1]
statb <- stat.test[, which(grepl("x", colnames(stat.test)))[1] : dim(stat.test)[2]]
reorder <- unlist(lapply(unique(stat.test$score), function(x) which(stat.test$score == x)))
statb <- statb[order(statb$x),]
stat.test <- cbind(stata, statb)

test_sigA <- bxp + stat_pvalue_manual(
    stat.test,  label = "p", tip.length = 0, size=4, step.increase = 0.1
)

test_sig3A <- test_sigA + xlab("") + ylab("UCell Score")
test_sig3A <- test_sig3A+theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 12, color="black"),
                               axis.text.y = element_text(size = 8, color="black"),
                               axis.title =element_text(size=12)) + theme(legend.position = "none")
test_sig3A 
ggsave("plots/Tcell_scores.pdf",test_sig3A, width = 4, height = 4.5)


#################
#################
#################

### T cell violins

orig <- c('#adb5bd','#7209b7')

# Load Data
setwd("~/Dropbox/Code_revision/")
d <- readRDS("data/seurat/RIPR_svzPlotDf.rda")
d <- d[sample(nrow(d)),]

d <- subset(d, d$orig.ident=="OldRIPR"|d$orig.ident=="OldControl")
d$orig.ident <- factor(d$orig.ident,  levels=c("OldControl","OldRIPR"), ordered=T)

d <- subset(d, d$Celltype=="T_cells")


a1 = 0.3
a2 = 0.7
s = 0.5


genes <- c("Mki67", "Pdcd1", "Tigit", "Tox", "Gzmk", "Gzmb","Gzma", "Cd69", "Ifng", "Tnf", "Il2", "Cd8a", "Cd8b1", "Prf1", "Cd44", "Klrg1", "Nkg7")


for (gene in genes){
    print(
        ggplot(d, aes(x=orig.ident, y= d[,gene])) +
            geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.2, h = 0)) +
            geom_violin(aes(fill=orig.ident), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
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
            stat_compare_means(method = "wilcox", comparisons = list(c("OldControl", "OldRIPR")), size=3, tip.length = 0.01, step.increase = 0.5)
        
    )
    
    ggsave(paste0("plots/",gene,"_Tcells_RIPR_old.pdf"), height=3, width=3)
}




###==========================================================================================================================
# DOTPLOT FC
##===================================================================================
df<- read.table("data/Tcell_Old_RIPRControl_DE.txt")
df <- tibble::rownames_to_column(df, "genes")
df<- df %>% mutate(celltype = "T_cells")

tcells <- c("Nr4a2","Tox","Tigit","Pdcd1","Gzmk","Il2","Ifng","Tnf","Prf1","Lamp1","Gzma","Gzmb","Cd44","Cd69","Mki67")

p <- df %>% subset(genes%in%tcells)
p$gene_factor <- factor(p$genes,  levels=tcells, ordered=T)


q <-ggplot(data = p, mapping = aes(x = celltype, y = gene_factor, color = avg_logFC)) +
    geom_point(size=6) +
    scale_color_distiller(palette="RdBu", limits = c(-1,1)*max(abs(p$avg_logFC))) +
    geom_point(shape = 1) +
    ggtitle("RIPR vs control Old") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
    theme(legend.position = "right")

q
ggsave("plots/Tcells_dotplot-new2.pdf", q, width = 3, height = 6.5)





q <- ggplot(data = p, mapping = aes(x = celltype, y = gene_factor, fill = avg_logFC)) +
    geom_tile() +  # Use geom_tile to create rectangles
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1) * max(abs(p$avg_logFC))) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 12, color = "black"), 
          axis.text.y = element_text(size = 12, color = "black")) +
    theme(legend.position = "right")+
    theme(legend.key.size = unit(0.4, 'cm'),legend.title = element_text(size=8))

q

ggsave("plots/Tcells_dotplot-July21.pdf", q, width = 2.5, height = 4.5)

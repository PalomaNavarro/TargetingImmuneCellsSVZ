## Violin plots for microglia 

rm(list=ls())
library(ggplot2)
library(ggpubr)
library(magrittr)

# Modify path below. All subsequent paths are relative.
setwd("Dropbox/Code")

# Load Data
load(file = "svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
d <- d[sample(nrow(d)),]


#Subset and plot
m <- subset(d, d$Celltype=="Microglia")
m$Age <- factor(m$Age, levels=c("y", "o"), ordered=T)


# Plot parameters
a1 = 0.25
a2 = 0.7
s = 0.4
ageColors <- c("deepskyblue", "firebrick")

genes <- c("Socs3","Il10ra", "Tnf","Il1b", "Lgals3", "Apoe", "Tmem119", "P2ry12", "Cd86", "Ifi27l2a", "H2-K1")

for (gene in genes){
    print(
        ggplot(m, aes(x=Age, y= m[,gene])) +
            geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.2, h = 0)) +
            geom_violin(aes(fill=Age), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
            scale_fill_manual(values=ageColors) +
            scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
            theme_classic()+
            theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
            theme(axis.text.y = element_text(size = 8,color="black")) +
            theme(strip.text.x = element_text(size = 10)) +
            labs(title=gene)+
            ylab(gene)+
            xlab(NULL) +
            theme(legend.position="none")+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.line = element_blank())+
            stat_compare_means(method = "wilcox", comparisons = list(c("o", "y")), size=3, tip.length = 0.01, step.increase = 0.5)
    )
    
    ggsave(paste0("plots/",gene,"_Microglia_violin.pdf"), height=3, width=3)
}


# Plot sum of genes in IL-10 pathway (Fig 1j)

        il10d <- c("Il10ra","Il10rb","Il10","Jak1","Tyk2","Socs3", "Stat3","Tgfb1","Tgfb2")
        
        # Subset data to signature genes
        sub_data <- m[, colnames(m) %in% il10d]
        meta <- m[, colnames(m) %in% c("Age", "Replicate", "Celltype")]
        sub_data$sub_response <- rowSums(sub_data)
        sub_data <- cbind(meta, sub_data)
        sub_data$Age <- factor(sub_data$Age, levels=c("y", "o"), ordered=T)
        
        
        # Plot parameters
        a1 = 0.25
        a2 = 0.7
        s = 0.4
        ageColors <- c("deepskyblue", "firebrick")
        symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                            symbols = c("****", "***", "**", "*", "ns"))
        
        gene <- "Sum of IL-10 genes"
        p <- ggplot(data=sub_data, aes(x=Age, y=sub_response)) +
            geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.2, h = 0)) +
            geom_violin(aes(fill=Age), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
            scale_fill_manual(values=ageColors) +
            #facet_wrap(~celltype_factor, nrow =2) +
            scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0.05, 0.065))) +
            theme_classic()+
            theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"))+
            theme(axis.text.y = element_text(size = 8,color="black"), axis.title=element_text(size=10) ) +
            labs(title=gene)+
            ylab("Sum of IL-10 pathway")+
            xlab(NULL) +
            theme(legend.position="none")+
            stat_compare_means(method = "wilcox", comparisons = list(c("o", "y")), size=3, tip.length = 0.01, step.increase = 0.5)
        p
        
        ggsave("plots/violins/sum_IL-10_genes_violin.pdf", p, height=3, width=3)
        

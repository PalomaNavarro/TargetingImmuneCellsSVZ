## Analyze T cells only


library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix)
library(scales)
library(ggthemes)
library(cowplot)
library(ggpubr)
sessionInfo()

setwd("~/Dropbox/Code/")
svz <- readRDS("data/seurat/Labeled_svz_seurat_IL10_conc_2023-06-13.rds")
tc <- subset(svz, Celltype =="T_cells")
rm(svz)


meta <- tc[[]]
pca <- Embeddings(object = tc, reduction = "pca")[, 1:20]
umap <- Embeddings(object = tc, reduction = "umap")
rna <- t(as.matrix(GetAssayData(object = tc[["RNA"]], slot = "data")))
d <- as.data.frame(cbind(meta, umap, pca, rna))



d <- readRDS("data/seurat/tcPlotDf_label_2023-06-13.rda")

## WT5 and Mutant 5

    d1<- subset(d,orig.ident=="Young"| orig.ident=="Control"|orig.ident=="Mut5"|orig.ident=="WT5")
    d1$orig_factor <- factor(d1$orig.ident,  levels=c("Young","Control","WT5", "Mut5"), ordered=T)
    
    orig.all <- c('#4cc9f0','#f84b35','#023047','#00afb9')
    
    
    a1 = 0.5
    a2 = 0.7
    s = 0.5
    my_comparisons <- list( c("Mut5", "WT5"), c("WT5", "Control"), c("Young", "Control"), c("Wt5","Mut5") )
    
    genes <- c("Pdcd1","Ifng","Gzmb","Socs3", "Stat3", "Mki67" )
    for (gene in genes){
        print(
            ggplot(d1, aes(x=orig_factor, y= d1[,gene])) +
                geom_point(size=s, alpha=a1, color="grey", position = position_jitter(w = 0.2, h = 0)) +
                geom_violin(aes(fill=orig_factor), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
                scale_fill_manual(values=orig.all) +
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
        
        ggsave(paste0("plots/",gene,"_Tcells_WT5_Mut5_Young.pdf"), height=3.5, width=4.5)
        
    }
    

## Plot only cells that express Gzmb

    dg <- filter(d1, Gzmb >0)
    dg<- subset(dg,orig.ident=="Young"| orig.ident=="Control"|orig.ident=="Mut5"|orig.ident=="WT5")
    dg$orig_factor <- factor(dg$orig.ident,  levels=c("Young","Control","WT5", "Mut5"), ordered=T)
    
    orig.all <- c('#f84b35','#023047','#00afb9')
    my_comparisons <- list( c("Mut5", "WT5"), c("WT5", "Control"), c("Young", "Control"),c("Mut5", "Control") )
    
    gene <- "Gzmb"
    p <- ggplot(dg, aes(x=orig_factor, y= Gzmb)) +
        geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.2, h = 0)) +
        geom_violin(aes(fill=orig_factor), alpha=a2, trim=T,scale="width", draw_quantiles = c(.5)) +
        scale_fill_manual(values=orig.all) +
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
    
    ggsave("plots/violins/Gzmb_Tcells_subset_nozero.pdf", p, height=3.5, width=4.5)


## WT1 and Mutant 1

    d2<- subset(d,orig.ident=="Young"| orig.ident=="Control"|orig.ident=="Mut1"|orig.ident=="WT1")
    d2$orig_factor <- factor(d2$orig.ident,  levels=c("Young","Control","WT1", "Mut1"), ordered=T)
    orig <- c('#4cc9f0','#f84b35','#065a60', '#90e0ef')

## Plot only cells that express Gzmb or Ifng

    dg <- filter(d2, Gzmb >0)
    
    gene <- "Gzmb"
    p <- ggplot(dg, aes(x=orig_factor, y= Gzmb)) +
        geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.2, h = 0)) +
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
    
    ggsave("plots/IL10_conc/Tcells/WT1_Mut1/Gzmb_Tcells_subset_nozero.pdf", p, height=3.5, width=4.5)



#=================================================================================================================================
## T cell scores
#=================================================================================================================================

        Control <- subset(tc, subset = orig.ident == "Control")
        WT_Il10 <- subset(tc, subset = orig.ident == "WT5")
        Mut_Il10 <- subset(tc, subset = orig.ident == "Mut5")
        
        library(SignatuR)
        
        programs <- GetSignature(SignatuR$Mm$Programs)
        names(programs)
        
        library(UCell)
        
        control.score <- AddModuleScore_UCell(Control, features = programs, assay = "RNA", name = NULL)
        wt.score <- AddModuleScore_UCell(WT_Il10, features = programs, assay = "RNA", name = NULL)
        mut.score <- AddModuleScore_UCell(Mut_Il10, features = programs, assay = "RNA", name = NULL)
        
        control_score_frame <- do.call("rbind", lapply(unique(names(programs)), 
                                                       function(x) data.frame(score = x, VALUE = as.numeric(control.score[[x]][,1]), TYPE = "CONTROL")))
        wt_score_frame <- do.call("rbind", lapply(unique(names(programs)), 
                                                  function(x) data.frame(score = x, VALUE = as.numeric(wt.score[[x]][,1]), TYPE = "WT")))
        mut_score_frame <- do.call("rbind", lapply(unique(names(programs)), 
                                                   function(x) data.frame(score = x, VALUE = as.numeric(mut.score[[x]][,1]), TYPE = "MUT")))
        
        combined <- do.call("rbind", list(control_score_frame,wt_score_frame, mut_score_frame))
        
        library(ggplot2)
        
        ##Clean program names
        combined <- combined[combined$score %in% c( "Tcell.stemness", "Tcell.cytotoxicity", "Tcell.exhaustion"),]
        combined$score[combined$score == "Tcell.stemness"] <- "Stemness"
        combined$score[combined$score == "Tcell.cytotoxicity"] <- "Cytotoxicity"
        combined$score[combined$score == "Tcell.exhaustion"] <- "Exhaustion"
        combined$TYPE <- factor(combined$TYPE, levels = c("CONTROL", "WT","MUT"))
        
        library(dplyr)
        library(ggpubr)
        library(rstatix)
        
        bxp <- ggboxplot(
            combined, x = "score", y = "VALUE", 
            color = "TYPE", palette = c('#ADB5BD','#023047','#00AFB9')
        )
        bxp + ylim(-0.1, NA)
        bxp
        
        ggsave("plots/Tcellscores.pdf", height=4.5, width=4)
        
        stat.test <- combined %>%
            group_by(score) %>%
            t_test(VALUE ~ TYPE) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj")
        stat.test
        
        
    
        ##Calculating significance for low-cell-counts
        
        ##Comparisons to control.
        stat.test <- stat.test[grepl("Control", stat.test$group1),]
        
        final_split <- split(combined, combined$score)
        
        do_test_empirical <- function(frame, var1 = "Wild-type Il10", var2 = "Control") {
            ##wherein I calculate the mean difference between control and a treatment in UCell score.
            ##then randomly sample from the pooled set of cells, matching # of cells, and calculate mean
            ##and ask how often our mean is greater than random expectation (label swapping).
            
            target_diff <- mean(frame$VALUE[frame$TYPE == var1]) - mean(frame$VALUE[frame$TYPE == var2])
            
            #randomly sampling from the entire frame.
            rand_N <- 10000
            random_differences <- list()
            for (x in 1:rand_N) {
                rand_var1 <- frame$VALUE[sample(1:dim(frame)[1], length(which(frame$TYPE == var1)))]
                rand_var2 <- frame$VALUE[sample(1:dim(frame)[1], length(which(frame$TYPE == var2)))]
                random_diff <- mean(rand_var1) - mean(rand_var2)
                random_differences[[toString(x)]] <- random_diff
            }
            
            curr_fdr <- length(which(unlist(random_differences) > target_diff)) / length(random_differences)
            
            out_frame <- data.frame(SCORE = frame$score[1], var1, var2, diff = target_diff, mean_random = mean(unlist(random_differences)),
                                    sd_random = sd(unlist(random_differences)), num_better = length(which(unlist(random_differences) > target_diff)),
                                    rand_N, curr_fdr)
            return(out_frame)
        }
        
        ##and then we'll just add another column to stat.test
        
        row_wise_res_list <- list()
        for (x in 1:dim(stat.test)[1]) {
            row_wise_res <- do_test_empirical(final_split[[stat.test$score[x]]], var1 = stat.test$group2[x], var2 = stat.test$group1[x])
            print(paste0(row_wise_res$SCORE, "-", row_wise_res$var1, "-", row_wise_res$var2))
            row_wise_res_list[[toString(x)]] <- row_wise_res
        }
        
        row_wise_res_frame <- do.call("rbind", row_wise_res_list)
        
        stat.test$empirical <- row_wise_res_frame$curr_fdr
        
        ##add significance codes
        stat.test <- stat.test %>% add_significance("empirical")
        
        stat.test <- stat.test %>%
            add_xy_position(x = "score", dodge = 0.8, scales = "fixed")
        
        stat.test$y.position[seq(1, dim(stat.test)[1], 2)] = 0.8
        stat.test$y.position[seq(2, dim(stat.test)[1], 2)] = 0.9
        
        stata <- stat.test[, 1:which(grepl("x", colnames(stat.test)))[1] - 1]
        statb <- stat.test[, which(grepl("x", colnames(stat.test)))[1] : dim(stat.test)[2]]
        reorder <- unlist(lapply(unique(combined$score), function(x) which(stat.test$score == x)))
        statb <- statb[order(statb$x),]
        stata <- stata[reorder,]
        stat.test <- cbind(stata, statb)
        
        test_sigA <- bxp + stat_pvalue_manual(  stat.test,  label = "empirical", tip.length = 0, size = rel(1.5), bracket.size = 0.1)
        
        test_sig3A <- test_sigA + xlab("") + ylab("UCell Score") + theme(text = element_text(size = rel(1.5)))
        test_sig3A <- test_sig3A+theme(    axis.text=element_text(size = rel(1.5)),  axis.title=element_text(size=rel(1.5),face="bold"))
        
        test_sig4A <- test_sig3A + theme(axis.line = element_line(colour = 'black', linewidth = rel(0.01))) + theme(axis.ticks = element_line(colour = 'black', linewidth = rel(0.01))) + theme(legend.position = "none")
        
        print(test_sig4A)
        
        
        
        
#=================================================================================================================================
# Differential expression with MAST
#=================================================================================================================================
        
        DefaultAssay(tc) <- "RNA"
        Idents(tc) <- "orig.ident"
        
        Microglia_de <-     FindMarkers(object = tc,
                                        ident.1 = "WT5",
                                        ident.2 = "Control",
                                        test.use = "MAST",
                                        max.cells.per.ident = 1000,
                                        min.pct = 0.0,
                                        logfc.threshold = log(0),
                                        random.seed = 3)
        
        Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
        view(Microglia_de)
        write.table(Microglia_de, file = "data/Tcells_WT5_vs_control_DE.txt", sep = "\t", quote = F)
        
        
        Microglia_de <-     FindMarkers(object = tc,
                                        ident.1 = "Mut5",
                                        ident.2 = "Control",
                                        test.use = "MAST",
                                        max.cells.per.ident = 1000,
                                        min.pct = 0.0,
                                        logfc.threshold = log(0),
                                        random.seed = 3)
        
        Microglia_de$fdr <- p.adjust(Microglia_de$p_val, method = "fdr", n = 31053)
        view(Microglia_de)
        write.table(Microglia_de, file = "data/Tcells_Mut5_vs_control_DE.txt", sep = "\t", quote = F)
        
        
#=================================================================================================================
## DF dot plot
#=================================================================================================================
        
        
        
        ## WT5 and Mut5
        tcells <- c("Nr4a2","Tox","Gzmk","Pdcd1","Ptprc","Klrg1","S1pr1", "Lef1", "Ccr7","Sell", "Il7r","Cd62l","Mki67","Tnf","Il2","Ifng","Cd69","Prf1","Gzma","Gzmb","Socs3")
        #Ptprc> CD45ra
    
        # Dotplot
        df1 <-read.table("data/DE/IL10/Tcell_WT5vsControl_DE.txt")
        df1 <- df1 %>% mutate(rep = "WT")
        df1$gene <- rownames(df1)
        df1 <- subset(df1, gene%in%tcells)
        
        df2 <-read.table("data/DE/IL10/Tcell_Mut5vsControl_DE.txt")
        df2 <- df2 %>% mutate(rep = "Mut")
        df2$gene <- rownames(df2)
        df2 <- subset(df2, gene%in%tcells)
        
        #combine
        df <- rbind (df1,df2)
        
        df$gene_factor <- factor(df$gene,  levels= tcells, ordered=T)
        df$rep_factor <- factor(df$rep,  levels= c("WT", "Mut"), ordered=T)
        

        ## FC plot with IL-10 pathway genes in T cells no p value
        
        df$gene_factor <- factor(df$gene,  levels= tcells, ordered=F)
        
        t <- ggplot(df, aes(x = rep_factor, y = fct_rev(gene_factor), fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_distiller(palette="RdBu", limits = c(-1,1.1)*max(abs(df$avg_log2FC))) +
            #scale_fill_distiller(palette="RdBu", limits = c(-2,2), oob=squish) +
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 10, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.key.size = unit(0.4, 'cm'),legend.title = element_text(size=8))
        
        t
        
        ggsave("plots/Tcells_dotplot_IL10_Wt_Mut5_Tcellsubtypes_square.pdf", width = 2.5, height =3.5)        
        
        ## WT1 and Mut1
        tcells <- c("Nr4a2","Tox","Gzmk","Pdcd1","Ptprc","Klrg1","S1pr1", "Lef1", "Ccr7","Sell", "Il7r","Cd62l","Mki67","Tnf","Il2","Ifng","Cd69","Prf1","Gzma","Gzmb","Socs3")
        #Ptprc> CD45ra
        
        # Dotplot
        df3 <-read.table("data/DE/IL10/Tcell_WT1vsControl_DE.txt")
        df3 <- df3 %>% mutate(rep = "WT1")
        df3$gene <- rownames(df3)
        df3 <- subset(df3, gene%in%tcells)
        
        df4 <-read.table("data/DE/IL10/Tcell_Mut1vsControl_DE.txt")
        df4 <- df4 %>% mutate(rep = "Mut1")
        df4$gene <- rownames(df4)
        df4 <- subset(df4, gene%in%tcells)
        
        #combine
        df <- rbind (df3,df4)
        
        df$gene_factor <- factor(df$gene,  levels= tcells, ordered=T)
        df$rep_factor <- factor(df$rep,  levels= c("WT1", "Mut1"), ordered=T)
        
        
        
        ## FC plot with IL-10 pathway genes in T cells no p value
        
        df$gene_factor <- factor(df$gene,  levels= tcells, ordered=F)
        
        t <- ggplot(df, aes(x = rep_factor, y = fct_rev(gene_factor), fill = avg_log2FC)) +
            geom_tile() +
            xlab(label = "Gene")+
            scale_fill_distiller(palette="RdBu", limits = c(-1,1.1)*max(abs(df$avg_log2FC))) +
            #scale_fill_distiller(palette="RdBu", limits = c(-2,2), oob=squish) +
            theme_classic() +
            theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_text(size = 10, color= "black")) +
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10,color= "black"))+
            scale_x_discrete(position = "top")+
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),  axis.line = element_blank())+
            scale_y_discrete(limits=rev)+
            theme(legend.key.size = unit(0.4, 'cm'),legend.title = element_text(size=8))
        
        t
        
        ggsave("plots/IL10_conc/Tcells/Tcells_dotplot_IL10_WT1_Mut1_square.pdf", width = 2.5, height =3.5)        
        
        


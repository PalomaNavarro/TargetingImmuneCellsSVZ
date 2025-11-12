
# Aging clock predictions IL-10 experiment with WT and engineered variant

library(Seurat)
library(tidyverse)
library(cowplot)
library(ggthemes)
library(scales)
library(glmnet)

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code")
svz <- readRDS("/data/Labeled_svz_seurat_IL10_conc_2023-06-13.rds")

# Inspect
    Convert_to_Dataframe <- function(svz) {
        DefaultAssay(svz) <- "RNA"
        svz[["SCT"]] <- NULL
        Celltypes <- c("Oligodendrocytes", "Microglia", "aNSCs_NPCs",
                       "Astrocytes_qNSCs", "Neuroblasts", "Endothelial")
        svz <- subset(svz, subset = Celltype %in% Celltypes)
        meta <- svz@meta.data
        meta <- meta[, c("orig.ident", "Celltype")]
        raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
        raw_counts <- raw_counts[, colSums(raw_counts) > 0]
        df <- as_tibble(cbind(meta, raw_counts))
        return(df)
    }
    
    df <- Convert_to_Dataframe(svz)
    df.meta <- df[,c(1:2)]
    df <- df[, -c(1:2)]
    
# Load Clock Training Object to match column structure
    svz_28 <- readRDS("data/bootstrap_pseudocell_15.rds")
    genes <- colnames(svz_28[1,6][[1]][[1]]) # use these genes to subset RIPR data. 20948 genes.
    
# Format new data
    df.genes <- colnames(df) # 21394
    df.missing <- setdiff(genes, df.genes)
    if (length(as.numeric(df.missing)) > 0) {
      print("There are some genes in the training data
            that are not in test data. Filling with Zeros.")
      missing_df <- matrix(0, nrow(df), length(df.missing))
      colnames(missing_df) <- df.missing
      rownames(missing_df) <- rownames(df)
      df <- cbind(df, missing_df)
    }
    df <- df[, genes] # Reorder and cut to size to match clock training data format
    
    
# Combine with and nest by metadata
    
    df <- cbind(df.meta, df)
    df <-  df %>% group_by(Celltype, orig.ident) %>% nest()
    
    
bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
        pseudocells <- c()
        # If dynamic then only sample with replacement if required due to shortage of cells.
        if (replace == "dynamic") {
            if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
        }
        for (i in c(1:n)) {
            batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
            pseudocells <- rbind(pseudocells, colSums(batch))
        }
        colnames(pseudocells) <- colnames(df)
        return(as_tibble(pseudocells))
    }
        
    # Apply boostrap.pseudocells using map()
    df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells))
    
    # Remove single cell data; keep just key metadata and pseudocells
    df2$data <- NULL
    saveRDS(df2, "data/bootstrap_pseudocell_IL10_conc_March24.rds")
    
#==================================================================================================
    ## PREDICTIONS
#==================================================================================================
    
    
    il10_data <- df2
    
    # Load IL10 data
    il10_data <- unnest(il10_data) # dim: 2400 21396
    meta <- as.data.frame(il10_data[, c(1:2)])
    umi <- il10_data[, -c(1:2)]
    normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
    logged <- log1p(normed * 10000)
    
    by_celltype <- as_tibble(cbind(meta, logged)) %>%
        group_by(Celltype) %>%
        nest()
    
    # Load models
    models <- readRDS("data/models_all_bootstrap.rds")
    colnames(models)[1] <- "Celltype"
    models$lognormalized <- NULL
    models$Celltype <- plyr::revalue(models$Celltype,
                                     c("Astrocyte_qNSC" = "Astrocytes_qNSCs",
                                       "Oligodendro" = "Oligodendrocytes",
                                       "aNSC_NPC" = "aNSCs_NPCs",
                                       "Neuroblast" = "Neuroblasts")
    )
    
    
    # combine
    by_celltype <-  dplyr::inner_join(by_celltype, models, by = "Celltype")
    custom_add_predictions <- function(data, model) {
        predictions <- predict(model, newx = as.matrix(data[,-1]), s = "lambda.min")
        p <- as_tibble(data.frame("Pred" = predictions[,1],
                                  "Sample" = data[,1]))
        return(p)
    }
    
    
    # Add predictions 
    by_celltype <- by_celltype %>%
        mutate(Predictions = map2(data, model, custom_add_predictions))
    
    output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)
    
    range(output$Pred) #0.8335413 34.0262098
    
    saveRDS(output, "data/IL10_conc_clock_predictions_March24.rds")


#==============================================================
    ## PLOTTING
#==============================================================

data<- output

data$Celltype <- factor(data$Celltype, levels = c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts",
                                                  "Microglia", "Oligodendrocytes",  "Endothelial" ))
data$Sample <- factor(data$orig.ident, levels = c("Young","Control", "WT1", "WT5", "WT50","Mut1", "Mut5", "Mut50" ))



# Dotplots

##**WT5 and Mut5 vs control**

        effect_table <- c()
        
        for (celltype in unique(data$Celltype)) {
            print(celltype)
            d1 <- filter(data, Celltype == celltype)
            print(head(d1))
            eff_Mut1 <- median(d1[d1$Sample == "Mut5", "Pred"]$Pred) - median(d1[d1$Sample == "Control", "Pred"]$Pred)
            eff_WT1 <- median(d1[d1$Sample == "WT5", "Pred"]$Pred) - median(d1[d1$Sample == "Control", "Pred"]$Pred)
            eff_WT <- median(d1[d1$Sample == "Mut5", "Pred"]$Pred) - median(d1[d1$Sample == "WT5", "Pred"]$Pred)
            
            effect_table <- rbind(effect_table, c(celltype,eff_Mut1,eff_WT1,eff_WT))
        }
        
        colnames(effect_table) <- c("Celltype","Mut5",  "WT5","ratio" )
        
        et <- as_tibble(effect_table)
        et_tidy <- pivot_longer(et, cols = c("Mut5", "WT5", "ratio"), names_to = "Age", values_to = "Effect_Size")
        et_tidy$Effect_Size <- as.numeric(et_tidy$Effect_Size)
        
        et_tidy$Celltype <- factor(et_tidy$Celltype, levels = c("Neuroblasts","aNSCs_NPCs","Astrocytes_qNSCs","Oligodendrocytes","Endothelial","Microglia"))
        et_tidy$Age <- factor(et_tidy$Age, levels = c("WT5","Mut5", "ratio"))
        
        filter(et_tidy, Celltype=="NA")
        
        ggplot(data = et_tidy, mapping = aes(x = Celltype, y = Age, size = abs(Effect_Size), color = Effect_Size)) +
            geom_point() +
            scale_color_distiller(palette="RdBu", limits = c(-0.8,1)*max(abs(et_tidy$Effect_Size))) +
            geom_point(shape = 1, color = "black") +
            theme_classic() +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
            theme(axis.text.x = element_text(size = 8, color = "black"), axis.text.y = element_text(size = 8, color = "black")) +
            theme(legend.position = "right") +
            theme(legend.key.size = unit(0.1, "cm"))+
            
            coord_flip()
        
        ggsave("plots/WT5_Mut5_vs_control_Il10_conc_dotplot_CHRONO_ratio.pdf", width = 3.5, height = 2.25, useDingbats=F)




##**WT1 and Mut1 vs control**

        effect_table <- c()
        
        for (celltype in unique(data$Celltype)) {
            print(celltype)
            d1 <- filter(data, Celltype == celltype)
            print(head(d1))
            eff_Mut1 <- median(d1[d1$Sample == "Mut1", "Pred"]$Pred) - median(d1[d1$Sample == "Control", "Pred"]$Pred)
            eff_WT1 <- median(d1[d1$Sample == "WT1", "Pred"]$Pred) - median(d1[d1$Sample == "Control", "Pred"]$Pred)
            eff_WT <- median(d1[d1$Sample == "Mut1", "Pred"]$Pred) - median(d1[d1$Sample == "WT1", "Pred"]$Pred)
            
            effect_table <- rbind(effect_table, c(celltype,eff_Mut1,eff_WT1,eff_WT))
        }
        
        colnames(effect_table) <- c("Celltype","Mut1",  "WT1","ratio" )
        
        et <- as_tibble(effect_table)
        et_tidy <- pivot_longer(et, cols = c("Mut1", "WT1", "ratio"), names_to = "Age", values_to = "Effect_Size")
        et_tidy$Effect_Size <- as.numeric(et_tidy$Effect_Size)
        
        et_tidy$Celltype <- factor(et_tidy$Celltype, levels = c("Neuroblasts","aNSCs_NPCs","Astrocytes_qNSCs","Oligodendrocytes","Endothelial","Microglia"))
        et_tidy$Age <- factor(et_tidy$Age, levels = c("WT1","Mut1", "ratio"))
        
        filter(et_tidy, Celltype=="NA")
        
        ggplot(data = et_tidy, mapping = aes(x = Celltype, y = Age, size = abs(Effect_Size), color = Effect_Size)) +
            geom_point() +
            scale_color_distiller(palette="RdBu", limits = c(-0.8,1)*max(abs(et_tidy$Effect_Size))) +
            geom_point(shape = 1, color = "black") +
            theme_classic() +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
            theme(axis.text.x = element_text(size = 8, color = "black"), axis.text.y = element_text(size = 8, color = "black")) +
            theme(legend.position = "right") +
            theme(legend.key.size = unit(0.15, "cm"))+
            
            coord_flip()
        
        ggsave("plots/WT1_Mut1_vs_control_Il10_conc_dotplot_CHRONO_ratio.pdf", width = 3.5, height = 2.25, useDingbats=F)  




### Histogram

    sub <- subset(data, Sample=="Young"|Sample=="Control"|Sample=="WT5"|Sample=="Mut5")
    
    sub$Sample <- factor(sub$Sample, levels = c("Young","Control", "WT5", "Mut5" ))
    sub$Celltype <- factor(sub$Celltype, levels = c("Microglia", "Oligodendrocytes",  "Endothelial","Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts" ))
    
    orig <- c('#4cc9f0','#f84b35','#065a60','#00afb9')
    
    
    
    ggplot(sub, aes(x=Pred, fill=Sample, color=Sample)) +
        geom_density(alpha = 0.7) +
        facet_wrap(. ~ Celltype, nrow=2) +
        scale_fill_manual(values = alpha(orig, 1)) +
        scale_color_manual(values = alpha(orig, 0.9)) +
        theme_classic() +
        xlim(c(0, 35)) + ylim(c(0, 0.3)) +
        ylab("Density") + xlab("Predicted Age")+
        theme(
            legend.text = element_text(size = 8),   # change legend label size
            legend.title = element_text(size = 8),
            strip.text.x = element_text(size = 8),
            strip.background = element_blank(),
            axis.text.x = element_text(color = "black", size=8),
            axis.title.x = element_text(color = "black", size=8),
            axis.title.y = element_text(color = "black", size=8),
            axis.text.y = element_text(color = "black", size=8),
            legend.key.size = unit(0.4, "cm")# change legend title size
        )
    
    ggsave("plots/IL10_conc_histo_chrono_all.pdf", width=6, height=4,useDingbats=F)


## subset WT5 and Mut 5 and make violins

sub <- subset(data,Sample=="Control"|Sample=="WT5"|Sample=="Mut5")
sub$Sample <- factor(sub$Sample, levels = c("Control", "WT5", "Mut5" ))


sub$Celltype <- factor(sub$Celltype, levels = c("Microglia","Endothelial","Oligodendrocytes","Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts" ))
sub$Sample <- factor(sub$Sample, levels = c("Control", "WT5", "Mut5" ))

orig <- c('#adb5bd','#023047','#00afb9')

    ggplot(data=sub, aes(x = Sample, y = Pred, fill =Sample)) +
        geom_violin(scale = "width", draw_quantiles = .5) +
        facet_grid(. ~ Celltype) +
        scale_fill_manual(values = alpha(orig, 0.7)) +
        ylab("Age Score") + 
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black", size=8),  # show x-axis text
            axis.text.y = element_text(color = "black", size=8),
            axis.title.y = element_text(color = "black", size=8),
            axis.title.x = element_blank(),
            strip.text.x = element_text(size = 9),
            strip.background = element_blank(),     
            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            legend.position = "none")
    )


ggsave(paste0("plots/Il10_conc_chrono_Violins_NSC", Sys.Date(), ".pdf"), width = 8, height = 2.5, useDingbats=F)



























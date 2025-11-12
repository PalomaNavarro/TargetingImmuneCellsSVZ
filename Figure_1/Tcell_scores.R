###Calculating T-cell state UCell scores

library(Seurat)

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/Code/")
t_ref_dataset <- readRDS("data/Brain_lung_Tcells_seurat.rds")

mega_big_query.lung <- subset(t_ref_dataset, subset = orig.ident == "LungTcell")
mega_big_query.brain <- subset(t_ref_dataset, subset = orig.ident == "BrainTcell")

#remotes::install_github("carmonalab/SignatuR")
library(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)

mega_big_query.lung$Celltype <- "T_cell"
mega_big_query.brain$Celltype <- "T_cell"

query.lung.score <- AddModuleScore_UCell(mega_big_query.lung, features = programs, assay = "RNA", name = NULL)
query.brain.score <- AddModuleScore_UCell(mega_big_query.brain, features = programs, assay = "RNA", name = NULL)

lung_score_frame <- do.call("rbind", lapply(unique(names(programs)), 
	function(x) data.frame(score = x, VALUE = as.numeric(query.lung.score[[x]][,1]), TYPE = "LUNG")))
brain_score_frame <- do.call("rbind", lapply(unique(names(programs)), 
	function(x) data.frame(score = x, VALUE = as.numeric(query.brain.score[[x]][,1]), TYPE = "BRAIN")))

combined <- do.call("rbind", list(lung_score_frame, brain_score_frame))
library(ggplot2)

##Clean program names
combined <- combined[combined$score %in% c( "Tcell.stemness", "Tcell.cytotoxicity", "Tcell.exhaustion"),]
combined$score[combined$score == "Tcell.stemness"] <- "Stemness"
combined$score[combined$score == "Tcell.cytotoxicity"] <- "Cytotoxicity"
combined$score[combined$score == "Tcell.exhaustion"] <- "Exhaustion"
combined$TYPE <- factor(combined$TYPE, levels = c("BRAIN", "LUNG"))

library(dplyr)
library(ggpubr)
library(rstatix)

bxp <- ggboxplot(combined, x = "score", y = "VALUE", 
  color = "TYPE", palette = c('#702670','#0e9594'))
  
bxp
stat.test <- combined %>%
  group_by(score) %>%
  wilcox_test(VALUE ~ TYPE)
  
stat.test
#score        .y.   group1 group2    n1    n2 statistic         p
#1 Cytotoxicity VALUE BRAIN  LUNG    3724  4976   9740081 2.90e-  5
#2 Exhaustion   VALUE BRAIN  LUNG    3724  4976  13229821 3.16e-280
#3 Stemness     VALUE BRAIN  LUNG    3724  4976   3798455 0        



##alphabetic-requirements.
stat.test <- stat.test[unlist(lapply(unique(combined$score), function(x) which(stat.test$score == x))),]

stat.test <- stat.test %>%
    add_xy_position(x = "score", dodge = 0.8)

stata <- stat.test[, 1:which(grepl("x", colnames(stat.test)))[1] - 1]
statb <- stat.test[, which(grepl("x", colnames(stat.test)))[1] : dim(stat.test)[2]]
reorder <- unlist(lapply(unique(stat.test$score), function(x) which(stat.test$score == x)))
statb <- statb[order(statb$x),]
stat.test <- cbind(stata, statb)

test_sigA <- bxp + stat_pvalue_manual(
    stat.test,  label = "p", tip.length = 0, size=3
)

test_sig3A <- test_sigA + xlab("") + ylab("UCell Score")
test_sig3A <- test_sig3A+theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 10, color="black"),
                               axis.text.y = element_text(size = 10, color="black"),
                               axis.title =element_text(size=10)) + theme(legend.position = "none")
test_sig3A 
ggsave("plots/brain_lung_T_cell_scores.pdf",test_sig3A, width = 4, height = 5)


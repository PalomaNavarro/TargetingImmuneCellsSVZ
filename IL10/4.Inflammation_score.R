###Generating an 'inflammatory score' signature for use with UCell
##DR - 2025
###

library(SignatuR)
data(SignatuR)
SignatuR <- AddNode(SignatuR,
  parent_node=SignatuR$Mm,
  name="IL10")

 positive_genes <<- c("H2-Aa", "H2-Ab1", "H2-DMa", "H2-Oa", "H2-Ob",
        "Stat1", "Bst2", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Ifitm3", "Ifi27", "Ifi27l2a", "B2m", "H2-K1", "H2-D1", "Ccl4", "Ccl3", "Ccl5")

    negative_genes <<- c("Il10ra", "Il10rb", "Jak1", "Tyk2", "Stat3", "Socs3", "Hmox1", "Blvra", "Etv3", "Cdkn2d")

negated_inflamm <- c("Il1a", "Il1b", "Il6", "Tnf", "Nfkbia")
##UCell can take into account the negative effects of genes by appending a '-' to the end of the gene name.
negated_inflamm <- unlist(lapply(negated_inflamm, function(x) paste0(x, '-')))

negative_genes <- c(negative_genes, negated_inflamm)

SignatuR <- AddSignature(SignatuR,
	node=SignatuR$Mm$IL10,
	name="MANUAL_PRO_INFLAMM",
	reference="A simple T cell signature",
	signature=unique(positive_genes))

SignatuR <- AddSignature(SignatuR,
	node=SignatuR$Mm$IL10,
	name="MANUAL_ANTI_INFLAMM",
	reference="A simple T cell signature",
	signature=unique(negative_genes))

saveRDS(SignatuR, "MANUAL_INFLAMM_SIGNATURE_FLIPPED_SIGN_JULY2025.rds")

calculate_UCELL_program <- function(cd8, outfix, programs, do_seurat = F, do_precalc = NULL, just_gen = F) {

library(UCell)

if(!is.null(do_precalc)) {

    print("doing precalc")
    if(!file.exists(do_precalc)) {

cd8_old <- cd8[, cd8$TREAT == "TREAT"]
cd8_young <- cd8[, cd8$TREAT != "TREAT"]

        print("need to generate ranks")
            cd8_old_dat <- GetAssayData(cd8_old, assay = "RNA", slot = "data")
    cd8_young_dat <- GetAssayData(cd8_young, assay = "RNA", slot = "data")
    OLD_ranks <- StoreRankings_UCell(cd8_old_dat)
    YOUNG_ranks <- StoreRankings_UCell(cd8_young_dat)
    out_list <- list()
    out_list[["TREAT"]] <- OLD_ranks
    out_list[["CONTROL"]] <- YOUNG_ranks
    saveRDS(out_list, do_precalc) 
    if(just_gen) {
        print("exit")
        return(NA)
    }
    }else{
        if(just_gen) {
            print("exit")
            return(NA)
        }
        all_ranks <- readRDS(do_precalc)
        OLD_ranks <- all_ranks[["TREAT"]]
        YOUNG_ranks <- all_ranks[["CONTROL"]]
    }

    query.old.score <- as.data.frame(ScoreSignatures_UCell(features = programs, precalc.rank = OLD_ranks))
    query.young.score <- as.data.frame(ScoreSignatures_UCell(features = programs, precalc.rank = YOUNG_ranks))

neo_programs <- paste0(names(programs), "_UCell")
old_score_frame <- do.call("rbind", lapply(neo_programs, 
	function(x) data.frame(score = x, VALUE = as.numeric(query.old.score[, x]), TYPE = "TREAT")))
young_score_frame <- do.call("rbind", lapply(unique(neo_programs), 
	function(x) data.frame(score = x, VALUE = as.numeric(query.young.score[, x]), TYPE = "CONTROL")))

combined <- do.call("rbind", list(young_score_frame, old_score_frame))

}else{
print("KNN method")

cd8_old <- cd8[, cd8$TREAT == "TREAT"]
cd8_young <- cd8[, cd8$TREAT != "TREAT"]

query.old.score <- AddModuleScore_UCell(cd8_old, features = programs, assay = "RNA", name = NULL)
query.old.score <- SmoothKNN(query.old.score, signature.names = names(programs), reduction="pca")
query.young.score <- AddModuleScore_UCell(cd8_young, features = programs, assay = "RNA", name = NULL)
query.young.score <- SmoothKNN(query.young.score, signature.names = names(programs), reduction="pca")
##smoothening...
neo_programs <- paste0(names(programs), "_kNN")

old_score_frame <- do.call("rbind", lapply(unique(neo_programs), 
	function(x) data.frame(score = x, VALUE = as.numeric(query.old.score[[x]][,1]), TYPE = "TREAT")))
young_score_frame <- do.call("rbind", lapply(unique(neo_programs), 
	function(x) data.frame(score = x, VALUE = as.numeric(query.young.score[[x]][,1]), TYPE = "CONTROL")))

combined <- do.call("rbind", list(young_score_frame, old_score_frame))
}

if(do_seurat) {
test <- NULL
try(test <- AddModuleScore(cd8, programs, name = names(programs), assay = "RNA"))

if(is.null(test)) {
    print("failed")
    #browser()
    try(test <- AddModuleScore(cd8, programs, name = names(programs), assay = "RNA", nbin = 12))
    if(is.null(test)) {
      try(test <- AddModuleScore(cd8, programs, name = names(programs), assay = "RNA", nbin = 2))
    }
}

if(is.null(test)) {
    print("seurat scoring totally fails")
    browser()
}
upreg <- data.frame(score = "MANUAL_PRO", VALUE = test$MANUAL_PRO_INFLAMM, TYPE = test$TREAT)
downreg <- data.frame(score = "MANUAL_ANTI", VALUE = test$MANUAL_ANTI_INFLAMM, TYPE = test$TREAT)
seurat_combined <- do.call("rbind", list(upreg, downreg))

seurat_combined$TYPE <- factor(seurat_combined$TYPE, levels = c("CONTROL", "TREAT"))
seurat_combined$CELL <- outfix
}

combined$TYPE <- factor(combined$TYPE, levels = c("CONTROL", "TREAT"))
orig_combined <- combined
orig_combined$CELL <- outfix
combined <- combined[combined$VALUE != 0,]

library(dplyr)
library(ggpubr)
library(rstatix)
bxp <- ggboxplot(
  combined, x = "score", y = "VALUE", 
  color = "TYPE", palette = c("#00AFBB", "#E7B800")
  )

stat.test <- NULL
try(stat.test <- combined %>%
  group_by(score) %>%
  t_test(VALUE ~ TYPE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
)
if(!is.null(stat.test)) {
stat.test <- stat.test %>%
  add_xy_position(x = "score", dodge = 0.8)

test_sig <- bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )

test_sig2 <- bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

test_sigA <- bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0, size = 12
  )

test_sig2A <- bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
##we'll clean up

test_sig3A <- test_sigA + xlab("") + ylab("UCell Score")
}else{
    test_sig3A <- bxp + xlab("") + ylab("UCell Score")
test_sig3A <- test_sig3A+theme(    axis.text=element_text(size = rel(1)),#(size=12), #text = element_text(size = 12),
        axis.title=element_text(size=rel(2.5),face="bold")) ##14
}
test_sig3A <- test_sig3A + theme(axis.text.x = element_text(angle = 90))

return_list <- list()
return_list[[1]] <- test_sig3A + ggtitle(outfix)
return_list[[2]] <- orig_combined

if(do_seurat) {
return_list[[3]] <- seurat_combined
}

return(return_list)
}

do_the_wrapper <- function (big_harmony_proc, target_signature, outfix, do_seurat = F, do_precalc = NULL) {

library(SignatuR)

library(UCell)

homo_names <- function(harmonized_seurat_filt) {
harmonized_seurat_filt@meta.data$Celltype <- as.character(harmonized_seurat_filt@meta.data$Celltype)
harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSC_NPC"] <- "aNSCs_NPC"
harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSCs_NPC"] <- "aNSCs_NPC"
harmonized_seurat_filt@meta.data$Celltype[harmonized_seurat_filt@meta.data$Celltype == "aNSCs_NPC"] <- "aNSCs_NPC"
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

cells_of_interest <- c("Microglia",  "Endothelial", "Astrocyte_qNSC", "aNSCs_NPC", "Neuroblast", "Oligodendro")

if(!is.null(do_precalc)) {
    print("attempting to pre-calculate ranks")
    if(!file.exists(paste0(do_precalc, "_Oligodendro.rds"))) {
        print("we're going to have to generate")
        big_harmony_proc <- homo_names(big_harmony_proc)
cells_of_interest <- intersect(cells_of_interest, unique(big_harmony_proc$Celltype))

    for (cell in cells_of_interest) {
    curr_set <- big_harmony_proc[, big_harmony_proc$Celltype == cell]
    err <- calculate_UCELL_program(curr_set, cell, programs, do_precalc = paste0(do_precalc, "_", cell, ".rds"), just_gen = T)
    print(err)
    print("done")
    print(cell)
    ##One runs the function once to generate pre-calculated ranks.
    ##then a second time to actually calculate UCell scores.
    }

    return(NA)
    }else{

spec <- readRDS(target_signature)
programs <- GetSignature(spec$Mm$IL10)
print(names(programs))

cells_files <- system(paste0("ls ", do_precalc, "*.rds"), intern = T)
cells_files <- gsub(paste0(do_precalc, "_"), "", cells_files)
cells_of_interest <- gsub(".rds", "", cells_files)

print("then can actually *use* precalc-ranks")
print("this number of cells")
print(length(cells_of_interest))
print(paste0(cells_of_interest, collapse = "@"))

cell_plots <- list()
for (cell in cells_of_interest) {
cell_plots[[cell]] <- do_the_thing_custom("filler", cell, programs, do_precalc = paste0(do_precalc, "_", cell, ".rds"))
}

}

}else{

spec <- readRDS(target_signature)
programs <- GetSignature(spec$Mm$IL10)
print(names(programs))

outfix <- paste0(outfix, "_MULTICELL_CUSTOM_MANUAL")
big_harmony_proc <- homo_names(big_harmony_proc)

cell_plots <- list()
for (cell in cells_of_interest) {
      print(cell)
    curr_set <- big_harmony_proc[, big_harmony_proc$Celltype == cell]
    if(dim(curr_set)[1] == 0) {
      print("...no?")
      browser()
      next()
    }
    if(length(unique(curr_set$TREAT)) != 2) {
      print("it's a binary comparison that fails")
      next()
    }

    cell_plots[[cell]] <- do_the_thing_custom(curr_set, cell, programs, do_seurat = do_seurat)
    print(paste0("finished ", cell))
}
}

pdf(paste0("SCORING/", outfix, "_custom_scoring_SMOOTHENED2.pdf"), width = 12, height = 12)
library(gridExtra)
for (cell in cells_of_interest) {
    print(cell_plots[[cell]][[1]])
}
##all together now.

just_cell_plots <- lapply(cell_plots, function(x) x[[1]])

grid.arrange(grobs = just_cell_plots, ncol = 3)
dev.off()

cellwise_supercollapse <- do.call("rbind", lapply(cell_plots, function(x) x[[2]]))
write.csv(cellwise_supercollapse, paste0("SCORING/", outfix, "_custom_scoring_MULTICELL.csv"), row.names = F)

}

####
##new treatments

library(Seurat)
library(MAST)
library(grDevices)
curr_study <- readRDS("data/IL10_svz_seurat.rds")
source("/extra/SVZ_ILUVATAR/custom_scoring_manual_MULTICELL.R")

target_signature <- "/extra/SVZ_ILUVATAR/MANUAL_INFLAMM_SIGNATURE_FLIPPED_SIGN_JULY2025.rds"
##keep in mind for the FUNCITONAL_NEW_EXP results, it's always compared to OLD.

for (treat in c("Mut1", "Mut5", "Mut50", "WT1", "WT5", "WT50", "Young")) {

  curr_sub <- curr_study[, curr_study$orig.ident %in% c("Control", treat)]
  curr_sub$TREAT <- "CONTROL"
  curr_sub$TREAT[curr_sub$orig.ident != "Control"] <- "TREAT"
  do_the_wrapper(curr_sub, "filler", "filler", do_precalc = paste0("SCORING/NEW_EXP_ranks_", treat, "_MULTI"))
  print(treat)
  
}

for (treat in c("Mut1", "Mut5", "Mut50", "WT1", "WT5", "WT50", "Young")) {

    if(file.exists(paste0("SCORING/NEW_EXP_ranks_", treat, "_MULTI_Oligodendro.rds"))) {
        print("done precalc")
        do_the_wrapper("filler", target_signature, paste0("FUNCTIONAL_NEW_EXP_", treat, "_INVERTED_MK2"), do_precalc = paste0("SCORING/NEW_EXP_ranks_", treat, "_MULTI"))
    }else{
        print("error")
    }
    }

rm(curr_study)
gc()
gc(environment())


###################################################################################
###################################################################################
###################################################################################
##Combining results into one large figure

score_set <- readLines("multi_cell.txt")

##Clean up names
score_set <- score_set[!grepl("RIPR_YOUNG", score_set)]
score_set <- score_set[!grepl("COMBO_ranks_Ripr.Il10", score_set)]
#score_set <- score_set[!grepl("NEW_EXP_ranks_Young.rds", score_set)]

score_names <- gsub(".rds", "", score_set)
score_names <- gsub("NEW_EXP_", "", score_names)
score_names <- gsub("CUST_UPD2_", "", score_names)
score_names <- gsub("_MICRO_CUSTOM", "", score_names)
score_names <- gsub("_custom_scoring_MULTICELL.csv", "", score_names)
score_names <- gsub("FUNCTIONAL_", "", score_names)
score_names <- gsub("_INVERTED_MK2", "", score_names)

##requested re-labels
score_names[score_names == "WT1"] <- "IL10 WT (1) on Old"
score_names[score_names == "WT_R1"] <- "IL10 WT (R1) on Old"
score_names[score_names == "WT5"] <- "IL10 WT (5) on Old"
score_names[score_names == "Mut1"] <- "IL10 Mutant (1) on Old"
score_names[score_names == "Mut5"] <- "IL10 Mutant (5) on Old"
score_names[score_names == "Mut1_R1"] <- "IL10 Mutant (R1) on Old"
score_names[score_names == "Young"] <- "IL10 Exp - Young Control"

score_frames <- lapply(score_set, read.csv)

for (x in 1:length(score_frames)) {
    score_frames[[x]]$TREAT <- score_names[x]
}

mega_cutdown <- c("IL10 Exp - Young Control", "IL10 WT (5) on Old", "IL10 Mutant (5) on Old")

score_frames <- score_frames[score_names %in% mega_cutdown]
names(score_frames) <- score_names[score_names %in% mega_cutdown]

YOUNG_EVERYTHING <- FALSE
young_flip <- TRUE
if(YOUNG_EVERYTHING) {

#make everything relative to YOUNG IL10 EXPERIMENT.
young_control <- score_frames[["IL10 Exp - Young Control"]]
young_control_control <- young_control[young_control$TYPE == "TREAT",] ##it was coded as young-vs-old. We want old-vs-young.
young_control_OLD <- young_control[young_control$TYPE == "CONTROL",] ##it was coded as young-vs-old. We want old-vs-young.

young_control_control$TYPE <- "CONTROL"
young_control_OLD$TYPE <- "TREAT"

restitch_age <- rbind(young_control_control, young_control_OLD)
##manually encoding...
a <- score_frames[["IL10 WT (5) on Old"]]
b <- score_frames[["IL10 Mutant (5) on Old"]]

##replace the old 'control' with the young 'control'
##first, remove the old control
a <- a[a$TYPE != "CONTROL",]
b <- b[b$TYPE != "CONTROL",]

#now, sub in the young control
a <- rbind(a, young_control_control)
b <- rbind(b, young_control_control)
a$TREAT <- "IL10 WT (5) on Old"
b$TREAT <- "IL10 Mutant (5) on Old"

score_frames <- list()
score_frames[[1]] <- restitch_age
score_frames[[2]] <- a
score_frames[[3]] <- b
##make sure they all have the same # of control cells
#lapply(score_frames, function(x) table(x$TYPE))
}

if(young_flip) {

##we can represent the young sample as old-vs-young,
  ##while the treatments are 'treatment - vs - old control'
  curr_frame <- score_frames[["IL10 Exp - Young Control"]]
  curr_frame$TYPE[curr_frame$TYPE == "CONTROL"] <- "OLD"
  curr_frame$TYPE[curr_frame$TYPE == "TREAT"] <- "CONTROL"
  curr_frame$TYPE[curr_frame$TYPE == "OLD"] <- "TREAT"
  ##sanity-check: the 'control' (young cells) should be much larger than the 'treatment'.
  score_frames[["IL10 Exp - Young Control"]] <- curr_frame
}

mega_score <- do.call("rbind", score_frames)

##let's just do one score at a time.
##save this for below, where we report absolute score values
saveRDS(mega_score, "mega_score_INVERTED_JULY2025_for_dotplot_gen.rds")

library(dplyr)
library(ggpubr)
library(rstatix)
library(dplyr)
if(young_flip) {
  mega_cutdown[1] <- "IL10 Old - vs - Young Controls"
 mega_score$TREAT[mega_score$TREAT == "IL10 Exp - Young Control"]  <-  "IL10 Old - vs - Young Controls"
}
mega_split <- split(mega_score, mega_score$score)

library(gridExtra)

ultra_split <- split(mega_score, mega_score$CELL)

cellwise_plots <- list()
cellwise_frames <- list()

for (cell in names(ultra_split)) {
    mega_frame <- ultra_split[[cell]]

    mega_split <- split(mega_frame, mega_frame$score)

treats_of_interest <- mega_cutdown

##massive_box calculates the fold-increase or fold-decrease in activity of a given signature in treated
##cells relative to their respective controls.
massive_box <- function(big_score_frame, score) {

treatwise_rows <- list()

treats_of_interest <- mega_cutdown

setdiff(big_score_frame$TREAT, treats_of_interest)
setdiff(treats_of_interest, big_score_frame$TREAT)

for (treat in treats_of_interest) {
    print(treat)
combined <- big_score_frame[big_score_frame$TREAT == treat,]

combined$TYPE <- factor(combined$TYPE, levels = c("TREAT", "CONTROL"))
orig_combined <- combined
combined$VALUE[combined$VALUE == 0] <- 0.01 ##Avoid undefined errors.

curr_diff_mean <- mean(combined$VALUE[combined$TYPE != "CONTROL"]) - mean(combined$VALUE[combined$TYPE == "CONTROL"])
curr_log_mean <- log2(mean(combined$VALUE[combined$TYPE != "CONTROL"])/mean(combined$VALUE[combined$TYPE == "CONTROL"]))
curr_log_median <- log2(median(combined$VALUE[combined$TYPE != "CONTROL"])/median(combined$VALUE[combined$TYPE == "CONTROL"]))

curr_row <- data.frame(TREAT = treat, score = score, mean_control = mean(combined$VALUE[combined$TYPE == "CONTROL"]),
  mean_treat = mean(combined$VALUE[combined$TYPE != "CONTROL"]), curr_diff_mean, curr_log_mean, T_test_p = NA, 
  T_test_stat = NA, curr_log_median)

treatwise_rows[[treat]] <- curr_row
}

row_frame <- do.call("rbind", treatwise_rows)

return(row_frame)
}

row_scores <- lapply(1:length(mega_split), function(x) massive_box(mega_split[[x]], names(mega_split)[x]))

treats_of_interest <- mega_cutdown

row_collapse <- do.call("rbind", row_scores)
row_collapse$score <- as.factor(row_collapse$score)
row_collapse$TREAT <- factor(row_collapse$TREAT, levels = treats_of_interest)

####
# do the 'ratio' of pro- and anti-inflammatory genes.
mini_split <- split(mega_frame, mega_frame$TREAT)

ratio_maker <- function(miniframe, score = "RATIO") {

control_frame <- miniframe[miniframe$TYPE == "CONTROL",]
treat_frame <- miniframe[miniframe$TYPE == "TREAT",]

treat_ratio <- log2(mean(treat_frame$VALUE[treat_frame$score == "MANUAL_ANTI_INFLAMM_UCell"])/mean(treat_frame$VALUE[treat_frame$score == "MANUAL_PRO_INFLAMM_UCell"]))
control_ratio <- log2(mean(control_frame$VALUE[control_frame$score == "MANUAL_ANTI_INFLAMM_UCell"])/mean(control_frame$VALUE[control_frame$score == "MANUAL_PRO_INFLAMM_UCell"]))

curr_frame <- rbind(data.frame(TREAT = miniframe$TREAT[1], score = "RATIO", RATIO = treat_ratio, 
  mean_PRO = mean(treat_frame$VALUE[treat_frame$score == "MANUAL_PRO_INFLAMM_UCell"]),
  mean_ANTI = mean(treat_frame$VALUE[treat_frame$score == "MANUAL_ANTI_INFLAMM_UCell"]),
  type = "TREAT"),
  data.frame(TREAT = miniframe$TREAT[1], score = "RATIO", RATIO = control_ratio, type = "CONTROL",
  mean_PRO = mean(control_frame$VALUE[control_frame$score == "MANUAL_PRO_INFLAMM_UCell"]),
  mean_ANTI = mean(control_frame$VALUE[control_frame$score == "MANUAL_ANTI_INFLAMM_UCell"]))
)

return(curr_frame)
}

ratio_frame <- do.call("rbind", lapply(mini_split, ratio_maker))

fake_generator <- function(POS, NEG, out_score) {
##Quantify differences.
fake_set_list <- list()
for(treat in unique(row_collapse$TREAT)) {

  curr_row <- data.frame(TREAT = treat, score = out_score,
  mean_POS = row_collapse$curr_log_mean[intersect(which(row_collapse$TREAT == treat), which(row_collapse$score == POS))],
  mean_NEG = row_collapse$curr_log_mean[intersect(which(row_collapse$TREAT == treat), which(row_collapse$score == NEG))])
    curr_row$ADD <- curr_row$mean_POS + -(curr_row$mean_NEG)
    curr_row$pro_score <- curr_row$mean_POS
    curr_row$anti_score <- curr_row$mean_NEG

  fake_set_list[[treat]] <- curr_row
}

out_fake <- do.call("rbind", fake_set_list)
return(out_fake)
}

simple_frame <- fake_generator("MANUAL_PRO_INFLAMM_UCell", "MANUAL_ANTI_INFLAMM_UCell", "SIMPLE")

gen_plot <- function(frame, title) {
library(ggplot2)
library(forcats)
curr_plot <- ggplot(frame, aes(x = score, y = curr_log_mean, fill = fct_inorder(TREAT))) + geom_bar(position="dodge", stat="identity")# + geom_text(aes(label = fct_inorder(TREAT)), position = position_dodge(width = 0.9), vjust = -0.25, angle = 45)#geom_boxplot(aes(fill = TREAT), position = position_dodge(0.9))
curr_plot <- curr_plot + theme_classic() + xlab("") + ylab("Log2 Ratio") + theme(legend.title=element_blank()) + ggtitle(title) + theme(axis.text.x = element_text(size = 12))
##NEXT
curr_plot <- curr_plot + geom_text(aes(label = TREAT, y = 0, vjust = 0, angle = 90), position = position_dodge(width = .9), size = 8)
return(curr_plot)
}

library(ggplot2)
simple_frame$score <- simple_frame$ADD
simple_plot <- gen_plot(simple_frame, "Pro- and flipped anti-inflamm") + theme(legend.position = "none") + ylab("Inflammatory Response")
simple_plot <- simple_plot + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

simple_frame_POS <- simple_frame[, c("TREAT", "mean_POS")]
colnames(simple_frame_POS)[2] <- "score"
simple_frame_POS$type = "PRO-Inflammatory"
simple_frame_NEG <- simple_frame[, c("TREAT", "mean_NEG")]
colnames(simple_frame_NEG)[2] <- "score"
simple_frame_NEG$type = "ANTI-Inflammatory"
simple_reverse_frame <- simple_frame_NEG
simple_reverse_frame$score <- -simple_reverse_frame$score
stack_frame <- rbind(simple_frame_POS, simple_frame_NEG)
stack_frame_rev <- rbind(simple_frame_POS, simple_reverse_frame)

stack_allscores <- ggplot(stack_frame, aes(x = fct_inorder(TREAT), y = score, fill = type)) + geom_col()
title <- "UCell Scores - log2(Treatment/Control)"
stack_allscores <- stack_allscores + theme_classic() + xlab("") + ylab("Log2 Ratio") + theme(legend.title=element_blank()) + ggtitle(title) + theme(axis.text.x = element_text(size = 12))
##NEXT
stack_allscores <- stack_allscores + geom_text(aes(label = fct_inorder(TREAT), y = 0, vjust = 0, angle = 90),  size = 6) #position = position_dodge(width = .9),

stack_allscores <- stack_allscores + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pair_all_scores <- ggplot(stack_frame, aes(x = fct_inorder(TREAT), y = score, fill = type)) + geom_bar(position = "dodge", stat = "identity")
pair_all_scores <- pair_all_scores + scale_fill_manual(values = c("#7BB856", "red"))
title <- "UCell Scores - log2(Treatment/Control)"
pair_all_scores <- pair_all_scores  + theme_classic() + xlab("") + ylab("Log2 Ratio") + theme(legend.title=element_blank()) + ggtitle(title) + theme(axis.text.x = element_text(size = 12))
##NEXT
pair_all_scores <- pair_all_scores + geom_text(aes(label = fct_inorder(TREAT), y = 0, hjust = 0, vjust = 2, angle = 90),  size = 6) + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

stack_allscores_REV <- ggplot(stack_frame_rev, aes(x = fct_inorder(TREAT), y = score, fill = type)) + geom_col()
title <- "UCell Scores - log2(Treatment/Control), with anti-inflammatory genes flipped"
stack_allscores_REV <- stack_allscores_REV + theme_classic() + xlab("") + ylab("Log2 Ratio") + theme(legend.title=element_blank()) + ggtitle(title) + theme(axis.text.x = element_text(size = 12))
##NEXT
stack_allscores_REV <- stack_allscores_REV + geom_text(aes(label = fct_inorder(TREAT), y = 0, vjust = 0, angle = 90),  size = 6) #position = position_dodge(width = .9),

stack_allscores_REV <- stack_allscores_REV + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pair_all_scores_REV <- ggplot(stack_frame_rev, aes(x = fct_inorder(TREAT), y = score, fill = type)) + geom_bar(position = "dodge", stat = "identity")
pair_all_scores_REV <- pair_all_scores_REV + scale_fill_manual(values = c("#7BB856", "red"))
title <- "UCell Scores - log2(Treatment/Control), with anti-inflammatory genes flipped"
pair_all_scores_REV <- pair_all_scores_REV  + theme_classic() + xlab("") + ylab("Log2 Ratio") + theme(legend.title=element_blank()) + ggtitle(title) + theme(axis.text.x = element_text(size = 12))
##NEXT
pair_all_scores_REV <- pair_all_scores_REV + geom_text(aes(label = fct_inorder(TREAT), y = 0, hjust = 0, vjust = 2, angle = 270),  size = 4) + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

cell <- unique(mega_frame$CELL)

simple_plot  <- simple_plot + ggtitle(cell)
pair_all_scores  <- pair_all_scores + ggtitle(cell)
pair_all_scores_REV  <- pair_all_scores_REV + ggtitle(cell)

ratio_frame$TREAT <- factor(ratio_frame$TREAT, levels = treats_of_interest)
ratio_frame <- ratio_frame[order(ratio_frame$TREAT),]
curr_ratio_plot <- ggplot(ratio_frame, aes(x = fct_inorder(TREAT), y = RATIO, fill = type)) + geom_bar(position="dodge", stat="identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_classic() + xlab("") + ylab("Log2 Ratio (Anti/Pro)") + ggtitle("Ratio of Anti- to Pro-Inflammatory Genes") +
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
geom_text(aes(label = fct_inorder(TREAT), y = 0, vjust = 0.5, hjust = 1, angle = 90), size = 6)
curr_ratio_plot <- curr_ratio_plot + theme(legend.title=element_blank()) + ggtitle(cell)

ratio_cont <- ratio_frame[ratio_frame$TREAT == "IL10 Old - vs - Young Controls",]
ratio_young <- ratio_cont[ratio_cont$type == "CONTROL",]
ratio_old <- ratio_cont[ratio_cont$type == "TREAT",]

ratio_young$type <- "Young Control"
ratio_old$type <- "Old Control"
ratio_frame <- ratio_frame[ratio_frame$type != "CONTROL",]
ratio_treat <- ratio_frame[ratio_frame$TREAT != "IL10 Old - vs - Young Controls",]
ratio_treat$type <- ratio_treat$TREAT

ratio_simple <- rbind(ratio_young, ratio_old, ratio_treat)
ratio_simple$type <- factor(ratio_simple$type)
curr_ratio_plot_simple <- ggplot(ratio_simple, aes(x = fct_inorder(type), y = RATIO, fill = type)) + geom_bar(stat="identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_classic() + xlab("") + ylab("Log2 Ratio (Anti/Pro)") + ggtitle("Ratio of Anti- to Pro-Inflammatory Genes") +
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
geom_text(aes(label = fct_inorder(type), y = 0, vjust = 0.5, hjust = 1, angle = 90), size = 6)
curr_ratio_plot_simple <- curr_ratio_plot_simple + theme(legend.title=element_blank()) + ggtitle(cell)

plot_set <- list(simple_plot, pair_all_scores, pair_all_scores_REV, curr_ratio_plot, curr_ratio_plot_simple)

cellwise_plots[[cell]] <- plot_set

out_frame <- data.frame(simple_frame, CELL = cell)

cellwise_frames[[cell]] <- out_frame
print('finished')
print(cell)
}

##cellwise plotting
simple_plot_res <- lapply(cellwise_plots, function(x) x[[1]])
cells_of_interest <- c("Microglia",  "Endothelial", "Astrocyte_qNSC", "aNSCs_NPC", "Neuroblast", "Oligodendro")
simple_plot_res <- lapply(cells_of_interest, function(x) simple_plot_res[[x]])

pair_plot_res <- lapply(cellwise_plots, function(x) x[[2]] + theme(legend.position = "none"))
pair_plot_res <- lapply(cells_of_interest,function(x) pair_plot_res[[x]])

pair_plot_resN <- lapply(cellwise_plots, function(x) x[[3]] + theme(legend.position = "none") )
pair_plot_resN <- lapply(cells_of_interest,function(x) pair_plot_resN[[x]])

stack_plot_res <- lapply(cellwise_plots, function(x) x[[4]] + theme(legend.position = "none"))
stack_plot_res <- lapply(cells_of_interest,function(x) stack_plot_res[[x]])

stack_plot_res_SIMP <- lapply(cellwise_plots, function(x) x[[5]] + theme(legend.position = "none"))
stack_plot_res_SIMP <- lapply(cells_of_interest,function(x) stack_plot_res_SIMP[[x]])

library(gridExtra)
if(YOUNG_EVERYTHING) {
pdf("MULTICELL_inflammation_signaling_results_relative_YOUNG_CONTROL.pdf", width = 22, height = 16)
}else{
pdf("MULTICELL_inflammation_signaling_results_relative_OLD_CONTROL_INVERTED_JULY2025.pdf", width = 22, height = 16)
}
print(grid.arrange(grobs = simple_plot_res, nrow = 2))
print(grid.arrange(grobs = pair_plot_res, nrow = 2))
print(grid.arrange(grobs = pair_plot_resN, nrow = 2))
print(grid.arrange(grobs = stack_plot_res, nrow = 2))
print(grid.arrange(grobs = stack_plot_res_SIMP, nrow = 2))
dev.off()

if(YOUNG_EVERYTHING) {
pdf("MICROGLIA_inflammation_signaling_results_relative_YOUNG_CONTROL.pdf", width = 8, height = 6)
}else{
pdf("MICROGLIA_inflammation_signaling_results_relative_OLD_CONTROL_INVERTED_JULY2025.pdf", width = 8, height = 8)
}
print(cellwise_plots$Microglia[[2]])
print(cellwise_plots$Microglia[[3]])
print(cellwise_plots$Microglia[[4]])
print(cellwise_plots$Microglia[[5]])
dev.off()


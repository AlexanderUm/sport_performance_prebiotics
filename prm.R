################################################################################
# List of parameter that will be used 
################################################################################
prm.ls <- list()

dir.create("out/supp", recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------------------------------------
# General 
#-------------------------------------------------------------------------------
prm.ls[["general"]] <- list("seed" = 4397547, 
                            "Group_cols" = c("Treatment", "Resp"), 
                            "Time_col_num" = "Time", 
                            "Time_col_fact" = "TestDay", 
                            "Prat_ID_col" = "PersonID", 
                            "Preb_col" = "Prebiotic", 
                            "n_core" = 4) 


#-------------------------------------------------------------------------------
# Data preparation (0_data_prep.R)
#-------------------------------------------------------------------------------
prm.ls[["data_prep"]] <- list("seq_path" = "data/seq",
                              "meta_path" = "data/sp_metadata_Britt.csv", 
                              "samp_to_remove" = c("R04", "R06", "R33"), # COVID
                              "remove_covid_from_phyloseq" = FALSE,
                              "min_reads_tax" = 10, 
                              "rare_depth" = 8450, 
                              "tax_king_keep" = c("d__Bacteria", "d__Archaea"), 
                              "tax_kick" = c("Genus" = "Mitochondria"), 
                              "tax_levels" = c("ASV", "Genus"))


#-------------------------------------------------------------------------------
# Primary outcomes 
#-------------------------------------------------------------------------------
prm.ls[["primary"]] <- list("lmm_form" = "TTE ~ Age + Sex + Time*Treatment + (1|PersonID)",
                            "emmeans_form" = "pairwise ~ Time|Treatment",
                            "samp_to_remove" = c("R04", "R06", "R33"), # COVID
                            "resp_var" = "TTE", 
                            "gr_var" = "Treatment", 
                            "out_dir_path" = "out/primary")

#-------------------------------------------------------------------------------
# Alpha 
#-------------------------------------------------------------------------------
prm.ls[["alpha"]] <- list("Tax_lvl" = "ASV",
                          "Norm" = "Rare",
                          "measures" = c("Observed", "Shannon",
                                         "InvSimpson", "PhyloDiverity"), 
                         "out_dir_path" = "out/alpha")


#-------------------------------------------------------------------------------
# Beta
#-------------------------------------------------------------------------------
prm.ls[["beta"]] <- list("Tax_lvl" = c("ASV", "Genus"),
                         "Norm" = "CSS",
                         "distances" = c("Unweighted UniFrac" = "unifrac", 
                                         "Weighted UniFrac" = "wunifrac", 
                                         "Jaccard" = "jaccard", 
                                         "Bray-Curtis" = "bray"),
                         "n_perm" = 999, 
                         "out_dir_path" = "out/beta")


#-------------------------------------------------------------------------------
# DA
#-------------------------------------------------------------------------------
prm.ls[["DA"]] <- list("Tax_lvl" = c("ASV", "Genus"),
                         "Norm" = c("CSS"),
                         "Pval_cutoff" = 0, # 0 means will not be used for filtering
                         "Qval_cutoff" = 0.1, 
                         "Min_prev" = 0.25, 
                         "Maas_fix_effect" = "Time", 
                         "Maas_rand_effect" = "PersonID", 
                         "Maas_method" = "LM", 
                         "Maas_pcorr_method" = "BH", 
                         "Maas_norm" = "NONE",
                         "Maas_tansform" = "NONE",
                         "heat_width" = 12,
                         "out_dir_path" = "out/DA")


#-------------------------------------------------------------------------------
# Correlations
#-------------------------------------------------------------------------------
# Correlation with taxa delta abundance
prm.ls[["Corr"]] <- list("Tax_lvl" = c("ASV", "Genus"),
                         "Norm" = c("CSS"), 
                         "Min_prev" = 0.35, 
                         "Corr_cols" = c("Delta_TTE"), 
                         "Samples_gr" =  list("Treatment" = list(c("CON", "PRE"), 
                                                                   "PRE", "CON")), 
                         "out_dir_path" = "out/Correlations", 
                         "Method" = "spearman", 
                         "Plot_top_n" = 10, 
                         "Plot_per_row" = 5)


#-------------------------------------------------------------------------------
lib.to.load <- c("phyloseq", "tidyverse", "qiime2R", 
                 "metagenomeSeq", "RColorBrewer", "broom", 
                 "ggsignif", "vegan", "usedist", "cowplot", "ape", 
                 "Maaslin2", "ComplexHeatmap", "circlize", "lmerTest", 
                 "emmeans")

for (i in lib.to.load) {library(i, character.only = TRUE)}



#-------------------------------------------------------------------------------
# Write out parameters file 
#-------------------------------------------------------------------------------
save(list = c("prm.ls"), file = "out/supp/prm.R")

rm(list = ls())

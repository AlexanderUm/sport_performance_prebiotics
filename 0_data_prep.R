#-------------------------------------------------------------------------------
# Parameters 
#-------------------------------------------------------------------------------
load("out/supp/prm.R")

seq.path <- prm.ls[["data_prep"]][["seq_path"]]

meta.path <- prm.ls[["data_prep"]][["meta_path"]]

samp.to.remove <- prm.ls[["data_prep"]][["samp_to_remove"]]

remove.covid.from.phyloseq <- prm.ls[["data_prep"]][["remove_covid_from_phyloseq"]]

min.reads <- prm.ls[["data_prep"]][["min_reads_tax"]]

rare.depth <- prm.ls[["data_prep"]][["rare_depth"]]

tax.king.keep <- prm.ls[["data_prep"]][["tax_king_keep"]]

tax.kick <- prm.ls[["data_prep"]][["tax_kick"]]

tax.lvls <- prm.ls[["data_prep"]][["tax_levels"]]

r.seed <- prm.ls[["general"]][["seed"]]


#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
source("R/phy_shorten_tax_names.R")
source("R/phy_norm.R")
source("R/phy_prevalence_calc.R")

set.seed(r.seed)


prev.tabs.ls <- list()
#-------------------------------------------------------------------------------
# Data import
#-------------------------------------------------------------------------------
# Import data into phyloseq 
ps.0 <- qza_to_phyloseq(features = paste0(seq.path, "/asv_table.qza"), 
                        tree = paste0(seq.path, "/tree/rooted-tree.qza"), 
                        taxonomy = paste0(seq.path, "/taxonomy_07.qza"))

meta.data <- read.csv(meta.path) %>% 
                  mutate(SeqID = as.character(X)) %>% 
                  rename(Sex = M_F, PersonID = Rnr) %>% 
                  mutate(Time = as.numeric(case_match(Timepoint, 
                                                      "TD1" ~ 0, 
                                                      "TD3" ~ 6)), 
                         Treatment = factor(case_match(Group, 
                                                       "GR1" ~ "CON", 
                                                       "GR2" ~ "PRE"), 
                                            levels = c("PRE", "CON")), 
                         Prebiotic = factor(Prebiotic, 
                                            levels = c("Placebo", "FOS", "RS", 
                                                       "GOS", "Inulin")), 
                         TestDay = factor(Timepoint, 
                                          levels = c("TD1", "TD3")), 
                         PersonID = as.factor(PersonID), 
                         Sex = as.factor(Sex), 
                         TTE_Change = ifelse(Delta_TTE == 0, 
                                             "No_Change", 
                                             ifelse(Delta_TTE > 0, 
                                                     "Increase", "Decrease"))) %>% 
                         mutate(TreatmentTestDay = as.factor(paste0(Treatment, 
                                                                ":", TestDay)),
                                TTE_log = log(TTE + 1), 
                                Resp = factor(ifelse(Delta_TTE >= 180, 
                                              "RESP", "NONE"), 
                                              levels = c("RESP", "NONE"))) %>% 
                         mutate(RespTestDay = as.factor(paste0(Resp, ":", TestDay)))


# Add metadata to phyloseq object
ps.meta <- meta.data %>% 
                filter(SeqID %in% intersect(meta.data$SeqID, 
                                  sample_names(ps.0))) %>% 
                        group_by(PersonID) %>%
                        filter(n() == 2) %>%
                        ungroup() %>% 
                        as.data.frame() %>% 
                        droplevels() %>% 
                        mutate(RowNAMES = SeqID) %>% 
                        column_to_rownames("RowNAMES")

if(remove.covid.from.phyloseq) {
  
  ps.meta <- meta.data %>% 
                filter(!PersonID %in% samp.to.remove) %>% 
                droplevels()
}

ps1 <- prune_samples(rownames(ps.meta), ps.0)

ps.meta <- ps.meta[sample_names(ps1), ]

sample_data(ps1) <- ps.meta
                

#-------------------------------------------------------------------------------
# Prune taxa 
#-------------------------------------------------------------------------------
# Taxa with less than 10 reads in total   
ps1 <- prune_taxa(taxa_sums(ps.0) > min.reads, ps1)

# Keep ASVs
ps1 <- prune_taxa(tax_table(ps1)[, "Kingdom"] %in% tax.king.keep, ps1)

# Kick ASVs
for(i.kick in 1:length(tax.kick)) {
  
  kick.inst <- tax.kick[i.kick]
  
  ps1 <- prune_taxa(!tax_table(ps1)[, names(kick.inst)] %in% kick.inst, ps1)
  
}


#-------------------------------------------------------------------------------
# Tax glom and count tranformation
#-------------------------------------------------------------------------------
ps.ls <- list()

for(i.lvl in tax.lvls)  {
  
  if(i.lvl == "ASV") {ps.inst <- ps1
  
  # Adjust taxa names
  taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% 
                              as.data.frame() %>% 
                              setNames("feature") %>% 
                              mutate(feature2 = if(n() > 1) {
                                 paste0(feature, "__asv", row_number())} else {
                                  paste0(feature, "__asv")}, .by = "feature") %>% 
                              pull(feature2)
  
  } else { 
    
    ps.inst <- tax_glom(ps1, i.lvl, NArm = FALSE) 
  
    # Adjust taxa names
    taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% 
                                as.data.frame() %>% 
                                setNames("feature") %>% 
                                mutate(feature2 = if(n() > 1) {
                                  paste0(feature, "__", 
                                         tolower(str_sub(i.lvl, 1, 1)), 
                                         row_number())} else {
                                  feature}, .by = "feature") %>% 
                                pull(feature2)
  }
  
 
  
  # Calculate prevalence per taxa 
  prev.tabs.ls[[i.lvl]] <- phy_prevalence_calc(ps.inst, 
                                    per_group_cols = c("Treatment", 
                                                       "TestDay", 
                                                       "Resp", 
                                                       "RespTestDay"))
  
  # Transform count 
  ps.rare <- rarefy_even_depth(ps.inst, 
                               rngseed = r.seed)
  
  ps.rare.log <- transform_sample_counts(ps.rare, function(x) {log2(x + 1)})
  
  ps.ls[[i.lvl]] <- list("Raw" = ps.inst, 
                         "Rare" = ps.rare, 
                         "Rare_log2" = ps.rare.log,
                         "CSS" = phy_css_norm(ps.inst), 
                         "CLR" = phy_clr_norm(ps.inst), 
                         "TSS" = phy_tss_norm(ps.inst, log2_transform = FALSE))
}


#-------------------------------------------------------------------------------
# Color schema 
#-------------------------------------------------------------------------------
# Groups color
color.sch.ls <- list("TestDay" = brewer.pal(7, "Dark2")[c(1,3)] %>% 
                                    setNames(levels(meta.data$TestDay)), 
                     "Treatment" = c(brewer.pal(7, "Set1")[c(2)], "gray35") %>% 
                                    setNames(levels(meta.data$Treatment)),
                     "Prebiotic" = c("gray35", brewer.pal(4, "Set1")) %>% 
                                    setNames(levels(meta.data$Prebiotic)),
                     "TreatmentTestDay" = brewer.pal(7, "Accent")[c(1,2,3,7)] %>% 
                                    setNames(levels(meta.data$TreatmentDay)), 
                     "Resp" = c("forestgreen","gray35") %>% 
                                    setNames(levels(meta.data$Resp)), 
                     "RespTestDay" = brewer.pal(7, "Accent")[c(1,2,3,7)] %>% 
                                        setNames(levels(meta.data$RespTestDay)))

#-------------------------------------------------------------------------------
# Write data
#-------------------------------------------------------------------------------
dir.create("out/supp/")

save(list = c("ps.ls", "ps.meta", "meta.data", 
              "color.sch.ls", "prev.tabs.ls"),
     file = "out/supp/0_data.RData")

rm(list = ls())

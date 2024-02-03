#-------------------------------------------------------------------------------
set.seed(549794)

lib.to.load <- c("phyloseq", "tidyverse", "qiime2R", 
                 "metagenomeSeq", "RColorBrewer")

for (i in lib.to.load) {library(i, character.only = TRUE)}

source("R/phy_shorten_tax_names.R")

#------------------------------------------------------------------------------
# Data import
#------------------------------------------------------------------------------
# Import data into phyloseq 
ps.0 <- qza_to_phyloseq(features = "data/asv_table.qza", 
                        tree = "data/tree/rooted-tree.qza", 
                        taxonomy = "data/taxonomy_07.qza")


meta.data <- read.csv("data/sp_metadata_Britt.csv") %>% 
                  mutate(SeqID = as.character(X)) %>% 
                  rename(Sex = M_F, PersonID = Rnr) %>% 
                  mutate(Time = as.numeric(case_match(Timepoint, 
                                                      "TD1" ~ 0, 
                                                      "TD3" ~ 6)), 
                         Treatment = factor(case_match(Group, 
                                                       "GR1" ~ "CON", 
                                                       "GR2" ~ "PRE"), 
                                            levels = c("CON", "PRE")), 
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
                         mutate(TreatmentDay = as.factor(paste0(Treatment, 
                                                                ":", TestDay)), 
                                TimeGroupInter = interaction(Time, Treatment),
                                TTE_log = log(TTE + 1), 
                                Resp = factor(ifelse(Delta_TTE >= 180, 
                                              "RESP", "NONE"), 
                                              levels = c("NONE", "RESP")))


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

ps1 <- prune_samples(rownames(ps.meta), ps.0)

ps.meta <- ps.meta[sample_names(ps1), ]

sample_data(ps1) <- ps.meta
                

#-------------------------------------------------------------------------------
# Filter out very rare (artifacts) and contaminating taxa
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Prune taxa 
#-------------------------------------------------------------------------------
# Taxa with less than 10 reads in total   
ps1 <- prune_taxa(taxa_sums(ps.0) > 10, ps1)

# Remove ASVs: 
# Kingdom: "d__Eukaryota", "Unassigned"
# Genus: "Mitochondria"
ps1 <- prune_taxa(!tax_table(ps1)[, "Genus"] %in% "Mitochondria", ps1)

ps1 <- prune_taxa(!tax_table(ps1)[, "Kingdom"] %in% c("d__Eukaryota", 
                                                      "Unassigned"), ps1)

#-------------------------------------------------------------------------------
# Make telling names 
#-------------------------------------------------------------------------------
taxa_names(ps1) <- phy_shorten_tax_names(ps1) %>% 
                          gsub("_group", "", .) %>% 
                          paste0(., "_asv", 1:length(taxa_names(ps1)))


# Genus level
ps1.genus <- tax_glom(ps1, "Genus")

taxa_names(ps1.genus) <- phy_shorten_tax_names(ps1.genus) %>% 
                            gsub("_group", "", .) %>% 
                            make.unique()

# Genus level
ps1.family <- tax_glom(ps1, "Family")

taxa_names(ps1.family) <- phy_shorten_tax_names(ps1.family) %>% 
                                                gsub("_group", "", .) %>% 
                                                make.unique()

#-------------------------------------------------------------------------------
# Phylseqs list
#-------------------------------------------------------------------------------
# CSS - ASV
ps1.css.asv <- ps1

otu_table(ps1.css.asv) <- ps1 %>% 
                            phyloseq_to_metagenomeSeq(.) %>% 
                            cumNorm(., p=cumNormStatFast(.)) %>% 
                            MRcounts(., norm=TRUE, log=TRUE) %>% 
                            as.data.frame() %>% 
                            otu_table(., taxa_are_rows = TRUE)


# CSS - Genus
ps1.css.genus <- ps1.genus

otu_table(ps1.css.genus) <- ps1.genus %>%
                                phyloseq_to_metagenomeSeq(.) %>%
                                cumNorm(., p=cumNormStatFast(.)) %>%
                                MRcounts(., norm=TRUE, log=TRUE) %>%
                                as.data.frame() %>%
                                otu_table(., taxa_are_rows = TRUE)

# Normalize count
ps1.rare <- rarefy_even_depth(ps1, rngseed = 0745048)

ps1.rare.genus <- rarefy_even_depth(ps1.genus, rngseed = 93475)

#-------------------------------------------------------------------------------
# Log transform
#-------------------------------------------------------------------------------
ps1.log.rare <- transform_sample_counts(ps1.rare, function(x) log(x + 1))

ps1.log.rare.genus <- transform_sample_counts(ps1.rare.genus, function(x) log(x + 1))


# Combine into a list
pss.ls <- list(raw_asv = ps1, 
               raw_genus = ps1.genus, 
               raw_family = ps1.family,
               css_asv = ps1.css.asv, 
               css_genus = ps1.css.genus, 
               rare_asv = ps1.rare, 
               rare_genus = ps1.rare.genus, 
               rare_log_asv = ps1.log.rare,
               rare_log_genus = ps1.log.rare.genus)

#-------------------------------------------------------------------------------
# Color schema 
#-------------------------------------------------------------------------------
# Groups color
color.sch.ls <- list("TestDay" = brewer.pal(7, "Dark2")[c(1,3)] %>% 
                                    setNames(levels(meta.data$TestDay)), 
                     "Treatment" = c("gray35", brewer.pal(7, "Set1")[c(2)]) %>% 
                                    setNames(levels(meta.data$Treatment)),
                     "Prebiotic" = c("gray35", brewer.pal(4, "Set1")) %>% 
                                    setNames(levels(meta.data$Prebiotic)),
                     "TreatmentDay" = brewer.pal(7, "Accent")[c(1,2,3,7)] %>% 
                                    setNames(levels(meta.data$TreatmentDay)), 
                     "Resp" = c("gray35", "forestgreen") %>% 
                                    setNames(levels(meta.data$Resp)))

#-------------------------------------------------------------------------------
# Write data
#-------------------------------------------------------------------------------
dir.create("out/supp/")

save(list = c("pss.ls", "ps.meta", "meta.data", "color.sch.ls"),
     file = "out/supp/0_data.RData")

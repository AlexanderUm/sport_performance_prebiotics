#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(549794)

lib.to.load <- c("phyloseq", "tidyverse", "ggsignif", "broom")

for (i in lib.to.load) {library(i, character.only = TRUE)}

rm(list = c("i", "lib.to.load"))

# Load data
load("out/supp/0_data.RData")


#-------------------------------------------------------------------------------
# Prepare data and calculate alpha diversity
#-------------------------------------------------------------------------------
ps.inst <- pss.ls$raw_asv

used.var <- c("TestDay", "Treatment", "PersonID", 
              "Prebiotic", "TTE_Change", "Resp")

div.measures <- c("Observed", "Shannon", "InvSimpson", "PhyloDiverity")

group.vars <- c("Treatment", "Resp")




# Calculate diversity indexes 
alpha.df <- estimate_richness(physeq = ps.inst, 
                              measures = div.measures)

if("PhyloDiverity" %in% div.measures) {
# Calculate Phylogeny Diversity (package "picante")
alpha.df$PhyloDiverity <- picante::pd(samp = t(otu_table(ps.inst)), 
                                      tree = phy_tree(ps.inst), 
                                      include.root = FALSE)  %>% 
                                      select("PD")  %>%  
                                      unlist()
}

# Add relevant columns 
alpha.df.long <- ps.meta  %>%  
                      select(all_of(used.var)) %>% 
                      bind_cols(alpha.df)  %>% 
                      gather(key = "Diversity_Index", 
                             value = "Value", 
                             -all_of(used.var)) 


#-------------------------------------------------------------------------------
# Pairwise comparison of Alpha diversity per group at TD1 vs TD3
#-------------------------------------------------------------------------------
# 1. Prepare data 
# Remove samples without a pairs > sort samples
alpha.df <- alpha.df.long %>% 
                  group_by(across(c("Diversity_Index", 
                                    "PersonID", 
                                    ))) %>% 
                  filter(n() == 2)  
                

# 2. Test (Paired Wilcoxon)
wil.res.ls <- list()

for(i.gr.col in group.vars) {
  
  alpha.wilc.res <- NULL
  
  for (i in unique(alpha.df$Diversity_Index)) {
  
  for (i1 in unique(alpha.df[[i.gr.col]])) {
    
    wil.df <- alpha.df %>% 
                filter(Diversity_Index == i, 
                       .data[[i.gr.col]] == i1) %>%
                group_by(TestDay) %>%
                arrange(PersonID, .by_group = TRUE) %>%
                ungroup()
    
    alpha.wilc.res <- wilcox.test(Value ~ TestDay, 
                                  paired = TRUE, 
                                  data = wil.df) %>%
              tidy() %>%
              select(c("p.value")) %>%
              mutate(Comparison = "TD1_vs_TD3",
                     !!i.gr.col := i1,
                     DivIndex = i) %>%
              rbind(alpha.wilc.res, .)
    }

  }
  
  wil.res.ls[[i.gr.col]] <- alpha.wilc.res
  
}


#-------------------------------------------------------------------------------
# Plot Alpha diversity 
#-------------------------------------------------------------------------------
# Make data frame for significance levels 
alpha.p.ls <- list()

for(i.gr.col in group.vars) { 
  
  wrap.form <- paste0("Diversity_Index ~ ", i.gr.col)
  
  max.div <- alpha.df %>% 
                  group_by(across(c(i.gr.col, "Diversity_Index"))) %>% 
                  slice(which.max(Value)) %>% 
                  mutate(Id = paste0(Diversity_Index, 
                                     "_", 
                                     .data[[i.gr.col]]), 
                         End = if_else(TestDay == "TD1", "TD3", "TD1"), 
                         y = Value*1.1)
  
  sig.df <- wil.res.ls[[i.gr.col]] %>% 
                mutate(Id = paste0(DivIndex, "_", .data[[i.gr.col]]), 
                       p.short = paste0("p=", round(p.value, 3))) %>% 
                filter(p.value <= 0.05) %>% 
                left_join(., max.div)
  
  # Plot alpha diversity 
  alpha.plot <- ggplot(alpha.df, aes(y = Value, x = TestDay)) + 
                    geom_point(aes(color = Prebiotic), 
                               size = 2, 
                               alpha = 1) +
                    geom_line(aes(group = PersonID, 
                                  color = Prebiotic)) +
                    geom_violin(fill = NA, alpha = 0.1) +
                    geom_signif(data = sig.df,
                                aes(xmin = TestDay,
                                    xmax = End,
                                    annotations = p.short,
                                    y_position = y),
                                textsize = 4, vjust = -0.1,
                                manual = TRUE, margin_top = 1) +
                    geom_point(data = sig.df,
                               aes(x = End, y = y*1.175), x=NA) +
                    facet_grid(as.formula(wrap.form), scales = "free") + 
                    theme_bw() + 
                    scale_color_manual(values = color.sch.ls$Prebiotic) + 
                    xlab("Test Day") +
                    theme(axis.title.x = element_text()) + 
                    guides(colour = guide_legend(order = 2), 
                           linetype = guide_legend(order = 1))

  alpha.p.ls[[i.gr.col]] <- alpha.plot
  
}
  

res.alpha <- list(Plot = alpha.p.ls, 
                  Statistics = wil.res.ls)

save(list = "res.alpha",
     file = "out/supp/2_alpha_res.RData")

#-------------------------------------------------------------------------------
# Load data and functions
#-------------------------------------------------------------------------------
load("out/supp/prm.R")
load("out/supp/0_data.RData")

source("R/phy_alpha.R")

#-------------------------------------------------------------------------------
# Parameters 
#-------------------------------------------------------------------------------
gr.cols <- prm.ls[["general"]][["Group_cols"]]

id.col <- prm.ls[["general"]][["Prat_ID_col"]]

preb.col <- prm.ls[["general"]][["Preb_col"]]

part.id.col <- prm.ls[["general"]][["Prat_ID_col"]]

seed <- prm.ls[["general"]][["seed"]]

time.col.fact <- prm.ls[["general"]][["Time_col_fact"]]

tax.lvl <- prm.ls[["alpha"]][["Tax_lvl"]]

norm <- prm.ls[["alpha"]][["Norm"]]

alpha.measur <- prm.ls[["alpha"]][["measures"]]

dir.out <- prm.ls[["alpha"]][["out_dir_path"]]

set.seed(seed)


################################################################################
#-------------------------------------------------------------------------------
# Calculate alpha diversity
#-------------------------------------------------------------------------------
alpha.res.ls <- list()

alpha.par.grid <- expand.grid("Tax_lvl" = tax.lvl, 
                              "Normalization" = norm, 
                              stringsAsFactors = FALSE)


for(i.prm in 1:nrow(alpha.par.grid)) {
  
  lvl.inst <- alpha.par.grid[i.prm, "Tax_lvl"]
  
  norm.inst <- alpha.par.grid[i.prm, "Normalization"]
  
  name.inst <- paste0(lvl.inst, "--", norm.inst)
  
  ps.inst <- ps.ls[[lvl.inst]][[norm.inst]]
  
  alpha.tab.inst <- phy_alpha(ps.inst, measures = alpha.measur)
  
  alpha.res.ls[["index"]][[name.inst]] <- alpha.tab.inst
  
  # Write results 
  dir.create(paste0(dir.out, "/tabs/"), recursive = TRUE, showWarnings = FALSE)
  
  write.csv(x = alpha.tab.inst, 
            file = paste0(dir.out, "/tabs/", name.inst, "--diversity.csv"))
  
}


#-------------------------------------------------------------------------------
# Pairwise comparison of Alpha diversity per group at TD1 vs TD3
#-------------------------------------------------------------------------------
comp.name <- paste(levels(ps.meta[[time.col.fact]]), collapse = "--vs--")

for(i.gr in gr.cols) {
  
  tests.grid <- expand.grid("Ind" = alpha.measur, 
                          "Group" = levels(ps.meta[[i.gr]]), 
                          stringsAsFactors = FALSE)

  for(i.prm.set in names(alpha.res.ls[["index"]])) {
    
    ind.tab.inst <- alpha.res.ls[["index"]][[i.prm.set]] %>% 
                      cbind(., ps.meta)
    
    wil.res.comb.df <- NULL
    
    for(i.test in 1:nrow(tests.grid)) {
      
      ind.inst <- tests.grid[[i.test, "Ind"]]
      
      gr.lvl.inst <- tests.grid[[i.test, "Group"]]
      
      form.inst <- paste0(ind.inst, "~", time.col.fact)
      
      # Prepare data 
      ind.tab.inst.f <- ind.tab.inst %>% 
                            dplyr::filter(.data[[i.gr]] == gr.lvl.inst) %>% 
                            filter(n() == 2, .by = all_of(id.col)) %>% 
                            arrange(across(all_of(c(time.col.fact, id.col)))) %>% 
                            droplevels()
      
      # Wilcoxon test 
      wil.res.comb.df <- wilcox.test(as.formula(form.inst), 
                                    paired = TRUE, 
                                    data = ind.tab.inst.f, 
                                    exact = FALSE) %>% 
                          tidy() %>%
                          select(c("p.value")) %>%
                          mutate(!!i.gr := gr.lvl.inst,
                                 DivIndex = ind.inst) %>%
                          rbind(wil.res.comb.df, .)
    }
    
    wil.res.comb.df <- wil.res.comb.df %>% 
                          mutate(Comparision = comp.name,
                                 Test = "Paired-Wilcoxon")
    
    alpha.res.ls[["Stat"]][[i.gr]][[i.prm.set]] <- wil.res.comb.df
    
    # Write results 
    dir.create(paste0(dir.out, "/tabs/"), recursive = TRUE, showWarnings = FALSE)
    
    write.csv(x = wil.res.comb.df, 
              file = paste0(dir.out, "/tabs/", i.gr, "--", 
                            i.prm.set, "--stat.csv"))
    
  }

}


#-------------------------------------------------------------------------------
# Plot Alpha diversity 
#-------------------------------------------------------------------------------
for(i.set in names(alpha.res.ls[["index"]])) {
  
  alpha.tab.l.inst <- alpha.res.ls[["index"]][[i.set]] %>% 
                        cbind(ps.meta[, c(gr.cols, 
                                          time.col.fact,
                                          preb.col, 
                                          part.id.col)]) %>% 
                        pivot_longer(cols = all_of(alpha.measur), 
                                     names_to = "Index")
  
  for(i.gr.col in gr.cols) {
    
     wrap.form <- paste0("Index ~ ", i.gr.col)
    
     max.div <- alpha.tab.l.inst %>% 
                    group_by(across(all_of(c(i.gr.col, "Index")))) %>% 
                    slice(which.max(value)) %>% 
                    mutate(Id = paste0(Index, 
                                       "_", 
                                       .data[[i.gr.col]]), 
                           Start = levels(ps.meta[[time.col.fact]])[1],
                           End = levels(ps.meta[[time.col.fact]])[2], 
                           y = value*1.1)
     
     sig.df <- alpha.res.ls[["Stat"]][[i.gr.col]][[i.set]] %>% 
                   mutate(Id = paste0(DivIndex, "_", .data[[i.gr.col]]), 
                          p.short = ifelse(round(p.value, 3) == 0, 
                                           "p<0.001", 
                                           paste0("p=", round(p.value, 3)))) %>% 
                   filter(p.value <= 0.05) %>% 
                   left_join(., max.div)
     
     # Plot alpha diversity 
     alpha.plot <- ggplot(alpha.tab.l.inst, 
                              aes(y = value, 
                                  x = .data[[time.col.fact]])) + 
                           geom_point(aes(color = .data[[preb.col]]), 
                                      size = 2, 
                                      alpha = 1) +
                           geom_line(aes(group = .data[[part.id.col]], 
                                         color = .data[[preb.col]])) +
                           geom_violin(fill = NA, alpha = 0.1) +
                           geom_signif(data = sig.df,
                                       aes(xmin = Start,
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
                           theme(axis.title.x = element_text(), 
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank()) + 
                           guides(colour = guide_legend(order = 2), 
                                  linetype = guide_legend(order = 1))
     
     
   alpha.res.ls[["Plots"]][[i.gr.col]][[i.set]] <- alpha.plot
   
   dir.create(paste0(dir.out, "/plots/"), recursive = TRUE, showWarnings = FALSE)
   
   ggsave(paste0(dir.out, "/plots/", i.gr.col, "--", i.set, ".png"), 
          plot = alpha.plot, 
          width = 4.25, 
          height = 1.25*length(alpha.measur) + 1, 
          dpi = 600)
    
  }
  
}

################################################################################
save(list = "alpha.res.ls",
     file = "out/supp/2_alpha_res.RData")

rm(list = ls())

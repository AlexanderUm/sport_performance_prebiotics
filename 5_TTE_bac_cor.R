#-------------------------------------------------------------------------------
# Load data and functions
#-------------------------------------------------------------------------------
load("out/supp/prm.R")
load("out/supp/0_data.RData")

#-------------------------------------------------------------------------------
# Parameters 
#-------------------------------------------------------------------------------
gr.cols <- prm.ls[["general"]][["Group_cols"]]

id.col <- prm.ls[["general"]][["Prat_ID_col"]]

preb.col <- prm.ls[["general"]][["Preb_col"]]

part.id.col <- prm.ls[["general"]][["Prat_ID_col"]]

seed <- prm.ls[["general"]][["seed"]]

n.core <- prm.ls[["general"]][["n_core"]]

time.col.fact <- prm.ls[["general"]][["Time_col_fact"]]

time.col.num <- prm.ls[["general"]][["Time_col_num"]]

tax.lvl <- prm.ls[["Corr"]][["Tax_lvl"]]

norm <- prm.ls[["Corr"]][["Norm"]]

sample.sets <- prm.ls[["Corr"]][["Samples_gr"]]

corr.cols <- prm.ls[["Corr"]][["Corr_cols"]]

min.prev <- prm.ls[["Corr"]][["Min_prev"]]

out.folder <- prm.ls[["Corr"]][["out_dir_path"]]

top.n.plot <- prm.ls[["Corr"]][["Plot_top_n"]]

plots.n.col <- prm.ls[["Corr"]][["Plot_per_row"]]


set.seed(seed)

dir.create(paste0(out.folder, "/plots/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(out.folder, "/tabs/"), showWarnings = FALSE, recursive = TRUE)

res.corr.comb <- list()

#-------------------------------------------------------------------------------
# Correlation between delta abundance and delta TTE
#-------------------------------------------------------------------------------
prm.cor.grid <- expand.grid("Taxa_lvl" = tax.lvl, 
                            "Normalization" = norm, 
                            stringsAsFactors = FALSE)

for(i.prm in 1:nrow(prm.cor.grid)) { 
  
  lvl.inst <- prm.cor.grid[i.prm, "Taxa_lvl"]
  
  norm.inst <- prm.cor.grid[i.prm, "Normalization"]
  
  # Extract otu table
  otu.inst <- ps.ls[[lvl.inst]][[norm.inst]] %>% 
                    otu_table() %>% 
                    as.matrix() %>% 
                    t() %>% 
                    as.data.frame() 
  
  # Add meta and calculate delta abundance between time point one and two
  otu.inst.delta <- otu.inst %>% 
                        bind_cols(., ps.meta[, c(id.col, 
                                                 time.col.fact, 
                                                 names(sample.sets), 
                                                 corr.cols)]) %>%
                        group_by(.data[[id.col]]) %>%
                        arrange(.data[[time.col.fact]], .by_group = TRUE) %>% 
                        mutate(across(all_of(colnames(otu.inst)), 
                                      function(x){x[2] - x[1]})) %>% 
                        ungroup() %>% 
                        filter(.data[[time.col.fact]] == levels(.data[[time.col.fact]])[1]) 
  
  prev.all.inst <- prev.tabs.ls[[lvl.inst]] %>%
                      mutate(feature = rownames(.))
  
  for(i.set.gr in names(sample.sets)) {
                  
    for(i.set in sample.sets[[i.set.gr]]) {
      
      # Prevalence data subset 
      if(length(i.set) == 1) {
        
        prev.inst <- prev.all.inst %>% 
                      select(all_of(c("feature", 
                                      paste0(i.set.gr, "--", i.set))))
      } else {
        
        prev.inst <- prev.all.inst %>% 
                        select(all_of(c("feature", "across_all")))
      }
        
    # Delta prevalence subset 
    otu.inst.delta.set <- otu.inst.delta %>% 
                              filter(.data[[i.set.gr]] %in% i.set) 
    
    # Correlation 
    cor.res.inst <- data.frame()
    
    for (i.tax in colnames(otu.inst)) {
      
      cor.res <- cor.test(otu.inst.delta.set[[i.tax]], 
                          otu.inst.delta.set[[corr.cols]], 
                          method = "spearman")
      
      res.inst <- c(feature = i.tax, 
                    rho = as.numeric(cor.res$estimate), 
                    pval = cor.res$p.value, 
                    Variable = i.set.gr,
                    Group = paste(i.set, collapse = "&"), 
                    "Taxa Level" = lvl.inst, 
                    method = cor.res$method) 
      
      cor.res.inst <- bind_rows(cor.res.inst, res.inst)
      
    }
    
    # Filter and add additional information 
    cor.res.f.inst <- cor.res.inst %>% 
                          mutate(across(all_of(c("rho", "pval")), as.numeric)) %>% 
                          left_join(., prev.inst, by = "feature") %>% 
                          rename("Prevalence" = 
                                   colnames(prev.inst)[colnames(prev.inst) != 
                                                        "feature"]) %>% 
                          filter(Prevalence >= min.prev) %>% 
                          arrange(desc(abs(rho))) %>% 
                          mutate(qval = p.adjust(pval, method = "BH"), 
                                 feature = factor(feature, 
                                                  levels = unique(feature))) 
    
    cor.res.f.out.inst <- cor.res.f.inst %>% 
                            select(feature, rho, pval, qval, Prevalence) %>% 
                            mutate(across(-feature, 
                                          function(x){round(x, 3)})) %>% 
                            mutate(across(-feature, 
                                          function(x){ifelse(x == 0, 
                                                             ">0.001", 
                                                             sprintf("%.3f", x))}))
                          
    # Collect data 
    res.corr.comb[[lvl.inst]][[norm.inst]][[i.set.gr]][[
                   paste(i.set, collapse = "_")]][["tab"]] <- cor.res.f.out.inst 
    
    # Write data
    write.csv(cor.res.f.out.inst, 
              paste0(out.folder, "/tabs/", 
                     "CorrRes--", i.set.gr, "(", 
                     paste(i.set, collapse = "_"), ")--",
                     lvl.inst, "--", 
                     norm.inst, "--",
                     min.prev, ".csv"))
    
    # Plot 
    #---------------------------------------------------------------------------
    p.cor.df <- cor.res.f.inst[1:top.n.plot, ] %>% 
                droplevels() %>% 
                mutate(across(all_of(c("pval", "qval")), 
                       function(x){ifelse(x < 0.01, 
                                          "<0.01", 
                                          paste0("=", sprintf("%.2f", 
                                                             round(x, 2))))})) %>% 
                mutate(lab_col = paste0("[", 
                       "\u03C1=", sprintf("%.2f", round(rho, 2)), "; ", 
                       "p", pval, "; ", 
                       "q", qval, "]"))
    
    p.df <- otu.inst.delta.set %>% 
              pivot_longer(cols = colnames(otu.inst), 
                           names_to = "feature", 
                           values_to = "Abundance") %>% 
              filter(feature %in% p.cor.df$feature) %>% 
              arrange(factor(feature, levels = levels(p.cor.df$feature))) %>% 
              left_join(., p.cor.df[, c("feature", "lab_col")], 
                        by = "feature") %>% 
              mutate(y_text = (max(Abundance) + 
                                (max(Abundance)-min(Abundance))*0.2), 
                     x_text = min(.data[[corr.cols]]), 
                     .by = "feature") %>% 
              mutate(feature = gsub("^(.{1})(.*)$", "\\1.\\2", feature)) %>% 
              mutate(feature = gsub("__", " ", feature)) %>% 
              mutate(feature = factor(feature, levels = unique(feature)))
   
    
    p <- ggplot(p.df, aes(x = .data[[corr.cols]], y = Abundance)) +
                geom_smooth(method = "lm", se = FALSE) +
                geom_point() +
                geom_text(aes(y = y_text, 
                               x = x_text, 
                               label = lab_col), 
                           inherit.aes = FALSE, hjust = 0, vjust = 0.75) +
                facet_wrap(~feature, 
                           scales = "free", 
                           ncol = plots.n.col) + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      strip.text = element_text(size = 6)) + 
                ylab(label = paste0("Abundance (", norm.inst, ")"))
    
    
    w.plot <- plots.n.col*2 + 1.5
      
    if(length(levels(p.df$feature)) < plots.n.col) {
      
      w.plot <- length(levels(p.df$feature))*2 + 1.5
      
    }
    
    h.plot <- ceiling(length(levels(p.df$feature))/plots.n.col*2) + 0.5
    
    res.corr.comb[[lvl.inst]][[norm.inst]][[i.set.gr]][[
      paste(i.set, collapse = "_")]][["plot"]] <- list(p = p, 
                                                       w = w.plot, 
                                                       h = h.plot)
    
    # Write plot 
    ggsave(paste0(out.folder, "/plots/", 
                   "Corr--", i.set.gr, "(", 
                   paste(i.set, collapse = "_"), ")--",
                   lvl.inst, "--", 
                   norm.inst, "--", 
                   min.prev, ".png"), 
           plot = p, width = w.plot, height = h.plot, dpi = 600)
    }
  }
}
  

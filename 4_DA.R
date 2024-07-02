#-------------------------------------------------------------------------------
# Load data and functions
#-------------------------------------------------------------------------------
load("out/supp/prm.R")
load("out/supp/0_data.RData")


source("R/filt_tax_phy.R")

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

tax.lvl <- prm.ls[["DA"]][["Tax_lvl"]]

norm <- prm.ls[["DA"]][["Norm"]]

qval.cut <- prm.ls[["DA"]][["Qval_cutoff"]]

pval.cut <- prm.ls[["DA"]][["Pval_cutoff"]]

min.prev <- prm.ls[["DA"]][["Min_prev"]]

rand.col <- prm.ls[["DA"]][["Maas_rand_effect"]]

fix.col <- prm.ls[["DA"]][["Maas_fix_effect"]]

maaslin.norm.method <- prm.ls[["DA"]][["Maas_norm"]]

maas.transform <- prm.ls[["DA"]][["Maas_tansform"]]

maas.cor.p.method <- prm.ls[["DA"]][["Maas_pcorr_method"]]

maas.method <- prm.ls[["DA"]][["Maas_method"]]

out.folder <- prm.ls[["DA"]][["out_dir_path"]]


set.seed(seed)

res.da.comb <- list()


################################################################################
# MAASLIN per group
################################################################################
dir.create(paste0(out.folder, "/tabs"), 
           recursive = TRUE, showWarnings = FALSE)

prm.da.grid <- expand.grid("Tax_lvl" = tax.lvl, 
                           "Normal" = norm, 
                           "Maas_normal" = maaslin.norm.method,
                           "Group_col" = gr.cols,
                           stringsAsFactors = FALSE)

for(i.prm in 1:nrow(prm.da.grid)) {
  
  lvl.inst <- prm.da.grid[i.prm, "Tax_lvl"]
  
  norm.inst <- prm.da.grid[i.prm, "Normal"]
  
  maas.norm.inst <- prm.da.grid[i.prm, "Maas_normal"]
  
  gr.col.inst <- prm.da.grid[i.prm, "Group_col"]
  
  tab.lbl <-  paste0(gr.col.inst, "--", 
                     lvl.inst, "--", 
                     norm.inst, "--", 
                     min.prev)
  
  tax.keep.inst <- prev.tabs.ls[[lvl.inst]] %>% 
                      select(starts_with(paste0(gr.col.inst, "--"))) %>% 
                      filter(if_any(everything(), ~. >= min.prev)) %>% 
                      rownames(.)
  
  ps.inst <- ps.ls[[lvl.inst]][[norm.inst]]  
  
  ps.inst <- prune_taxa(intersect(taxa_names(ps.inst), tax.keep.inst), ps.inst)
  
  res.maas.long <- NULL
  
  summary.long <- NULL
  
  for(i.gr.lvl in levels(ps.meta[[gr.col.inst]])) {
    
    meta.gr.inst <- ps.meta %>% 
                      filter(.data[[gr.col.inst]] == i.gr.lvl)
    
    asv.gr.inst <- ps.inst %>% 
                      prune_samples(rownames(meta.gr.inst), .) %>% 
                      otu_table()  %>% 
                      as.matrix()  %>% 
                      as.data.frame() 
    
    # Make a directory for the maaslin output 
    out.dir.inst <- paste0(out.folder, "/maaslin/", 
                           gr.col.inst, "--",
                           lvl.inst, "--", 
                           norm.inst, "--",
                           min.prev, "--",
                           i.gr.lvl)
    
    dir.create(out.dir.inst, showWarnings = FALSE, recursive = TRUE)
    
    maas.out <- Maaslin2(
                  input_data =  asv.gr.inst, 
                  input_metadata = meta.gr.inst, 
                  output =  out.dir.inst, 
                  fixed_effects = fix.col, 
                  random_effects = rand.col, 
                  correction = maas.cor.p.method, 
                  cores = n.core, 
                  min_abundance = 0, 
                  min_prevalence = 0, 
                  min_variance = 0, 
                  normalization = maaslin.norm.method, 
                  transform = maas.transform, 
                  analysis_method = maas.method, 
                  max_significance = qval.cut)
    
    
    res.da.comb[[lvl.inst]][[norm.inst]][[gr.col.inst]][[
                  "tabs"]][["ind"]][[
                    i.gr.lvl]] <- maas.out
    
    res.maas.long <- maas.out$results %>% 
                          mutate(Taxa_lvl = lvl.inst, 
                                 Norm_method = norm.inst,
                                 Test_group = gr.col.inst, 
                                 Strata = i.gr.lvl)  %>% 
                          bind_rows(res.maas.long, .)
    
    
    # Summary table that will be printed in results 
    descr.df <- asv.gr.inst %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      mutate(SeqID = rownames(.)) %>% 
                      pivot_longer(-SeqID, 
                                   names_to = "feature", 
                                   values_to = "Abundance") %>% 
                      left_join(., meta.gr.inst, by = "SeqID") %>% 
                      reframe(Mean_abund = mean(Abundance), 
                              SD_abund = sd(Abundance),
                              Median_abund = median(Abundance), 
                              .by = all_of(c("feature", time.col.fact))) %>% 
                      pivot_wider(values_from = -c("feature", time.col.fact), 
                                  names_from = .data[[time.col.fact]])
    
    summary.long <- maas.out$results %>% 
                          select(feature, coef, stderr, 
                                 pval, qval, N, N.not.zero) %>% 
                          left_join(., descr.df, by = "feature") %>% 
                           mutate(across(-feature, function(x){round(x, 3)})) %>% 
                           mutate(Strata = i.gr.lvl) %>% 
                           rbind(summary.long, .)
  }
  
  res.da.comb[[lvl.inst]][[norm.inst]][[gr.col.inst]][[
                 "tabs"]][["long"]] <- res.maas.long
  
  summary.wide <- summary.long %>% 
                    pivot_wider(values_from = -c("Strata", "feature"), 
                                names_from = Strata) %>% 
                    mutate(Group_var = gr.col.inst, 
                           Fix_effect = fix.col, 
                           Rand_effect = rand.col, 
                           Normalization = norm.inst) %>% 
                    mutate(feature = gsub("^(.{1})(.*)$", "\\1.\\2", feature)) %>% 
                    mutate(feature = gsub("__", " ", feature))
  
  res.da.comb[[lvl.inst]][[norm.inst]][[gr.col.inst]][[
                 "tabs"]][["summary"]] <- summary.wide
  
  write.csv(summary.wide, 
            paste0(out.folder,"/tabs/maas_summary--", tab.lbl, ".csv"), 
            row.names = FALSE)
}


################################################################################
# Visualization 
################################################################################
dir.create(paste0(out.folder, "/plots"), 
           recursive = TRUE, showWarnings = FALSE)

prm.heat.grid <- expand.grid("Norm" = norm,
                             "Group_col" = gr.cols,
                             stringsAsFactors = FALSE)

for(i.prm in 1:nrow(prm.heat.grid)) {
  
  # Variables
  norm.inst <- prm.heat.grid[i.prm, "Norm"]
  
  gr.col.inst <- prm.heat.grid[i.prm, "Group_col"]
  
  
  abund.lbl <- paste0("Abundance (", norm.inst, ")")
  
  plot.lbl <- paste0(gr.col.inst, "--", 
                     norm.inst, "--", 
                     min.prev)
  
  # Order metadata
  ps.meta.ord <- ps.meta %>% 
                    arrange(.data[[gr.col.inst]], 
                            .data[[time.col.fact]], 
                            .data[[preb.col]], 
                            .data[[id.col]])
  
  #-----------------------------------------------------------------------------
  # Extract data for heat maps - combined taxonomic levels 
  maas.res.inst.sig <- NULL
  
  asv.gr.inst <- NULL
  
  maas.res.violin.inst <- NULL
  
  for(lvl.inst in tax.lvl) {
    
    # Extract maaslin results 
    maas.res.inst <- res.da.comb[[lvl.inst]][[norm.inst]][[gr.col.inst]][[
                                      "tabs"]][["long"]] 
    
    # Significant taxa in first level of the group 
    sig.tax.maas <- maas.res.inst %>% 
                      filter(pval <= pval.cut | qval <= qval.cut) %>% 
                      filter(Strata == unique(Strata)[1])
    
    # Wide table with maaslin results 
    if(nrow(sig.tax.maas) > 0) {
      
      maas.res.violin.inst <- maas.res.inst %>% 
                                filter(feature %in% sig.tax.maas$feature, 
                                       .by = Strata) %>% 
                                rbind(maas.res.violin.inst, .)
      
      maas.res.inst.sig0 <- maas.res.inst %>% 
                            filter(feature %in% sig.tax.maas$feature, 
                                   .by = Strata) %>% 
                            select(feature, coef, pval, qval, Strata) %>% 
                            pivot_wider(names_from = "Strata", 
                                        values_from = c("coef", "pval", "qval")) 
      
      # Add a column to split heat map rows 
      gr.lvls <- c(levels(ps.meta.ord[[gr.col.inst]]), 
                      paste(levels(ps.meta.ord[[gr.col.inst]]), collapse = "&"))
      
      maas.res.inst.sig <- sig.tax.maas %>% 
                              select(feature, Strata) %>% 
                              pivot_wider(values_from = Strata, 
                                          names_from = Strata) %>% 
                              unite(., Group_lvl, -feature, na.rm = TRUE, sep = "&") %>% 
                              left_join(maas.res.inst.sig0, ., by = "feature") %>% 
                              mutate(Group_lvl = factor(Group_lvl, levels = gr.lvls)) %>% 
                              arrange(Group_lvl) %>% 
                              arrange(Group_lvl, 
                                      desc(abs(across(all_of(grep("coef_", 
                                                                  names(.), 
                                                                  value = TRUE)))))) %>% 
                              mutate(split = lvl.inst) %>% 
                              rbind(maas.res.inst.sig, .)
      
      # Extract asv 
      asv.gr.inst <- ps.ls[[lvl.inst]][[norm.inst]] %>% 
                        prune_taxa(taxa_names(.) %in% sig.tax.maas$feature, .) %>% 
                        otu_table()  %>% 
                        as.matrix()  %>% 
                        as.data.frame() %>% 
                        rbind(asv.gr.inst, .) 
                        
      
      }
  }
  
  asv.gr.inst <- asv.gr.inst %>%                         
                    .[maas.res.inst.sig[["feature"]], rownames(ps.meta.ord)]
  
  # Adjust taxa names in used objects 
  rownames(asv.gr.inst) <- rownames(asv.gr.inst) %>% 
                                gsub("^(.{1})(.*)$", "\\1.\\2", .) %>% 
                                gsub("__", " ", .)
  
  maas.res.inst.sig <- maas.res.inst.sig %>% 
                          mutate(feature = gsub("^(.{1})(.*)$", "\\1.\\2", feature)) %>% 
                          mutate(feature = gsub("__", " ", feature))
  
  maas.res.violin.inst <- maas.res.violin.inst %>% 
                              mutate(feature = gsub("^(.{1})(.*)$", "\\1.\\2", feature)) %>% 
                              mutate(feature = gsub("__", " ", feature))
    
  
  ##############################################################################
  # Heat map
  #-----------------------------------------------------------------------------
  ht_opt$TITLE_PADDING = unit(c(2.5, 2.5), "points")
  
  # Color scheme
  #-----------------------------------------------------------------------------
  col_fun = colorRamp2(c(min(asv.gr.inst), 
                         0, 
                         max(asv.gr.inst)), 
                       c("blue", "white", "gold"))
  
  # Top annotation 
  #-----------------------------------------------------------------------------
  top.anot.df <- ps.meta.ord %>% 
                  select(all_of(c(gr.col.inst, time.col.fact, preb.col))) 
  
  top.anot <- HeatmapAnnotation(df = top.anot.df, 
                                col = setNames(list(color.sch.ls[[gr.col.inst]],
                                                    color.sch.ls[[time.col.fact]], 
                                                    color.sch.ls[[preb.col]]), 
                                               c(gr.col.inst, 
                                                 time.col.fact, 
                                                 preb.col)),
                                gp = gpar(col = "gray10", lwd = 0.001), 
                                simple_anno_size = unit(0.35, "cm"), 
                                annotation_name_gp = gpar(fontsize = 8), 
                                annotation_name_side = "left")
  
  ps.meta.ord$Split_col <- paste0(ps.meta.ord[[gr.col.inst]], ":", 
                                   ps.meta.ord[[time.col.fact]]) %>% 
                            factor(., levels = unique(.))
  
  
  # Row annotation
  #-----------------------------------------------------------------------------
  # Limits for bar annotation 
  coef.both <- maas.res.inst.sig %>% 
                    select(starts_with("coef_")) 
  
  bar.anot.lim <- c(min(coef.both)*1.5, 
                    max(coef.both)*1.2)
  
  # Annotation object to add to heat maps
  row.anot.plot <- NULL
  
  for(i.gr in levels(ps.meta.ord[[gr.col.inst]])) {
    
    color.inst <- color.sch.ls[[gr.col.inst]][i.gr]
    
    p.inst <- maas.res.inst.sig %>% 
                pull(paste0("pval_", i.gr)) 
    
    q.inst <- maas.res.inst.sig %>% 
                pull(paste0("qval_", i.gr)) %>% 
                round(., 2) %>% 
                sprintf("%.2f", .) %>% 
                ifelse(. == "0.00", "<0.01", .)
    
    coef.inst <- maas.res.inst.sig %>% 
                    pull(paste0("coef_", i.gr))
    
    anot.inst <- ifelse(p.inst < 0.01, 
                        paste0(" p<0.01[", q.inst, "]  "), 
                        paste0(" p=", sprintf("%.2f", p.inst), 
                               "[", q.inst, "]  "))
    
    row.anot.plot <- row.anot.plot + 
                        rowAnnotation(
                         "pval" = anno_text(anot.inst,  
                                               gp = gpar(fontsize = 6)), 
                          "Coef" = anno_barplot(coef.inst, 
                                             width = unit(1.5, "cm"), 
                                             baseline = 0, 
                                             ylim = bar.anot.lim, 
                                             gp = gpar(fontsize = 6, 
                                                       fill = color.inst)), 
                         annotation_name_side = "top", 
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 8))
  }
  
  # Heat map plot size
  w.heat <- ncol(asv.gr.inst)/10
  
  h.heat <- nrow(asv.gr.inst)/6
  
  # Plot heatmap 
  #-----------------------------------------------------------------------------
  heat.map.sig <- Heatmap(asv.gr.inst, 
                          col = col_fun ,
                          name = abund.lbl, 
                          width = unit(w.heat, "in"), 
                          height = unit(h.heat, "in"), 
                          column_split = ps.meta.ord$Split_col, 
                          column_title_gp = gpar(fontsize = 10),
                          row_title_gp = gpar(fontsize = 7, 
                                              fontface = "bold", 
                                              fill = "white", 
                                              border = "gray"), 
                          row_split = maas.res.inst.sig$split, 
                          # row_title_rot = 0,
                          cluster_column_slices = FALSE,
                          cluster_rows = FALSE,
                          top_annotation = top.anot, 
                          cluster_columns = FALSE, 
                          show_column_names = FALSE, 
                          show_row_dend = FALSE, 
                          row_names_side = "left", 
                          rect_gp = gpar(col = "gray40", lwd = 0.001),
                          row_names_gp = gpar(fontsize = 6.5, 
                                              fontface = "italic"), 
                          column_title = paste0("Significantly different taxa (", 
                                                gr.col.inst, ")")) + 
                  row.anot.plot
  
  
  res.da.comb[["comb"]][[norm.inst]][[gr.col.inst]][[
                 "plots"]][["heat"]] <- list("p" = heat.map.sig, 
                                             "w" = w.heat + 8,
                                             "h" = h.heat + 1)
  
  # Save the plot 
  png(filename = paste0(out.folder,"/plots/heat--", plot.lbl, ".png"), 
      width = w.heat + 8, 
      height = h.heat + 1, 
      units = "in", res = 600)
  draw(heat.map.sig)
  dev.off()
  
  
  ##############################################################################
  # Violin plot 
  
  # Data prep
  #-----------------------------------------------------------------------------
  # Main data frame 
  violin.main.df <- asv.gr.inst %>% 
                        t() %>% 
                        cbind(., ps.meta.ord) %>% 
                        pivot_longer(cols = rownames(asv.gr.inst), 
                                     names_to = "feature", 
                                     values_to = abund.lbl) %>% 
                        left_join(., maas.res.inst.sig, by = "feature") %>% 
                        mutate(feature = factor(feature, 
                                                levels = maas.res.inst.sig$feature))
  
  # Max y per facet 
  max.y.df <- violin.main.df %>% 
                  slice(which.max(.data[[abund.lbl]]), 
                        .by = c("feature", "split")) %>% 
                  mutate(y_lab = .data[[abund.lbl]] * 1.15, 
                         inv_point = .data[[abund.lbl]] * 1.35) %>% 
                  select(feature, split, y_lab)
    
  
  # Annotation df
  sig.anot.df <- maas.res.violin.inst %>% 
                  left_join(., maas.res.inst.sig[c("feature", "split")], 
                            by = "feature") %>% 
                  mutate(across(all_of(c("coef", "pval", "qval")), 
                                function(x){round(x, digits = 2)})) %>% 
                  mutate(x_start = paste0(Strata, ":", 
                                          levels(ps.meta.ord[[time.col.fact]])[1]), 
                         x_end = paste0(Strata, ":", 
                                        levels(ps.meta.ord[[time.col.fact]])[2]), 
                         lab = ifelse(qval < 0.01, 
                                      paste0("q<0.01 [e=", 
                                             sprintf("%.2f", coef), "]"), 
                                      paste0("q=", sprintf("%.2f", qval), 
                                             " [e=", sprintf("%.2f", coef), "]"))) %>% 
                  mutate(feature = factor(feature, 
                                          levels = levels(violin.main.df$feature)), 
                         split = factor(split)) %>% 
                  left_join(., max.y.df, by = c("feature", "split"))
  
  # Plot dimensions 
  n.tax <- length(levels(violin.main.df$feature))
  
  w.violin <- 4*2.5 + 1
  
  if(n.tax < 4) {w.violin <- n.tax*2.5 + 1}
  
  h.violin <- ceiling(n.tax/4)*2 + 0.5
  
  # Plot data violin plot
  violin.plot <- ggplot(violin.main.df, aes(x = Split_col, 
                                           y = .data[[abund.lbl]])) + 
                        geom_violin() + 
                        geom_point(aes(colour = .data[[preb.col]])) +
                        geom_line(aes(group = PersonID, 
                                      color = Prebiotic)) +
                        geom_signif(data = sig.anot.df,
                                    aes(xmin = x_start,
                                        xmax = x_end,
                                        annotations = lab,
                                        y_position = y_lab),
                                    textsize = 2, vjust = -0.1,
                                    manual = TRUE, margin_top = 1) +
                        geom_point(data = sig.anot.df,
                                   aes(x = x_start, y =  y_lab*1.05), x=NA) +
                        facet_wrap(c("feature", "split"), 
                                   scales = "free_y", ncol = 4) + 
                        theme_bw() +
                        theme(axis.title.x = element_text()) + 
                        guides(colour = guide_legend(order = 2), 
                               linetype = guide_legend(order = 1)) +
                        theme(axis.text.x = element_text(angle=45, 
                                                         vjust=1.1, 
                                                         hjust=1), 
                              strip.text = element_text(size = 6), 
                              axis.title.x = element_blank(), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank()) + 
                        scale_color_manual(values = color.sch.ls[[preb.col]])
        
  res.da.comb[["comb"]][[norm.inst]][[gr.col.inst]][[
                 "plots"]][["violin"]] <- list("p" = violin.plot, 
                                               "w" = w.violin, 
                                               "h" = h.violin)
  
  
  
  # Save violin plot 
  ggsave(filename = paste0(out.folder,"/plots/violin--", plot.lbl, ".png"), 
         plot = violin.plot, 
         width = w.violin, 
         height = h.violin, dpi = 600)
}


save(list = c("res.da.comb"),
     file = "out/supp/4_DA_res.RData")

# rm(list = ls())

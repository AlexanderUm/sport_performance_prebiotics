#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(3490759)

lib.to.load <- c("phyloseq", "tidyverse", 
                 "ggsignif", "broom", "cowplot", 
                 "Maaslin2", "MicrobiomeStat", 
                 "ComplexHeatmap", "circlize")

for (i in lib.to.load) {library(i, character.only = TRUE)}

rm(list = c("i", "lib.to.load"))


source("R/filt_tax_phy.R")

# Load data
load("out/supp/0_data.RData")



#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------

group.vars <- c("Treatment", "Resp")

qval.cut <- 0.1

pval.cut <- 0.01

min.prev <- 0.30

maaslin.norm.method <- "TSS"

used.ps <- c("raw_asv", "raw_genus")

out.folder <- "out/maaslin/"

#-------------------------------------------------------------------------------
# MAASLIN per group
#-------------------------------------------------------------------------------
maas.res.full <- list()

maas.res.long.df <- list()

for (i.ps in used.ps) {
  
  ps.inst <- pss.ls[[i.ps]] %>% 
              tax_filt_phy(., prev = min.prev, 
                           group_col = "Treatment")
              
  for (i.gr in group.vars) {
    
     res.maas.long <-NULL
  
    for (i in unique(ps.meta[[i.gr]]))  { 
    
      dir.create(paste0(out.folder, "gr_", i), recursive = TRUE)
    
      ps.gr <- prune_samples(ps.meta[[i.gr]] == i, ps.inst)
    
      ps.gr.asv <- ps.gr %>% 
                    otu_table()  %>% 
                    as.matrix()  %>% 
                    as.data.frame() 
    
      ps.gr.meta <- ps.meta[ps.meta[[i.gr]] == i, ]
    
      maas.out <- Maaslin2(
                    input_data =  ps.gr.asv, 
                    input_metadata = ps.gr.meta, 
                    output =  paste0(out.folder, "gr_", i), 
                    fixed_effects = c("Time"), 
                    random_effects = "PersonID", 
                    correction = "BH", 
                    cores = 4, 
                    min_abundance = 0, 
                    min_prevalence = 0, 
                    min_variance = 0, 
                    normalization = maaslin.norm.method, 
                    transform = "LOG", 
                    analysis_method = "LM", 
                    max_significance = 0.1)
    
    
    maas.res.full[[i.ps]][[i]] <- maas.out
    
    res.maas.long <- maas.out$results %>% 
                          mutate(Strata = i, 
                                 Level = i.ps, 
                                 Group = i.gr)  %>% 
                          bind_rows(res.maas.long, .)
    }
     
     maas.res.long.df[[i.gr]][[i.ps]] <- res.maas.long
    
  }

}



################################################################################
# Visualization 
################################################################################
lab.tax.lvl <- c(raw_asv = "ASVs", 
                  raw_genus = "Genera")



heat.all.ls <- list()

sig.heat.ls <- list()

viol.plot.ls <- list()

for(i.ps in used.ps) {
  
  for(i.gr in group.vars) {
  #----------------------------------------------------------------------------
  # Abundance data preparation 
  #----------------------------------------------------------------------------
  ps.meta.ord <- ps.meta %>% 
                 arrange(all_of(i.gr), TestDay, Prebiotic, PersonID)
    
    
  strata.lvl <- levels(ps.meta[[i.gr]])
    
  sig.tax.maas <- maas.res.long.df[[i.gr]][[i.ps]] %>% 
                        filter(pval <= pval.cut) 
  
  sig.tax.PRE <- maas.res.long.df[[i.gr]][[i.ps]] %>% 
                            filter(Strata == strata.lvl[2], 
                                   feature %in% sig.tax.maas$feature) %>% 
                            select(feature, coef, pval, qval)
  
  sig.tax.CON <- maas.res.long.df[[i.gr]][[i.ps]] %>% 
                    filter(Strata == strata.lvl[1], 
                           feature %in% sig.tax.maas$feature) %>% 
                    arrange(match(feature, sig.tax.PRE$feature)) %>%
                    select(feature, coef, pval, qval) %>% 
                    rename(coef.CON = coef, 
                           pval.CON = pval, 
                           qval.CON = qval)      
  
  sig.tax.comb <- full_join(sig.tax.PRE, sig.tax.CON, by = "feature") %>% 
                    arrange(desc(coef))
  
  row.split.sig <- ifelse(sig.tax.comb$qval <= qval.cut, 
                           paste0("q < ", qval.cut), 
                           ifelse(sig.tax.comb$pval <= pval.cut, 
                                  paste0("p < ", pval.cut), 
                                  paste0("p < ", pval.cut))) %>% 
                          factor(., levels = c(paste0("q < ", qval.cut), 
                                            paste0("p < ", pval.cut), 
                                            "p > 0.05"))
  
  otu.inst.rel <- pss.ls[[i.ps]] %>% 
                    transform_sample_counts(function(x) (x/sum(x))*100) %>% 
                    otu_table() %>% 
                    as.matrix() %>% 
                    as.data.frame() %>%   
                    as.matrix() %>% 
                    .[, as.character(ps.meta.ord$SeqID)]
  
  otu.inst <- otu.inst.rel
  
  otu.inst[otu.inst == 0] <- min(otu.inst[otu.inst != 0])/2
  
  otu.inst <- log(otu.inst)
  
  
  otu.inst.f <- otu.inst[sig.tax.comb$feature, ]
  
  #----------------------------------------------------------------------------
  # Heatmap top Annotation 
  #----------------------------------------------------------------------------
  top.anot.df <- ps.meta.ord %>% 
                  select(TestDay, all_of(i.gr), Prebiotic) %>% 
                  rename("Test Day" = TestDay)
  
  top.anot <- HeatmapAnnotation(df = top.anot.df, 
                                col = setNames(list(color.sch.ls$TestDay,
                                                    color.sch.ls[[i.gr]], 
                                                    color.sch.ls$Prebiotic), 
                                           c("Test Day", i.gr, "Prebiotic")),
                                gp = gpar(col = "gray10", lwd = 0.001), 
                                simple_anno_size = unit(0.35, "cm"), 
                                annotation_name_gp = gpar(fontsize = 8), 
                                annotation_name_side = "left")
  
  #----------------------------------------------------------------------------
  # Heat map row annotation 
  #----------------------------------------------------------------------------
  # Data for annotation 
  q.pre.anot <- round(sig.tax.comb[, c("pval", 
                                        "qval")], 2)  
  
  q.pre.anot.chr <- ifelse(q.pre.anot$pval < 0.01, 
                           paste0(" p<0.01"), 
                           paste0(" p=", sprintf("%.2f", q.pre.anot$pval))) %>% 
                    paste0(., "(", sprintf("%.2f", q.pre.anot$qval), ") ") %>% 
                    gsub("\\NA\\(\\NA\\) ", " ", .)
  
  
  q.con.anot <- round(sig.tax.comb[, c("pval.CON", 
                                       "qval.CON")], 2)  
  
  q.con.anot.chr <- ifelse(q.con.anot$pval.CON < 0.01, 
                           paste0(" p<0.01"), 
                           paste0(" p=", sprintf("%.2f",q.con.anot$pval.CON))) %>% 
                      paste0(., "(", sprintf("%.2f", q.con.anot$qval.CON), ") ") %>% 
                      gsub("\\NA\\(\\NA\\) ", " ", .)
  
  
  bar.anot.lim <- c((min(sig.tax.maas$coef)*1.5), 
                    (max(sig.tax.maas$coef)*1.1))
  

  # Row annotation
  if(i.gr == "Treatment") {
    
    row.anot.plot <- rowAnnotation(
                    "q-PRE" =  anno_text(q.pre.anot.chr,  
                                         gp = gpar(fontsize = 6)), 
                     PRE = anno_barplot(sig.tax.comb$coef, 
                                        width = unit(1.5, "cm"), 
                                        baseline = 0, 
                                        ylim = bar.anot.lim, 
                                        gp = gpar(fontsize = 6, 
                                              fill = color.sch.ls[[i.gr]][2])), 
                    "q-CON" =  anno_text(q.con.anot.chr, 
                                         gp = gpar(fontsize = 6)),
                    CON = anno_barplot(sig.tax.comb$coef.CON, 
                                       width = unit(1.5, "cm"), 
                                       baseline = 0, 
                                       ylim = bar.anot.lim, 
                                       gp = gpar(fontsize = 6, 
                                              fill = color.sch.ls[[i.gr]][1])), 
                                   annotation_name_side = "top", 
                                   annotation_name_rot = 0,
                                   annotation_name_gp = gpar(fontsize = 8))
  } else {
  
    row.anot.plot <- rowAnnotation(
                    "q-PRE" =  anno_text(q.pre.anot.chr,  
                                         gp = gpar(fontsize = 6)), 
                     RESP = anno_barplot(sig.tax.comb$coef, 
                                         width = unit(1.5, "cm"), 
                                         baseline = 0, 
                                         ylim = bar.anot.lim, 
                                         gp = gpar(fontsize = 6, 
                                                   fill = color.sch.ls[[i.gr]][2])), 
                    "q-CON" =  anno_text(q.con.anot.chr, 
                                         gp = gpar(fontsize = 6)),
                    NONE = anno_barplot(sig.tax.comb$coef.CON, 
                                        width = unit(1.5, "cm"), 
                                        baseline = 0, 
                                        ylim = bar.anot.lim, 
                                        gp = gpar(fontsize = 6, 
                                                  fill = color.sch.ls[[i.gr]][1])), 
                                   annotation_name_side = "top", 
                                   annotation_name_rot = 0,
                                   annotation_name_gp = gpar(fontsize = 8))
  }
 
    # Color scheme
    col_fun = colorRamp2(c(min(otu.inst.f), 
                           0, 
                           max(otu.inst.f)), 
                         c("blue", "white", "gold"))
    
  #----------------------------------------------------------------------------
  # Heat map
  #----------------------------------------------------------------------------
   
   cols.split <- paste0(ps.meta.ord[[i.gr]], ":", 
                        ps.meta.ord$TestDay)
    
   heat.map.sig <- Heatmap(otu.inst.f, 
                             col = col_fun ,
                             name = "Abundance(log)",
                             column_split = cols.split, 
                             column_title_gp = gpar(fontsize = 10),
                             row_title_gp = gpar(fontsize = 10, 
                                                 fill = "gray80"), 
                             row_split = row.split.sig,
                             cluster_column_slices = FALSE,
                             cluster_rows = FALSE,
                             top_annotation = top.anot, 
                             right_annotation = row.anot.plot,
                             cluster_columns = FALSE, 
                             show_column_names = FALSE, 
                             show_row_dend = FALSE, 
                             row_names_side = "left", 
                             rect_gp = gpar(col = "gray40", lwd = 0.001),
                             row_names_gp = gpar(fontsize = 8, 
                                                  fontface = "italic"), 
                             column_title = paste0("Significant different ", 
                                                  lab.tax.lvl[i.ps]))
   
   sig.heat.ls[[i.gr]][[i.ps]] <- heat.map.sig
   
   
   #----------------------------------------------------------------------------
   # Free clustering heatmaps 
   #----------------------------------------------------------------------------
   otu.inst.f.all <- otu.inst[unique(maas.res.long.df[[i.gr]][[i.ps]]$feature), ]
   
   heat.all.ls[[i.gr]][[i.ps]][["all"]] <- Heatmap(otu.inst.f.all, 
                                      col = col_fun ,
                                      name = "Abundance(log)",
                                      top_annotation = top.anot, 
                                      cluster_columns = TRUE, 
                                      show_column_names = FALSE, 
                                      show_row_names = FALSE, 
                                      show_row_dend = FALSE)
   
   heat.all.ls[[i.gr]][[i.ps]][["sig"]] <- Heatmap(otu.inst.f, 
                                      col = col_fun ,
                                      name = "Abundance(log)",
                                      top_annotation = top.anot, 
                                      show_column_names = FALSE, 
                                      show_row_dend = FALSE, 
                                      row_names_side = "left")
   
   #-------------------------------------------------------------------------------
   # Violin plots 
   #-------------------------------------------------------------------------------
   viol.sig.tax <- maas.res.long.df[[i.gr]][[i.ps]] %>% 
                      filter(qval <= qval.cut,
                             Strata == levels(ps.meta.ord[[i.gr]])[2]) %>% 
                      arrange(-abs(coef))
   
   viol.sig.stat <- maas.res.long.df[[i.gr]][[i.ps]] %>% 
                      filter(feature %in% viol.sig.tax$feature) 
   
   viol.sig.stat$anot <- ifelse(viol.sig.stat$qval < qval.cut, 
                                paste0("q<0.01"), 
                                paste0("q=", 
                                       sprintf("%.2f", 
                                               viol.sig.stat$qval))) %>% 
                                paste0(., "[", 
                                       sprintf("%.2f", 
                                               viol.sig.stat$coef), "]") 
   
   viol.df <- otu.inst[viol.sig.tax$feature, ] %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    bind_cols(., ps.meta.ord) %>% 
                    select(c(viol.sig.tax$feature, 
                             "Treatment", "PersonID", "Time", 
                             "Prebiotic", "TestDay", 
                             "TreatmentDay", "Resp")) %>% 
                    pivot_longer(cols = c(viol.sig.tax$feature), 
                                 names_to = "feature", 
                                 values_to = "Value") %>%
                    mutate(x.axis.var = paste0(.data[[i.gr]], ":" ,TestDay), 
                           feature = factor(feature, 
                                            levels = viol.sig.tax$feature))

   
   # Make data frame for significance levels 
   max.div <- viol.df %>% 
                  group_by(across(c(i.gr, 
                                    "feature", 
                                    "x.axis.var"))) %>% 
                   slice(which.max(Value)) %>%
                   ungroup() %>%
                   mutate(y = max(Value)*1.3)
   
   sig.df <- viol.sig.stat %>% 
                 mutate(!!i.gr := Strata) %>% 
                 left_join(max.div, ., by = c(i.gr, 
                                              "feature")) %>%
                 select(feature, x.axis.var, 
                        y, anot, TestDay) %>%
                 pivot_wider(names_from = TestDay, 
                             values_from = x.axis.var) %>% 
                 mutate(feature = factor(feature, 
                               levels = viol.sig.tax$feature))
   
   
   viol.plot <- ggplot(viol.df, aes(y = Value, x = x.axis.var)) + 
                         geom_point(aes(color = Prebiotic), 
                                    size = 2, 
                                    alpha = 1) +
                         geom_line(aes(group = PersonID, 
                                       color = Prebiotic)) +
                         geom_violin(fill = NA, alpha = 0.1) +
                         geom_signif(data = sig.df,
                                     aes(xmin = TD1,
                                         xmax = TD3,
                                         annotations = anot,
                                         y_position = y),
                                     textsize = 2, vjust = -0.1,
                                     manual = TRUE, margin_top = 1) +
                         geom_point(data = sig.df,
                                    aes(x = TD1, y = y*1.3), x=NA) +
                         facet_grid(.~ feature) + 
                         theme_bw() + 
                         scale_color_manual(values = color.sch.ls$Prebiotic) + 
                         xlab("Test Day") +
                         ylab("Abundance(log)") +
                         theme(axis.title.x = element_text()) + 
                         guides(colour = guide_legend(order = 2), 
                                linetype = guide_legend(order = 1)) +
                         theme(axis.text.x = element_text(angle=45, 
                                                          vjust=1.1, 
                                                          hjust=1), 
                               strip.text = element_text(size = 6)) 
   
   viol.plot.ls[[i.gr]][[i.ps]] <- viol.plot
   
  }}


#-------------------------------------------------------------------------------
# Save DA plot 
#-------------------------------------------------------------------------------
dir.create("out/plot/DA", recursive = TRUE)

for(i in names(sig.heat.ls)) {
  
  for(j in names(sig.heat.ls[[i]])) {
    
    png(paste0("out/plot/DA/DA_heatmap_", i, "_", j, ".png"), 
    width = 12, height = 2.75, units = "in", res = 600)
    
    print(sig.heat.ls[[i]][[j]])
    
    dev.off()
  }
}

comb.heat.ls <- list()

for(i in names(sig.heat.ls)) {
  # Combined heatmap
  comb.heat = sig.heat.ls[[i]][[2]] %v% sig.heat.ls[[i]][[1]]
  
  comb.heat.ls[[i]] <- comb.heat
  
  png(paste0("out/plot/DA/DA_heat_comb_", i, ".png"), 
      width = 12.5, height = 5, units = "in", res = 600)
  
  draw(comb.heat, ht_gap = unit(0.75, "cm"))
  
  dev.off()

}


# Violin plots
for(i in names(viol.plot.ls)) {
  
  for(j in names(viol.plot.ls[[i]])) {
    
    # Violin plots
    ggsave(filename = paste0("out/plot/DA/DA_violin_", i, "_", j, ".png"), 
           plot = viol.plot.ls[[i]][[j]], width = 12, 
           height = 2.75, dpi = 600)
    
  }
}


#-------------------------------------------------------------------------------
# Save data
#-------------------------------------------------------------------------------
res.tab <- list()

for(i1 in names(maas.res.long.df)) {
  
  for(i2 in names(maas.res.long.df[[i1]])) {
    
   i.maas <- maas.res.long.df[[i1]][[i2]] 
   
   for(i3 in unique(i.maas$Strata)) {
     
     i.maas.f <- i.maas %>% 
       filter(Strata == i3) %>% 
       arrange(pval)
     
     res.tab[[i1]][[gsub(".*_", "", i2)]][[i3]] <- i.maas.f
     
   }
   
  }

} 


res.da <- list(Plots = list(Violin = viol.plot.ls, 
                            Heatmap = comb.heat.ls), 
               Results = res.tab, 
               Par = list(qval.cut = qval.cut, 
                          pval.cut = pval.cut, 
                          min.prev = min.prev))


save(list = c("res.da"),
     file = "out/supp/4_DA_res.RData")

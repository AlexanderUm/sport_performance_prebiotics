#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(549794)

lib.to.load <- c("phyloseq", "tidyverse", 
                 "ggsignif", "broom", "vegan", "cowplot")

for (i in lib.to.load) {library(i, character.only = TRUE)}

rm(list = c("i", "lib.to.load"))

# Load data
load("out/supp/0_data.RData")


#-------------------------------------------------------------------------------
# Variables for PERMANOVA
#-------------------------------------------------------------------------------
used.ps <- c("css_asv", "css_genus")

# Calculate distances and statistically test
used.dist <- c("unifrac", "wunifrac", "jaccard", "bray")

used.var <- c("TestDay", "Time", "Treatment", "Sex", "Age",
              "PersonID", "Prebiotic", "TreatmentDay", "Resp")

group.vars <- c("Treatment", "Resp")

nperm <- 999


################################################################################
# Calculate distances for the full data set 
#-------------------------------------------------------------------------------
ls.dist <- list()

for(i.ps in used.ps) {
  
  ps.inst <- pss.ls[[i.ps]]
  
  for (i.dist in used.dist) {
    
    if(i.dist == "jaccard") {
      
      dist.inst <- phyloseq::distance(ps.inst, i.dist, binary = TRUE)
      
    } else {
      
    dist.inst <- phyloseq::distance(ps.inst, i.dist)
      
    }
    
    ls.dist[[i.ps]][[i.dist]] <- dist.inst
  }
  
}


#-------------------------------------------------------------------------------
# PERMANOVA for the full data set
#-------------------------------------------------------------------------------
# Formula list for PERMANOVA testing

ls.form.full <- list("dist ~ PersonID + Time*Treatment", 
                     "dist ~ PersonID + Time*Resp")

res.adon.df <- NULL

res.adon.ls <- list()

for(i.ps in used.ps) {
  
  for (i.dist in used.dist) {
    
    dist <- ls.dist[[i.ps]][[i.dist]]
    
    for (i.form in ls.form.full) {
      
      res.adon <- adonis2(formula = as.formula(i.form), 
                               data = ps.meta, 
                               by = "terms", 
                               permutations = nperm, 
                               parallel = 4) %>% 
                      tidy() %>% 
                      mutate(Distance = i.dist,
                             PS = i.ps,
                             Formula = as.character(i.form)) 
      
      res.adon.df <- bind_rows(res.adon.df, res.adon)
      
      i.form.name <- gsub(".*\\*", "", i.form)
      
      res.adon.ls[[i.form.name]][[i.ps]][[i.dist]] <- res.adon
        
}}}# End of for the loop


################################################################################
# PERMANOVA subset by Time point and Treatment
#-------------------------------------------------------------------------------
used.starata.col <- c("TestDay")
                      
test.col <- c("Treatment", "Resp")


res.adon.str.df <- NULL

res.adon.str.ls <- NULL

for(i.str.col in used.starata.col) {

  for(i.str.lvl in unique(ps.meta[[i.str.col]])) {
    
    # Subset metadata
    ps.meta.inst <- ps.meta[ps.meta[[i.str.col]] == i.str.lvl, ]
    
    for(i.ps in used.ps) {
      
      # Subset phyloseq object
      ps.inst <- pss.ls[[i.ps]] %>% 
                    prune_samples(rownames(ps.meta.inst), .)
      
      for (i.dist in used.dist) {
        
        # Calculate distance
        if(i.dist == "jaccard") {
          
          dist <- phyloseq::distance(ps.inst, i.dist, binary = TRUE)
          
        } else {
          
          dist <- phyloseq::distance(ps.inst, i.dist)
          
        }
        
        for(i.gr in test.col) {
          
        i.form <- paste0("dist ~ ", i.gr)

        res.adon.str <- adonis2(formula = as.formula(i.form), 
                                         data = ps.meta.inst, 
                                         by = "terms", 
                                         permutations = nperm, 
                                         parallel = 4)  %>% 
                                tidy() %>% 
                                mutate(Distance = i.dist,
                                       PS = i.ps,
                                       Formula = as.character(i.form), 
                                       Strata = paste0(i.str.col, "--", 
                                                       i.str.lvl))  
        
        res.adon.str.df <- bind_rows(res.adon.str.df, res.adon.str)
        
        res.adon.str.ls[[i.gr]][[i.ps]][[i.dist]][[i.str.lvl]] <- res.adon.str
        
}}}}} # End of for the loop
  
  
################################################################################
# Plot ordination
#-------------------------------------------------------------------------------
dists.names <- c("Unweighted UniFrac", "Weighted UniFrac", 
                 "Jaccard", "Bray-Curtis") %>% 
                 setNames(used.dist)


pcoa.plots.ls <- list() 

for(i.ps in used.ps) {
  
  for (i.dist in used.dist) { 
    
    dist <- ls.dist[[i.ps]][[i.dist]]
    
    # Make PCoA object
    a.pcoa <- ape::pcoa(dist)
    
 
    axis.pcoa <- a.pcoa$vectors[, 1:2] %>%
                        as.data.frame() %>%
                        bind_cols(., ps.meta[, used.var])
    
    # Extract percentage
    var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
    
    pcoa.p <- ggplot(axis.pcoa,
                         aes(x = Axis.1, 
                             y = Axis.2)) +
                    stat_ellipse(geom = "polygon", 
                                 aes(fill = .data[["TreatmentDay"]]), 
                                 alpha = 0.25, level = 0.90, size =2) +
                    geom_path(aes(color = .data[["Prebiotic"]],
                                  group = .data[["PersonID"]]),
                              alpha = 0.5) +
                    geom_point(aes(color = .data[["Prebiotic"]],
                                   shape = .data[["TestDay"]]),
                               alpha = 0.75, size = 2) +
                    ggtitle(paste0("Distance: ", dists.names[i.dist])) +
                    scale_color_manual(values = color.sch.ls$Prebiotic) +
                    scale_fill_manual(values = color.sch.ls$TreatmentDay) +
                    coord_fixed() +
                    theme_bw() +
                    xlab(paste0(colnames(axis.pcoa)[1], " [", var.c[1], "%]")) +
                    ylab(paste0(colnames(axis.pcoa)[2], " [", var.c[2], "%]")) +
                    guides(color = guide_legend(order=1),
                           fill = guide_legend(order=3),
                           shape = guide_legend(order=2)) + 
                    labs(fill = "Treatment: Test Day",
                         color = "Prebiotic",
                         shape = "Test Day") + 
                    theme(plot.title = element_text(size=11, face="italic"))
    
    
    pcoa.plots.ls[[i.ps]][[i.dist]] <- pcoa.p
    
}}

     
#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
pcoa.plots.ls.comb <- list()

for (i.p.ls in names(pcoa.plots.ls)) {
  
  p.ls.inst <- pcoa.plots.ls[[i.p.ls]]

  p.ls.inst.f <- lapply(p.ls.inst, function(x){x+theme(legend.position="none")})
  
  legd.inst <- get_legend(p.ls.inst[[1]])

  p.grid.1 <- plot_grid(plotlist = p.ls.inst.f, 
                        ncol = 2)
  
  p.grid.2 <- plot_grid(p.grid.1, 
                        legd.inst, 
                        ncol = 2, 
                        rel_widths = c(0.8, .2))
  
  pcoa.plots.ls.comb[[i.p.ls]] <- p.grid.2
}


#-------------------------------------------------------------------------------
# Combine Alpha Diversity and Beta diversity plots
#-------------------------------------------------------------------------------
load("out/supp/2_alpha_res.RData")

alpha.gr <- plot_grid(res.alpha$Plot$Treatment + theme(legend.position="none"),
                      NULL, nrow = 2, rel_heights = c(0.99, 0.01))

div.plot.comb <- plot_grid(alpha.gr,
                      NULL,
                      pcoa.plots.ls.comb$css_asv, 
                      nrow = 1, 
                      rel_widths = c(0.275, 0.02, 0.725),
                      labels = c("A", "", "B"))

save_plot("out/plot/fig_3_diversity_beta.png", 
       plot = plot(div.plot.comb), 
       base_width = 12, base_height = 7, units = "in", dpi = 600)


#-------------------------------------------------------------------------------
# Write results
#-------------------------------------------------------------------------------
res.beta <- list(Plot = list("Combined diversity" = div.plot.comb, 
                             "Combined beta diversity" = pcoa.plots.ls.comb, 
                             "Individual beta diversity" = pcoa.plots.ls), 
                 Results = list("Full" = res.adon.ls,
                                "Stratified" = res.adon.str.ls), 
                 Par = list("nperm" = nperm))

save(list = c("res.beta"),
     file = "out/supp/3_beta_res.RData")



 
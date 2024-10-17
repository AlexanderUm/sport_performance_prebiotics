#-------------------------------------------------------------------------------
# Load data and functions
#-------------------------------------------------------------------------------
load("out/supp/prm.R")
load("out/supp/0_data.RData")
load("out/supp/2_alpha_res.RData")

source("R/phy_dists_ls.R")

#-------------------------------------------------------------------------------
# Parameters 
#-------------------------------------------------------------------------------
gr.cols <- prm.ls[["general"]][["Group_cols"]]

id.col <- prm.ls[["general"]][["Prat_ID_col"]]

preb.col <- prm.ls[["general"]][["Preb_col"]]

part.id.col <- prm.ls[["general"]][["Prat_ID_col"]]

seed <- prm.ls[["general"]][["seed"]]

time.col.fact <- prm.ls[["general"]][["Time_col_fact"]]

time.col.num <- prm.ls[["general"]][["Time_col_num"]]

tax.lvl <- prm.ls[["beta"]][["Tax_lvl"]]

norm <- prm.ls[["beta"]][["Norm"]]

dists.vec <- prm.ls[["beta"]][["distances"]]

nperm <- prm.ls[["beta"]][["n_perm"]]

dir.out <- prm.ls[["beta"]][["out_dir_path"]]

set.seed(seed)

res.beta.comb <- list()

################################################################################
# Calculate distances for the full data set 
#-------------------------------------------------------------------------------
par.gird <- expand.grid("Tax_lvl" = tax.lvl, 
                        "Normal" = norm, 
                        stringsAsFactors = FALSE)

ls.dist <- list()

for(i.prm in 1:nrow(par.gird)) {
  
  lvl.inst <- par.gird[i.prm, "Tax_lvl"]
  
  norm.inst <- par.gird[i.prm, "Normal"] 
  
  ps.inst <- ps.ls[[lvl.inst]][[norm.inst]]
  
  dist.inst <- phy_dist_ls(ps.inst, dists = dists.vec) %>% 
                           setNames(names(dists.vec))
  
  ls.dist[[lvl.inst]][[norm.inst]] <- dist.inst
  
  rm(list = c("lvl.inst", "norm.inst", "ps.inst", "dist.inst"))
}
  

################################################################################
# PERMANOVA
#-------------------------------------------------------------------------------
# Define parameters 
def.prm.ls <- list("Tax_lvl" = tax.lvl, 
                   "Normal" = norm,
                   "Distances" = names(dists.vec))

prm.gird.perm <- NULL

for(i.gr in gr.cols) {
  
  def.prm.ls[["Group"]] <- i.gr
  
  # All samples 
  def.prm.ls[["Strata"]] <- "All" 
  
  def.prm.ls[["Strata_col"]] <- NA
  
  prm.grid.inst <- expand.grid(def.prm.ls, stringsAsFactors = FALSE) 
  
  prm.grid.inst$Formula <- paste0("dist ~ ", 
                                  id.col, " + ", 
                                  time.col.num, "*", 
                                  i.gr)
  
  prm.gird.perm <- rbind(prm.gird.perm, prm.grid.inst)  
  
  # Strata per group 
  def.prm.ls[["Strata"]] <- levels(ps.meta[[i.gr]])
  
  def.prm.ls[["Strata_col"]] <- i.gr
  
  prm.grid.inst <- expand.grid(def.prm.ls, stringsAsFactors = FALSE) 
  
  prm.grid.inst$Formula <- paste0("dist ~ ", 
                                  id.col, " + ", 
                                  time.col.num)
  
  prm.gird.perm <- rbind(prm.gird.perm, prm.grid.inst)
  
  # Time 
  def.prm.ls[["Strata"]] <- levels(ps.meta[[time.col.fact]])
  
  def.prm.ls[["Strata_col"]] <- time.col.fact

  prm.grid.inst <- expand.grid(def.prm.ls, stringsAsFactors = FALSE) 
  
  prm.grid.inst$Formula <- paste0("dist ~ ", i.gr)
  
  prm.gird.perm <- rbind(prm.gird.perm, prm.grid.inst)
  
  rm(list = c("prm.grid.inst", "i.gr"))
  
}


# Run PERMANOVA
res.adon.ls <- list()

for(i.prm in 1:nrow(prm.gird.perm)) {
  
  # Parameters
  lvl.inst <- prm.gird.perm[i.prm, "Tax_lvl"]
  
  norm.inst <- prm.gird.perm[i.prm, "Normal"] 
  
  dist.mesure.inst <- prm.gird.perm[i.prm, "Distances"]
  
  strata.inst <- prm.gird.perm[i.prm, "Strata"]
  
  strata.col.inst <- prm.gird.perm[i.prm, "Strata_col"]
  
  gr.inst <- prm.gird.perm[i.prm, "Group"]
  
  form.inst <- prm.gird.perm[i.prm, "Formula"]
  
  dist.inst <- prm.gird.perm[i.prm, "Distances"]
  
  # Extract distance and metadata 
  meta.inst <- ps.meta 
  
  dist <- ls.dist[[lvl.inst]][[norm.inst]][[dist.inst]]
  
  if(!is.na(strata.col.inst)) {
    
    meta.inst <- meta.inst %>% 
                  filter(.data[[strata.col.inst]] == strata.inst)
    
    dist <- dist_subset(dist, rownames(meta.inst))
    
  }
  
  # Run adonis2
  res.adon <- adonis2(formula = as.formula(form.inst), 
                      data = meta.inst, 
                      by = "terms", 
                      permutations = nperm, 
                      parallel = 4) %>% 
                tidy() %>% 
                suppressWarnings() %>% 
                bind_cols(., prm.gird.perm[i.prm, ]) %>% 
                add_row()
  
  res.adon.ls[[lvl.inst]][[norm.inst]][[gr.inst]][[strata.inst]][[dist.inst]] <- 
    res.adon
}


# Combine and Write out results 
comb.prm.grid <- prm.gird.perm %>% 
                    select(-all_of(c("Strata_col", "Formula", "Distances"))) %>% 
                    unique()

for(i.prm in 1:nrow(comb.prm.grid)) {
  
  lvl.inst <- comb.prm.grid[i.prm, "Tax_lvl"]
  
  norm.inst <- comb.prm.grid[i.prm, "Normal"]
  
  gr.inst <- comb.prm.grid[i.prm, "Group"]
  
  strata.inst <- comb.prm.grid[i.prm, "Strata"]
  
  comb.res <- res.adon.ls[[lvl.inst]][[norm.inst]][[gr.inst]][[strata.inst]] %>% 
                bind_rows() %>% 
                as.data.frame()
  
  comb.res[is.na(comb.res)] <- ""
  
  dir.create(paste0(dir.out, "/tabs/"), recursive = TRUE, showWarnings = FALSE)
  
  write.csv(comb.res, 
            paste0(dir.out, "/tabs/adonis_", 
                   paste(comb.prm.grid[i.prm, ], collapse = "--"), 
                   ".csv"), 
            row.names = FALSE)
  
  res.beta.comb[[lvl.inst]][[norm.inst]][[gr.inst]][["Stat"]][[strata.inst]] <- comb.res
  
}


################################################################################
# Plot ordination
#-------------------------------------------------------------------------------
plot.prm.grid <- expand.grid("Taxa_lvl" = tax.lvl, 
                             "Normalization" = norm, 
                             "Distances" = names(dists.vec), 
                             "Group" = gr.cols,
                             stringsAsFactors = FALSE)

pcoa.plots.ls <- list() 

for(i.prm in 1:nrow(plot.prm.grid)) {

lvl.inst <- plot.prm.grid[i.prm, "Taxa_lvl"]

norm.inst <- plot.prm.grid[i.prm, "Normalization"]

dist.inst <- plot.prm.grid[i.prm, "Distances"]

gr.inst <- plot.prm.grid[i.prm, "Group"]

dist <- ls.dist[[lvl.inst]][[norm.inst]][[dist.inst]]
    
# Make PCoA object
a.pcoa <- ape::pcoa(dist)
    
axis.pcoa <- a.pcoa$vectors[, 1:2] %>%
                        as.data.frame() %>%
                        bind_cols(., ps.meta)
    
# Extract percentage
var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
    
pcoa.p <- ggplot(axis.pcoa,
                         aes(x = Axis.1, 
                             y = Axis.2)) +
                    stat_ellipse(geom = "polygon", 
                                 aes(fill = .data[[paste0(gr.inst, time.col.fact)]]), 
                                 alpha = 0.25, level = 0.90, size =2) +
                    geom_path(aes(color = .data[[preb.col]],
                                  group = .data[[id.col]]),
                              alpha = 0.5) +
                    geom_point(aes(color = .data[[preb.col]],
                                   shape = .data[[time.col.fact]]),
                               alpha = 0.75, size = 2) +
                    ggtitle(paste0("Distance: ", dist.inst)) +
                    scale_color_manual(values = color.sch.ls[[preb.col]]) +
                    scale_fill_manual(values = color.sch.ls[[paste0(gr.inst, time.col.fact)]]) +
                    coord_fixed() +
                    theme_bw() +
                    xlab(paste0(colnames(axis.pcoa)[1], " [", var.c[1], "%]")) +
                    ylab(paste0(colnames(axis.pcoa)[2], " [", var.c[2], "%]")) +
                    guides(color = guide_legend(order=1),
                           fill = guide_legend(order=3),
                           shape = guide_legend(order=2)) + 
                    labs(fill = paste0(gr.inst, ": Test Day"),
                         color = "Prebiotic",
                         shape = "Test Day") + 
                    theme(plot.title = element_text(size=14, face="italic"),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), 
                          axis.title = element_text(size = 12), 
                          axis.text = element_text(size = 12), 
                          legend.text = element_text(size = 12), 
                          legend.title = element_text(size = 12, face = "bold"))
    
    
    pcoa.plots.ls[[lvl.inst]][[norm.inst]][[gr.inst]][[dist.inst]] <- pcoa.p
    
}

     
#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
plot.prm.grid <- expand.grid("Taxa_lvl" = tax.lvl, 
                             "Normalization" = norm, 
                             "Group" = gr.cols,
                             stringsAsFactors = FALSE)

pcoa.plots.ls.comb <- list()

for (i.prm in 1:nrow(plot.prm.grid)) {
  
  # Parameters 
  lvl.inst <- plot.prm.grid[i.prm, "Taxa_lvl"]
  
  norm.inst <- plot.prm.grid[i.prm, "Normalization"]
  
  gr.inst <- plot.prm.grid[i.prm, "Group"]
  
  # Beta plots 
  p.ls.inst <- pcoa.plots.ls[[lvl.inst]][[norm.inst]][[gr.inst]]

  p.ls.inst.f <- lapply(p.ls.inst, function(x){x+theme(legend.position="none")})
  
  legd.inst <- suppressWarnings(get_legend(p.ls.inst[[1]]))

  p.grid.1 <- plot_grid(plotlist = p.ls.inst.f, 
                        ncol = 2)
  
  p.grid.2 <- plot_grid(p.grid.1, 
                        legd.inst, 
                        ncol = 2, 
                        rel_widths = c(0.8, .2))
  
  # Combine with alpha 
  alpha.p.ls.inst <- alpha.res.ls[["Plots"]][[gr.inst]][[1]]
  
  alpha.gr <- plot_grid(alpha.p.ls.inst+
                       theme(legend.position="none"),
                       NULL, nrow = 2, rel_heights = c(0.99, 0.01))
  
  div.plot.comb <- plot_grid(alpha.gr,
                             NULL,
                             p.grid.2, 
                             nrow = 1, 
                             rel_widths = c(0.275, 0.02, 0.725),
                             labels = c("A", "", "B"))
  
  
  dir.create(paste0(dir.out, "/plots/"), recursive = TRUE, showWarnings = FALSE)
  
  save_plot(paste0(dir.out, "/plots/divers_comb--", 
                   gr.inst, "--", lvl.inst, ".png"),
            plot = plot(div.plot.comb), 
            base_width = 12, base_height = 7, units = "in", dpi = 600)
  
  
  res.beta.comb[[lvl.inst]][[norm.inst]][[gr.inst]][["Plots"]][["only_beta"]] <- 
      p.grid.2
  
  res.beta.comb[[lvl.inst]][[norm.inst]][[gr.inst]][["Plots"]][["with_alpha"]] <- 
      div.plot.comb
}


#-------------------------------------------------------------------------------
# Write out results 
save(list = c("res.beta.comb"),
     file = "out/supp/3_beta_res.RData")

rm(list = ls())
 
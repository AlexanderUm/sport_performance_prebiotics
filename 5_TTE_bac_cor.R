#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(549794)

lib.to.load <- c("phyloseq", "tidyverse", "ggsignif", "broom", "lmerTest")

for (i in lib.to.load) {library(i, character.only = TRUE)}

rm(list = c("i", "lib.to.load"))

# Load data
load("out/supp/0_data.RData")

source("R/filt_tax_phy.R")

#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------
min.prev <- 0.25

min.prev2 <- 0.5

sample.sets <- list(c("CON", "PRE"), "CON", "PRE")

tax.lvl.delta <- c("css_genus", "css_asv")

n.plot <- 10

#-------------------------------------------------------------------------------
# Correlation between delta abundance and delta TTE
#-------------------------------------------------------------------------------
delta.cor.res.ls <- list()

delta.cor.plot.ls <- list()

for(i.lvl in tax.lvl.delta) {
  
  # Extract otu table
  otu.inst <- pss.ls[[i.lvl]] %>% 
                  tax_filt_phy(., prev = min.prev) %>%
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame() 
  
  # Add meta and calculate delta abundance between time point one and two
  otu.inst.delta <- otu.inst %>% 
                      bind_cols(., ps.meta[, c("PersonID", 
                                               "TestDay", 
                                               "Treatment", 
                                               "TTE")]) %>%
                      group_by(PersonID) %>%
                      arrange(TestDay, .by_group = TRUE) %>% 
                      mutate(across(where(is.numeric), 
                                       function(x){x[2] - x[1]})) %>% 
                      ungroup() %>% 
                      filter(TestDay == "TD1") 
  
  for(i.set in sample.sets) {
    
    otu.inst.delta.set <- otu.inst.delta %>% 
                            filter(Treatment %in% i.set) 
    
    otu.inst.set <- otu.inst %>%
                      filter(ps.meta$Treatment %in% i.set) 
    
    # Correlation 
    cor.res.inst <- data.frame()
    
    for (i.tax in colnames(otu.inst)) {
      
      cor.res <- cor.test(otu.inst.delta.set[[i.tax]], 
                             otu.inst.delta.set$TTE, 
                             method = "spearman")
      
      # Collect data 
      res.inst <- c(Taxa = i.tax, 
                    rho = as.numeric(cor.res$estimate), 
                    pval = cor.res$p.value, 
                    Variable = "delta TTE",
                    Group = paste(i.set, collapse = "&"), 
                    "Taxa Prev" = sum(otu.inst.set[, i.tax] != 0)/nrow(otu.inst.set),
                    "Taxa Level" = i.lvl, 
                    method = cor.res$method) 
      
      cor.res.inst <- bind_rows(cor.res.inst, res.inst)
      
    }
    
    cor.res.inst <- cor.res.inst %>% 
                          mutate(across(c("rho", "pval", "Taxa Prev"), 
                                        as.numeric)) %>%
                          filter(.data[["Taxa Prev"]] >= min.prev2) %>%
                          arrange(desc(abs(rho))) 
    
    delta.cor.res.ls[[i.lvl]][[paste(i.set, collapse = "_")]] <- cor.res.inst
    
    
    # Plot scatterplot 
    p.tax <- cor.res.inst[1:n.plot, "Taxa"]
    
    p.df <- otu.inst.delta.set[, c("TTE", p.tax)] %>% 
              pivot_longer(-TTE, names_to = "Taxa", values_to = "Abundance") %>% 
              mutate(Taxa = factor(Taxa, levels = p.tax))
    
    p <- ggplot(p.df, aes(x = TTE, y = Abundance)) +
                                geom_smooth(method = "lm") +
                                geom_point() +
                                facet_wrap(~Taxa, scales = "free", ncol = 5) + 
                                theme_bw()
    
    delta.cor.plot.ls[[i.lvl]][[paste(i.set, collapse = "_")]] <- p
    
  }
}


#-------------------------------------------------------------------------------
# Correlation between abundance and TTE at test day
#-------------------------------------------------------------------------------
cor.res.ls.tp <- list()

cor.plot.ls.tp <- list()

tax.lvl.tp <- c("css_genus", "css_asv", "rare_log_asv", "rare_log_genus")

test.var <- "TTE"


for(i.lvl in tax.lvl.tp) {
  
  # Extract otu table
  otu.inst <- pss.ls[[i.lvl]] %>% 
                  tax_filt_phy(., prev = min.prev) %>%
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame() 
  
  
  for (i.tp in unique(ps.meta$Timepoint)) {
    
    otu.inst.tp <- otu.inst %>% 
                    bind_cols(ps.meta[test.var]) %>%
                    filter(ps.meta$TestDay == i.tp)
    
    
    cor.res.inst <- data.frame()
    
    for (i.tax in colnames(otu.inst)) {
      
      cor.res <- cor.test(otu.inst.tp[[i.tax]], 
                          otu.inst.tp[[test.var]], 
                          method = "spearman")
      
      # Collect data 
      res.inst <- c(Taxa = i.tax, 
                    rho = as.numeric(cor.res$estimate), 
                    pval = cor.res$p.value, 
                    Variable = test.var,
                    TestDay = i.tp, 
                    "Taxa Prev" = sum(otu.inst.tp[, i.tax] != 0)/nrow(otu.inst.tp),
                    "Taxa Level" = i.lvl, 
                    method = cor.res$method) 
      
      cor.res.inst <- bind_rows(cor.res.inst, res.inst)
      
    }
    
    cor.res.inst <- cor.res.inst %>% 
                          mutate(across(c("rho", "pval", "Taxa Prev"), 
                                                      as.numeric)) %>%
                          filter(.data[["Taxa Prev"]] >= min.prev2) %>%
                          arrange(desc(abs(rho))) 
    
    cor.res.ls.tp[[i.lvl]][[i.tp]] <- cor.res.inst
    
    # Plot scatter plot 
    p.tax <- cor.res.inst[1:n.plot, "Taxa"]
    
    p.df <- otu.inst.tp[, c(test.var, p.tax)] %>% 
                  pivot_longer(-test.var, 
                               names_to = "Taxa", 
                               values_to = "Abundance") %>% 
                  mutate(Taxa = factor(Taxa, levels = p.tax))
    
    p <- ggplot(p.df, aes(x = .data[[test.var]], 
                          y = Abundance)) +
                  geom_smooth(method = "lm") +
                  geom_point() +
                  facet_wrap(~Taxa, scales = "free", ncol = 5) + 
                  theme_bw()
    
    cor.plot.ls.tp[[i.lvl]][[i.tp]] <- p
    
  }
}


res.cor <- list(Delta = list(Plot = delta.cor.plot.ls, 
                             Results = delta.cor.res.ls), 
                "Time Point" = list(Plot = cor.plot.ls.tp, 
                                    Results = cor.res.ls.tp), 
                Par = list("Minimum prevalence" = min.prev2))

save(list = c("res.cor"), file = "out/supp/5_cor.Rdata")

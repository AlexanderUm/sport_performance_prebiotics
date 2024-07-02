#-------------------------------------------------------------------------------
# Calculate prevalece of taxa across a phyloseq
#-------------------------------------------------------------------------------

phy_prevalence_calc <- function(phy_raw_count, 
                                per_group_cols = NULL) {
  
  otu.tab <- otu_table(phy_raw_count) %>% 
              as.matrix() %>% 
              as.data.frame()
  
  if(taxa_are_rows(phy_raw_count)) {
    
    otu.tab <- as.data.frame(t(otu.tab))
    
  }
  
  prev.tab <- data.frame(across_all = colSums(otu.tab != 0)/nrow(otu.tab))
  
  if(!is.null(per_group_cols)) {
    
    meta.tab <- sample_data(phy_raw_count) %>% 
                  as.matrix() %>% 
                  as.data.frame()
    
    for(i.col in per_group_cols) {
      
      for(i.col.lvl in unique(meta.tab[[i.col]])) {
        
        otu.tab.f <- otu.tab[meta.tab[[i.col]] == i.col.lvl, ]
        
        prev.tab[, paste0(i.col, "--", i.col.lvl)] <- 
          colSums(otu.tab.f != 0)/nrow(otu.tab.f)
        
      }
      
    }
    
  }
  
 return(prev.tab)
  
}
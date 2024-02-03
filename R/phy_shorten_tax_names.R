#-----------------------------------------------------------------------------
# Function to make short names. 
#-----------------------------------------------------------------------------

phy_shorten_tax_names <- function(phy_obj) {
  
  tax.t <- as(tax_table(phy_obj), "matrix")
  
  tax.t[tax.t == ""] <- NA
  
  tax.out <- c()
  
  excle.names <- paste(c("_NA", "uncultured", "gut metagenome",
                        "metagenome", "human_gut", "uncultured_archaeon", 
                         "<empty>", "Incertae_Sedis"), 
                       collapse = "|")
    
  tax.lvl.pref <- substr(colnames(tax.t), 1, 1)
  
  ucg.pref <- paste(paste0(tax.lvl.pref, "_UCG"), collapse = "|")
  
  
  for (r in 1:nrow(tax.t)) {
    
    ind.tax <- paste0(tax.lvl.pref, "_", tax.t[r, ])
    
    l.tax <-  tail(ind.tax[!grepl(excle.names, 
                                  ind.tax)], 1)
    
    if(length(l.tax) == 0) {
      
      l.tax <- NA
      
    }
    
    # Unidentified genus fix 
    
    if(grepl(ucg.pref, l.tax)) {
      
      no.ucg.tax <- ind.tax[!grepl(ucg.pref, ind.tax)] 
      
      l.tax.p <- tail(no.ucg.tax[!grepl(excle.names, 
                                       no.ucg.tax)], 1)
      
      l.tax <- paste0(l.tax.p, "_", gsub(".*_", "", l.tax))
      
    }
    
    
    tax.out <- c(tax.out, gsub(" ", ".", l.tax))
    
  }
  
  tax.out <- gsub("-", "_", tax.out)
  
  tax.out <- gsub("\\[|\\]", "", tax.out)
  
  return(tax.out)
}


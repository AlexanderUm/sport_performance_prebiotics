
#-------------------------------------------------------------------------------
# Function: Filter phyloseqs based on the prevalence and 
#           root the tree at midpoint.  
#-------------------------------------------------------------------------------
tax_filt_phy <- function(phy_obj, 
                          prev = 0, 
                          abund = 0, 
                          prev_is_fraction = TRUE, 
                          abund_is_fraction = TRUE,
                          group_col = NULL) {
  
  # Required libraries 
  require(phyloseq)
  require(ape)
  require(phangorn)
  require(dplyr)
  
  #-----------------------------------------------------------------------------
  # Messages 
  #-----------------------------------------------------------------------------
  
  if (!is.null(group_col)) {
    
    message("Groups column is selected. Minimal prevalece and abundance are
              applyed per group!")
  }
  
  
  #-----------------------------------------------------------------------------
  # Extract objects 
  #-----------------------------------------------------------------------------
  
  # OTU table 
  otu <- as(otu_table(phy_obj), "matrix")
  
  # Transpose if necessary
  if(taxa_are_rows(phy_obj)){otu <- t(otu)}
  
  # Make binary OTU table
  otu.b <- otu
  
  otu.b[otu.b > 0] <- 1
  
  # If abundance is fraction convert to proportion 
  if (abund_is_fraction) {   
    
    otu <- t(apply(otu, 1, function(x){x/sum(x)})) 
    
    }
  
  # Metadata column 
  if(!is.null(group_col)) {
    
    meta.d <- phy_obj |> 
              sample_data() |> 
              as.matrix() |>
              as.data.frame()
    
    meta.col <- meta.d[, group_col]
  
  }
  
  #-----------------------------------------------------------------------------
  # Filter based on prevalence 
  #-----------------------------------------------------------------------------
  
  # Start: If group is NOT provided 
  if (is.null(group_col)) {
    
    # For prevalence - count per OTU
    if (prev_is_fraction) { 
      
      prev.count <- colSums(otu.b) / nrow(otu.b)
    
    } else { prev.count <- colSums(otu.b) }
    
    tax.to.keep.prev <- names(prev.count)[prev.count > prev]
    
    
    # For abundance - count per OTU
    abund.count <- colSums(otu)
    
    tax.to.keep.abund <- names(abund.count)[abund.count > abund]
    
    
    tax.to.keep <- intersect(tax.to.keep.abund, 
                             tax.to.keep.prev)
    
    # Finish: If group is NOT provided
    # Start : If group is provided 
  } else {
    
    # Empty object for prevalence 
    tax.to.keep.prev <- c()
    
    # Empty object for count 
    tax.to.keep.abund <- c()
    
    # Start: Loop through levels in the provided metadata column 
    for (lvl in unique(meta.col)) {
      
      # For prevalence 
      otu.prev.lvl <- otu.b[lvl %in% meta.col, ]
      
      if(prev_is_fraction) {
        
        prev.count.lvl <- colSums(otu.prev.lvl) / nrow(otu.prev.lvl)
        
      } else {
        
        prev.count.lvl <- colSums(otu.prev.lvl)
        
      }
      
      tax.to.keep.prev.lvl <- names(prev.count.lvl)[prev.count.lvl > prev]
      
      tax.to.keep.prev <- c(tax.to.keep.prev, tax.to.keep.prev.lvl)
      
      
      # For abundance 
      
      otu.abund.lvl <- otu[lvl %in% meta.col, ]
      
      abund.count.lvl <- colSums(otu.abund.lvl)
      
      tax.to.keep.abund.lvl <- names(abund.count.lvl)[abund.count.lvl > abund]
      
      tax.to.keep.abund <- c(tax.to.keep.abund, tax.to.keep.abund.lvl)
      
      
    } # Finish: Loop through levels in the provided metadata column 
    
    tax.to.keep.prev <- unique(tax.to.keep.prev)
    
    tax.to.keep.abund <- unique(tax.to.keep.abund)
    
    tax.to.keep <- intersect(tax.to.keep.abund, 
                             tax.to.keep.prev)
    
  } # Finish : If group is provided 
  
  phy.obj.filt <- prune_taxa(tax.to.keep, phy_obj)
  
  return(phy.obj.filt)
  
} # END OF FUNCTION


  


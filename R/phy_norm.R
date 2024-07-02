#-------------------------------------------------------------------------------
# Normalize ASV count with CSS
#-------------------------------------------------------------------------------
phy_css_norm <- function(phyloseq, 
                         log2_transform = TRUE) {
  
  require("phyloseq")
  require("tidyverse")
  require("metagenomeSeq")
  
  otu_table(phyloseq) <- phyloseq %>% 
                        phyloseq_to_metagenomeSeq(.) %>% 
                        cumNorm(., p=cumNormStatFast(.)) %>% 
                        MRcounts(., norm=TRUE, log=log2_transform) %>% 
                        as.data.frame() %>% 
                        otu_table(., taxa_are_rows = TRUE)
  
  return(phyloseq)
  
}


#-------------------------------------------------------------------------------
# Normalize ASV count with TSS
#-------------------------------------------------------------------------------
phy_tss_norm <- function(phyloseq, 
                         log2_transform = TRUE) {
  
  require("phyloseq")
  require("tidyverse")
  require("vegan")
  
  otu.tab <- phyloseq %>% 
                otu_table() %>% 
                as.matrix()
  
  if(taxa_are_rows(phyloseq)) { MARGIN = 2
    
  } else { MARGIN = 1}
  
  otu.tab.trans <- decostand(otu.tab, MARGIN, method = "total") %>% 
                      as.data.frame()
  
  if(log2_transform) {
    
    otu.tab.trans <-  replace(otu.tab.trans, 
                                otu.tab.trans==0, 
                                min(otu.tab.trans[otu.tab.trans>0])/2) %>% 
                      log2(.)
  }
  
  otu_table(phyloseq) <- otu_table(otu.tab.trans, 
                                   taxa_are_rows = taxa_are_rows(phyloseq))
  
  return(phyloseq)
  
}


phy_clr_norm <- function(phyloseq) {
  
  require("phyloseq")
  require("tidyverse")
  require("vegan")
  
  otu.tab <- phyloseq %>% 
    otu_table() %>% 
    as.matrix()
  
  if(taxa_are_rows(phyloseq)) { MARGIN = 2
  
  } else { MARGIN = 1}
  
  otu.tab.trans <- decostand(otu.tab, 
                             MARGIN, 
                             method = "clr", 
                             pseudocount = 1) %>% 
                    as.data.frame()
  
  otu_table(phyloseq) <- otu_table(otu.tab.trans, 
                                   taxa_are_rows = taxa_are_rows(phyloseq))
  
  return(phyloseq)
}
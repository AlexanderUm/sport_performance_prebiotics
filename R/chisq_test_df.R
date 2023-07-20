
chisq_test_df <- function(your_df) {
  
  # Make a list with unique combination of column names 
  col.name.comb <- expand.grid(colnames(your_df), 
                               colnames(your_df)) %>% 
                  filter(Var1 != Var2) %>% 
                  mutate(Count = 1:nrow(.)) %>% 
                  group_by(Count) %>%
                  group_split()
  
  col.name.comb <- lapply(col.name.comb, 
                          function(x){select(x, -"Count") %>% 
                              mutate_all(as.character) %>% 
                              unlist(., use.names = FALSE) %>% 
                              sort()})
  
  col.name.comb <- unique(col.name.comb)          
  
  #-----------------------------------------------------------------------------
  # Run test 
  #-----------------------------------------------------------------------------
  
  test.res <- NULL
  
  for(i in col.name.comb) {
    
    test.res <- chisq.test(your_df[,i[1]], your_df[,i[2]]) %>% 
      tidy() %>% 
      mutate(Var1 = i[1], Var2 = i[2], .) %>% 
      rbind(test.res, .)
  }
  
  test.res <- select(test.res, c("Var1", "Var2", "p.value", 
                                 "statistic", "method"))
  
  return(test.res)
  
}
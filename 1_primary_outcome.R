#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(549794)

lib.to.load <- c("phyloseq", "tidyverse", "ggsignif", "broom")

for (i in lib.to.load) {library(i, character.only = TRUE)}

rm(list = c("i", "lib.to.load"))

# Load data
load("out/supp/0_data.RData")


#-------------------------------------------------------------------------------
# LMM for data
#-------------------------------------------------------------------------------
samp.to.remove <- c("R04", "R06", "R33")

prime.data <- meta.data %>% 
                  filter(!PersonID %in% samp.to.remove)

mod.rand.full <- "TTE ~ Age + Sex + Time*Treatment + (1|PersonID)"

# Fit model 
r.mod.full <- lmerTest::lmer(mod.rand.full, data = prime.data)

mod.sum <- summary(r.mod.full)

# Pairwise comparisons
paired.mod.sum <- emmeans::emmeans(r.mod.full, 
                                   pairwise ~ Time|Treatment, 
                                   type = "response")


#-------------------------------------------------------------------------------
# Model assessment  
#-------------------------------------------------------------------------------
# Homogeneity of variance plot.
homo.var <- plot(r.mod.full)

# QQ plot of the residuals.
tdat <- data.frame(predicted=predict(r.mod.full), 
                   residual = residuals(r.mod.full))

res.qq <- ggplot(tdat,aes(sample=residual)) + 
              stat_qq() + 
              stat_qq_line() + 
              theme_bw()

#-------------------------------------------------------------------------------
# Visualization TEE
#-------------------------------------------------------------------------------
# Data for significance 
max.div <- prime.data %>% 
                dplyr::select(c("TTE", "Treatment")) %>% 
                group_by(Treatment) %>% 
                slice(which.max(TTE)) %>% 
                mutate(y.adj = (TTE*1.1))

sig.df <- paired.mod.sum$contrasts %>% 
                as.data.frame() %>% 
                mutate(p.short = paste0("p=", round(p.value, 4)), 
                       Start = "TD1", 
                       End = "TD3") %>% 
                left_join(., max.div, by = "Treatment")


# Plot data
perf.dif <- ggplot(prime.data, aes(y = TTE, x = TestDay)) + 
                      geom_point(aes(color = Prebiotic), 
                                 size = 2, 
                                 alpha = 1) +
                      geom_line(aes(group = PersonID, 
                                    color = Prebiotic, 
                                    linetype = TTE_Change), 
                                alpha = 0.5, 
                                size = 0.6) + 
                      geom_violin(fill = NA, 
                                  alpha = 0.1, 
                                  fatten = 0.5) +
                      stat_summary(fun.y=mean, 
                                   geom="point", 
                                   shape=18, 
                                   size=4, 
                                   color="black", 
                                   fill="black") +
                      geom_signif(data = sig.df,
                                  aes(xmin = Start,
                                      xmax = End,
                                      annotations = p.short,
                                      y_position = y.adj),
                                  textsize = 4, vjust = -0.1,
                                  manual = TRUE, margin_top = 1) +
                      geom_point(data = sig.df,
                                 aes(x = End, y = y.adj*1.1), x=NA) +
                      facet_grid(~ Treatment) + 
                      theme_bw() +
                      scale_linetype_manual("TTE Shift", 
                                            values = c(Increase = 1, 
                                                       Decrease = 6)) + 
                      scale_color_manual(values = color.sch.ls$Prebiotic) + 
                      xlab("Test Day") +
                      ylab("TTE (sec)")
                      theme(axis.title.x = element_text()) + 
                      guides(colour = guide_legend(order = 2), 
                             linetype = guide_legend(order = 1))

# Homogenuaty
#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------
ggsave(filename = "out/plot/primary_out.png", plot = perf.dif, 
       width = 5, height = 4)

prim.res <- list(Model = list("Full Model" = r.mod.full, 
                               "Model Summary" = mod.sum, 
                               "Paired Comparisons" = paired.mod.sum,
                               "Formula" = mod.rand.full),
                 Plots = list(Main = perf.dif, 
                                Homogeneity = homo.var, 
                                QQplot = res.qq))

save(list = "prim.res", 
     file = "out/supp/1_primary.RData")

#-------------------------------------------------------------------------------
# Load data and functions
#-------------------------------------------------------------------------------
load("out/supp/prm.R")
load("out/supp/0_data.RData")

#-------------------------------------------------------------------------------
# Parameters 
#-------------------------------------------------------------------------------
seed <- prm.ls[["general"]][["seed"]]

lmm.form <- prm.ls[["primary"]][["lmm_form"]]

emmeans.form <- prm.ls[["primary"]][["emmeans_form"]]

samp.to.remove <- prm.ls[["primary"]][["samp_to_remove"]]

resp.var <- prm.ls[["primary"]][["resp_var"]]

gr.var <- prm.ls[["primary"]][["gr_var"]]

out.dir <- prm.ls[["primary"]][["out_dir_path"]]


set.seed(seed)

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# LMM for data
#-------------------------------------------------------------------------------
prime.data <- meta.data %>% 
                  filter(!PersonID %in% samp.to.remove)


# Fit model 
r.mod.full <- lmerTest::lmer(as.formula(lmm.form), 
                             data = prime.data)

mod.sum <- summary(r.mod.full)
  
# Pairwise comparisons
paired.mod.sum <- emmeans::emmeans(r.mod.full, 
                                   as.formula(emmeans.form), 
                                   type = "response")


as.data.frame(mod.sum$coefficients)

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
                dplyr::select(all_of(c(resp.var, gr.var))) %>% 
                slice(which.max(.data[[resp.var]]), .by = gr.var) %>% 
                mutate(y.adj = (.data[[resp.var]]*1.1))

sig.df <- paired.mod.sum$contrasts %>% 
                as.data.frame() %>% 
                mutate(p.short = paste0("p=", round(p.value, 3)), 
                       Start = "TD1", 
                       End = "TD3") %>% 
                left_join(., max.div, by = gr.var)


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
                      ylab("TTE (sec)") +
                      theme(axis.title.x = element_text(), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank()) + 
                      guides(colour = guide_legend(order = 2), 
                             linetype = guide_legend(order = 1))

#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------
ggsave(filename = paste0(out.dir, "/primary_out.png"), 
       plot = perf.dif, 
       width = 5, height = 4)

mod.sum$coefficients %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){round(x, 4)})) %>% 
  write.csv(., file = paste0(out.dir, "/lmm_coeficients.csv"), 
            row.names = FALSE)

paired.mod.sum$contrasts %>% 
  as.data.frame() %>% 
  write.csv(., file = paste0(out.dir, "/emmens_contrasts.csv"), 
            row.names = FALSE)

prim.res <- list(Model = list("Full Model" = r.mod.full, 
                               "Model Summary" = mod.sum, 
                               "Paired Comparisons" = paired.mod.sum,
                               "Formula" = lmm.form),
                 Plots = list(Main = perf.dif, 
                                Homogeneity = homo.var, 
                                QQplot = res.qq))

save(list = "prim.res", 
     file = "out/supp/1_primary.RData")

rm(list = ls())
# Set the main theme for ggplot2
main_theme = theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.5),
        strip.text = element_text(colour = 'black', size = 6),
        strip.background = element_rect(colour = 'black', fill = 'grey'),
        axis.title = element_text(color = 'black',size = 6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6),
        legend.position = "none")
# PCoA plot
PCoA_plot_fun <- function(dist) {
  ord.fun <-  cmdscale(dist,  k = 2, eig = T, add = T)
  pcoa.plot <- data.frame(Group = metadata$Group, scores(ord.fun)) %>%
    mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed'))) %>%
    ggplot(aes(x = Dim1, y = Dim2)) + 
    geom_point(size = 1, alpha = 0.8, shape = 21, colour = "black", aes(fill = Group)) + 
    stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1, 
                 show.legend = FALSE, level = 0.95) +
    scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    scale_color_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
    # scale_x_continuous(expand = c(0.03, 0.03)) +
    # scale_y_continuous(expand = c(0.03, 0.03)) +
    labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), digits = 3), "%)", sep = "")) +
    main_theme +
    theme(legend.background = element_blank(),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key = element_blank(),
          legend.position = c(0.85, 0.85),
          legend.key.size = unit(0.4, 'cm'))
  return(pcoa.plot)
}

# Assuming vars is defined somewhere earlier in your code
similar_deter_fun <- function(dist) {
  vars <- c('G1_C', 'G1_T', 'G2_C', 'G2_T', 'G3_C', 'G3_T', 'G4_C', 
            'G4_T', 'G5_C', 'G5_T', 'G6_C', 'G6_T')
  similar_data <- lapply(vars, function(x) 
    usedist::dist_subset(dist, 
                         grep(x, metadata$Sample_name, value = TRUE))) %>%
    do.call(cbind, .) %>%
    data.frame() %>%
    gather("tem_group", "distance") %>%
    cbind(Gully_id = rep(c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH'), each = 20),
          Group = rep(c('Un-collapsed', 'Collapsed'), each = 10, times = 6)) %>%
    select(-tem_group) %>%
    mutate(Gully_id = factor(Gully_id, levels = c('EB', 'ML', 'RS', 'SLH', 'HSX', 'HH')),
           Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))
}

## Linear mixed models test the effect of permafrost thawing on microbial diversity
lmm_fun <- function(vars, df) {
  library(lme4)
  library(lmerTest)
  lmm_dist_modes <- lapply(vars, function(x) {
    lmer(substitute(i ~ Group + Time + Slope + MAP + (1 | Gully_id), list(i = as.name(x))), 
         data = df)})
  summary.model <- function(model){
    F.value <- anova(model)$'F value'
    p.value <- anova(model)$'Pr(>F)'
    p.stars <- function(p.values) {
      unclass(symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "")))}
    sig <- p.stars(p.value)
    results<-data.frame(F.value, p.value, sig)
    return(results)
  }
  df <- NULL
  for(i in 1:length(vars)) {
    tmp <- summary.model(lmm_dist_modes[[i]])
    if (is.null(df)){
      df <- tmp
    } else {
      df <- rbind(df, tmp)
    }
  }
  result_lmm <-data.frame(dist_index = rep(vars, each = 4), 
                          variables = rep(c("Group", "Time", "Slope", "MAP"), 
                                          length(vars)), df)
  return(result_lmm)
}


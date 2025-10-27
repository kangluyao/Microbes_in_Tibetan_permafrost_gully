#' Plot scatter plot with correlation annotation
#'
#' @param data Data frame or tibble
#' @param x_var Column name for x-axis (string)
#' @param y_var Column name for y-axis (string)
#' @param method Correlation method ("pearson", "spearman", or "kendall")
#' @return ggplot object
plot_correlation <- function(data, x_var, y_var, method = "pearson") {
  if (!all(c(x_var, y_var) %in% names(data))) {
    stop("x_var or y_var not found in data.")
  }
  
  # Create scatter plot with regression line
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.6, size = 3, color = "#0072B2") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", se = TRUE) +
    theme_minimal() +
    labs(
      title = paste0(y_var, " vs. ", x_var, " (", method, " Correlation)"),
      x = x_var,
      y = y_var
    ) +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Add p-value annotation
  p <- p + stat_cor(
    method = method,
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5,
    aes(label = ..p.label..)
  )
  
  # Add R² from linear model
  model <- lm(as.formula(paste(y_var, "~", x_var)), data = data)
  r_squared <- summary(model)$r.squared
  
  x_range <- range(data[[x_var]])
  y_range <- range(data[[y_var]])
  
  p <- p + annotate(
    "text",
    x = x_range[1] + 0.1 * diff(x_range),
    y = y_range[2] - 0.15 * diff(y_range),
    label = paste0("R^2 = ", format(r_squared, digits = 3)),
    hjust = 0,
    size = 5,
    color = "darkgreen"
  )
  
  return(p)
}
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
  ord.fun <-  cmdscale(tax_singlem_dist,  k = 2, eig = T, add = T)
  # Permanova test the difference in compositional variance
  perm_test <- adonis2(tax_singlem_dist ~ Group, data = metadata)
  r2 = perm_test$R2[1]
  p_value = perm_test$`Pr(>F)`[1]
  # Calculate Centroids
  scores <- data.frame(Group = metadata$Group, scores(ord.fun))
  centroids <- aggregate(cbind(Dim1, Dim2) ~ Group, 
                         data = scores, FUN = mean)
  # Combine sample coordinates (scores) with their group information
  plot_data <- merge(scores, centroids, by = "Group", 
                     suffixes = c("", ".centroid")) %>%
    mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')))
  # Plot
  pcoa.plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2, 
                                     color = Group, shape = Group)) +
    # stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1, color = NA, 
    #              type = "t", level = 0.95, linetype = 1, size = 1) +
    # stat_ellipse(aes(colour = Group), alpha = 0.2, size = 1, 
    #              show.legend = FALSE, level = 0.95) +
    geom_point(size = 1.5) +  
    geom_point(data = centroids, 
               aes(x = Dim1, y = Dim2, color = Group), 
               size = 2, shape = 18, show.legend = F) + 
    geom_segment(aes(xend = Dim1.centroid, yend = Dim2.centroid), 
                 linetype = 1, size = 0.25) +
    scale_color_manual(name ="Group", 
                       values = c("Un-collapsed" = "#79ceb8", 
                                  "Collapsed" = "#e95f5c"),
                       breaks = c("Un-collapsed", "Collapsed")) + 
    scale_shape_manual(values = c("Un-collapsed" = 16, "Collapsed" = 17)) +
    labs(x = paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), 
                                     digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), 
                                     digits = 3), "%)", sep = "")) +
    annotate("text", x = Inf, y = Inf, 
             label = paste("PERMANOVA R2 =", round(r2, 2), "p =", round(p_value, 3)), 
             hjust = 1.05, vjust = 1.1, size = 2) +
    main_theme
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


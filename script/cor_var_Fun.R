## rename the colnames
colnames(tp_alpha_div) <- c('group', "Chao1", 'Shannon', 'Simpson', 'Site', 'Sitegroup1', 'latitude', 'longitude')
tp_alpha_div[1:5, 1:5]

## test the relationships between alpha diversity, net propoties and environmental variables
y <- cbind(tp_total_alpha_div, tp_net_indexes)
env.df <- data.frame(sample_data(tp_physeq))
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a300", "FI", "FrI", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
x <- env.df[, vars]

#linner mixed model for matrix to matrix
lmm.mat.cal <- function(y, x){
  y <- as.matrix(y)
  x <- as.matrix(x)
  df<-NULL
  for(i in colnames(y)){
    for(j in colnames(x)){
      a <- y[, i, drop = F]
      b <- x[, j, drop = F]
      mode <- lmer(a ~ b + (1|Sitegroup1), data = env.df, na.action=na.omit)
      coeff <- summary(mode)$coefficients[2,1]
      r.square <- MuMIn::r.squaredGLMM(mode)[1]
      if (coeff>0) r = sqrt(r.square)
      else r = (-1) * sqrt(r.square)
      tmp <- c(i, j, r, r.square, anova(mode)$Pr)
      if(is.null(df)){
        df <- tmp  
      }
      else{
        df <- rbind(df, tmp)
      }    
    }
  }
  df<-data.frame(row.names=NULL,df)
  colnames(df)<-c("Diversity","Env","Correlation","r.square", "Pvalue")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$AdjPvalue<-rep(0,dim(df)[1])
  df$Correlation<-round(as.numeric(as.character(df$Correlation)),3)
  #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
  # 1 -> donot adjust
  # 2 -> adjust Env + Type (column on the correlation plot)
  # 3 -> adjust Diversity + Type (row on the correlation plot for each type)
  # 4 -> adjust Diversity (row on the correlation plot)
  # 5 -> adjust Env (panel on the correlation plot)
  adjustment_label<-c("NoAdj","AdjEnvAndType","AdjDiversityAndType","AdjDiversity","AdjEnv")
  adjustment<-5
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$Diversity)){
      for(j in unique(df$Type)){
        sel<-df$Diversity==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Diversity)){
      sel<-df$Diversity==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  df$Diversity <-factor(df$Diversity, ordered = T, levels = rev(colnames(y)))
  df$Env <-factor(df$Env, ordered = T, levels = colnames(x))
  return(df)
}

lmm.matrix <- lmm.mat.cal(y, x)
# write.csv(lmm.matrix, file = './tibet_dada2_asv/results/tables/env_div_net_propotie.csv')




### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a300", "BIX", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results

#model selection
library(MASS)
library(glmulti)
A1 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a300 + BIX + HIX +
                TN + Conductivity + Salinity + Mg + K + Na, data=env_div_agg,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top <- weightable(A1)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_1 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_1

modes_inf <- NULL
for(i in 1:nrow(top_1)){
  rse_sum <- summary(A1@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(car::vif(A1@objects[[i]]) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
  } else {
    modes_inf <- rbind(modes_inf,tmp)
  } 
}
modes_inf <- cbind(top_1, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)


### heatmap using standardized regression coefficients to explore the relationship between the LCBD and environment factors
results_plot_data <- data.frame(group = c(rep('Total community', nrow(results)), rep('Carbon cycling', nrow(results_carbon))),
                                rbind(results, results_carbon))

results_plot_data$vars <- factor(results_plot_data$vars,levels = rev(vars))
results_plot_data$group <- factor(results_plot_data$group,
                                  levels=c('Total community', 'Carbon cycling'))
p_env_div <- ggplot(aes(x=LCBD, y=vars, fill=sd.coeff), data=results_plot_data) +
  geom_tile() +
  scale_fill_gradient2(low='#1b9e77', mid='white', high='#d95f02') +
  geom_text(aes(label=sig), color="black", size=6) +
  labs(y=NULL, x=NULL, fill='Standardized regression coefficients') +
  facet_wrap( .~ group, ncol = 2) +
  theme_bw()+
  theme(legend.position="bottom", 
        panel.border = element_blank(),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
###  plot the linear regression relationships between the LCBD and best explained variables
#### total community
p_linear <- env_div_agg %>%
  dplyr::select(LCBD, MAP, HIX) %>%
  tidyr::gather(varibales, value, MAP:HIX, factor_key=TRUE) %>%
  ggplot(aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c('#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')
#### carbon cycling community
plot_carbon_data <- env_div_agg_carbon %>%
  dplyr::select(LCBD, MAP, SUVA254) %>%
  tidyr::gather(varibales, value, MAP:SUVA254, factor_key=TRUE)
p_carbon_linear <- ggplot(plot_carbon_data, aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c( '#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')
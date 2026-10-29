# Multivariables analysis
## variation partitioning with two-way PERMANOVA
# Regression approach with distance as response and raw predictors
# reference:https://media.nature.com/original/nature-assets/ismej/journal/v12/n2/extref/ismej2017183x2.txt
# Set work directory
setwd('e:/thermokarst_gully/')
wd_16s <- file.path(getwd(),"data/16S/rdp")
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")

# Load required packages
library(vegan)
library(tidyverse)
library(ggplot2)

# Data input
source("script/0.read_data.R")

# Write a function for distance-based regression
DistResponseReg <-function (DistMatrix,X) {
  DistMatrix <- as.matrix(DistMatrix)
  X <- as.matrix(X)
  n <- nrow(DistMatrix)
  p <- ncol(X)
  row.wt = rep(1, nrow(DistMatrix))
  col.wt = rep(1, ncol(DistMatrix))
  st <- sum(col.wt)
  sr <- sum(row.wt)
  row.wt <- row.wt/sr
  col.wt <- col.wt/st
  DistMatrix <- -0.5*(DistMatrix*DistMatrix)
  row.mean <- apply(row.wt * DistMatrix, 2, sum)
  col.mean <- apply(col.wt *t(DistMatrix), 2, sum)
  col.mean <- col.mean - sum(row.mean * col.wt)
  DistMatrix <- sweep(DistMatrix, 2, row.mean)
  G <- t(sweep(t(DistMatrix), 2, col.mean))
  H<-X%*%solve(t(X)%*%X)%*%t(X)
  I<-diag(n)
  predicted <- H%*%G%*%H
  residuals <- (I-H)%*%G%*%(I-H)
  MS_regression<-sum(diag(predicted))/p
  MS_residual<-sum(diag(residuals))/(n-p)
  F<-MS_regression/MS_residual
  MS_Total=sum(diag(G))/n;
  RsqAdj=1-MS_residual/MS_Total;
  
  result <- list(Rsq_Adj=RsqAdj,F_value=F,res_matrix=residuals,pred_matrix=predicted)
  return(result)
  
}

# variation partitioning with distance matrix as response
# and X (e.g., environment) and W (e.g., space or connectivity)
VarPartDistResponse<-function (Dist_Matrix,X,W,Number_Permutations) {
  X <- as.matrix(X)
  X <- apply(X, 2, scale)
  W <- as.matrix(W)
  W <- apply(W, 2, scale)
  Dist_Matrix <- as.matrix(Dist_Matrix)
  Number_Predictors_X <- ncol(X)
  Number_Predictors_W <- ncol(W)
  n<-nrow(X)
  
  XW <- as.matrix(cbind(X,W));
  result <- DistResponseReg(Dist_Matrix,XW)
  abc <- result$Rsq_Adj
  Fabc <- result$F_value
  
  result <- DistResponseReg(Dist_Matrix,X)
  ab <- result$Rsq_Adj
  Fab <- result$F_value
  residuals_X <- result$res_matrix 
  predicted_X <- result$pred_matrix
  
  result <- DistResponseReg(Dist_Matrix,W)
  bc <- result$Rsq_Adj
  Fbc <- result$F_value
  residuals_W <- result$res_matrix 
  predicted_W <- result$pred_matrix
  
  # unique fraction of contribution related to X
  a <- abc - bc
  # unique fraction of contribution related to W
  c <- abc- ab
  # common fraction of contribution between X and W
  b <- abc - a - c
  # residual fraction
  d <- 1-abc
  
  Fa=(a/Number_Predictors_X)/(d/(n-Number_Predictors_X-Number_Predictors_W));
  Fc=(c/Number_Predictors_W)/(d/(n-Number_Predictors_X-Number_Predictors_W));
  
  Prob_abc=1/Number_Permutations; Prob_ab=1/Number_Permutations; Prob_bc=1/Number_Permutations; Prob_a=1/Number_Permutations; Prob_c=1/Number_Permutations;
  
  # permutations test
  for (i in 1:(Number_Permutations-1)) {
    # testing fraction a; notice that we permute the residual values in W and not in X
    permuted_rows=sample(n,replace=FALSE)
    # permuting the residual matrix, which is from the distance, and hence the need to permute
    # rows and columns in the same way, hence the use of permuted_rows for columns and rows below
    # testing fraction a
    # Yperm=predicted_W+residuals_W[permuted_rows,permuted_rows] # implement permutation of residuals in the future
    
    # testing fraction a
    result <- DistResponseReg(Dist_Matrix,XW[permuted_rows,])
    abcRnd <- result$Rsq_Adj
    FabcRnd <- result$F_value
    result <- DistResponseReg(Dist_Matrix,W[permuted_rows,])
    bcRnd <- result$Rsq_Adj
    FbcRnd <- result$F_value
    aRnd=abcRnd-bcRnd;
    dRnd=1-abcRnd;
    FaRnd=(aRnd/Number_Predictors_X)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
    if (FaRnd >= Fa) {Prob_a<-Prob_a+1/Number_Permutations}
    
    # testing fraction c
    result <- DistResponseReg(Dist_Matrix,X[permuted_rows,])
    abRnd <- result$Rsq_Adj
    FabRnd <- result$F_value
    cRnd=abcRnd-abRnd;
    FcRnd=(cRnd/Number_Predictors_W)/(dRnd/(n-Number_Predictors_X-Number_Predictors_W));
    if (FcRnd >= Fc) {Prob_c<-Prob_c+1/Number_Permutations}
    
    # testing abc
    if (FabcRnd >= Fabc) {Prob_abc<-Prob_abc+1/Number_Permutations}
    # testing ab
    if (FabRnd >= Fab) {Prob_ab<-Prob_ab+1/Number_Permutations}
    # testing bc
    if (FbcRnd >= Fbc) {Prob_bc<-Prob_bc+1/Number_Permutations}
    
  }
  
  result <- mat.or.vec(7,2)
  
  result[1,1] <- abc
  result[2,1] <- ab
  result[3,1] <- bc
  result[4,1] <- a
  result[5,1] <- c
  result[6,1] <- b
  result[7,1] <- d
  
  result[1,2] <- Prob_abc
  result[2,2] <- Prob_ab
  result[3,2] <- Prob_bc
  result[4,2] <- Prob_a
  result[5,2] <- Prob_c
  result[6,2] <- NA
  result[7,2] <- NA
  
  colnames(result) <- c("Estimate","p-value")
  rownames(result) <- c("abc","ab","bc","a","c","b","d")
  result<-as.data.frame(result)
  return(result)
}

# model selection for linear model with distance as response
#reference:
SelectionDistResponseReg<-function (Dist_Matrix,X,Number_Permutations,alpha) {
  X <- as.matrix(X)
  Number_Predictors<-ncol(X)
  n<-nrow(X)
  Dist_Matrix <- as.matrix(Dist_Matrix)
  X <- apply(X, 2, scale)
  result<-DistResponseReg(Dist_Matrix,X);
  F_Observed<-result$F_value
  # global test with all predictors
  Prob_Global<-1/Number_Permutations
  for (i in 1:(Number_Permutations-1)) {
    X_permuted<-X[sample(n,replace=FALSE),]
    result<-DistResponseReg(Dist_Matrix,X_permuted)
    F_Random<-result$F_value
    if (F_Random >= F_Observed) {Prob_Global<-Prob_Global+1/Number_Permutations}
  }
  
  # set of these as NA in the global model is not significant
  Variables_In_Model <- NA
  Rsq_Final_Model <- NA
  if (Prob_Global < alpha) {
    # calculate contribution of each predictor separately as they're orthgonal
    F_Ind_X<-mat.or.vec(Number_Predictors,1)
    for (i in 1:Number_Predictors) {
      F_Ind_X[i]<-DistResponseReg(Dist_Matrix,X[,i])[2]
    }
    # start selection
    Variables_In_Model<-as.matrix(which.max(F_Ind_X))
    F_Ind_X[Variables_In_Model[1]]<-NA
    #Original_Columns=as.matrix(1:n)
    #Original_Columns[Variable_In_Model]<-NA
    found <- FALSE
    while (found == FALSE) { 
      # this could be made faster by not calculating again the F for the first variable entered, but this version is more general
      # contrast current model with the largest variable contribution not entered in the model
      candidate_model <- cbind(X[,Variables_In_Model],X[,which.max(F_Ind_X)])
      result_candidate <- DistResponseReg(Dist_Matrix,candidate_model)
      F_candidate_Obs <- result_candidate[2]
      # test wheter the entered variable improves fit
      Prob_F <- 1/Number_Permutations
      for (i in 1:(Number_Permutations-1)) {
        candidate_predictor <- which.max(F_Ind_X)
        candidate_model <- cbind(X[,Variables_In_Model],X[sample(n,replace=FALSE),candidate_predictor])
        result_candidate <- DistResponseReg(Dist_Matrix,candidate_model);
        F_candidate_Rnd <- result_candidate$F_value
        if (F_candidate_Rnd >= F_candidate_Obs) {Prob_F<-Prob_F+1/Number_Permutations}
      }
      if (Prob_F > alpha) {found <- TRUE} else {
        F_Ind_X[candidate_predictor] <- NA
        Variables_In_Model <- append(Variables_In_Model,candidate_predictor)
      }
      
      # nrow(na.omit(OriginalColumns))
    }
    result<-DistResponseReg(Dist_Matrix,X[,Variables_In_Model])
    Rsq_Final_Model<-result$Rsq_Adj
  }
  
  result <- list(Global_P=Prob_Global,Rsq_Selected_Model=Rsq_Final_Model,Selected_Variables=Variables_In_Model)
  return(result)
}

# Compute the distance matrix of microbial community structure
siteXspe <- data.frame(t(mags_abun_tab))
env_vars <- metadata[, -c(1:4)]
trait_df <- cwm_results
mag_tab <- decostand(siteXspe, 'hellinger')
commu.dist<-vegdist(mag_tab, 'bray', uppe = F)

# check collinearities among all variables with a threshold of VIF > 10
env_df <- as.matrix(env_vars)
env_df <- decostand(env_df, "standardize")

trait_df <- as.matrix(trait_df)
trait_df <- decostand(trait_df, "standardize")


all.vars <- data.frame(env_df, trait_df)
ord <- dbrda(mag_tab ~., all.vars)
ord.vif <- vif.cca(ord)
ord.vif
# Remove MAT (vif: 110.4) and rerun the model
ord <- dbrda(mag_tab ~., all.vars[,-which(colnames(all.vars) == "MAT")])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing simple.compound.degradation (vif: 47.8) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation"))])
ord.vif <- vif.cca(ord)
ord.vif
#Continue removing N.compound.transport (vif: 42.9) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport"))])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing Genome_length (vif: 30.5) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "genome_length"))])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing C1.compounds.assimilation (vif: 23.3) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "Genome_length",
                                    "Resource.Acquisition.Substrate.assimilation.C1.compounds"))])
ord.vif <- vif.cca(ord)
# Continue removing envelope.stress (vif: 17.0) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "genome_length",
                                    "Resource.Acquisition.Substrate.assimilation.C1.compounds",
                                    "Stress.Tolerance.Specific.envelope.stress"))])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing desiccation.osmotic.salt.stress (vif: 15.5) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "genome_length",
                                    "Resource.Acquisition.Substrate.assimilation.C1.compounds",
                                    "Stress.Tolerance.Specific.envelope.stress",
                                    "Stress.Tolerance.Specific.desiccation.osmotic.salt.stress"))])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing organophosphorus.transport (vif: 13.4) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "genome_length",
                                    "Resource.Acquisition.Substrate.assimilation.C1.compounds",
                                    "Stress.Tolerance.Specific.envelope.stress",
                                    "Stress.Tolerance.Specific.desiccation.osmotic.salt.stress",
                                    "Resource.Acquisition.Substrate.uptake.organophosphorus.transport"))])
ord.vif <- vif.cca(ord)
ord.vif
# Continue removing NO3_N (vif: 11.1) and rerun the model
ord <- dbrda(mag_tab ~., 
             all.vars[ , -which(colnames(all.vars) %in% 
                                  c("MAT", 
                                    "Resource.Acquisition.Substrate.degradation.simple.compound.degradation",
                                    "Resource.Acquisition.Substrate.uptake.N.compound.transport",
                                    "genome_length",
                                    "Resource.Acquisition.Substrate.assimilation.C1.compounds",
                                    "Stress.Tolerance.Specific.envelope.stress",
                                    "Stress.Tolerance.Specific.desiccation.osmotic.salt.stress",
                                    "Resource.Acquisition.Substrate.uptake.organophosphorus.transport",
                                    "NO3_N"))])
ord.vif <- vif.cca(ord)
ord.vif
# Construct the environmental and trait martrix
sel_trais_final <- c("Max growth rate",
                     "Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization",
                     "Resource Acquisition:Substrate uptake:carbohydrate transport",
                     "Resource Acquisition:Substrate assimilation:N compounds",
                     "Resource Acquisition:Substrate assimilation:P compounds",
                     "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration",
                     "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration",
                     "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation",
                     "Resource Use:Chemotrophy:chemolithoautotrophy",
                     "Resource Use:Phototrophy:photosystem",
                     "Stress Tolerance:Specific:temperature:low",
                     "Stress Tolerance:Specific:pH stress",
                     "Stress Tolerance:Specific:oxidative stress")

                                     
env_mat <- env_df[, c(intersect(colnames(env_df), names(ord.vif)))]
trait_mat <- trait_df[, sel_trais_final]
#variables selection
env_sel_vars <- env_mat[ , SelectionDistResponseReg(commu.dist, env_mat, 999, 0.05)$Selected_Variables]
trait_sel_vars <- trait_mat[ , SelectionDistResponseReg(commu.dist, trait_mat, 999, 0.05)$Selected_Variables]
#conduct VPA
vp <- VarPartDistResponse(commu.dist, env_sel_vars, trait_sel_vars, 999)

#plot
#first, variable partition plot
theme_set(theme_grey())
env_effect <- vp[4, 1] * 100
trait_effect <- vp[5, 1] * 100
shared <- vp[6, 1] * 100

# Load required libraries
library(ggplot2)
library(ggforce)

# Round for cleaner labels
env_lab   <- paste0("Environment\n", round(env_effect, 1), "%")
trait_lab <- paste0("Traits\n", round(trait_effect, 1), "%")
shared_lab <- paste0(round(shared, 1), "%")

# Data frame defining the two circles
# Adjust x-distance between centers to control overlap size visually
circles <- data.frame(
  x0 = c(-0.7, 0.7),
  y0 = c(0, 0),
  r  = c(1.5, 1.5),
  group = c("Environment", "Traits")
)
# Set color palette for the circles
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 2
cols = gg_color_hue(4)
# Plot the Venn diagram with ggplot2 and ggforce
p <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = group),
    alpha = 0.5, color = NA, linewidth = 0.8
  ) +
  coord_fixed() +
  scale_fill_manual(values = c("Environment" = cols[1], "Traits" = cols[2])) +
  annotate("text", x = -1.35, y = 0, label = env_lab, size = 5, fontface = "bold") +
  annotate("text", x = 1.35, y = 0, label = trait_lab, size = 5, fontface = "bold") +
  annotate("text", x = 0, y = 0, label = shared_lab, size = 5, fontface = "bold") +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Variance Partitioning: Environment vs. Traits")

print(p)

# Save the plot
ggsave(file.path("E:/thermokarst_gully/revision/result/Venn_plot.pdf"),
       p, width = 6, height = 5, dpi = 600)


# ========================================
# RLQ分析完整流程
# ========================================

# 1. 加载必要的包
library(vegan)      # 用于Hellinger转换
library(ade4)       # 用于RLQ分析
library(factoextra) # 用于聚类可视化（可选）

# 2. 读取数据（根据你的实际文件格式调整）
# Rdata: 环境变量表（行=样本，列=环境变量，类别型）
# Ldata: 物种丰度表（行=样本，列=物种）
# Qdata: 性状表（行=物种，列=性状，连续型）

# 示例读取方式（根据实际情况修改）
Rdata <- metadata[, -c(1:3)]
Ldata <- t(mags_abun_tab)
Qdata <- trait_data

# ========================================
# 3. 数据检查和对齐
# ========================================

# 检查数据维度
dim(Rdata)
dim(Ldata)
dim(Qdata)

# 确保三个表的样本和物种名称一致
# L表（物种丰度）的行名应与R表（环境）的行名一致
# L表（物种丰度）的列名应与Q表（性状）的行名一致

# 对齐样本顺序
common_samples <- intersect(rownames(Rdata), rownames(Ldata))
Rdata <- Rdata[common_samples, ]
Ldata <- Ldata[common_samples, ]

# 对齐物种顺序
common_species <- intersect(colnames(Ldata), rownames(Qdata))
Ldata <- Ldata[, common_species]
Qdata <- Qdata[common_species, ]

length(common_samples)
length(common_species)

# ========================================
# 4. Hellinger转换（标准化物种丰度数据）
# ========================================

# 应用Hellinger转换到物种丰度表
mags_abun_hell <- decostand(Ldata, method = "hellinger")

# ========================================
# 5. 准备RLQ分析的三个表
# ========================================

# R表：环境变量（类别型）
R_table <- Rdata

# L表：物种丰度（已经Hellinger转换）
L_table <- mags_abun_hell

# Q表：物种性状（连续型）
Q_table <- Qdata

# ========================================
# 6. 对三个表进行初步排序分析
# ========================================

# 重要：为了确保权重一致，先进行L表（物种丰度）的对应分析
# 然后用L表的权重来约束R表和Q表的分析

# L表：对物种丰度进行对应分析（这是核心）
dudiL <- dudi.coa(L_table, scannf = FALSE, nf = 3)

# 提取L表的行权重（样本权重）和列权重（物种权重）
row_weights <- dudiL$lw  # 样本权重
col_weights <- dudiL$cw  # 物种权重

range(row_weights)
range(col_weights)

# R表：对环境变量进行分析，使用L表的行权重
# 使用row.w参数指定行权重
dudiR <- dudi.hillsmith(R_table, row.w = row_weights, scannf = FALSE, nf = 3)
# 如果metadata是类别型，使用dudi.hillsmith
# 如果metadata全是连续型，使用：
# dudiR <- dudi.pca(R_table, row.w = row_weights, scannf = FALSE, nf = 3)

# Q表：对性状进行PCA，使用L表的列权重
dudiQ <- dudi.pca(Q_table, row.w = col_weights, scannf = FALSE, nf = 3, 
                  scale = TRUE, center = TRUE)
# row.w: 指定行权重（这里是物种权重）
# scale = TRUE: 标准化变量
# center = TRUE: 中心化

# ========================================
# 7. 执行RLQ分析
# ========================================
rlq_result <- rlq(dudiR, dudiL, dudiQ, scannf = FALSE, nf = 2)
# nf = 2: 保留前2个RLQ轴

# 查看RLQ分析结果摘要
summary(rlq_result)

# 查看特征值（解释的变异）
print(rlq_result$eig)

# 查看累积解释百分比
print(cumsum(rlq_result$eig) / sum(rlq_result$eig) * 100)

# ========================================
# 8. 可视化RLQ结果
# ========================================

# 基础图形
plot(rlq_result)

# 或者更详细的可视化
par(mfrow = c(2, 2))
s.arrow(rlq_result$l1, main = "环境变量")
s.label(rlq_result$lQ, main = "物种性状")
s.label(rlq_result$mR, main = "样本得分")
s.label(rlq_result$mQ, main = "物种得分")
par(mfrow = c(1, 1))

# ========================================
# 9. 提取物种在RLQ轴上的坐标用于聚类
# ========================================

# 提取物种在前两个RLQ轴上的坐标
species_scores <- rlq_result$lQ[, 1:2]

print(head(species_scores))

# ========================================
# 10. 基于RLQ轴进行层次聚类
# ========================================

# 计算欧氏距离
dist_matrix <- dist(species_scores, method = "euclidean")

# Ward层次聚类
ward_clust <- hclust(dist_matrix, method = "ward.D2")

# 使用Calinski-Harabasz准则确定最优聚类数
# 计算不同聚类数的CH指数
library(fpc)
CH_index <- numeric(10)
for (k in 2:10) {
  clusters <- cutree(ward_clust, k = k)
  CH_index[k] <- calinhara(species_scores, clusters)
}

# 找出最优聚类数
optimal_k <- which.max(CH_index)
optimal_k

# 使用最优聚类数进行分组
final_clusters <- cutree(ward_clust, k = optimal_k)

# 查看每个聚类的物种数
cat("\n各聚类的物种数:\n")
print(table(final_clusters))

# ========================================
# 11. 可视化聚类结果
# ========================================

# 树状图
plot(ward_clust, main = "物种聚类树状图", xlab = "", sub = "")
rect.hclust(ward_clust, k = optimal_k, border = 2:6)

# RLQ空间中的聚类散点图
plot(species_scores[, 1], species_scores[, 2],
     col = final_clusters, pch = 19,
     xlab = "RLQ axis 1", ylab = "RLQ axis 2",
     main = "物种在RLQ空间中的聚类")
text(species_scores[, 1], species_scores[, 2],
     labels = rownames(species_scores), cex = 0.6, pos = 3)
legend("topright", legend = paste("Cluster", 1:optimal_k),
       col = 1:optimal_k, pch = 19)

# ========================================
# 12. 计算相关比（Correlation Ratios）
# ========================================

# 相关比表示性状与聚类分组的关联强度
correlation_ratios <- function(traits, clusters) {
  n <- length(clusters)
  k <- length(unique(clusters))
  ratios <- numeric(ncol(traits))
  names(ratios) <- colnames(traits)
  
  for (i in 1:ncol(traits)) {
    # 总平方和
    grand_mean <- mean(traits[, i], na.rm = TRUE)
    SS_total <- sum((traits[, i] - grand_mean)^2, na.rm = TRUE)
    
    # 组间平方和
    SS_between <- 0
    for (j in unique(clusters)) {
      group_data <- traits[clusters == j, i]
      group_mean <- mean(group_data, na.rm = TRUE)
      SS_between <- SS_between + length(group_data) * (group_mean - grand_mean)^2
    }
    
    # 相关比 = 组间方差 / 总方差
    ratios[i] <- SS_between / SS_total
  }
  
  return(ratios)
}

corr_ratios <- correlation_ratios(Q_table, final_clusters)
print(sort(corr_ratios, decreasing = TRUE))

# ========================================
# 13. 空间自相关检验
# ========================================

# 如果有空间坐标信息，可以进行空间自相关检验
# 这里假设Rdata中有经纬度信息（lon, lat）

if (all(c("lon", "lat") %in% colnames(Rdata))) {
  # 创建空间权重矩阵（基于距离）
  coords <- Rdata[, c("lon", "lat")]
  
  # 计算地理距离
  geo_dist <- dist(coords)
  
  # 创建邻近矩阵（例如：最近的5个邻居）
  # 这里使用简单的距离反比权重
  nb <- dnearneigh(as.matrix(coords), d1 = 0, d2 = max(geo_dist)/10)
  lw <- nb2listw(nb, style = "W")
  
  # 对前两个RLQ轴的残差进行空间自相关检验
  # 残差 = 各轴得分 - 总得分
  residuals_axis1 <- rlq_result$mR[, 1] - rowSums(rlq_result$mR)
  residuals_axis2 <- rlq_result$mR[, 2] - rowSums(rlq_result$mR)
  
  # Moran's I 检验
  cat("\n=== 空间自相关检验 ===\n")
  cat("RLQ轴1残差:\n")
  print(moran.test(residuals_axis1, lw))
  
  cat("\nRLQ轴2残差:\n")
  print(moran.test(residuals_axis2, lw))
} else {
  cat("\n未找到空间坐标信息（lon, lat），跳过空间自相关检验\n")
}

# ========================================
# 14. 保存结果
# ========================================

# 保存聚类结果
results_df <- data.frame(
  Species = rownames(species_scores),
  Cluster = final_clusters,
  RLQ_axis1 = species_scores[, 1],
  RLQ_axis2 = species_scores[, 2]
)

# write.csv(results_df, "rlq_clustering_results.csv", row.names = FALSE)

# 保存相关比
# write.csv(data.frame(Trait = names(corr_ratios), 
#                      Correlation_Ratio = corr_ratios),
#           "correlation_ratios.csv", row.names = FALSE)

cat("\n=== 分析完成 ===\n")
cat("聚类结果已保存到 results_df\n")
cat("相关比已保存到 corr_ratios\n")
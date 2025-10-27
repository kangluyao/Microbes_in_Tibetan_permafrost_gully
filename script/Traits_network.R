# Reading data
library(microtrait)
microtrait_results <- readRDS(file.path(wd_fun, "/MAGs/microtraits/thermokarst_gully.microtraitresults.rds"))

# Normalizing the traits matrices by genome length
microtrait_results_metadata_norm = microtrait_results %>% trait.normalize(normby = "genome_length")

spec_names <- data.frame(microtrait_results_metadata_norm$trait_matrixatgranularity3)$id
trait_data <- data.frame(microtrait_results_metadata_norm$trait_matrixatgranularity3)[,-1]
rownames(trait_data) <- spec_names

# Calculate CWM for each trait and each community
trait_data <- trait_data[rownames(mags_abun_tab), ]
relative_abundance <- mags_abun_tab/100
cwm_results <- data.frame(matrix(NA, nrow = ncol(relative_abundance), ncol = ncol(trait_data)))
row.names(cwm_results) <- colnames(relative_abundance)
colnames(cwm_results) <- colnames(trait_data)

for (i in 1:ncol(trait_data)) {
  for (j in 1:ncol(relative_abundance)) {
    cwm_results[j, i] <- sum(relative_abundance[, j] * trait_data[, i], na.rm = TRUE)/sum(relative_abundance[, j])
  }
}

print(cwm_results)

## Test the overall difference in the community weighted mean traits.
library(vegan)
#determine the dissimilarity matrix based on the bray-curties distance
cwm_results_final <- cwm_results[ ,colSums(cwm_results[])>0]
cwm_results_final <- cwm_results_final[metadata$Sample_name, ]

dim(cwm_results_final)
# cwm_results_trans <- scale(cwm_results_final, center = TRUE, scale = T)
cwm_results_trans <- decostand(cwm_results_final, method = "log")

#############基于丰度相关性的微生物共发生网络
net.cal <- function(x){
  library(WGCNA)
  library(impute)
  library(preprocessCore)
  library(igraph)
  
  x <- as.data.frame(x)
  cp <- corAndPvalue(x,
                     use = "pairwise.complete.obs",
                     method = "spearman",
                     alternative = c("two.sided")
  )
  
  cp.r <- cp$cor
  
  # 处理NA值
  cp.r[is.na(cp.r)] <- 0
  cp$p[is.na(cp$p)] <- 1
  
  # FDR校正
  cp.p <- matrix(p.adjust(cp$p, method = "fdr"), 
                 nrow = nrow(cp$p), 
                 ncol = ncol(cp$p))
  rownames(cp.p) <- rownames(cp$p)
  colnames(cp.p) <- colnames(cp$p)
  
  # 应用阈值
  cp.r[cp.p > 0.05 | abs(cp.r) < 0.7] <- 0
  
  # 强制对称化（处理数值精度问题）
  cp.r <- as.matrix(cp.r)  # 确保是矩阵格式
  cp.r[lower.tri(cp.r)] <- t(cp.r)[lower.tri(cp.r)]  # 用上三角覆盖下三角
  
  # 再次检查并处理NA
  cp.r[is.na(cp.r)] <- 0
  
  # 设置对角线为0
  diag(cp.r) <- 0
  
  # 构建网络
  igraph <- graph_from_adjacency_matrix(cp.r, 
                                        mode = "undirected", 
                                        weighted = TRUE, 
                                        diag = FALSE)
  
  igraph <- simplify(igraph)
  bad.vs <- V(igraph)[degree(igraph) == 0] 
  igraph <- delete_vertices(igraph, bad.vs)
  
  E(igraph)$correlation <- E(igraph)$weight
  E(igraph)$weight <- abs(E(igraph)$weight)
  
  results <- list(cp.r, igraph)
  return(results)
}

# Construct traits igraph project
CK_cwm_results_trans <- cwm_results_trans[grep("_C", rownames(cwm_results_trans)), ]
COL_cwm_results_trans <- cwm_results_trans[grep("_T", rownames(cwm_results_trans)), ]
CK_trait_net <- net.cal(CK_cwm_results_trans)
COL_trait_net <- net.cal(COL_cwm_results_trans)
# Extract the igraph object
CK_trait_igraph <- CK_trait_net[[2]]
CK_trait_igraph 
plot(CK_trait_igraph, layout = layout.circle)

COL_trait_igraph <- COL_trait_net[[2]]
COL_trait_igraph 
plot(COL_trait_igraph, layout = layout.circle)

#igraph提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(CK_trait_igraph, 'E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/CK_trait_igraph.graphml', format = 'graphml')
write_graph(COL_trait_igraph, 'E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/COL_trait_igraph.graphml', format = 'graphml')

##节点特征
#节点数量
g <- COL_trait_igraph
vcount(g)

#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(g)$degree <- degree(g)
V(g)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree_distribution(g)[-1]
degree_num <- 1:max(V(g)$degree)

par(mfrow = c(1, 2))
hist(V(g)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- knn(g, V(g))$knn
plot(V(g)$degree, neighbor_degree, log = 'xy', 
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

#加权度（Weighted degree）
V(g)$weight_degree <- strength(g)
V(g)$weight_degree

#接近中心性（Closeness centrality）
V(g)$closeness_centrality <- closeness(g)
V(g)$closeness_centrality

#介数中心性（Betweenness centrality）
V(g)$betweenness_centrality <- betweenness(g)
V(g)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(g)$eigenvector_centrality <- evcent(g)$vector
V(g)$eigenvector_centrality

#探索三种描述节点中心性的特征的关系
#探索三种描述节点中心性的特征的关系
library(car)

scatter3d(V(g)$closeness_centrality, V(g)$betweenness_centrality, V(g)$eigenvector_centrality, 
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality', 
          surface = FALSE)

#探索节点度和节点中心性的关系，如与特征向量中心性的关系
plot(V(g)$degree, V(g)$eigenvector_centrality, 
     xlab = 'Degree', ylab = 'Eigenvector centrality')

#输出列表
node_list <- data.frame(
  node_id = V(g)$name, 
  degree = V(g)$degree, 
  weight_degree = V(g)$weight_degree, 
  closeness_centrality = V(g)$closeness_centrality, 
  betweenness_centrality = V(g)$betweenness_centrality, 
  eigenvector_centrality = V(g)$eigenvector_centrality)

head(node_list)
write.table(node_list, 'E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/node_list_collapsed.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##边特征
#边的数量
ecount(g)

#权重（Weighted），已在数据读入时转化获得
E(g)$weight

#边介数中心性（Edge betweenness centrality）
E(g)$betweenness_centrality <- edge.betweenness(g)
E(g)$betweenness_centrality

#输出列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  # correlation = E(g)$corr, 
  betweenness_centrality = E(g)$betweenness_centrality
)
head(edge_list)

write.table(edge_list, 'E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/edge_list_collapsed.txt', sep = '\t', row.names = FALSE, quote = FALSE)


# Test the difference in network attributes between un-collapsed and collapsed 
# soil microbial traits
ck_node_attr <- fread("E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/G_control_node.csv")
col_node_attr <- fread("E:/thermokarst_gully/data/metagenome/MAGs/microtraits/network/G_collapsed_node.csv")

node_attr <- data.frame(Group = c(rep("Un-collapsed", nrow(ck_node_attr)),
                                  rep("Collapsed", nrow(col_node_attr))),
                        rbind(ck_node_attr, col_node_attr)) %>%
  select(name, type, Group, betweenness_centrality, 
         closeness_centrality, degree, eigenvector_centrality)

# Prepare the data for plot
node_attr_plot_data <- rbind(node_attr %>% 
                               filter(!type == "Others") %>%
                               select(Group, degree, eigenvector_centrality, 
                                      betweenness_centrality, 
                                      closeness_centrality) %>%
                               pivot_longer(-Group, names_to = "node_attr", 
                                            values_to = "value") %>%
                               mutate(type = rep("Total", nrow(.))) %>%
                               select(Group, type, everything()),
                             node_attr %>% 
                               filter(!type == "Others") %>%
                               select(Group, type, degree, 
                                      eigenvector_centrality, 
                                      betweenness_centrality, 
                                      closeness_centrality) %>%
                               pivot_longer(-c(Group, type), 
                                            names_to = "node_attr", 
                                            values_to = "value"))

# Comparison using wilcox test
facet_labeller <- as_labeller(c(
  `betweenness_centrality` = "Betweenness Centrality",
  `closeness_centrality` = "Closeness Centrality",
  `degree` = "Degree",
  `eigenvector_centrality` = "Eigenvector Centrality",
  `Total` = "Total",
  `Resource Acquisition` = "Resource Acquisition",
  `Resource Use` = "Resource Use",
  `Stress Tolerance` = "Stress Tolerance"
))

library(ggpubr)
group_comparisons <- list(c('Un-collapsed', 'Collapsed'))
net_attrs_comparison <- node_attr_plot_data  %>% 
  mutate(Group = factor(Group, levels = c('Un-collapsed', 'Collapsed')),
         type = factor(type, levels = c("Total", 'Resource Acquisition', 
                                        'Resource Use', 'Stress Tolerance')),
         node_attr = factor(node_attr, 
                            levels = c("betweenness_centrality", 
                                       "closeness_centrality", "degree", 
                                       "eigenvector_centrality"))) %>%
  ggplot(aes(Group, value, fill = Group)) +
  geom_half_violin(position = position_nudge(x = 0.25), 
                   side = "r", width = 0.6, color = NA, alpha = 0.75) +
  # geom_boxplot(width = 0.4, size = 0.75, outlier.color = NA) +
  geom_jitter(aes(fill = Group), size = 1.5, 
              width = 0.15, stroke = 0, pch = 21, alpha = 0.75) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.25) +
  stat_compare_means(comparisons = group_comparisons, paired = F,
                     p.adjust.method = "BH", label = "p.signif", 
                     bracket.size = 0.5, bracket.width = 0.1,  size = 2.5,
                     tip.length = 0.00, method = "wilcox.test") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("#79ceb8", "#e95f5c", "#5cc3e8", "#ffdb00")) +
  facet_grid(node_attr ~ type, scales = "free", 
             labeller = facet_labeller) + #
  main_theme +
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")))

net_attrs_comparison

# Save plot
ggsave(file.path("E:/thermokarst_gully/revision/result/net_attrs_comparison.pdf"),
       net_attrs_comparison, width = 183, height = 120, units = "mm")

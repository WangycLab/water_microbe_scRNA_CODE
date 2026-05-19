library(plot1cell)
library(circlize)
library(dplyr)

############################################################
# 1. read data
############################################################
setwd("D:/ran/噬菌体-单菌结合项目-20250525/阶段性进展/海水污水处理/20250829/")
# seu <- readRDS("filtered_seurat_rmrRNA_0.1_10.rds")
seu <- readRDS("filtered_seurat_0.1_dim12res0.2-fromlinux-20260326.rds")
color_list <- c(
  "#9DD0C7","#A0CBE8","#F9A825","#499894","#B15928","#ED7A90",
  "#EEC79F","#B07AA1","#6DA3E5","#EDC948","#9C755F",
  "#797BB7","#9BC985","#6BB7CA","#4E79A7","#AFB42B"
  
)
Cluster_levels <- levels(seu$seurat_clusters)
Cluster_colors <- colorRampPalette(color_list)(length(Cluster_levels))
names(Cluster_colors) <- Cluster_levels 

#############################################################
# 去除全是rRNA的cluster
#############################################################
Idents(seu) <- seu$seurat_clusters
seu <- subset(seu, idents = "3", invert = TRUE)
old_levels <- levels(seu$seurat_clusters)
new_levels <- as.character(0:(length(old_levels) - 1))

seu$seurat_clusters <- factor(
  seu$seurat_clusters,
  levels = old_levels,
  labels = new_levels
)

Idents(seu) <- seu$seurat_clusters
table(seu$seurat_clusters)
seu$celltype <- seu$seurat_clusters
table(seu$Species)
dim(seu)

#############################################################
# 0、9、15 to 0
#############################################################
# 1. cluster
seu$seurat_clusters_old <- as.character(Idents(seu))

# 2. cluster label
cur_cluster <- as.character(Idents(seu))

# 3.  Cluster 0、9、15 to cluster 0
cur_cluster[cur_cluster %in% c("0", "9", "15")] <- "0"

# 4. 
old_levels <- sort(unique(as.numeric(cur_cluster)))
new_levels <- as.character(0:(length(old_levels) - 1))
map_vec <- setNames(new_levels, as.character(old_levels))

new_cluster <- unname(map_vec[cur_cluster])
names(new_cluster) <- colnames(seu)

# 5. 
seu <- AddMetaData(seu, metadata = new_cluster, col.name = "seurat_clusters")
seu$seurat_clusters <- factor(seu$seurat_clusters, levels = new_levels)

# 6. 
Idents(seu) <- "seurat_clusters"

# 7. results
table(seu$seurat_clusters)

# plot
DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.1
) +
  scale_color_manual(values = color_list) +  # 使用上面生成的颜色
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )





circ_data <- prepare_circlize_data(seu, scale = 0.8)

############################################################
# 2 Cluster 修复
############################################################
circ_data <- prepare_circlize_data(seu, scale = 0.65)
cluster_levels <- as.character(0:15)  # 13是根据Cluster的数量进行更改的

circ_data$Cluster <- factor(
  circ_data$seurat_clusters,
  levels = cluster_levels
)

############################################################
# 3 颜色定义
############################################################

# cluster_colors <- c(
#   "#9DD0C7","#D1392B","#A0CBE8","#499894","#B15928","#ED7A90",
#   "#EEC79F","#B07AA1","#6DA3E5","#EDC948","#9C755F",
#   "#797BB7","#B6992D","#9BC985","#6BB7CA","#4E79A7","#8A89A6",
#   "#818181"
# )
cluster_colors <- c(
  "#9DD0C7","#A0CBE8","#F9A825","#499894","#B15928","#ED7A90",
  "#EEC79F","#B07AA1","#6DA3E5","#EDC948","#9C755F",
  "#797BB7","#9BC985","#6BB7CA","#4E79A7","#AFB42B"
  
)

cluster_colors <- setNames(cluster_colors, cluster_levels)

group_colors <- c(
  SP="#A2BAD9",
  AP="#E2B1A9",
  EP="#ACD1B1"
)

species_levels <- names(sort(table(circ_data$Species), decreasing=TRUE))

color_list_species <- c(
  "#1B9E77","#D1F3FC","#FFD92F","#7570B3","#B2DF8A","#CD5C5C",
  "#F2CCA0","#CAB2D6","#A6761D","#FC8D62","#6A9171",
  "#FF9DA7","#8CD17D","#C0F4CA","#FEFCBF","#FABFD2","#79706E",
  "#D9D9D9"
)

species_colors <- setNames(
  color_list_species[seq_along(species_levels)],
  species_levels
)

############################################################
# 4 排序
############################################################

cluster_levels <- levels(circ_data$Cluster)

cluster_colors <- setNames(
  cluster_colors[seq_along(cluster_levels)],
  cluster_levels
)

circ_data <- circ_data %>%
  arrange(Cluster, x_polar)

############################################################
# 5 主 circlize 图
############################################################

circos.clear()

par(mar=c(1,1,1,7))

plot_circlize(
  circ_data,
  do.label=TRUE,
  pt.size=0.05,
  kde2d.n = 150,
  bg.color="white",
  col.use=cluster_colors,
  repel=TRUE,
  label.cex=0.8
)




############################################################
# 6 Cluster × Group
############################################################

group_df <- circ_data %>%
  group_by(Cluster, Group) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(Cluster) %>%
  mutate(
    start=cumsum(lag(n, default=0))/sum(n),
    end=cumsum(n)/sum(n)
  )

############################################################
# spacer
############################################################

circos.track(
  ylim=c(0,1),
  track.height=0.02,
  bg.border=NA
)

############################################################
# 7 Group 
############################################################

circos.track(
  ylim=c(0,1),
  track.height=0.015,
  bg.border=NA,
  panel.fun=function(x,y){
    
    sector <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    
    df <- group_df[group_df$Cluster==sector,]
    
    if(nrow(df)==0) return()
    
    par(lend="round")
    
    circos.segments(
      df$start * xlim[2],
      0.5,
      df$end * xlim[2],
      0.5,
      col = group_colors[df$Group],
      lwd = 3.5
    )

  }
)

############################################################
# spacer
############################################################

circos.track(
  ylim=c(0,1),
  track.height=0.02,
  bg.border=NA
)

############################################################
# 8 Cluster × Species
############################################################

species_df <- circ_data %>%
  group_by(Cluster, Species) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(Cluster) %>%
  mutate(
    start=cumsum(lag(n, default=0))/sum(n),
    end=cumsum(n)/sum(n)
  )

############################################################
# 9 Species 外圈（圆头线段）
############################################################

circos.track(
  ylim=c(0,1),
  track.height=0.015,
  bg.border=NA,
  panel.fun=function(x,y){
    
    sector <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    
    df <- species_df[species_df$Cluster==sector,]
    
    if(nrow(df)==0) return()
    
    par(lend="round")
    
    circos.segments(
      df$start * xlim[2],
      0.5,
      df$end * xlim[2],
      0.5,
      col = species_colors[df$Species],
      lwd = 3.5
    )
    
  }
)

############################################################
# 10 Gap Label
############################################################

############################################################
# 11 Legend
############################################################

par(xpd=TRUE)

legend(
  x=1.45,
  y=1.2,
  legend=names(cluster_colors),
  fill=cluster_colors,
  title="Cluster",
  bty="n",
  cex=0.7
)

legend(
  x=2,
  y=1.2,
  legend=names(group_colors),
  fill=group_colors,
  title="Group",
  bty="n",
  cex=0.7
)

legend(
  x=1.45,
  y=0.18,
  legend=names(species_colors),
  fill=species_colors,
  title="Species",
  bty="n",
  cex=0.7
)



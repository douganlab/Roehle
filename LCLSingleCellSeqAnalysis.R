# LCLSingleCellSeqAnalysis.R
# Copyright (C) 2020 Lestat Ali
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

source("LCLPlotHelpers.R")

library(Seurat)

#### Data Preparation ####

# Import 10x data for vehicle-treated cells
veh_raw_data <- Read10X(data.dir = "data/10x/veh/raw_feature_bc_matrix")
veh_raw_seurat <- CreateSeuratObject(counts = veh_raw_data, project = "veh", min.cells = 3, min.features = 100)
veh_raw_seurat$percent_mt <- PercentageFeatureSet(veh_raw_seurat, pattern = "^mt-")

# Import 10x data for LCL-treated cells
lcl_raw_data <- Read10X(data.dir = "data/10x/lcl/raw_feature_bc_matrix")
lcl_raw_seurat <- CreateSeuratObject(counts = lcl_raw_data, project = "lcl", min.cells = 3, min.features = 100)
lcl_raw_seurat$percent_mt <- PercentageFeatureSet(lcl_raw_seurat, pattern = "^mt-")

#### Merge & Cluster ####

vl_merged_raw_seurat <- merge(veh_raw_seurat, lcl_raw_seurat, add.cell.ids = c("veh", "lcl"))
count_min <- 500
count_max <-  25000 # Less than 3 SDs of mean 
mt_cutoff <- mean(vl_merged_raw_seurat$percent_mt) + 2 * sd(vl_merged_raw_seurat$percent_mt)
cell_filter <- (vl_merged_raw_seurat$nCount_RNA < count_max 
                & vl_merged_raw_seurat$nCount_RNA > count_min
                & vl_merged_raw_seurat$percent_mt < mt_cutoff)
vl_merged_seurat <- vl_merged_raw_seurat[, cell_filter]
vl_merged_seurat <- NormalizeData(vl_merged_seurat)
vl_merged_seurat <- FindVariableFeatures(vl_merged_seurat, nfeatures = 2500)
vl_merged_seurat <- ScaleData(vl_merged_seurat, features = rownames(vl_merged_seurat))
vl_merged_seurat <- RunPCA(vl_merged_seurat, features = VariableFeatures(object = vl_merged_seurat))
vl_merged_seurat <- RunTSNE(vl_merged_seurat, dims = 1:19, perplexity = 60)

# First clustering round: unsupervised Louvain
vl_merged_seurat <- FindNeighbors(vl_merged_seurat, dims = 1:19, k.param = 20)
vl_merged_seurat <- FindClusters(vl_merged_seurat, resolution = 0.6)

# Sub-cluster B cells
vl_cluster6_seurat <- subset(vl_merged_seurat, seurat_clusters == "6")
vl_cluster6_seurat <- RunTSNE(vl_cluster6_seurat, dims = 1:19)
vl_cluster6_seurat <- FindNeighbors(vl_cluster6_seurat, dim = 1:2, reduction = "tsne", k.param = 20)
vl_cluster6_seurat <- FindClusters(vl_cluster6_seurat, resolution = 0.5)
# Subcluster 0 is enriched for macrophage markers (e.g. Lyz2 and Mafb)

# Sub-cluster T cells
vl_cluster8_seurat <- subset(vl_merged_seurat, seurat_clusters == "8")
vl_cluster8_seurat <- RunTSNE(vl_cluster8_seurat, dims = 1:19)
vl_cluster8_seurat <- FindNeighbors(vl_cluster8_seurat, dim = 1:2, reduction = "tsne", k.param = 20)
vl_cluster8_seurat <- FindClusters(vl_cluster8_seurat, resolution = 0.5)
draw.cluster.tsne(vl_cluster8_seurat) 
# Subclusters 1 & 3 are enriched for macrophage markers (e.g. Lyz2 and Mafb)

# Collect actively phagocytosing marophages into their own cluster
phagocytosed_t_cells <- WhichCells(vl_cluster6_seurat, idents = "0")
phagocytosed_b_cells <- WhichCells(vl_cluster8_seurat, idents = c("1", "3"))
active_phagocytes <- c(phagocytosed_t_cells, phagocytosed_b_cells)
vl_merged_seurat <- SetIdent(vl_merged_seurat, cells = active_phagocytes, "15")
levels(vl_merged_seurat) <- 0:15 # Just for ordering sanity

#### Plot tSNE with final clustering (Fig 5A) ####

draw.cluster.tsne(vl_merged_seurat)

#### Cluster Marker Gene Violin Plots (Fig 5B) #### 

all_cluster_ids <- levels(vl_merged_seurat)
for (cluster_id in all_cluster_ids) {
  markers <- FindMarkers(vl_merged_seurat, ident.1 = cluster_id, test.use = "MAST")
}

interesting_markers <- c("Csf3r", "Siglecf", "C1qa", "Mafb", "Mki67", "Cd3d",
                         "Ccr7", "Cd79a", "Cd209a", "Xcr1", "Siglech", "Cpa3")

draw.violin.sandwich(vl_merged_seurat, interesting_markers)

#### Macrophager Marker Gene Heatmap (Fig 5C) #### 

heatmap_markers <- read.csv("data/marker_genes_for_heatmap.csv", header = F, stringsAsFactors = F)$V1
heatmap_clusters <- as.character(c(1, 2, 5, 9, 11, 15, 6, 8))
all_avg_exp <- AverageExpression(vl_merged_seurat)$RNA
marker_avg_exp <- all_avg_exp[heatmap_markers, heatmap_clusters]
marker_scaled_exp <- t(scale(t(marker_avg_exp)))
marker_scaled_exp_df <- reshape2::melt(marker_scaled_exp)
colnames(marker_scaled_exp_df) <- c("Gene", "ClusterID", "Expression")
marker_scaled_exp_df$ClusterID <- factor(marker_scaled_exp_df$ClusterID, levels = heatmap_clusters)
marker_scaled_exp_df$Gene <- factor(marker_scaled_exp_df$Gene, levels = rev(heatmap_markers))

heatmap <- ggplot(data = marker_scaled_exp_df, aes(y = Gene, x = ClusterID, fill = Expression)) +
  geom_tile(colour = "white") + 
  scale_fill_gradient2(guide = "none", high = "red", mid = "white", low = "blue", midpoint = 0, name = "") +
  scale_x_discrete(expand = expand_scale(add = 0), position = "top",
                   labels = c("Macro-\nMafb", "Macro-\nMaf", "Macro-\nVegfa", "Macro-\nMgl2",
                              "Macro-\nCycling", "Macro-\nPhagocytic", "T Cell", "B Cell")) + 
  scale_y_discrete(expand = expand_scale(add = 0)) +
  theme_classic() + xlab("") + ylab("") + 
  theme(axis.line.x = element_line(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Arial", size = 9, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, family = "Arial", size = 9))
heatmap

#### RNA Reads By Macrophage Cluster (Fig 5D) #### 

all_macro_cells <- WhichCells(vl_merged_seurat, idents = c(1, 2, 5, 9, 11, 15))
macro_rna_df <- data.frame("Cell" = all_macro_cells)
macro_rna_df$ClusterID <- Idents(vl_merged_seurat)[all_macro_cells]
macro_rna_df$RNACount <- vl_merged_seurat$nCount_RNA[all_macro_cells]
macro_cluster_labels <- c("Mafb", "Maf", "Vegfa", "Mgl2", "cycling", "phagocytic")
rna_count_plot <- ggplot(data = macro_rna_df, aes(x = ClusterID, y = RNACount)) +
  geom_violin(alpha = 0.3, aes(fill = ClusterID), colour = "white", scale = "width") +
  geom_jitter(width = 0.25, aes(colour = ClusterID), size = 0.75) +
  scale_colour_manual(values = lcl_cluster_colors, guide = "none") +
  scale_fill_manual(values = lcl_cluster_colors, guide = "none") +
  scale_x_discrete(labels = macro_cluster_labels) +
  scale_y_continuous(labels = seq(0, 25, 5)) +
  ylab("RNA Read Count (x1000)") + xlab("Macrophage Cluster") +
  theme.lcl2() + theme(axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(size = 10),
                       axis.text.y = element_text(size = 10),
                       axis.title.y = element_text(size = 12)) 
rna_count_plot

#### Cluster Size As Percent of Each Sample (Fig 5E) #### 

sample_composition_df <- data.frame("ClusterID" = rep(all_cluster_ids, 2),
                                    "SampleID" = c(rep("veh", 16), rep("lcl", 16)),
                                    "Percent" = 0)

cluster_visualization_order <- c("0", "3", "4", "1", "2", "5", "9", "11", "15",
                                 "6", "8", "7", "12", "10", "13", "14")
sample_composition_df$ClusterID <- factor(sample_composition_df$ClusterID, levels = cluster_visualization_order)

sample_composition_df$SampleID <- factor(sample_composition_df$SampleID, levels = c("veh", "lcl"))
for (i in 1:nrow(sample_composition_df)) {
  current_cluster <- sample_composition_df$ClusterID[i]
  current_sample <- sample_composition_df$SampleID[i]
  cluster_size <- length(which(vl_merged_seurat$seurat_clusters == current_cluster 
                               & vl_merged_seurat$orig.ident == current_sample))
  sample_size <- length(which(vl_merged_seurat$orig.ident == current_sample))
  sample_composition_df$Percent[i] <- 100 * cluster_size / sample_size
}

sample_composition_plot <- ggplot(data = sample_composition_df, aes(x = ClusterID, y = Percent)) +
  geom_col(aes(fill = SampleID), colour = "white", size = 0.25, position = "dodge") +
  scale_y_continuous(expand = expand_scale(0), limits = c(0, 25)) +
  scale_fill_manual(values = c("black", "#43B5B9"), name = "", labels = c("Vehicle", "LCL161")) +
  scale_x_discrete(labels = lcl_cluster_names) +
  ylab("% of Sample") +
  theme.lcl2()
sample_composition_plot

#### MHC-II+ Ly6c2+ Highlighted Cells (Fig 5F) #### 

veh_cells <- Cells(vl_merged_seurat)[startsWith(Cells(vl_merged_seurat), "veh")]
lcl_cells <- Cells(vl_merged_seurat)[startsWith(Cells(vl_merged_seurat), "lcl")]

mhc2_ly6c2_cells <- Cells(vl_merged_seurat)[((as.vector(vl_merged_seurat$RNA["H2-DMb1",]) > 1)
                                             | (as.vector(vl_merged_seurat$RNA["H2-DMb2",]) > 0.25)
                                             | (as.vector(vl_merged_seurat$RNA["H2-DMa",]) > 1))
                                            & (as.vector(vl_merged_seurat$RNA["Ly6c2",]) > 0)]

draw.highlighted.tsne(vl_merged_seurat[,veh_cells], mhc2_ly6c2_cells)
draw.highlighted.tsne(vl_merged_seurat[,lcl_cells], mhc2_ly6c2_cells)


# LCLPlotHelpers.R
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

library(ggplot2)

lcl_cluster_colors <- c("0" = "#F5CB1A",
                        "1" = "palevioletred",
                        "2" = "#F9989A", 
                        "3" = "#6CC6C3", 
                        "4" ="#CBC8BB",
                        "5" = "darkseagreen", 
                        "6" = "#CB449F",
                        "7" = "#3D5EA9", 
                        "8" = "#68397A",
                        "9" = "burlywood2",
                        "10" = "#ADACCB",
                        "11" = "#B98E92",
                        "12" = "steelblue4",
                        "13" = "#DF0000", 
                        "14" = "lightskyblue",
                        "15" = "green4")

lcl_cluster_names <- c("0" ="Granulocyte 1",
                       "1" = "Macro-Mafb",
                       "2" = "Macro-Maf",
                       "3" = "Granulocyte 2",
                       "4" = "Granulocyte 3",
                       "5" = "Macro-Vegfa",
                       "6" = "T cells",
                       "7" = "cDC1-Ccl22",
                       "8" = "B cells",
                       "9" = "Macro-Mgl2",
                       "10" = "cDC2-Cd209a",
                       "11" = "Macro-Cycling",
                       "12" = "cDC1-Clec9a",
                       "13" = "pDC",
                       "14" = "Mast cells",
                       "15" = "Macro-Phagocytic")

theme.lcl <- function() {
  return(theme_bw() + theme(text = element_text(family = "Arial"),
                            panel.grid = element_blank(),
                            panel.border = element_blank(),
                            legend.text = element_text(size = 16),
                            legend.title = element_text(size = 16),
                            axis.text.x = element_text(size = 13),
                            axis.text.y = element_text(size = 13),
                            axis.title = element_text(size = 16),
                            axis.line = element_line(size = 0.5, arrow = arrow(length = unit(0.1, "inches")))))
}

theme.lcl2 <- function() {
  return(theme.lcl() + 
           theme(axis.title.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.text = element_text(size = 13),
                 legend.position = "top",
                 axis.line = element_line(arrow = NULL)))
}

draw.cluster.tsne <- function(seurat_object, save = FALSE) {
  embeddings_df <- as.data.frame(Reductions(seurat_object, slot = "tsne")@cell.embeddings)
  embeddings_df$Cluster <- Idents(seurat_object)
  plot <- ggplot(data = embeddings_df, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(shape = 19, aes(colour = Cluster), alpha = 0.9, size = 0.8) +
    scale_colour_manual(values = lcl_cluster_colors) +
    guides(colour = guide_legend(override.aes = list(size = 7))) + 
    xlab("tSNE Dimension 1") + ylab("tSNE Dimension 2") +
    theme.lcl()
  if (save) {
    ggsave("tsne_clusters.pdf", plot, width = 12, height = 9.5, device = cairo_pdf)
  } else {
    plot
  }
}

draw.highlighted.tsne <- function(seurat_object, highlight_cells, title = "") {
  embeddings_df <- as.data.frame(Reductions(seurat_object, slot = "tsne")@cell.embeddings)
  fg_embeddings_df <- embeddings_df[rownames(embeddings_df) %in% highlight_cells,]
  bg_embeddings_df <- embeddings_df[!(rownames(embeddings_df) %in% highlight_cells),]
  plot <- ggplot(data = NULL, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(data = bg_embeddings_df, shape = 19, colour = "grey80", alpha = 0.75, size = 0.5) +
    geom_point(data = fg_embeddings_df, shape = 19, colour = "green4", alpha = 0.9, size = 0.75) +
    xlab("tSNE Dimension 1") + ylab("tSNE Dimension 2") + ggtitle(title) +
    theme.lcl() + theme(axis.title = element_text(size = 12),
                        axis.text.x = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text.y = element_blank())
  return(plot)
}

draw.violin.sandwich <- function(seurat_object, genes) {
  count_data <- GetAssayData(seurat_object$RNA)
  cluster_ids <- Idents(seurat_object)
  gene_expr <- count_data[genes,]
  plottable_expr <- reshape2::melt(t(as.matrix(gene_expr)))
  plottable_expr$Cluster <- rep(cluster_ids, length(genes))
  colnames(plottable_expr) <- c("CellID", "Gene", "Expression", "ClusterID")
  
  plottable_expr$ClusterID <- factor(plottable_expr$ClusterID, levels = viz_cluster_levels)
  expr_violinp <- ggplot(data = plottable_expr, aes(x = ClusterID, y = Expression)) +
    geom_violin(aes(fill = ClusterID, colour = ClusterID), size = 0.8, alpha = 1, scale = "width") + 
    facet_grid(rows = vars(Gene)) + 
    scale_colour_manual(values = lcl_cluster_colors) +
    scale_fill_manual(values = lcl_cluster_colors) +
    xlab("") + ylab("") +
    guides(fill = "none", colour = "none") +
    theme.lcl2() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "cm"),
          plot.title = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  return(expr_violinp)
}

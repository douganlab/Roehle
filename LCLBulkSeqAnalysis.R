# LCLBulkSeqAnalysis.R
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
source("LCLAnalysisHelpers.R")

library(DESeq2)
library(ggplot2)
library(Hmisc)

#### Load In Vivo Data ####

vivo_metadata <- read.csv("data/bulk_vivo_metadata.csv", stringsAsFactors = F)
rownames(vivo_metadata) <- vivo_metadata$ID
vivo_metadata$Treatment <- factor(vivo_metadata$Treatment, levels = c("Vehicle", "LCL161"))

vivo_raw_counts <- read.csv("data/bulk_vivo_counts.csv", row.names = 1)
vivo_raw_counts <- vivo_raw_counts[, vivo_metadata$ID]

vivo_dds <- DESeqDataSetFromMatrix(vivo_raw_counts, vivo_metadata, ~ Treatment)
vivo_analyzed_dds <- DESeq(vivo_dds)

vivo_res <- lfcShrink(vivo_analyzed_dds, "Treatment_LCL161_vs_Vehicle", type = "apeglm")

#### Load In Vitro Data #### 

vitro_metadata <- read.csv("data/bulk_vitro_metadata.csv", stringsAsFactors = F)
rownames(vitro_metadata) <- vitro_metadata$RID

vitro_raw_counts <- read.csv("data/bulk_vitro_counts.csv", row.names = 1)
vitro_raw_counts <- vitro_raw_counts[, vitro_metadata$RID]

vitro_dds <- DESeqDataSetFromMatrix(vitro_raw_counts, vitro_metadata, ~ Treatment)
vitro_dds$Treatment <- relevel(vitro_dds$Treatment, ref = "LCL")
vitro_analyzed_dds <- DESeq(vitro_dds)

vitro_res_lt <- lfcShrink(vitro_analyzed_dds, "Treatment_LCL.LT_vs_LCL", type = "apeglm")
vitro_res_lt_sig <- subset(vitro_res_lt, padj < 0.05)

#### Determine Overlap Between In Vitro & In Vivo #### 

matched_vivo_res <- vivo_res[rownames(vitro_res_lt_sig),]

vitro_vivo_combo_df <- data.frame("Gene" = rownames(vitro_res_lt_sig),
                               "VitroAvgExpression" = vitro_res_lt_sig$baseMean,
                               "VivoAvgExpression" = matched_vivo_res$baseMean,
                               "VitroLFC" = vitro_res_lt_sig$log2FoldChange,
                               "VivoLFC" = matched_vivo_res$log2FoldChange,
                               "VitroRawP" = vitro_res_lt_sig$pvalue,
                               "VivoRawP" = matched_vivo_res$pvalue,
                               "VitroAdjP" = vitro_res_lt_sig$padj,
                               "VivoAdjP" = matched_vivo_res$padj,
                               stringsAsFactors = F)

# Must have concordant LFCs
vitro_vivo_combo_df <- subset(vitro_vivo_combo_df, VitroLFC * VivoLFC > 0)
vitro_vivo_combo_df <- subset(vitro_vivo_combo_df, !is.na(VivoLFC))

# Combine P values
meta_pval_df <- apply.fisher.method(vitro_vivo_combo_df$VitroAdjP, vitro_vivo_combo_df$VivoAdjP)
vitro_vivo_combo_df$MetaP <- meta_pval_df$p.adj

# Combine expression data
vitro_vivo_combo_df$MetaAvgExpression <- rowMeans(cbind(vitro_vivo_combo_df$VitroAvgExpression, vitro_vivo_combo_df$VivoAvgExpression))
vitro_vivo_combo_df$MetaLFC <- rowMeans(cbind(vitro_vivo_combo_df$VitroLFC, vitro_vivo_combo_df$VivoLFC))

# For readability
vitro_vivo_combo_df <- vitro_vivo_combo_df[order(vitro_vivo_combo_df$MetaLFC, decreasing = T),]

#### Violin Plot of Common Genes (Fig 6F) ####

vitro_vivo_combo_sub_df <- subset(vitro_vivo_combo_df, (MetaAvgExpression > 100
                                                        & MetaP < 0.05
                                                        & VitroRawP < 0.05
                                                        & VivoRawP < 0.05))

vitro_vivo_volcano <- ggplot(data = vitro_vivo_combo_sub_df, aes(x = MetaLFC, y = -log10(MetaP))) +
  geom_point(data = subset(vitro_vivo_combo_sub_df, MetaLFC > 0), size = 2, shape = 19, alpha = 0.9, colour = "red2") +
  geom_point(data = subset(vitro_vivo_combo_sub_df, MetaLFC < 0), size = 2, shape = 19, alpha = 0.9, colour = "royalblue3") +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 2), expand = expand_scale(0)) +
  xlab("Fold Change (log)") + ylab("FDR-Adjusted P (-log)") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme.lcl() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.line = element_line(arrow = NULL))
vitro_vivo_volcano

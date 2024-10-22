# ==============================================================================
# Script Name: GeoMx TBL PCA + UMAP from DESeq2
# Description: Using normalzied counts from DESeq2 to project counts into space.
# Unlike in PTB, we need to apply batch correction to get a strong separation. 
# Author: Jonathan Perrie
# Date: 2024-10-21
# ==============================================================================

setwd("C:/Users/Jonathan/Documents/UCLA/Pelligrini/projects/geomx/")

library(DESeq2)
library(dplyr)
library(ggrepel)
library(openxlsx)
library(sva)
library(circlize)
library(ggplot2)
library(ggsignif)
library(umap)

# ==============================================================================
# Input and output parameters
# ==============================================================================

raw_count_input <- "current_data/geomx_integrated_raw_data_2024.csv"
pca_output <- "plots/TBL_PCA.png"
umap_output <- "plots/TBL_UMAP.png"
pca_bc_output <- "plots/TBL_PCA_corrected.png"
umap_bc_output <- "plots/TBL_UMAP_corrected.png"
legend_output <- "plots/TBL_PCA_UMAP_legend.png"

# ==============================================================================
# DESeq2
# ==============================================================================

mat <- read.table(raw_count_input,sep=",",row.names=1,header=T,check.names=F)

# Define boolean flags for filtering data
tbl_bool <- grepl("\\| TBL ",colnames(mat))

tbl_mat <- mat[,tbl_bool]

treatment <- trimws(sapply(strsplit(colnames(tbl_mat),"\\|"),"[",1))
batch <- trimws(sapply(strsplit(colnames(tbl_mat),"\\|"),function(x) paste(x[2:(length(x)-1)], collapse = " ")))

meta <- data.frame(treatment=trimws(treatment),batch=trimws(batch))
rownames(meta) <- colnames(tbl_mat)
meta$treatment <- factor(meta$treatment, levels=c("Core", "Giant", "Mantle", "Infiltrate"))
meta$batch <- as.factor(meta$batch)

dds <- DESeqDataSetFromMatrix(countData = tbl_mat,
                              colData = meta,
                              design = ~ treatment+batch)
dds <- DESeq(dds) 

# ==============================================================================
# DESeq2 counts PCA/UMAP
# ==============================================================================

# Extract normalized counts and filter out specific conditions
counts <- assay(vst(dds))
meta <- meta

# Calculate gene variability and perform PCA on the most variable genes
geneVars <- rowVars(as.matrix(counts))
select <- rownames(counts)[order(geneVars, decreasing = TRUE)[1:500]]
pca <- prcomp(t(counts[select,]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
group <- meta$treatment

# Prepare data for plotting
intgroup.df <- data.frame("label" = group, "shape" = sub("^\\D+", "", meta$batch))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))

# Generate PCA plot without a legend
g1 <- ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke = 1.33) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed() +
  xlim(-pca_range, pca_range) +
  ylim(-pca_range, pca_range) +
  theme_classic() +
  labs(color = "Region", shape = "Patient") +
  scale_shape_manual(values = c(1, 0, 6, 5)) +
  scale_color_manual(values = c("#00BA38", "#4F7942", "#619CFF", "#000080")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20, family = 'sans'),
    legend.position = "none",  # Remove legend
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    legend.spacing.y = unit(2, "mm")
  )

ggsave(file = pca_output, plot = g1, device = "png", width=7.5,height=6, dpi = 400)

set.seed(0)

# Run UMAP using PCA results as initialization
umap_result <- umap(pca$x)

# Prepare data for plotting
d <- data.frame(PC1 = umap_result$layout[,1], PC2 = umap_result$layout[,2], group = group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))

# Generate UMAP plot
g2 <- ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke = 1.33) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  coord_fixed() +
  xlim(-pca_range, pca_range) +  # Set x and y axis limits
  ylim(-pca_range, pca_range) +
  theme_classic() +
  labs(color = "Region", shape = "Patient") +
  scale_shape_manual(values = c(1, 0, 6, 5)) +
  scale_color_manual(values = c("#00BA38", "#4F7942", "#619CFF", "#000080")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20, family = 'sans'),
    legend.position = "none",  # Remove legend
    legend.box = "horizontal",
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    legend.spacing.y = unit(2, "mm")
  ) 

ggsave(file = umap_output, plot = g2, device = "png", width=7.5,height=6, dpi = 400)

# ==============================================================================
# Batch corrected PCA/UMAP
# ==============================================================================

# Extract normalized counts and filter out specific conditions
counts <- assay(vst(dds))
batch <- meta$batch
group <- meta$treatment

# Batch correction applied to data, control for all four groups 
counts <- ComBat(dat = as.matrix(counts), batch = batch, mod = model.matrix(~ group))

# Calculate gene variability and perform PCA on the most variable genes
geneVars <- rowVars(as.matrix(counts))
select <- rownames(counts)[order(geneVars, decreasing = TRUE)[1:500]]
pca <- prcomp(t(counts[select,]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
group <- meta$treatment

# Prepare data for plotting
intgroup.df <- data.frame("label" = group, "shape" = sub("^\\D+", "", meta$batch))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))

# Generate PCA plot without a legend
g1 <- ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke = 1.33) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed() +
  xlim(-pca_range, pca_range) +
  ylim(-pca_range, pca_range) +
  theme_classic() +
  labs(color = "Region", shape = "Patient") +
  scale_shape_manual(values = c(1, 0, 6, 5)) +
  scale_color_manual(values = c("#00BA38", "#4F7942", "#619CFF", "#000080")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20, family = 'sans'),
    legend.position = "none",  # Remove legend
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    legend.spacing.y = unit(2, "mm")
  )

ggsave(file = pca_bc_output, plot = g1, device = "png", width=7.5,height=6, dpi = 400)

set.seed(0)

# Run UMAP using PCA results as initialization
umap_result <- umap(pca$x)

# Prepare data for plotting
d <- data.frame(PC1 = umap_result$layout[,1], PC2 = umap_result$layout[,2], group = group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))

# Generate UMAP plot
g2 <- ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke = 1.33) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  coord_fixed() +
  xlim(-pca_range, pca_range) +  # Set x and y axis limits
  ylim(-pca_range, pca_range) +
  theme_classic() +
  labs(color = "Region", shape = "Patient") +
  scale_shape_manual(values = c(1, 0, 6, 5)) +
  scale_color_manual(values = c("#00BA38", "#4F7942", "#619CFF", "#000080")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20, family = 'sans'),
    legend.position = "none",  # Remove legend
    legend.box = "horizontal",
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    legend.spacing.y = unit(2, "mm")
  ) 

ggsave(file = umap_bc_output, plot = g2, device = "png", width=7.5,height=6, dpi = 400)

# ==============================================================================
# Legend
# ==============================================================================

# Set up the ggplot object
g <- ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke = 1.33) +
  theme_classic() +
  labs(color = "Region", shape = "Patient") +
  scale_shape_manual(values = c(1, 0, 6, 5)) +
  scale_color_manual(values = c("#00BA38", "#4F7942", "#619CFF", "#000080")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20, family = 'sans'),
    legend.position = c(0.5, 0.5),  # Central position for legend
    legend.box = "horizontal",
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    legend.spacing.y = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 3))
  )

# Build the ggplot object to extract the legend
legend <- ggplot_gtable(ggplot_build(g))
legend_grob <- legend$grobs[[which(sapply(legend$grobs, function(x) x$name) == "guide-box")]]

# Turn off the default device to avoid empty file saving
dev.off()

# Draw the legend on a new page using grid system
grid::grid.newpage()
grid::grid.draw(legend_grob)

# Save the extracted legend to a PNG file
ggsave(file = legend_output, plot = grid::grid.grabExpr(grid::grid.draw(legend_grob)), device = "png", width = 6, height = 6, dpi = 400)

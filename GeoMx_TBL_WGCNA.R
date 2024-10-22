# ==============================================================================
# Script Name: GeoMx TBL WGCNA
# Description: Computing gene modules for TBL data and identifying most
# correlated genes in each module. 
# Author: Jonathan Perrie
# Date: 2024-10-21
# ==============================================================================

setwd("C:/Users/Jonathan/Documents/UCLA/Pelligrini/projects/geomx/")

library(DESeq2)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(sva)
library(ComplexHeatmap)
library(WGCNA)
library(circlize)

# ==============================================================================
# Input and output parameters
# ==============================================================================

raw_count_input <- "current_data/geomx_integrated_raw_data_2024.csv"
thrsh_output <- "plots/TBL_threshold.png"
me_hm_output <- "plots/TBL_me_hm.png"
me_hm_tx_output <- "plots/TBL_me_hm_tx.png"
me_tx_corr_output <- "current_data/TBL_me_hm_tx.xlsx"

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
# Select gene threshold for WGCNA - # genes to include for module analysis 
# ==============================================================================

counts <- assay(vst(dds))
batch <- meta$batch
group <- meta$treatment

# Batch correction applied to data, control for all four groups 
counts <- ComBat(dat = as.matrix(counts), batch = batch, mod = model.matrix(~ group))

geneVars <- rowVars(as.matrix(counts))
sorted_geneVars <- sort(geneVars, decreasing = TRUE)

# ==============================================================================
# Plot gene threshold 
# ==============================================================================

# Select the top 50 and every 100th gene after that
top_genes <- c(seq(1,25),unique(floor(seq(25, length(sorted_geneVars)^0.9, length.out = 50)^1.05)))

# Prepare a data frame for plotting
gene_labels <- names(sorted_geneVars)[top_genes]
gene_values <- sorted_geneVars[top_genes]
plot_df <- data.frame(
  Gene = gene_labels,
  Variance = gene_values,
  Rank = top_genes
)

# Add some random jitter to avoid overlap
plot_df$Variance_jitter <- plot_df$Variance + runif(nrow(plot_df), -0.05, 0.05)

# Prepare a full dataset for plotting all points
all_points_df <- data.frame(
  Rank = seq_along(sorted_geneVars),
  Variance = sorted_geneVars
)

gg_thresh <- ggplot() +
  # Plot all points with a faint line
  geom_line(data = all_points_df, aes(x = Rank, y = Variance), color = "gray80") +  # Faint gray line for all points
  # Red dots for selected points with text
  geom_point(data = plot_df, aes(x = Rank, y = Variance_jitter), color = "red", size = 2) +  # Red dots for labeled points
  # Text labels with ggrepel
  geom_text_repel(data = plot_df, aes(x = Rank, y = Variance_jitter, label = Gene), 
                  hjust = 0.25, vjust = 0.25, size = 3, alpha = 0.6, max.overlaps = Inf) +  # Use ggrepel for text labels
  # Add a vertical line at gene 3000
  geom_vline(xintercept = 3000, color = "blue", linetype = "dashed", linewidth = 1) +  # Dashed vertical line at Rank 3000
  theme_minimal() +
  labs(x = "Genes",
       y = "Variance") +
  theme(
    axis.line.x = element_line(color = "black"),  # Add bottom spine (x-axis line)
    axis.line.y = element_line(color = "black"),  # Add left spine (y-axis line)
    axis.title.x = element_text(size = 24),       # Make x-axis label ("Genes") big
    axis.title.y = element_text(size = 24),       # Make y-axis label ("Variance") big
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size=16), # Remove x-axis tick values if desired
    axis.ticks.x = element_blank(),               # Remove x-axis ticks if desired
    panel.grid = element_blank(),                  # Remove grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set white background
    plot.background = element_rect(fill = "white", color = NA)    # Set plot background to white
  )

ggsave(file = thrsh_output, plot = gg_thresh, device = "png", width=9, height=6, dpi = 400)

# ==============================================================================
# WGCNA
# ==============================================================================

# params to adjust: 
# deepSplit = 0-4 (lower yields fewer modules)
# power: 8-12 or use estimate pickSoftThreshold(x)$powerEstimate
# gene count: 3k 
# clustering: Ward's linkage 
# minClusterSize: IGH genes tend to form a small module 

select <- names(sorted_geneVars)[1:3000]
datExpr <- as.data.frame(t(counts[select,]))
adjacency <- adjacency(datExpr, power = 12, type = "signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
geneTree <- hclust(as.dist(1 - TOM), method = "ward.D2")

dynamicMods <- cutreeDynamic(
  dendro = geneTree, 
  distM = TOM, 
  deepSplit = 1, 
  pamRespectsDendro = TRUE, 
  respectSmallClusters = TRUE, 
  minClusterSize = 10
)

moduleEigengenes <- moduleEigengenes(datExpr, colors = labels2colors(dynamicMods))$eigengenes

# ==============================================================================
# WGCNA module eigengene heatmap 
# ==============================================================================

df <- as.matrix(moduleEigengenes)

df_rownames <- sapply(strsplit(rownames(moduleEigengenes), " \\| "), function(x) paste(x[2], as.integer(x[3]), sep = "-"))
rownames(df) <- df_rownames

module_cols <- gsub("ME", "", colnames(df))
colnames(df) <- gsub("ME", "", colnames(df))

me_count <- table(labels2colors(dynamicMods))
colnames(df) <- paste0(colnames(df), " (", me_count[colnames(df)], ")")

ha <- HeatmapAnnotation(
  df = data.frame(
    Patient = meta$batch,
    Region = meta$treatment
  ),
  col = list(
    Patient = c("TBL 16" = "orange", "TBL 42" = "cyan", "TBL 62" = "rosybrown2"),
    Region = c(
      "Core" = "#00BA38",
      "Giant" = "#4F7942",
      "Mantle" = "#619CFF",
      "Infiltrate" = "#000080"
    )
  ),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontface = "bold")
)

breaks <- c(min(df, 0), 0, max(df))
colors <- c("blue", "white", "red")
color_mapping <- colorRamp2(breaks, colors)

names(module_cols) <- colnames(df)
me_ha <- rowAnnotation(
  df = data.frame("ME" = colnames(df)),
  col = list("ME" = module_cols),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontface = "bold"),
  show_legend = FALSE
)

hm_me <- Heatmap(
  t(df), 
  name = "ME score", 
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  right_annotation = me_ha,
  bottom_annotation = ha,
  clustering_method_columns = "ward.D2",
  clustering_method_rows = "ward.D2",
  col = color_mapping,
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 16),
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean"
)

cell_size_px <- 500  # Size of each cell in pixels, adjust as needed
extra_height_px <- 200  # Extra height for annotations, adjust as needed
extra_width_px <- 1000 

n_rows <- nrow(df)
n_cols <- ncol(df)
dpi <- 400
png_width_in <- (n_rows * cell_size_px/4 + extra_height_px) / dpi
png_height_in <- (n_cols * cell_size_px/6 + extra_width_px) / dpi

png(me_hm_output, width = png_width_in, height = png_height_in, units = "in", res = dpi)
draw(hm_me)
dev.off()

# ==============================================================================
# WGCNA module gene correlations 
# ==============================================================================

geneModuleCorrelation <- cor(datExpr, moduleEigengenes, use = "p")
threshold <- 0.6
highlyCorrelatedGenes <- list()

for (module in colnames(moduleEigengenes)) {
  geneNames <- select[labels2colors(dynamicMods) == gsub("ME", "", module)]
  geneIndices <- which(geneModuleCorrelation[geneNames, module] > threshold)
  geneNames <- geneNames[geneIndices]
  correlations <- geneModuleCorrelation[geneNames, module]
  
  highlyCorrelatedGenes[[module]] <- data.frame(Gene = geneNames, Correlation = correlations)
}

top3CorrelatedGenes <- list()

for (module in names(highlyCorrelatedGenes)) {
  orderedGenes <- highlyCorrelatedGenes[[module]][order(highlyCorrelatedGenes[[module]]$Correlation, decreasing = TRUE), ]
  topGenes <- head(orderedGenes, n = 5)
  
  top3CorrelatedGenes[[module]] <- topGenes$Gene
}

gmc_df <- geneModuleCorrelation[unique(unname(unlist(top3CorrelatedGenes))), paste0("ME", module_cols)[row_order(hm_me)]]
colnames(gmc_df) <- paste0(gsub("ME", "", colnames(gmc_df)), " (", me_count[gsub("ME", "", colnames(gmc_df))], ")")

me_ha <- HeatmapAnnotation(
  df = data.frame("ME" = colnames(gmc_df)),
  col = list("ME" = module_cols),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontface = "bold"),
  show_legend = FALSE
)

breaks <- c(-1, 0, 1)
colors <- c("blue", "white", "red")
color_mapping <- colorRamp2(breaks, colors)

hm_tx <- Heatmap(
  gmc_df,
  name = "Pearson", 
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "euclidean",
  bottom_annotation = me_ha,
  col = color_mapping,
  column_names_gp = grid::gpar(fontsize = 16),
  column_names_rot = 90,
  row_names_gp = grid::gpar(fontsize = 8),
  show_heatmap_legend = FALSE
)

# hm_tx specs 
cell_size_px <- 500  # Size of each cell in pixels, adjust as needed
extra_height_px <- 200  # Extra height for annotations, adjust as needed
extra_width_px <- 900 

n_rows <- nrow(gmc_df)
n_cols <- ncol(gmc_df)
png_height_in <- (n_rows * cell_size_px/6 + extra_height_px) / dpi
png_width_in <- (n_cols * cell_size_px/4 + extra_width_px) / dpi
png(me_hm_tx_output, width = png_width_in, height = png_height_in, units = "in", res = dpi)
draw(hm_tx)
dev.off()

# ==============================================================================
# WGCNA module gene correlation to excel 
# ==============================================================================

threshold <- 0.3  # Threshold to only display genes with correlation > 0.3 
me_color_list <- list()  # To store gene names per module
me_corr_list <- list()   # To store gene correlations per module
maxlen <- 0              # Initialize max length for padding NA values

# Iterate over each unique module color
for (me_color in unique(labels2colors(dynamicMods))) {
  # Select genes in the current module
  curr_genes <- select[labels2colors(dynamicMods) == me_color]
  
  # Filter genes by correlation threshold
  corr_values <- geneModuleCorrelation[curr_genes, paste0("ME", me_color)]
  filtered_genes <- curr_genes[corr_values >= threshold]
  
  # Update maxlen to reflect the largest number of filtered genes in any module
  if (length(filtered_genes) > maxlen) {
    maxlen <- length(filtered_genes)
  }
  
  # Store the filtered genes for the current module
  me_color_list[[me_color]] <- filtered_genes
}

# Iterate again for sorting and padding
for (me_color in unique(labels2colors(dynamicMods))) {
  curr_genes <- me_color_list[[me_color]] 
  
  # Sort the genes by their correlation values, in decreasing order
  sorted_genes <- curr_genes[order(geneModuleCorrelation[curr_genes, paste0("ME", me_color)], decreasing = TRUE)]
  sorted_corr <- geneModuleCorrelation[curr_genes, paste0("ME", me_color)][order(geneModuleCorrelation[curr_genes, paste0("ME", me_color)], decreasing = TRUE)]
  
  # Pad with NA values if the number of genes is less than maxlen
  me_color_list[[me_color]] <- c(sorted_genes, rep(NA, maxlen - length(sorted_genes)))
  me_corr_list[[me_color]] <- c(sorted_corr, rep(NA, maxlen - length(sorted_corr)))
}

# Convert the lists to data frames
me_color_df <- as.data.frame(me_color_list)
me_corr_df <- as.data.frame(me_corr_list)

# Create an Excel workbook and add sheets
wb <- createWorkbook()

addWorksheet(wb, "Gene Names")
addWorksheet(wb, "Correlations")

# Write the data frames to the Excel sheets
writeData(wb, sheet = "Gene Names", x = me_color_df)
writeData(wb, sheet = "Correlations", x = me_corr_df)

# Save the workbook
saveWorkbook(wb, file = me_tx_corr_output, overwrite = TRUE)
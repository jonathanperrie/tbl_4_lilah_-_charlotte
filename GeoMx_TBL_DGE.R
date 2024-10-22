# ==============================================================================
# Script Name: GeoMx TBL DGE using DESeq2
# Description: Computing DGE for mantle vs core and then generating pathway 
# enrichment plot based on Enrichr hits. 
# Author: Jonathan Perrie
# Date: 2024-10-21
# ==============================================================================

setwd("C:/Users/Jonathan/Documents/UCLA/Pelligrini/projects/geomx/")

library(DESeq2)
library(dplyr)
library(openxlsx)
library(GO.db)
library(enrichR)
library(ggplot2)

# ==============================================================================
# Input and output parameters
# ==============================================================================

raw_count_input <- "current_data/geomx_integrated_raw_data_2024.csv"
dge_output <- "current_data/DGE_TBL.xlsx"
path_sig_output <- "plots/TBL_pathways.png"
path_excel_output <- "current_data/TBL_pathways.xlsx"

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
meta$treatment[meta$treatment == "Giant"] = "Core"
meta$treatment <- factor(meta$treatment, levels=c("Core", "Mantle", "Infiltrate"))
meta$batch <- as.factor(meta$batch)

dds <- DESeqDataSetFromMatrix(countData = tbl_mat,
                              colData = meta,
                              design = ~ treatment+batch)
dds <- DESeq(dds) 

# ==============================================================================
# DGE for TBL MvC
# ==============================================================================

pad_list <- function(x,maxlen) {
  c(x, rep(NA, maxlen - length(x)))
}

process_res <- function(res,thr,fc){
  res<-res[!is.na(res$padj),]
  res<-res[res$padj<=thr,]
  
  neg_res = res[res$log2FoldChange<(-fc),]
  plus_res = res[res$log2FoldChange>fc,]
  
  a<-rownames(neg_res[order(abs(neg_res$log2FoldChange),decreasing=TRUE),])
  b<-rownames(plus_res[order(abs(plus_res$log2FoldChange),decreasing=FALSE),])
  
  res<-rbind(res[a,],res[b,])
  return(res)
}

tbl_res_mc <- results(dds,contrast=c("treatment","Mantle","Core"))
tbl_cm = rownames(filter(as.data.frame(process_res(tbl_res_mc,5E-2,0)),log2FoldChange<0))
tbl_mc = rev(rownames(filter(as.data.frame(process_res(tbl_res_mc,5E-2,0)),log2FoldChange>0)))

# ==============================================================================
# Write DGE to file 
# ==============================================================================


maxlen <- max(c(length(tbl_cm),length(tbl_mc)))

genelists_df <- data.frame("TBL CvM"=pad_list(tbl_cm,maxlen), "TBL MvC"=pad_list(tbl_mc,maxlen), check.names=FALSE)

wb <- createWorkbook()

results_list <- list(
  "gene list" = genelists_df,
  "TBL MvC" = process_res(tbl_res_mc,5E-2,0)
)

lapply(names(results_list), function(sheet_name) {
  res <- results_list[[sheet_name]]
  addWorksheet(wb, sheet_name)
  if (sheet_name == "gene list") {
    writeData(wb, sheet_name, as.data.frame(res), rowNames = FALSE)
  } else {
    writeData(wb, sheet_name, as.data.frame(res), rowNames = TRUE)
  }
})

saveWorkbook(wb, dge_output, overwrite = TRUE)

# ==============================================================================
# Pathway query with Enrichr
# ==============================================================================

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr")  # Human genes
}
if (websiteLive) dbs <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2023")

geneLists <- list()
geneLists[["core>mantle"]] <- tbl_cm
geneLists[["mantle>core"]] <- tbl_mc

pathways <- list()

# Function to get the ancestors of a GO term and format as a semicolon-separated list
get_GO_ancestors <- function(go_id) {
  ancestors <- as.list(GOBPANCESTOR)[[go_id]]
  if (is.null(ancestors) || !length(ancestors)) return("")  # Return empty string if no ancestors found
  paste(ancestors, collapse = ";")
}

# Enrichment analysis for each module and database
for (module in names(geneLists)) {
  fm_dbs <- enrichr(genes = geneLists[[module]], dbs)
  
  # Store enrichment results with filtering for significant pathways and selecting top 40
  pathways[[module]] <- setNames(
    lapply(dbs, function(enrichr_db) {
      enriched_results <- fm_dbs[[enrichr_db]]
      
      # Sort by Adjusted.P.value and select top 40
      enriched_results <- enriched_results %>%
        arrange(Adjusted.P.value) %>%
        head(n = 40)
      
      # Add GO_IDs and formatted list of ancestors for GO terms, if applicable
      if (grepl("GO", enrichr_db)) {
        enriched_results$GO_ID <- gsub(".*\\((GO:\\d+)\\)", "\\1", enriched_results$Term)
        enriched_results$Ancestors <- sapply(enriched_results$GO_ID, get_GO_ancestors)
      }
      
      enriched_results$region <- module
      enriched_results
    }),
    dbs
  )
  Sys.sleep(3)  # Pause to prevent API rate limits or similar issues
}


# Post-process to check for sub-hits
for (module in names(pathways)) {
  for (db in names(pathways[[module]])) {
    results <- pathways[[module]][[db]]
    
    if ("Ancestors" %in% colnames(results)) {
      # Check if any GO ID is in another's ancestor list
      results$Is_Sub_Hit <- sapply(results$GO_ID, function(id) {
        any(sapply(results$Ancestors, function(anc_list) grepl(id, anc_list)))
      })
      pathways[[module]][[db]] <- results
    }
  }
}

shorthand_map <- list(
  "core_mantle_GO_Biological_Process_2023" = "CvMGOBP2023",
  "mantle_core_GO_Biological_Process_2023" = "MvCGOBP2023"
)

# Create an Excel workbook
wb <- createWorkbook()

# Loop over each module and database, and write the results to the workbook
for (module in names(pathways)) {
  for (db in names(pathways[[module]])) {
    # Construct the long sheet name
    sheet_name_long <- paste(module, db, sep = "_")
    sheet_name_long <- gsub("[^[:alnum:]_]", "_", sheet_name_long)  # Replace non-alphanumeric characters
    
    # Check if there's a predefined shorthand for the sheet name
    sheet_name <- ifelse(sheet_name_long %in% names(shorthand_map),
                         shorthand_map[[sheet_name_long]],
                         substr(sheet_name_long, 1, 31))  # Default to first 31 characters if no shorthand
    
    # Add a new sheet for this module and database
    addWorksheet(wb, sheet_name)
    
    # Get the results for the module and database
    results <- pathways[[module]][[db]]
    
    # Write the data to the sheet
    writeData(wb, sheet_name, results)
  }
}

# Save the workbook to a file
saveWorkbook(wb, path_excel_output, overwrite = TRUE)


# ==============================================================================
# Select specific hits from the pathways 
# ==============================================================================

# GO terms of interest
# finding deepest level GO terms: filter(pathways$`core>mantle`$GO_Biological_Process_2023, !Is_Sub_Hit)
go_terms <- c(
  "Regulation Of Inflammatory Response",
  "Lysosomal Lumen Acidification (GO:0007042)",
  "Regulation Of T Cell Proliferation (GO:0042129)",
  "Receptor Internalization (GO:0031623)",
  "Regulation Of T Cell Receptor Signaling Pathway (GO:0050862)",
  "B Cell Receptor Signaling Pathway (GO:0050853)",
  "Regulation Of B Cell Differentiation (GO:0045577)",
  "T Cell Activation (GO:0042110)"
)

# Function to escape parentheses in GO terms
escape_parentheses <- function(input_string) {
  escaped_string <- gsub("(", "\\(", input_string, fixed = TRUE)
  escaped_string <- gsub(")", "\\)", escaped_string, fixed = TRUE)
  return(escaped_string)
}

go_terms <- unname(sapply(unique(go_terms), escape_parentheses))

# Function to find pathway hits based on significant terms
find_pathways <- function(paths, terms, n = 80) {
  term_bool <- as.logical(rowSums(sapply(terms, function(x) { grepl(x, paths$Term[1:n]) })))
  return(paths[term_bool, ])
}
# Finding pathway hits for different regions and databases
path_hits <- list()
regions <- names(pathways)
databases <- c("GO_Biological_Process_2023")

for (region in regions) {
  for (db in dbs) {
    path_hits[[paste(region)]] <- find_pathways(pathways[[region]][[db]], go_terms)
  }
}

# ==============================================================================
# Prepare pathway hits for plotting 
# ==============================================================================

# Combine pathway results into a single dataframe and assign unique row names
path_hit_names <- names(path_hits)
path_df <- do.call(rbind, lapply(path_hit_names, function(name) {
  df <- path_hits[[name]]
  if (nrow(df) > 0) {
    rownames(df) <- paste(name, seq_len(nrow(df)), sep = "_")
  }
  return(df)
}))

# Function to capitalize the first letter of a string
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Function to customize the capitalization and format of pathway names
process_pathway_name <- function(name) {
  name <- gsub("\\s*\\(GO:[^)]+\\)", "", name)  # Removes GO terms
  name <- gsub("\\s*WP[^)]+", "", name)         # Removes WikiPathway terms
  name <- gsub("\\s*R-HSA[^)]+", "", name)      # Removes Reactome terms
  
  name <- tolower(name)  # Normalize text to lowercase
  name <- gsub("b cell", "B cell", name, perl = TRUE)
  name <- gsub("^t cell", "T cell", name, perl = TRUE)
  name <- gsub(" t cell", " T cell", name, perl = TRUE)
  name <- gsub(" ii", " II", name, perl = TRUE)
  name <- gsub("i cell", "I cell", name, perl = TRUE)
  name <- gsub("tcr", "TCR", name, perl = TRUE)
  name <- gsub("mhc", "MHC", name, perl = TRUE)
  
  name <- firstup(name)  # Capitalize the first letter of the processed string
  paste(strwrap(name, width = 48), collapse = "\n")  # Wrap text to 48 characters
}

formatted_pathNames <- sapply(path_df$Term, process_pathway_name)
path_df$term <- factor(formatted_pathNames)

# Calculate fraction of overlap in pathways
path_df$frac <- sapply(strsplit(path_df$Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
path_df$region <- factor(firstup(as.character(path_df$region)))

# Random order for pathways
set.seed(123)  # Set seed for reproducibility
path_df$unique_term <- paste(path_df$term, path_df$region, sep = "_")

# Custom order for pathways
path_df$random_order <- order(10^(9 - as.numeric(path_df$region) + 1) * -log10(path_df$Adjusted.P.value), decreasing = TRUE)
unique_terms_ordered <- path_df$unique_term[path_df$random_order]
terms_ordered <- unique(path_df$term[path_df$random_order])
path_df$term <- factor(path_df$term, levels = terms_ordered)

# Relabel regions for visualization
path_df$region <- factor(path_df$region, levels = levels(path_df$region), labels = c("Core > Mantle", "Mantle > Core"))

# Create a ggplot visualization of pathways
gg_path <- ggplot(path_df, aes(x = region, y = term, size = frac, color = -log10(Adjusted.P.value))) +
  geom_point() +
  scale_color_gradientn(
    colors = c("blue", "red"),
    limits = c(-log10(0.05), 10),
    breaks = c(-log10(0.05), -log10(1e-5), -log10(1e-10)),
    labels = c("0.05", "1e-5", "1e-10"),
    oob = scales::squish
  ) +
  scale_y_discrete(position = "right") +
  labs(
    size = "Gene ratio",
    color = expression(paste("Adjusted ", italic("P"), " value")),
    x = "",  # Removing x-axis label
    y = ""   # Removing y-axis label
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    text = element_text(size = 20, family = "sans"),
    axis.text.y.right = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
    legend.position = "left",
    panel.background = element_rect(fill = "white", color = NA),  # Set white background
    plot.background = element_rect(fill = "white", color = NA)    # Set plot background to white
  )

ggsave(file = path_sig_output, plot = gg_path, device = "png", width=7.5,height=6, dpi = 400)
# ==============================================================================
# SECTION 1: Setup & Helper Functions
# ==============================================================================

# Load required libraries
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(ggplot2)

# Helper function: Biased color palette (skews towards higher values)
biased_palette <- function(pal = "Reds", n = 100, power = 2) {
  base_col <- colorRampPalette(RColorBrewer::brewer.pal(9, pal))(n)
  idx <- ceiling((1:n)^power / n^power * n)
  base_col[idx]
}

# Helper function: 'Staircase' sorting within groups
# Sorts columns to creating a visually organized heatmap diagonal structure
within_group_order_stair <- function(M, annot, group_col = "motif_annotation",
                                     group_levels,
                                     row_priority = c("Mlig", "Smed", "Sman"),
                                     decreasing = c(TRUE, TRUE, TRUE),
                                     na_fill = 0) {
  # Validation
  stopifnot(all(colnames(M) %in% rownames(annot)))
  stopifnot(all(row_priority %in% rownames(M)))
  if (length(decreasing) == 1) decreasing <- rep(decreasing, length(row_priority))
  
  # Prepare annotation subset matching matrix columns
  annot <- annot[colnames(M), , drop = FALSE]
  split_cols <- split(colnames(M), annot[[group_col]])
  ord <- character(0)
  
  # Iterate through grouped levels
  for (g in group_levels[group_levels %in% names(split_cols)]) {
    cols <- split_cols[[g]]
    if (length(cols) <= 1) { ord <- c(ord, cols); next }
    
    # Extract sub-matrix and fill NAs
    Ms <- M[row_priority, cols, drop = FALSE]
    X  <- Ms
    X[is.na(X)] <- na_fill
    
    # Create sort keys based on row priority values
    keys <- lapply(seq_along(row_priority), function(i) {
      v <- as.numeric(X[i, ])
      if (decreasing[i]) -v else v
    })
    keys <- c(keys, list(cols)) # Tie-breaker
    
    # Execute sort
    o <- do.call(order, keys)
    ord <- c(ord, cols[o])
  }
  ord
}

# Helper function: Process data into Matrix and Annotation DF
# This standardizes the matrix creation logic used in both branches below
create_matrix_and_annotation <- function(data_df) {
  # Create wide matrix: Species x Cluster_Annotation
  mat <- data_df %>%
    group_by(species, cluster_ann) %>%
    summarise(value = max(max_norm, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = cluster_ann, values_from = value) %>%
    arrange(species)
  
  mat_m <- as.matrix(column_to_rownames(mat, "species"))
  
  # Filter rows to specific species order
  mat_m <- mat_m[c("Mlig", "Smed", "Sman"), , drop = FALSE]
  
  # Create column annotations
  # Lookup 'final_cluster' from original data to ensure accuracy
  lookup_fc <- data_df %>% distinct(cluster_ann, final_cluster)
  
  ann_col <- tibble(cluster_ann = colnames(mat_m)) %>%
    left_join(lookup_fc, by = "cluster_ann") %>%
    separate(
      cluster_ann,
      into   = c("final_cluster_v2", "motif_annotation"),
      sep    = " \\| ",
      remove = FALSE
    )
  
  ann_col_df <- as.data.frame(
    ann_col[, c("motif_annotation", "final_cluster_v2", "final_cluster")]
  )
  rownames(ann_col_df) <- ann_col$cluster_ann
  
  # Add non-conserved flag
  ann_col_df$nc_flag <- ifelse(
    ann_col_df$final_cluster == "non-conserved",
    "non-conserved",
    "other"
  )
  ann_col_df$nc_flag <- factor(ann_col_df$nc_flag, levels = c("non-conserved", "other"))
  
  # Set Group factor levels
  group_levels <- c("DLX", "ETS", "ZIP6", "ZNT2", "ZNT2-21") # Added ZNT2-21 for safety
  ann_col_df$motif_annotation <- factor(ann_col_df$motif_annotation, levels = group_levels)
  
  return(list(mat = mat_m, ann = ann_col_df))
}


# ==============================================================================
# SECTION 2: Initial Data Cleaning & Preprocessing
# ==============================================================================

# Read raw data
df <- read_csv("/Users/gyang/Downloads/Guang_Figure5_table.csv", col_types = cols())

# Identify rows containing 'POU' in any string column
string_cols <- df %>% select(where(is.character))
mask <- apply(string_cols, 1, function(row) any(str_detect(row, "POU")))
cat("Removed", sum(mask), "rows containing 'POU'.\n")

# Filter data
df <- df[!mask, ]

# Parse key_unique into species and motif_annotation
df <- df %>%
  mutate(
    species = str_split(key_unique, "_", simplify = TRUE)[, 1],
    motif_annotation = str_split(key_unique, "_", simplify = TRUE)[, 2]
  )

# Handle final_cluster logic: resolve 'non-conserved' to original cluster ID
df$final_cluster_v2 <- df$final_cluster
df <- df %>%
  mutate(final_cluster_v2 = if_else(
    final_cluster_v2 == "non-conserved",
    cluster,
    final_cluster_v2
  ))

# Select columns and save intermediate file
df_out <- df %>%
  select(max_norm, species, motif_annotation, final_cluster, final_cluster_v2)
write_csv(df_out, "/Users/gyang/Downloads/Guang_Figure5_table_plot.csv")


# ==============================================================================
# SECTION 3: Data Processing - Branch A (Standard ZNT2 Logic)
# ==============================================================================

# Reload data for plotting phase
df_plot <- data.frame(fread("/Users/gyang/Downloads/Guang_Figure5_table_plot.csv"))

# Filter unwanted annotations
df_plot <- df_plot[!df_plot$motif_annotation %in% c("promoter", "distal"), ]

# --- Logic Branch A ---
# Remove ZNT2-1, merge ZNT2-2 into ZNT2
df_A <- df_plot %>%
  filter(!grepl("ZNT2-1", motif_annotation)) %>%
  mutate(motif_annotation = gsub("ZNT2-2", "ZNT2", motif_annotation)) %>%
  mutate(cluster_ann = paste0(final_cluster_v2, " | ", motif_annotation))

# Generate Matrix and Annotation
res_A <- create_matrix_and_annotation(df_A)
mat_m_A    <- res_A$mat
ann_col_A  <- res_A$ann

# Perform staircase sorting
group_levels_A <- c("DLX", "ETS", "ZIP6", "ZNT2")
ord_cols_A <- within_group_order_stair(
  mat_m_A, ann_col_A,
  group_levels = group_levels_A,
  row_priority = c("Mlig", "Smed", "Sman")
)

# Ordered objects for Branch A
mat_ord_A     <- mat_m_A[, ord_cols_A, drop = FALSE]
ann_col_ord_A <- ann_col_A[ord_cols_A, , drop = FALSE]


# ==============================================================================
# SECTION 4: Data Processing - Branch B (Alternative ZNT2 Logic)
# ==============================================================================

# --- Logic Branch B ---
# Remove ZNT2-2, map ZNT2-1 -> ZNT2 -> ZNT2-21
df_B <- df_plot %>%
  filter(!grepl("ZNT2-2", motif_annotation))

# Sequential text replacement for annotation
df_B[] <- lapply(df_B, function(x) gsub("ZNT2-1", "ZNT2", x))
df_B[] <- lapply(df_B, function(x) gsub("ZNT2", "ZNT2-21", x))

df_B <- df_B %>%
  mutate(
    max_norm = as.numeric(max_norm),
    cluster_ann = paste0(final_cluster_v2, " | ", motif_annotation)
  )

# Generate Matrix and Annotation
res_B <- create_matrix_and_annotation(df_B)
mat_m_B    <- res_B$mat
ann_col_B  <- res_B$ann

# Perform staircase sorting (Note updated group level ZNT2-21)
group_levels_B <- c("DLX", "ETS", "ZIP6", "ZNT2-21")
ord_cols_B <- within_group_order_stair(
  mat_m_B, ann_col_B,
  group_levels = group_levels_B,
  row_priority = c("Mlig", "Smed", "Sman")
)

# Ordered objects for Branch B
mat_ord_B     <- mat_m_B[, ord_cols_B, drop = FALSE]
ann_col_ord_B <- ann_col_B[ord_cols_B, , drop = FALSE]


# ==============================================================================
# SECTION 5: Merging & Final Visualization
# ==============================================================================

# --- 1. Combine Results ---
# Merge unique annotations from both branches
ann_col_combined <- unique(rbind(ann_col_ord_A, ann_col_ord_B))

# Merge matrices: Branch A full + Specific columns from Branch B (ZNT2-21)
# Instead of hardcoded indices, dynamically select columns containing "ZNT2-21"
znt21_cols <- grep("ZNT2-21", colnames(mat_ord_B))

if(length(znt21_cols) == 0) {
  warning("No columns with 'ZNT2-21' found in Matrix B. Merging only Matrix A.")
  mat_combined <- mat_ord_A
} else {
  mat_combined <- cbind(mat_ord_A, mat_ord_B[, znt21_cols, drop = FALSE])
}

# --- 2. Final Sorting Strategy ---
# Prepare annotation for sorting
ann_for_sort <- ann_col_combined
ann_for_sort$col_id <- rownames(ann_for_sort)

# Logic: Sort final_cluster numerically (e.g. cluster1, cluster2, cluster10)
# and place 'non-conserved' at the end.
all_clusters <- unique(ann_for_sort$final_cluster)
cluster_names <- all_clusters[all_clusters != "non-conserved"]
sorted_cluster_names <- cluster_names[order(as.numeric(str_remove(cluster_names, "cluster")))]
final_cluster_levels <- c(sorted_cluster_names, "non-conserved")

ann_for_sort$final_cluster <- factor(ann_for_sort$final_cluster, levels = final_cluster_levels)

# Logic: Sort by Motif Annotation Group
display_group_levels <- c("DLX", "ETS", "ZIP6", "ZNT2")
ann_for_sort$motif_annotation <- factor(ann_for_sort$motif_annotation, levels = display_group_levels)

# Execute Sort: Group -> Final Cluster
ann_sorted <- ann_for_sort %>%
  arrange(motif_annotation, final_cluster)

# Reorder Matrix and Annotation to match sorted order
new_col_order   <- ann_sorted$col_id
mat_final       <- mat_combined[, new_col_order]
ann_final       <- ann_col_combined[new_col_order, , drop = FALSE]

# --- 3. Plotting ---
ann_colors <- list(
  nc_flag = c("non-conserved" = "black", "other" = "white")
)

p <- pheatmap(
  mat_final,
  color = colorRampPalette(brewer.pal(9, "Reds"))(100),
  breaks = seq(0, 1, length.out = 101),
  legend_breaks = c(0, 0.5, 1),
  annotation_col = ann_final,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "grey85",
  border_color = NA,
  main = "Sorted by Group -> Final Cluster"
)

# Save plot
ggsave(p$gtable, filename = "0219-revision.pdf", width = 10, height = 6)


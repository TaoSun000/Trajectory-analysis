#run test for a work at Dr. Wang Lab.
#task 2 on April 17
# trajectory analysis
# continue from Ubuntu results of all_samples.counts.txt and SraRunTable.csv.

# -------------------------------------------------------------
# Step 1. Setup and Package Installation
# -------------------------------------------------------------

# Install Bioconductor version 3.18
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Install monocle3 from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("cole-trapnell-lab/monocle3")

# Load required packages
library(monocle3)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)
library(Matrix)

# Set working directory
setwd("E:/UBC_wang_qn2")

# -------------------------------------------------------------
# Step 2. Load Count Matrix and Metadata
# -------------------------------------------------------------

counts <- read.delim("all_samples.counts.txt", comment.char = "#", row.names = 1, check.names = FALSE)
colnames(counts) <- gsub("\\.sorted\\.bam$", "", colnames(counts))
counts <- counts[, !colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length")]

meta <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE) %>%
  filter(!duplicated(Run)) %>%
  filter(Run %in% colnames(counts))

meta <- meta[match(colnames(counts), meta$Run), ]
rownames(meta) <- meta$Run
stopifnot(all(colnames(counts) == rownames(meta)))

# -------------------------------------------------------------
# Step 3. Create cell_data_set Object
# -------------------------------------------------------------

gene_metadata <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts)
)

cds <- new_cell_data_set(
  expression_data = as.matrix(counts),
  cell_metadata = meta,
  gene_metadata = gene_metadata
)

# -------------------------------------------------------------
# Step 4. Preprocessing and Dimension Reduction
# -------------------------------------------------------------

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)

# -------------------------------------------------------------
# Step 5. Clustering and Trajectory Inference
# -------------------------------------------------------------

cds <- cluster_cells(cds, resolution = 1e-2)
cds <- learn_graph(cds)
cds <- order_cells(cds)

# -------------------------------------------------------------
# Step 6. Plot Cells
# -------------------------------------------------------------

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE,
           label_leaves = TRUE, label_branch_points = TRUE) %>%
  ggsave("plot_by_pseudotime.png", path = "E:/UBC_wang_qn2", width = 6, height = 6)

plot_cells(cds, color_cells_by = "tissue") %>%
  ggsave("plot_by_tissue.png", path = "E:/UBC_wang_qn2")

plot_cells(cds, color_cells_by = "AGE") %>%
  ggsave("plot_by_age.png", path = "E:/UBC_wang_qn2")

plot_cells(cds, color_cells_by = "sex") %>%
  ggsave("plot_by_sex.png", path = "E:/UBC_wang_qn2")

plot_cells(cds, color_cells_by = "cluster") %>%
  ggsave("plot_by_cluster.png", path = "E:/UBC_wang_qn2")

# -------------------------------------------------------------
# Step 7. Differential Expression Analysis Over Pseudotime
# -------------------------------------------------------------

deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph") %>%
  arrange(q_value)

write.csv(deg_pseudotime, "E:/UBC_wang_qn2/deg_pseudotime.csv")

# -------------------------------------------------------------
# Step 8. Visualize Top Genes and Heatmap
# -------------------------------------------------------------

sig_genes <- subset(deg_pseudotime, q_value < 0.01)
top_genes <- rownames(sig_genes[order(sig_genes$morans_I, decreasing = TRUE)[1:50], ])

expr_matrix <- t(exprs(cds[top_genes, ]))
expr_matrix <- expr_matrix[order(pseudotime(cds)), ]
expr_scaled <- t(scale(t(expr_matrix)))

png("E:/UBC_wang_qn2/top50_pseudotime_heatmap.png", width = 1200, height = 800)
heatmap_result <- pheatmap(expr_scaled, cluster_cols = FALSE, show_rownames = TRUE)
dev.off()

clusters <- cutree(heatmap_result$tree_row, k = 4)
gene_cluster_df <- data.frame(gene_id = names(clusters), cluster = clusters)
write.csv(gene_cluster_df, "E:/UBC_wang_qn2/top50_gene_clusters.csv", row.names = FALSE)

# -------------------------------------------------------------
# Step 9. GO Enrichment Analysis
# -------------------------------------------------------------

all_sig_genes <- rownames(sig_genes)
bitr_result <- bitr(all_sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_results <- enrichGO(
  gene = bitr_result$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(as.data.frame(go_results), "E:/UBC_wang_qn2/go_enrichment_results.csv", row.names = FALSE)

# GO Summary
library(dplyr)
go_top_summary <- go_results %>%
  as.data.frame() %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(GeneRatio = Count / as.numeric(sub("/.*", "", BgRatio))) %>%
  select(Description, p.adjust, Count, GeneRatio, geneID)

write.csv(go_top_summary, "E:/UBC_wang_qn2/go_top20_summary_table.csv", row.names = FALSE)

message("Trajectory analysis completed and key results saved.")

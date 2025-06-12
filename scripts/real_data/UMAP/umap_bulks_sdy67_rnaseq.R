# umap for the predicted bulks from harp and compare it to real and predicted from sortedrna-seq
devtools::document()
library(tidyverse)
library(ggplot2)
library(SRSU)
library(DTD)
library(stringr)
zimmermann_data <- readRDS("data/source/sdy67_rnaseq_ref.rds")

bulk.pheno.all <- zimmermann_data$bulk_pheno
bulk.counts.all <- zimmermann_data$bulk_counts

# Define output directory and data name
output_dir <- "data/generated/real_data/sdy67_rnaseq_reference_150train_100test_5run"

if (!dir.exists(output_dir)) {
  stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}
# Create a directory for statistical analysis
stats_dir <- file.path(output_dir, "statistical_analysis")

if (!dir.exists(output_dir)) {
  stop("Error: Directory does not exist: ", output_dir, ". You need to run the benchmark first.")
}

irun <- 5 # this is the closest run to the average — it is the most representative run
harp <- readRDS(file.path(output_dir, paste0("output_harp_run_", irun, ".rds")))
bulk_count_test <- readRDS(file.path(output_dir, paste0("bulk_counts_test_run_", irun, ".rds")))
harp_ref <- harp$reference_profiles$estimated_reference_second
test.samples <- colnames(bulk_count_test)
all.samples <- colnames(bulk.counts.all)
train.samples <- setdiff(all.samples, test.samples) # this is only when we use 150 samples, and 100 test
bulk_pheno_test <- bulk.pheno.all[, test.samples]

celltype <- colnames(harp_ref)
gene_names <- rownames(harp_ref)

true_bulk_expression_test <- bulk.counts.all[gene_names, test.samples]

colnames(true_bulk_expression_test) <- paste0("test_", colnames(true_bulk_expression_test))
true_bulk_expression_train <- bulk.counts.all[gene_names, train.samples]
colnames(true_bulk_expression_train) <- paste0("train_", colnames(true_bulk_expression_train))
# true bulk_rna expresiion of test samples
predicted_bulk_expression_harp_test <- harp_ref %*% bulk_pheno_test[celltype, test.samples]

predicted_bulk_expression_harp_train <- harp_ref %*% bulk.pheno.all[celltype, train.samples]
colnames(predicted_bulk_expression_harp_test) <- paste0("predicted_test", colnames(predicted_bulk_expression_harp_test))

colnames(predicted_bulk_expression_harp_train) <- paste0("predicted_train", colnames(predicted_bulk_expression_harp_train))


reference.profile <- compute_reference_harp(sc_library = zimmermann_data$sc_library)
reference.profile <- reference.profile[gene_names, celltype]
predicted_bulks_sc <- reference.profile %*% bulk.pheno.all[celltype, colnames(bulk.pheno.all)]
colnames(predicted_bulks_sc) <- paste0("sc_", colnames(predicted_bulks_sc))

bulk.df <- cbind(
  true_bulk_expression_test, true_bulk_expression_train,
  predicted_bulk_expression_harp_test, predicted_bulk_expression_harp_train, predicted_bulks_sc
)


# UMAP Embedding
set.seed(12)
message("Starting UMAP embedding...")
umap.embedded.bulks <- uwot::umap(
  X = t(bulk.df),
  metric = "correlation",
  ret_model = TRUE,
  n_neighbors = 20,
  min_dist = 2,
  spread = 2
)
umap.frame <- data.frame(
  "umap1" = umap.embedded.bulks$embedding[, 1],
  "umap2" = umap.embedded.bulks$embedding[, 2],
  "samples" = colnames(bulk.df)
)

umap.frame$origin <- ifelse(grepl("predicted_", rownames(umap.frame)), "Predicted Bulk Expression (Harp)",
  ifelse(grepl("sc", rownames(umap.frame)), "Predicted Bulk Expression (sorted RNA-seq)",
    "PBMC (Zimmermann et al.)"
  )
)

umap.frame$group <- ifelse(grepl("test", rownames(umap.frame)), "test", "train")

umap_bulk <- ggplot(data = umap.frame, aes(
  x = umap1, y = umap2, col = origin,
  alpha = ifelse(origin == "Predicted Bulk Expression (sorted RNA-seq)", 1, group)
)) +
  geom_point(size = 2) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.position = c(0.7, 0.95)
  ) +
  labs(
    title = "",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Bulk samples",
    alpha = "Group samples"
  ) +
  scale_color_manual(
    values = alpha(c("#D16103", "#52854C", "#882255"))
  ) +
  scale_alpha_manual(
    values = c("test" = 1, "train" = 0.5),
    labels = c("Test Set", "Train Set")
  ) +
  guides(alpha = "none") +
  guides(color = guide_legend(override.aes = list(size = 2)))

# save the plot
ggsave(file.path(stats_dir, "umap_zimmerman_monaco_harp_predicted.pdf"),
  plot = umap_bulk, width = 11, height = 11
)
